// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "tests/Unit/TestingFramework.hpp"

#include <array>
#include <cstddef>
#include <functional>
#include <initializer_list>  // IWYU pragma: keep
#include <memory>
#include <pup.h>
#include <string>
#include <unordered_map>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"  // IWYU pragma: keep
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Direction.hpp"
#include "Domain/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Actions/ImposeBoundaryConditions.hpp"
#include "Parallel/PhaseDependentActionList.hpp"  // IWYU pragma: keep
#include "Time/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "tests/Unit/ActionTesting.hpp"
#include "tests/Unit/TestHelpers.hpp"

// IWYU pragma: no_include <boost/functional/hash/extensions.hpp>

// IWYU pragma: no_include "DataStructures/VariablesHelpers.hpp"  // for Variables
// IWYU pragma: no_include "NumericalAlgorithms/DiscontinuousGalerkin/SimpleMortarData.hpp"
// IWYU pragma: no_include "Parallel/PupStlCpp11.hpp"

// IWYU pragma: no_forward_declare ActionTesting::InitializeDataBox
// IWYU pragma: no_forward_declare Tensor
// IWYU pragma: no_forward_declare Variables
// IWYU pragma: no_forward_declare dg::Actions::ImposeDirichletBoundaryConditions

namespace {
constexpr size_t Dim = 2;

struct Var : db::SimpleTag {
  static std::string name() noexcept { return "Var"; }
  using type = Scalar<DataVector>;
};

struct PrimitiveVar : db::SimpleTag {
  static std::string name() noexcept { return "PrimitiveVar"; }
  using type = Scalar<DataVector>;
};

struct BoundaryCondition {
  tuples::TaggedTuple<Var> variables(const tnsr::I<DataVector, Dim>& /*x*/,
                                     double /*t*/,
                                     tmpl::list<Var> /*meta*/) const noexcept {
    return tuples::TaggedTuple<Var>{Scalar<DataVector>{{{{30., 40., 50.}}}}};
  }

  tuples::TaggedTuple<PrimitiveVar> variables(
      const tnsr::I<DataVector, Dim>& /*x*/, double /*t*/,
      tmpl::list<PrimitiveVar> /*meta*/) const noexcept {
    return tuples::TaggedTuple<PrimitiveVar>{
        Scalar<DataVector>{{{{15., 20., 25.}}}}};
  }
  // clang-tidy: do not use references
  void pup(PUP::er& /*p*/) noexcept {}  // NOLINT
};

struct BoundaryConditionTag {
  using type = BoundaryCondition;
};

template <bool HasPrimitiveAndConservativeVars>
struct System {
  static constexpr const size_t volume_dim = Dim;
  static constexpr bool is_in_flux_conservative_form = true;
  static constexpr bool has_primitive_and_conservative_vars =
      HasPrimitiveAndConservativeVars;

  using variables_tag = Tags::Variables<tmpl::list<Var>>;

  struct conservative_from_primitive {
    using return_tags = tmpl::list<Var>;
    using argument_tags = tmpl::list<PrimitiveVar>;

    static void apply(const gsl::not_null<Scalar<DataVector>*> var,
                      const Scalar<DataVector>& primitive_var) {
      get(*var) = 2.0 * get(primitive_var);
    }
  };
};

using exterior_bdry_vars_tag =
    Tags::Interface<Tags::BoundaryDirectionsExterior<Dim>,
                    Tags::Variables<tmpl::list<Var>>>;

template <typename Metavariables>
struct component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using const_global_cache_tags = tmpl::list<BoundaryConditionTag>;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<db::AddSimpleTags<
              Tags::Time,
              Tags::Interface<Tags::BoundaryDirectionsExterior<Dim>,
                              Tags::Coordinates<Dim, Frame::Inertial>>,
              exterior_bdry_vars_tag>>>>,
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Testing,
          tmpl::list<
              dg::Actions::ImposeDirichletBoundaryConditions<Metavariables>>>>;
};

template <bool HasPrimitiveAndConservativeVars>
struct Metavariables {
  using system = System<HasPrimitiveAndConservativeVars>;
  using component_list = tmpl::list<component<Metavariables>>;

  using boundary_condition_tag = BoundaryConditionTag;
  using initial_data_tag = boundary_condition_tag;
  enum class Phase { Initialization, Testing, Exit };
};

template <bool HasConservativeAndPrimitiveVars>
void run_test() {
  using metavariables = Metavariables<HasConservativeAndPrimitiveVars>;
  using my_component = component<metavariables>;
  // Just making up two "external" directions
  const auto external_directions = {Direction<2>::lower_eta(),
                                    Direction<2>::upper_xi()};

  ActionTesting::MockRuntimeSystem<metavariables> runner{{BoundaryCondition{}}};
  {
    tnsr::I<DataVector, Dim> arbitrary_coords{
        DataVector{3, std::numeric_limits<double>::signaling_NaN()}};
    db::item_type<Tags::Interface<Tags::BoundaryDirectionsExterior<Dim>,
                                  Tags::Coordinates<Dim, Frame::Inertial>>>
        external_bdry_coords{{{Direction<2>::lower_eta(), arbitrary_coords},
                              {Direction<2>::upper_xi(), arbitrary_coords}}};
    db::item_type<exterior_bdry_vars_tag> exterior_bdry_vars;
    for (const auto& direction : external_directions) {
      exterior_bdry_vars[direction].initialize(3);
    }

    ActionTesting::emplace_component_and_initialize<my_component>(
        &runner, 0,
        {1.2, std::move(external_bdry_coords),
         std::move(exterior_bdry_vars)});
  }
  runner.set_phase(metavariables::Phase::Testing);

  ActionTesting::next_action<my_component>(make_not_null(&runner), 0);

  // Check that BC's were indeed applied.
  const auto& external_vars =
      ActionTesting::get_databox_tag<my_component, exterior_bdry_vars_tag>(
          runner, 0);
  db::item_type<exterior_bdry_vars_tag> expected_vars{};
  for (const auto& direction : external_directions) {
    expected_vars[direction].initialize(3);
  }
  get<Var>(expected_vars[Direction<2>::lower_eta()]) =
      Scalar<DataVector>({{{30., 40., 50.}}});
  get<Var>(expected_vars[Direction<2>::upper_xi()]) =
      Scalar<DataVector>({{{30., 40., 50.}}});
  CHECK(external_vars == expected_vars);
}

}  // namespace

SPECTRE_TEST_CASE("Unit.DiscontinuousGalerkin.Actions.BoundaryConditions",
                  "[Unit][NumericalAlgorithms][Actions]") {
  run_test<false>();
  run_test<true>();
}
