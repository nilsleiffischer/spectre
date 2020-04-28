// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>
#include <unordered_map>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"  // IWYU pragma: keep
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DenseVector.hpp"
#include "Framework/ActionTesting.hpp"
#include "NumericalAlgorithms/Convergence/HasConverged.hpp"
#include "ParallelAlgorithms/Actions/SetData.hpp"
#include "ParallelAlgorithms/LinearSolver/ConjugateGradient/ElementActions.hpp"  // IWYU pragma: keep
#include "ParallelAlgorithms/LinearSolver/ConjugateGradient/InitializeElement.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"  // IWYU pragma: keep
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

// IWYU pragma: no_include <boost/variant/get.hpp>
// IWYU pragma: no_forward_declare db::DataBox

namespace {

struct DummyOptionsGroup {};

struct VectorTag : db::SimpleTag {
  using type = DenseVector<double>;
};

using fields_tag = VectorTag;
using operator_applied_to_fields_tag =
    db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, fields_tag>;
using operand_tag = LinearSolver::Tags::Operand<fields_tag>;
using operator_applied_to_operand_tag =
    db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, operand_tag>;
using residual_tag = LinearSolver::Tags::Residual<fields_tag>;

template <typename Metavariables>
struct ElementArray {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<tmpl::list<VectorTag>>,
                     LinearSolver::cg_detail::InitializeElement<
                         fields_tag, DummyOptionsGroup>>>,
      Parallel::PhaseActions<typename Metavariables::Phase,
                             Metavariables::Phase::Testing,
                             tmpl::list<LinearSolver::cg_detail::PrepareStep<
                                 fields_tag, DummyOptionsGroup>>>>;
};

struct Metavariables {
  using component_list = tmpl::list<ElementArray<Metavariables>>;
  enum class Phase { Initialization, Testing, Exit };
};

}  // namespace

SPECTRE_TEST_CASE(
    "Unit.ParallelAlgorithms.LinearSolver.ConjugateGradient.ElementActions",
    "[Unit][ParallelAlgorithms][LinearSolver][Actions]") {
  using element_array = ElementArray<Metavariables>;

  ActionTesting::MockRuntimeSystem<Metavariables> runner{{}};

  // Setup mock element array
  ActionTesting::emplace_component_and_initialize<element_array>(
      make_not_null(&runner), 0, {DenseVector<double>(3, 0.)});
  ActionTesting::next_action<element_array>(make_not_null(&runner), 0);

  // DataBox shortcuts
  const auto get_tag = [&runner](auto tag_v) -> decltype(auto) {
    using tag = std::decay_t<decltype(tag_v)>;
    return ActionTesting::get_databox_tag<element_array, tag>(runner, 0);
  };
  const auto tag_is_retrievable = [&runner](auto tag_v) {
    using tag = std::decay_t<decltype(tag_v)>;
    return ActionTesting::tag_is_retrievable<element_array, tag>(runner, 0);
  };
  const auto set_tag = [&runner](auto tag_v, const auto& value) {
    using tag = std::decay_t<decltype(tag_v)>;
    ActionTesting::simple_action<element_array,
                                 ::Actions::SetData<tmpl::list<tag>>>(
        make_not_null(&runner), 0, value);
  };

  ActionTesting::set_phase(make_not_null(&runner),
                           Metavariables::Phase::Testing);

  // Can't test the other element actions because reductions are not yet
  // supported. The full algorithm is tested in
  // `Test_ConjugateGradientAlgorithm.cpp` and
  // `Test_DistributedConjugateGradientAlgorithm.cpp`.

  {
    INFO("InitializeElement");
    CHECK(get_tag(LinearSolver::Tags::IterationId<DummyOptionsGroup>{}) ==
          std::numeric_limits<size_t>::max());
    tmpl::for_each<tmpl::list<operator_applied_to_fields_tag, operand_tag,
                              operator_applied_to_operand_tag, residual_tag>>(
        [&tag_is_retrievable](auto tag_v) {
          using tag = tmpl::type_from<decltype(tag_v)>;
          CAPTURE(db::tag_name<tag>());
          CHECK(tag_is_retrievable(tag{}));
        });
    CHECK_FALSE(get_tag(LinearSolver::Tags::HasConverged<DummyOptionsGroup>{}));
  }

  SECTION("InitializeHasConverged") {
    ActionTesting::simple_action<
        element_array, LinearSolver::cg_detail::InitializeHasConverged<
                           fields_tag, DummyOptionsGroup>>(
        make_not_null(&runner), 0,
        Convergence::HasConverged{{1, 0., 0.}, 1, 0., 0.});
    CHECK(get_tag(LinearSolver::Tags::HasConverged<DummyOptionsGroup>{}));
  }
  SECTION("UpdateOperand") {
    set_tag(LinearSolver::Tags::IterationId<DummyOptionsGroup>{}, size_t{0});
    set_tag(operand_tag{}, DenseVector<double>(3, 2.));
    set_tag(residual_tag{}, DenseVector<double>(3, 1.));
    ActionTesting::next_action<element_array>(make_not_null(&runner), 0);
    ActionTesting::simple_action<
        element_array,
        LinearSolver::cg_detail::UpdateOperand<fields_tag, DummyOptionsGroup>>(
        make_not_null(&runner), 0, 2.,
        Convergence::HasConverged{{1, 0., 0.}, 1, 0., 0.});
    CHECK(get_tag(LinearSolver::Tags::Operand<VectorTag>{}) ==
          DenseVector<double>(3, 5.));
    CHECK(get_tag(LinearSolver::Tags::IterationId<DummyOptionsGroup>{}) == 1);
    CHECK(get_tag(LinearSolver::Tags::HasConverged<DummyOptionsGroup>{}));
  }
}
