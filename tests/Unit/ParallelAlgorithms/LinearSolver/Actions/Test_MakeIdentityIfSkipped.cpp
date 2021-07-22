// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DenseVector.hpp"
#include "Framework/ActionTesting.hpp"
#include "NumericalAlgorithms/Convergence/Tags.hpp"
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/LinearSolver/Actions/MakeIdentityIfSkipped.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Literals.hpp"
#include "Utilities/TMPL.hpp"

namespace {

struct TestOptionsGroup {};

struct FieldsTag : db::SimpleTag {
  using type = DenseVector<double>;
};

struct SourceTag : db::SimpleTag {
  using type = DenseVector<double>;
};

struct CheckRunIfSkippedTag : db::SimpleTag {
  using type = bool;
};

struct TestLinearSolver {
  using options_group = TestOptionsGroup;
  using fields_tag = FieldsTag;
  using source_tag = SourceTag;
};

struct RunIfSkipped {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    db::mutate<CheckRunIfSkippedTag>(
        make_not_null(&box),
        [](const gsl::not_null<bool*> flag) { *flag = true; });
    return {std::move(box)};
  }
};

template <typename Metavariables>
struct ElementArray {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<
              tmpl::list<Convergence::Tags::HasConverged<TestOptionsGroup>,
                         FieldsTag, SourceTag, CheckRunIfSkippedTag>>>>,
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Testing,
          tmpl::list<LinearSolver::Actions::make_identity_if_skipped<
                         TestLinearSolver, RunIfSkipped>,
                     Parallel::Actions::TerminatePhase>>>;
};

struct Metavariables {
  using element_array = ElementArray<Metavariables>;
  using component_list = tmpl::list<element_array>;
  enum class Phase { Initialization, Testing, Exit };
};

void test_make_identity_if_skipped(const bool skipped) {
  CAPTURE(skipped);
  Convergence::HasConverged has_converged{skipped ? 0_st : 1_st, 0};
  using element_array = typename Metavariables::element_array;
  ActionTesting::MockRuntimeSystem<Metavariables> runner{{}};
  ActionTesting::emplace_component_and_initialize<element_array>(
      make_not_null(&runner), 0,
      {has_converged, DenseVector<double>{3, 1.}, DenseVector<double>{3, 2.},
       false});
  ActionTesting::set_phase(make_not_null(&runner),
                           Metavariables::Phase::Testing);
  while (not ActionTesting::get_terminate<element_array>(runner, 0)) {
    ActionTesting::next_action<element_array>(make_not_null(&runner), 0);
  }
  const auto get_tag = [&runner](auto tag_v) -> decltype(auto) {
    using tag = std::decay_t<decltype(tag_v)>;
    return ActionTesting::get_databox_tag<element_array, tag>(runner, 0);
  };
  CHECK(get_tag(SourceTag{}) == DenseVector<double>{3, 2.});
  CHECK(get_tag(FieldsTag{}) ==
        (skipped ? DenseVector<double>{3, 2.} : DenseVector<double>{3, 1.}));
  CHECK(get_tag(CheckRunIfSkippedTag{}) == skipped);
}

}  // namespace

SPECTRE_TEST_CASE("Unit.ParallelLinearSolver.Actions.MakeIdentityIfSkipped",
                  "[Unit][ParallelAlgorithms][LinearSolver][Actions]") {
  test_make_identity_if_skipped(false);
  test_make_identity_if_skipped(true);
}
