// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>

#include "DataStructures/DataBox/Tag.hpp"
#include "Framework/ActionTesting.hpp"
#include "Parallel/Actions/Goto.hpp"  // IWYU pragma: keep
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"  // IWYU pragma: keep
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace {
struct Label1;
struct Label2;
struct Label3;

template <typename Label>
struct Counter : db::SimpleTag {
  using type = size_t;
};

template <typename Label>
struct HasConverged : db::ComputeTag {
  using argument_tags = tmpl::list<Counter<Label>>;
  static bool function(const size_t counter) noexcept { return counter >= 2; }
  static std::string name() noexcept { return "HasConverged"; }
};

template <typename Label>
struct Increment {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    db::mutate<Counter<Label>>(
        make_not_null(&box),
        [](const gsl::not_null<size_t*> counter) noexcept { (*counter)++; });
    return {std::move(box)};
  }
};

template <typename Label>
struct Reset {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    db::mutate<Counter<Label>>(
        make_not_null(&box),
        [](const gsl::not_null<size_t*> counter) noexcept { *counter = 0; });
    return {std::move(box)};
  }
};

/// [component]
template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;

  using repeat_until_phase_action_list = tmpl::flatten<tmpl::list<
      Actions::RepeatUntil<
          HasConverged<Label1>,
          tmpl::list<
              Increment<Label1>, Reset<Label2>,
              Actions::RepeatUntil<
                  HasConverged<Label2>,
                  tmpl::list<Increment<Label2>, Increment<Label3>>, Label2>>,
          Label1>,
      Parallel::Actions::TerminatePhase>>;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::Initialization,
          tmpl::list<ActionTesting::InitializeDataBox<
              db::AddSimpleTags<Counter<Label1>, Counter<Label2>,
                                Counter<Label3>>,
              db::AddComputeTags<HasConverged<Label1>, HasConverged<Label2>>>>>,

      Parallel::PhaseActions<
          typename Metavariables::Phase, Metavariables::Phase::TestGoto,
          tmpl::list<Actions::Goto<Label1>, Actions::Label<Label2>,
                     Actions::Label<Label1>, Actions::Goto<Label2>>>,

      Parallel::PhaseActions<typename Metavariables::Phase,
                             Metavariables::Phase::TestRepeatUntil,
                             repeat_until_phase_action_list>>;
};
/// [component]

/// [metavariables]
struct Metavariables {
  using component_list = tmpl::list<Component<Metavariables>>;

  enum class Phase { Initialization, TestGoto, TestRepeatUntil, Exit };
};
/// [metavariables]
}  // namespace

/// [test case]
SPECTRE_TEST_CASE("Unit.Parallel.GotoAction", "[Unit][Parallel][Actions]") {
  using component = Component<Metavariables>;
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<Metavariables>;
  MockRuntimeSystem runner{{}};
  ActionTesting::emplace_component_and_initialize<component>(
      &runner, 0, {size_t{0}, size_t{0}, size_t{0}});

  ActionTesting::set_phase(make_not_null(&runner),
                           Metavariables::Phase::TestGoto);
  runner.force_next_action_to_be<component, Actions::Label<Label1>>(0);
  runner.next_action<component>(0);
  CHECK(runner.get_next_action_index<component>(0) == 3);

  runner.force_next_action_to_be<component, Actions::Goto<Label1>>(0);
  runner.next_action<component>(0);
  CHECK(runner.get_next_action_index<component>(0) == 2);

  runner.force_next_action_to_be<component, Actions::Goto<Label2>>(0);
  runner.next_action<component>(0);
  CHECK(runner.get_next_action_index<component>(0) == 1);

  ActionTesting::set_phase(make_not_null(&runner),
                           Metavariables::Phase::TestRepeatUntil);
  while (not ActionTesting::get_terminate<component>(runner, 0)) {
    runner.next_action<component>(0);
  }
  CHECK(ActionTesting::get_databox_tag<component, HasConverged<Label1>>(runner,
                                                                        0));
  CHECK(ActionTesting::get_databox_tag<component, HasConverged<Label2>>(runner,
                                                                        0));
  CHECK(ActionTesting::get_databox_tag<component, Counter<Label1>>(runner, 0) ==
        2);
  CHECK(ActionTesting::get_databox_tag<component, Counter<Label2>>(runner, 0) ==
        2);
  CHECK(ActionTesting::get_databox_tag<component, Counter<Label3>>(runner, 0) ==
        4);

  // Test zero iterations of the `RepeatUntil` loop
  ActionTesting::set_phase(make_not_null(&runner),
                           Metavariables::Phase::TestRepeatUntil);
  runner.next_action<component>(0);
  CHECK(runner.get_next_action_index<component>(0) ==
        tmpl::index_of<typename component::repeat_until_phase_action_list,
                       Parallel::Actions::TerminatePhase>::value);
  // Make sure `Increment` was not called for this situation where the
  // condition is already fulfilled at the start.
  CHECK(ActionTesting::get_databox_tag<component, Counter<Label1>>(runner, 0) ==
        2);
  CHECK(ActionTesting::get_databox_tag<component, Counter<Label2>>(runner, 0) ==
        2);
  CHECK(ActionTesting::get_databox_tag<component, Counter<Label3>>(runner, 0) ==
        4);
}
/// [test case]
