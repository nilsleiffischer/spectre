// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Evolution/Systems/Cce/OptionTags.hpp"
#include "ParallelAlgorithms/Initialization/MergeIntoDataBox.hpp"
#include "Time/Tags.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/Rational.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Cce {
namespace Actions {

/*!
 * \ingroup ActionsGroup
 * \brief Initializes the contents of the `CharacteristicEvolution` component
 * for performing the time evolution of the system, which is the singleton that
 * handles the main evolution system for CCE computations.
 *
 * \details Sets up the \ref DataBoxGroup to be ready to perform the
 * time-stepping associated with the CCE system.
 *
 * \ref DataBoxGroup changes:
 * - Modifies: nothing
 * - Adds:
 *  - `Tags::TimeStepId`
 *  - `Tags::Next<Tags::TimeStepId>`
 *  - `Tags::TimeStep`
 *  - `Tags::Time`
 *  -
 * ```
 * Tags::HistoryEvolvedVariables<
 * metavariables::evolved_coordinates_variables_tag,
 * db::add_tag_prefix<Tags::dt,
 * metavariables::evolved_coordinates_variables_tag>>
 * ```
 *  -
 * ```
 * Tags::HistoryEvolvedVariables<
 * ::Tags::Variables<metavariables::evolved_swsh_tag>,
 * ::Tags::Variables<metavariables::evolved_swsh_dt_tag>>
 * ```
 * - Removes: nothing
 *
 * \note This action relies on the `SetupDataBox` aggregated initialization
 * mechanism, so `Actions::SetupDataBox` must be present in the `Initialization`
 * phase action list prior to this action.
 */
template <typename Metavariables>
struct InitializeCharacteristicEvolutionTime {
  using initialization_tags = tmpl::list<InitializationTags::TargetStepSize>;
  using const_global_cache_tags =
      tmpl::list<::Tags::TimeStepper<TimeStepper>>;

  using coordinate_variables_tag =
      typename Metavariables::evolved_coordinates_variables_tag;
  using evolved_swsh_variables_tag =
      ::Tags::Variables<tmpl::list<typename Metavariables::evolved_swsh_tag>>;
  using simple_tags =
      tmpl::list<::Tags::TimeStepId, ::Tags::Next<::Tags::TimeStepId>,
                 ::Tags::TimeStep, ::Tags::Time,
                 ::Tags::HistoryEvolvedVariables<coordinate_variables_tag>,
                 ::Tags::HistoryEvolvedVariables<evolved_swsh_variables_tag>>;
  using compute_tags = tmpl::list<::Tags::SubstepTimeCompute>;

  template <typename DbTags, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTags>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    const double initial_time_value = db::get<Tags::StartTime>(box);
    const double step_size = db::get<InitializationTags::TargetStepSize>(box);

    const Slab single_step_slab{initial_time_value,
                                initial_time_value + step_size};
    const Time initial_time = single_step_slab.start();
    const TimeDelta fixed_time_step =
        TimeDelta{single_step_slab, Rational{1, 1}};
    TimeStepId initial_time_id{true, 0, initial_time};
    const auto& time_stepper = db::get<::Tags::TimeStepper<TimeStepper>>(box);
    TimeStepId second_time_id =
        time_stepper.next_time_id(initial_time_id, fixed_time_step);

    typename ::Tags::HistoryEvolvedVariables<coordinate_variables_tag>::type
        coordinate_history;

    typename ::Tags::HistoryEvolvedVariables<evolved_swsh_variables_tag>::type
        swsh_history;
    db::mutate_assign(make_not_null(&box), simple_tags{},
                      std::move(initial_time_id), std::move(second_time_id),
                      fixed_time_step, initial_time_value,
                      std::move(coordinate_history), std::move(swsh_history));
    return std::make_tuple(std::move(box));
  }
};

}  // namespace Actions
}  // namespace Cce
