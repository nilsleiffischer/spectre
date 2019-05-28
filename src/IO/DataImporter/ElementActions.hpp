// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/ElementId.hpp"
#include "IO/DataImporter/DataFileReader.hpp"
#include "IO/Observer/ArrayComponentId.hpp"
#include "Parallel/ConstGlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Utilities/MakeString.hpp"
#include "Utilities/Requires.hpp"

namespace importer {
namespace Actions {

struct RegisterElementWithSelf;

/*!
 * \brief Register an element with the `importer::DataFileReader` component.
 *
 * Invoke this action on each element of an array parallel component to register
 * them with the `importer::DataFileReader`.
 */
struct RegisterWithImporter {
  template <typename... DbTags, typename... InboxTags, typename Metavariables,
            size_t Dim, typename ActionList, typename ParallelComponent,
            Requires<sizeof...(DbTags) != 0> = nullptr>
  static auto apply(db::DataBox<tmpl::list<DbTags...>>& /*box*/,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    Parallel::ConstGlobalCache<Metavariables>& cache,
                    const ElementIndex<Dim>& array_index,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    const std::string element_name = MakeString{}
                                     << ElementId<Dim>(array_index);
    auto& local_reader_component =
        *Parallel::get_parallel_component<
             importer::DataFileReader<Metavariables>>(cache)
             .ckLocalBranch();
    Parallel::simple_action<importer::Actions::RegisterElementWithSelf>(
        local_reader_component,
        observers::ArrayComponentId(
            std::add_pointer_t<ParallelComponent>{nullptr},
            Parallel::ArrayIndex<ElementIndex<Dim>>(array_index)),
        element_name);
  }
};

}  // namespace Actions
}  // namespace importer
