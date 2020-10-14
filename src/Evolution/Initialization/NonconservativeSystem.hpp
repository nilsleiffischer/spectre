// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>  // IWYU pragma: keep  // for move

#include "DataStructures/DataBox/DataBox.hpp"
#include "Evolution/Initialization/InitialData.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Initialization/MergeIntoDataBox.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace domain {
namespace Tags {
template <size_t VolumeDim>
struct Mesh;
}  // namespace Tags
}  // namespace domain
// IWYU pragma: no_forward_declare db::DataBox
/// \endcond

namespace Initialization {
namespace Actions {
/// \ingroup InitializationGroup
/// \brief Allocate variables needed for evolution of nonconservative systems
///
/// Uses:
/// - DataBox:
///   * `Tags::Mesh<Dim>`
///
/// DataBox changes:
/// - Adds:
///   * System::variables_tag
///
/// - Removes: nothing
/// - Modifies: nothing
///
/// \note This action relies on the `SetupDataBox` aggregated initialization
/// mechanism, so `Actions::SetupDataBox` must be present in the
/// `Initialization` phase action list prior to this action.
template <typename Metavariables>
struct NonconservativeSystem {
  using system = typename Metavariables::system;
  static_assert(not system::is_in_flux_conservative_form,
                "System is in flux conservative form");
  static constexpr size_t dim = system::volume_dim;
  using variables_tag = typename system::variables_tag;
  using simple_tags = db::AddSimpleTags<variables_tag>;
  using compute_tags = db::AddComputeTags<>;

  template <typename DbTagsList, typename... InboxTags, typename ArrayIndex,
            typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/, ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    using Vars = typename variables_tag::type;
    db::mutate_assign(
        make_not_null(&box), simple_tags{},
        Vars{db::get<domain::Tags::Mesh<dim>>(box).number_of_grid_points()});

    return std::make_tuple(std::move(box));
  }
};
}  // namespace Actions
}  // namespace Initialization
