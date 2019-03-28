// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <tuple>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "DataStructures/Variables.hpp"
#include "Domain/Tags.hpp"
#include "Parallel/ConstGlobalCache.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "Utilities/TypeTraits.hpp"

namespace dg {
namespace Actions {

/// \ingroup ActionsGroup
/// \ingroup DiscontinuousGalerkinGroup
/// \brief Compute \f$\hat{n} \cdot F\f$ on the boundaries for a
/// non-conservative system.
///
/// Uses:
/// - System:
///   - `variables_tag`
/// - DataBox:
///   - `DirectionsTag`,
///   - Interface items as required by
///   `Metavariables::system::normal_dot_fluxes`
///
/// DataBox changes:
/// - Adds: nothing
/// - Removes: nothing
/// - Modifies:
///   - `Tags::Interface<
///   DirectionsTag, db::add_tag_prefix<Tags::NormalDotFlux, variables_tag>>`
template <typename DirectionsTag, typename F>
struct ComputeFluxes {
 private:
  template <typename ReturnTagsList, typename ArgumentTagsList>
  struct apply_flux_impl;
  template <typename... ReturnTags, typename... ArgumentTags>
  struct apply_flux_impl<tmpl::list<ReturnTags...>,
                         tmpl::list<ArgumentTags...>> {
    static void apply(
        const gsl::not_null<db::item_type<ReturnTags>*>... return_tags,
        const db::item_type<ArgumentTags>&... argument_tags,
        const db::item_type<DirectionsTag>& directions) {
      for (const auto& direction : directions) {
        F::apply(make_not_null(&return_tags->at(direction))...,
                 argument_tags.at(direction)...);
      }
    }
  };

 public:
  template <typename DbTags, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTags>&&> apply(
      db::DataBox<DbTags>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    using interface_argument_tags =
        tmpl::transform<typename F::argument_tags,
                        tmpl::bind<Tags::Interface, DirectionsTag, tmpl::_1>>;
    using interface_return_tags =
        tmpl::transform<typename F::return_tags,
                        tmpl::bind<Tags::Interface, DirectionsTag, tmpl::_1>>;
    db::mutate_apply<interface_return_tags, interface_argument_tags>(
        apply_flux_impl<interface_return_tags, interface_argument_tags>{},
        make_not_null(&box), get<DirectionsTag>(box));
    return std::forward_as_tuple(std::move(box));
  }
};
}  // namespace Actions
}  // namespace dg
