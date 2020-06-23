// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <vector>

#include "Domain/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/MeshHierarchy.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/Tags.hpp"
#include "Utilities/TaggedTuple.hpp"

#include "Parallel/Printf.hpp"

namespace LinearSolver::multigrid {

template <size_t Dim, typename OptionsGroup>
struct ElementsAllocator {
  template <typename ElementArray>
  using array_allocation_tags =
      tmpl::list<domain::Tags::InitialRefinementLevels<Dim>,
                 Tags::BaseRefinementLevels<Dim>,
                 Tags::ParentRefinementLevels<Dim>,
                 domain::Tags::InitialExtents<Dim>, Tags::ParentExtents<Dim>,
                 Tags::ArraySection<Tags::MultigridLevel, ElementArray>,
                 Tags::MultigridLevel, Tags::CoarsestGridPoints<OptionsGroup>>;

  template <typename ElementArray, typename Metavariables,
            typename... InitializationTags>
  static void apply(
      Parallel::CProxy_ConstGlobalCache<Metavariables>& global_cache,
      tuples::TaggedTuple<InitializationTags...>
          initialization_items) noexcept {
    auto& local_cache = *(global_cache.ckLocalBranch());
    auto& element_array =
        Parallel::get_parallel_component<ElementArray>(local_cache);
    const auto& domain = Parallel::get<domain::Tags::Domain<Dim>>(local_cache);
    auto& initial_refinement_levels =
        get<domain::Tags::InitialRefinementLevels<Dim>>(initialization_items);
    auto& initial_extents =
        get<domain::Tags::InitialExtents<Dim>>(initialization_items);
    auto& base_refinement_levels =
        get<Tags::BaseRefinementLevels<Dim>>(initialization_items);
    auto& parent_refinement_levels =
        get<Tags::ParentRefinementLevels<Dim>>(initialization_items);
    auto& parent_extents =
        get<Tags::ParentExtents<Dim>>(initialization_items);
    auto& multigrid_level = get<Tags::MultigridLevel>(initialization_items);
    do {
      // Store the grid as base before coarsening it
      base_refinement_levels = initial_refinement_levels;
      initial_refinement_levels = parent_refinement_levels;
      initial_extents = parent_extents;
      // Construct coarsened (parent) grid
      std::tie(parent_refinement_levels, parent_extents) =
          LinearSolver::multigrid::coarsen(
              initial_refinement_levels, initial_extents,
              get<Tags::CoarsestGridPoints<OptionsGroup>>(
                  initialization_items));
      // Create element IDs for all elements on this level
      std::vector<ElementId<Dim>> all_element_ids{};
      for (const auto& block : domain.blocks()) {
        const std::vector<ElementId<Dim>> element_ids = initial_element_ids(
            block.id(), initial_refinement_levels[block.id()], multigrid_level);
        all_element_ids.insert(all_element_ids.begin(), element_ids.begin(),
                               element_ids.end());
      }
      // Create an array section for this refinement level
      std::vector<CkArrayIndex> all_array_indices(all_element_ids.size());
      std::transform(all_element_ids.begin(), all_element_ids.end(),
                     all_array_indices.begin(),
                     [](const ElementId<Dim>& element_id) noexcept {
                       return Parallel::ArrayIndex<ElementId<Dim>>(element_id);
                     });
      get<Tags::ArraySection<Tags::MultigridLevel, ElementArray>>(
          initialization_items) =
          CProxySection_AlgorithmArray<ElementArray,
                                       typename ElementArray::array_index>::
              ckNew(element_array.ckGetArrayID(), all_array_indices.data(),
                    all_array_indices.size());
      // Create the elements for this refinement level and distribute them among
      // processors
      const int number_of_procs = Parallel::number_of_procs();
      for (size_t i = 0; i < all_element_ids.size(); ++i) {
        element_array(all_element_ids[i])
            .insert(global_cache, initialization_items,
                    static_cast<int>(i) % number_of_procs);
      }
    //   Parallel::printf("Elements on level %zu: %s\n", multigrid_level,
    //                    all_element_ids);
    //   Parallel::printf("Grid points on this level: %s\n",
    //                    initial_extents);
      ++multigrid_level;
    } while (initial_refinement_levels != parent_refinement_levels or
             initial_extents != parent_extents);
    element_array.doneInserting();
  }
};  // namespace LinearSolver::multigrid

}  // namespace LinearSolver::multigrid