// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>
#include <string>
#include <tuple>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Index.hpp"
#include "Domain/ElementId.hpp"
#include "Domain/InterfaceComputeTags.hpp"
#include "Domain/InterfaceHelpers.hpp"
#include "Domain/Mesh.hpp"
#include "Domain/OrientationMapHelpers.hpp"
#include "Domain/Tags.hpp"
#include "Informer/Verbosity.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Parallel/ConstGlobalCache.hpp"
#include "Parallel/InboxInserters.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/ParallelComponentHelpers.hpp"
#include "Parallel/Printf.hpp"
#include "ParallelAlgorithms/DiscontinuousGalerkin/HasReceivedFromAllMortars.hpp"
#include "ParallelAlgorithms/Initialization/MergeIntoDataBox.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/ElementCenteredSubdomainData.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/OverlapHelpers.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Tags.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace LinearSolver::Schwarz::detail {

template <typename SubdomainDataType, typename OptionsGroup>
struct SubdomainDataBufferTag : db::SimpleTag {
  static std::string name() noexcept {
    return "SubdomainData(" + option_name<OptionsGroup>() + ")";
  }
  using type = SubdomainDataType;
};

template <typename FieldsTag, typename OptionsGroup, typename SubdomainOperator>
struct InitializeElement {
 private:
  using fields_tag = FieldsTag;
  using residual_tag =
      db::add_tag_prefix<LinearSolver::Tags::Residual, fields_tag>;
  static constexpr size_t Dim = SubdomainOperator::volume_dim;
  using SubdomainData =
      ElementCenteredSubdomainData<Dim,
                                   db::get_variables_tags_list<residual_tag>>;
  using subdomain_solver_tag =
      Tags::SubdomainSolver<LinearSolver::Serial::Gmres<SubdomainData>,
                            OptionsGroup>;
  template <typename Tag>
  using overlaps_tag = Tags::Overlaps<Tag, Dim, OptionsGroup>;

 public:
  using initialization_tags =
      tmpl::list<domain::Tags::InitialExtents<Dim>, subdomain_solver_tag>;
  using initialization_tags_to_keep = tmpl::list<subdomain_solver_tag>;
  using const_global_cache_tags = tmpl::list<Tags::MaxOverlap<OptionsGroup>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
                    const ElementId<Dim>& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    const auto& element_mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    const size_t element_num_points = element_mesh.number_of_grid_points();

    return std::make_tuple(
        ::Initialization::merge_into_databox<
            InitializeElement,
            db::AddSimpleTags<
                SubdomainDataBufferTag<SubdomainData, OptionsGroup>>,
            db::AddComputeTags<
                domain::Tags::InternalDirections<Dim>,
                domain::Tags::BoundaryDirectionsInterior<Dim>,
                domain::Tags::InterfaceCompute<
                    domain::Tags::InternalDirections<Dim>,
                    domain::Tags::Direction<Dim>>,
                domain::Tags::InterfaceCompute<
                    domain::Tags::BoundaryDirectionsInterior<Dim>,
                    domain::Tags::Direction<Dim>>,
                domain::Tags::InterfaceCompute<
                    domain::Tags::InternalDirections<Dim>,
                    Tags::IntrudingExtentsCompute<Dim, OptionsGroup>>,
                domain::Tags::InterfaceCompute<
                    domain::Tags::InternalDirections<Dim>,
                    Tags::IntrudingOverlapWidthCompute<Dim, OptionsGroup>>,
                domain::Tags::LogicalCoordinates<Dim>,
                Tags::ElementWeightCompute<Dim, OptionsGroup>,
                domain::Tags::InterfaceCompute<
                    domain::Tags::InternalDirections<Dim>,
                    Tags::IntrudingOverlapWeightCompute<Dim, OptionsGroup>>>>(
            std::move(box), SubdomainData{element_num_points}));
  }
};

template <size_t Dim, typename OptionsGroup, typename OverlapData>
struct OverlapDataInboxTag
    : public Parallel::InboxInserters::Map<
          OverlapDataInboxTag<Dim, OptionsGroup, OverlapData>> {
  using temporal_id = size_t;
  using type = std::unordered_map<temporal_id, OverlapMap<Dim, OverlapData>>;
};

template <typename FieldsTag, typename OptionsGroup, typename SubdomainOperator>
struct SendOverlapData {
 private:
  using fields_tag = FieldsTag;
  using residual_tag =
      db::add_tag_prefix<LinearSolver::Tags::Residual, fields_tag>;
  static constexpr size_t Dim = SubdomainOperator::volume_dim;
  using SubdomainData =
      ElementCenteredSubdomainData<Dim,
                                   db::get_variables_tags_list<residual_tag>>;
  using OverlapData = typename SubdomainData::OverlapData;
  using overlap_data_inbox_tag =
      OverlapDataInboxTag<Dim, OptionsGroup, OverlapData>;

 public:
  using const_global_cache_tags = tmpl::list<Tags::MaxOverlap<OptionsGroup>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::ConstGlobalCache<Metavariables>& cache,
      const ElementId<Dim>& /*element_id*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    // Skip communicating if the overlap is empty
    if (UNLIKELY(db::get<Tags::MaxOverlap<OptionsGroup>>(box) == 0)) {
      return {std::move(box)};
    }

    const auto& residual = get<residual_tag>(box);
    const auto& element = get<domain::Tags::Element<Dim>>(box);
    const auto& full_extents = get<domain::Tags::Mesh<Dim>>(box).extents();
    const auto& all_intruding_extents =
        get<domain::Tags::Interface<domain::Tags::InternalDirections<Dim>,
                                    Tags::IntrudingExtents<Dim, OptionsGroup>>>(
            box);
    const auto& temporal_id =
        get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);

    auto& receiver_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);

    // Send data on intruding overlaps to the corresponding neighbors
    for (const auto& direction_and_neighbors : element.neighbors()) {
      const auto& direction = direction_and_neighbors.first;
      const auto& orientation = direction_and_neighbors.second.orientation();
      const auto direction_from_neighbor = orientation(direction.opposite());
      const auto& intruding_extents = all_intruding_extents.at(direction);
      for (const auto& neighbor : direction_and_neighbors.second) {
        auto overlap_data = LinearSolver::Schwarz::data_on_overlap(
            residual, full_extents, intruding_extents, direction);
        Parallel::receive_data<overlap_data_inbox_tag>(
            receiver_proxy[neighbor], temporal_id,
            std::make_pair(
                std::make_pair(direction_from_neighbor, element.id()),
                std::move(overlap_data)));
      }
    }
    return {std::move(box)};
  }
};

template <typename FieldsTag, typename OptionsGroup, typename SubdomainOperator>
struct ReceiveOverlapData {
 private:
  using fields_tag = FieldsTag;
  using residual_tag =
      db::add_tag_prefix<LinearSolver::Tags::Residual, fields_tag>;
  static constexpr size_t Dim = SubdomainOperator::volume_dim;
  using SubdomainData =
      ElementCenteredSubdomainData<Dim,
                                   db::get_variables_tags_list<residual_tag>>;
  using OverlapData = typename SubdomainData::OverlapData;
  using overlap_data_inbox_tag =
      OverlapDataInboxTag<Dim, OptionsGroup, OverlapData>;

 public:
  using const_global_cache_tags = tmpl::list<Tags::MaxOverlap<OptionsGroup>>;
  using inbox_tags = tmpl::list<overlap_data_inbox_tag>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& /*element_id*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    if (UNLIKELY(db::get<Tags::MaxOverlap<OptionsGroup>>(box) == 0)) {
      return {std::move(box)};
    }
    auto& inbox = tuples::get<overlap_data_inbox_tag>(inboxes);
    const auto& temporal_id =
        get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
    const auto temporal_received = inbox.find(temporal_id);
    if (temporal_received != inbox.end()) {
      db::mutate<SubdomainDataBufferTag<SubdomainData, OptionsGroup>>(
          make_not_null(&box),
          [&temporal_received](
              const gsl::not_null<SubdomainData*> subdomain_data) noexcept {
            for (auto& overlap_id_and_data : temporal_received->second) {
              const auto& overlap_id = overlap_id_and_data.first;
              auto& received_overlap_data = overlap_id_and_data.second;
              subdomain_data->overlap_data.insert_or_assign(
                  overlap_id, std::move(received_overlap_data));
            }
          });
      inbox.erase(temporal_received);
    }
    return {std::move(box)};
  }

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            size_t Dim>
  static bool is_ready(
      const db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& /*element_id*/) noexcept {
    if (UNLIKELY(db::get<Tags::MaxOverlap<OptionsGroup>>(box) == 0)) {
      return true;
    }
    return dg::has_received_from_all_mortars<overlap_data_inbox_tag>(
        get<LinearSolver::Tags::IterationId<OptionsGroup>>(box),
        get<domain::Tags::Element<Dim>>(box), inboxes);
  }
};

template <size_t Dim, typename OptionsGroup, typename OverlapSolution>
struct OverlapSolutionInboxTag
    : public Parallel::InboxInserters::Map<
          OverlapSolutionInboxTag<Dim, OptionsGroup, OverlapSolution>> {
  using temporal_id = size_t;
  using type =
      std::unordered_map<temporal_id, OverlapMap<Dim, OverlapSolution>>;
};

template <typename FieldsTag, typename OptionsGroup, typename SubdomainOperator>
struct SolveSubdomain {
 private:
  using fields_tag = FieldsTag;
  using residual_tag =
      db::add_tag_prefix<LinearSolver::Tags::Residual, fields_tag>;
  static constexpr size_t Dim = SubdomainOperator::volume_dim;
  using SubdomainData =
      ElementCenteredSubdomainData<Dim,
                                   db::get_variables_tags_list<residual_tag>>;
  using OverlapSolution = typename SubdomainData::OverlapData;
  using overlap_solution_inbox_tag =
      OverlapSolutionInboxTag<Dim, OptionsGroup, OverlapSolution>;

 public:
  using const_global_cache_tags =
      tmpl::list<Tags::MaxOverlap<OptionsGroup>,
                 LinearSolver::Tags::Verbosity<OptionsGroup>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::ConstGlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    // Copy the central element residual into the subdomain data
    db::mutate<SubdomainDataBufferTag<SubdomainData, OptionsGroup>>(
        make_not_null(&box),
        [](const gsl::not_null<SubdomainData*> subdomain_data,
           const db::const_item_type<residual_tag>& residual) noexcept {
          subdomain_data->element_data = residual;
        },
        db::get<residual_tag>(box));

    // Retrieve the subdomain residual that we assembled by communicating
    // overlap data across neighbors
    const auto& subdomain_residual =
        db::get<SubdomainDataBufferTag<SubdomainData, OptionsGroup>>(box);

    // Allocate workspace memory for applying the subdomain operator
    const size_t num_points =
        db::get<domain::Tags::Mesh<Dim>>(box).number_of_grid_points();
    SubdomainOperator subdomain_operator{num_points};
    SubdomainData subdomain_result_buffer{num_points};
    // Solve the subdomain problem
    const auto& subdomain_solver =
        get<Tags::SubdomainSolverBase<OptionsGroup>>(box);
    auto subdomain_solve_result = subdomain_solver(
        [&box, &subdomain_result_buffer,
         &subdomain_operator](const SubdomainData& arg) noexcept {
          // The subdomain operator can retrieve any information on the
          // subdomain geometry that is available through the DataBox. The user
          // is responsible for communicating this information across neighbors
          // if necessary.
          db::apply<typename SubdomainOperator::element_operator>(
              box, arg, make_not_null(&subdomain_result_buffer),
              make_not_null(&subdomain_operator));
          tmpl::for_each<tmpl::list<domain::Tags::InternalDirections<Dim>,
                                    domain::Tags::BoundaryDirectionsInterior<
                                        Dim>>>([&box, &arg,
                                                &subdomain_result_buffer,
                                                &subdomain_operator](
                                                   auto directions_v) noexcept {
            using directions = tmpl::type_from<decltype(directions_v)>;
            using face_operator =
                typename SubdomainOperator::template face_operator<directions>;
            interface_apply<directions, typename face_operator::argument_tags,
                            get_volume_tags<face_operator>>(
                face_operator{}, box, arg,
                make_not_null(&subdomain_result_buffer),
                make_not_null(&subdomain_operator));
          });
          return subdomain_result_buffer;
        },
        subdomain_residual,
        make_with_value<SubdomainData>(subdomain_residual, 0.));
    const auto& subdomain_solve_has_converged = subdomain_solve_result.first;
    auto& subdomain_solution = subdomain_solve_result.second;
    if (not subdomain_solve_has_converged or
        subdomain_solve_has_converged.reason() ==
            Convergence::Reason::MaxIterations) {
      Parallel::printf(
          "WARNING: Subdomain solver on element %s did not converge in %zu "
          "iterations.\n",
          element_id, subdomain_solve_has_converged.num_iterations());
    } else if (UNLIKELY(
                   static_cast<int>(
                       get<LinearSolver::Tags::Verbosity<OptionsGroup>>(box)) >=
                   static_cast<int>(::Verbosity::Verbose))) {
      Parallel::printf(
          "%s Subdomain solver converged in %zu iterations (%s): %e -> %e \n",
          element_id, subdomain_solve_has_converged.num_iterations(),
          subdomain_solve_has_converged.reason(),
          subdomain_solve_has_converged.initial_residual_magnitude(),
          subdomain_solve_has_converged.residual_magnitude());
    }

    // Apply weighting
    if (LIKELY(db::get<Tags::MaxOverlap<OptionsGroup>>(box) > 0)) {
      subdomain_solution.element_data *=
          get(db::get<Tags::Weight<OptionsGroup>>(box));
    }

    // Apply solution to central element
    db::mutate<fields_tag>(
        make_not_null(&box),
        [&subdomain_solution](
            const gsl::not_null<db::item_type<fields_tag>*> fields) noexcept {
          *fields += subdomain_solution.element_data;
        });

    // Send overlap solutions back to the neighbors that they are on
    if (LIKELY(db::get<Tags::MaxOverlap<OptionsGroup>>(box) > 0)) {
      auto& receiver_proxy =
          Parallel::get_parallel_component<ParallelComponent>(cache);
      const auto& element = db::get<domain::Tags::Element<Dim>>(box);
      const auto& temporal_id =
          db::get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
      for (auto& overlap_id_and_solution : subdomain_solution.overlap_data) {
        const auto& overlap_id = overlap_id_and_solution.first;
        const auto& direction = overlap_id.first;
        const auto& neighbor_id = overlap_id.second;
        const auto& orientation =
            element.neighbors().at(direction).orientation();
        const auto direction_from_neighbor = orientation(direction.opposite());
        auto& overlap_solution = overlap_id_and_solution.second;
        Parallel::receive_data<overlap_solution_inbox_tag>(
            receiver_proxy[neighbor_id], temporal_id,
            std::make_pair(
                std::make_pair(direction_from_neighbor, element.id()),
                std::move(overlap_solution)));
      }
    }
    return {std::move(box)};
  }
};

template <typename FieldsTag, typename OptionsGroup, typename SubdomainOperator>
struct ReceiveOverlapSolution {
 private:
  using fields_tag = FieldsTag;
  using residual_tag =
      db::add_tag_prefix<LinearSolver::Tags::Residual, fields_tag>;
  static constexpr size_t Dim = SubdomainOperator::volume_dim;
  using SubdomainData =
      ElementCenteredSubdomainData<Dim,
                                   db::get_variables_tags_list<residual_tag>>;
  using OverlapSolution = typename SubdomainData::OverlapData;
  using overlap_solution_inbox_tag =
      OverlapSolutionInboxTag<Dim, OptionsGroup, OverlapSolution>;

 public:
  using const_global_cache_tags = tmpl::list<Tags::MaxOverlap<OptionsGroup>>;
  using inbox_tags = tmpl::list<overlap_solution_inbox_tag>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& /*element_id*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    if (UNLIKELY(db::get<Tags::MaxOverlap<OptionsGroup>>(box) == 0)) {
      return {std::move(box)};
    }
    auto& inbox = tuples::get<overlap_solution_inbox_tag>(inboxes);
    const auto& temporal_id =
        get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
    const auto temporal_received = inbox.find(temporal_id);
    if (temporal_received != inbox.end()) {
      db::mutate<fields_tag>(
          make_not_null(&box),
          [&temporal_received](
              const gsl::not_null<db::item_type<fields_tag>*> fields,
              const Index<Dim>& full_extents, const auto& all_intruding_extents,
              const std::unordered_map<Direction<Dim>, Scalar<DataVector>>&
                  all_intruding_overlap_weights) noexcept {
            for (const auto& overlap_id_and_solution :
                 temporal_received->second) {
              const auto& overlap_id = overlap_id_and_solution.first;
              const auto& direction = overlap_id.first;
              const auto& overlap_solution = overlap_id_and_solution.second;
              const auto& intruding_extents =
                  all_intruding_extents.at(direction);
              LinearSolver::Schwarz::add_overlap_data(
                  fields,
                  overlap_solution *
                      get(all_intruding_overlap_weights.at(direction)),
                  full_extents, intruding_extents, direction);
            }
          },
          db::get<domain::Tags::Mesh<Dim>>(box).extents(),
          db::get<domain::Tags::Interface<
              domain::Tags::InternalDirections<Dim>,
              Tags::IntrudingExtents<Dim, OptionsGroup>>>(box),
          db::get<domain::Tags::Interface<domain::Tags::InternalDirections<Dim>,
                                          Tags::Weight<OptionsGroup>>>(box));
      inbox.erase(temporal_received);
    }
    return {std::move(box)};
  }

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            size_t Dim>
  static bool is_ready(
      const db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& /*element_id*/) noexcept {
    if (UNLIKELY(db::get<Tags::MaxOverlap<OptionsGroup>>(box) == 0)) {
      return true;
    }
    return dg::has_received_from_all_mortars<overlap_solution_inbox_tag>(
        get<LinearSolver::Tags::IterationId<OptionsGroup>>(box),
        get<domain::Tags::Element<Dim>>(box), inboxes);
  }
};

}  // namespace LinearSolver::Schwarz::detail
