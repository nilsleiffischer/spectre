// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "IO/Observer/Tags.hpp"
#include "Informer/LogActions.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Tags.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/InterMeshOperators.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/MeshHierarchy.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/Tags.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace logging {
template <typename Metavariables>
struct Logger;
}  // namespace logging
/// \endcond

namespace LinearSolver::multigrid::detail {

template <size_t Dim, typename FieldsTag, typename OptionsGroup,
          typename SourceTag>
struct InitializeElement {
 private:
  using fields_tag = FieldsTag;
  using operator_applied_to_fields_tag =
      db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, fields_tag>;
  using source_tag = SourceTag;

 public:
  using initialization_tags =
      tmpl::list<Tags::BaseRefinementLevels<Dim>,
                 Tags::ParentRefinementLevels<Dim>,
                 domain::Tags::InitialExtents<Dim>, Tags::MultigridLevel>;
  using initialization_tags_to_keep = tmpl::list<Tags::MultigridLevel>;
  using const_global_cache_tags =
      tmpl::list<LinearSolver::Tags::Iterations<OptionsGroup>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ElementId<Dim>& element_id, const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    const bool is_coarsest_grid =
        get<domain::Tags::InitialRefinementLevels<Dim>>(box) ==
            get<Tags::ParentRefinementLevels<Dim>>(box) and
        get<domain::Tags::InitialExtents<Dim>>(box) ==
            get<Tags::ParentExtents<Dim>>(box);
    std::optional<ElementId<Dim>> parent_id =
        is_coarsest_grid ? std::nullopt
                         : std::make_optional(
                               LinearSolver::multigrid::parent_id(element_id));
    auto child_ids = LinearSolver::multigrid::child_ids(
        element_id,
        get<Tags::BaseRefinementLevels<Dim>>(box)[element_id.block_id()],
        get<domain::Tags::InitialRefinementLevels<Dim>>(
            box)[element_id.block_id()]);

    // Parallel::printf(
    //     "Multigrid level %zu element %s has parent %s and children %s.\n",
    //     get<Tags::MultigridLevel>(box), element_id,
    //     (parent_id ? get_output(*parent_id) : "none"), child_ids);

    auto observation_key_suffix =
        std::make_optional(get<Tags::MultigridLevel>(box) == 0
                               ? std::string{""}
                               : (std::string{"Level"} +
                                  get_output(get<Tags::MultigridLevel>(box))));

    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    // Assuming all multigrid levels have the same extents
    auto parent_mesh =
        parent_id ? domain::Initialization::create_initial_mesh(
                        db::get<Tags::ParentExtents<Dim>>(box), *parent_id)
                  : Mesh<Dim>{};
    // Parallel::printf("Mesh: %s with parent mesh %s\n", mesh, parent_mesh);

    auto restriction_operator =
        parent_id ? LinearSolver::multigrid::restriction_operator<
                        Metavariables::massive_operator>(
                        mesh, parent_mesh, element_id.segment_ids(),
                        parent_id->segment_ids())
                  : std::array<Matrix, Dim>{};

    auto prolongation_operator =
        parent_id ? LinearSolver::multigrid::prolongation_operator(
                        parent_mesh, mesh, element_id.segment_ids(),
                        parent_id->segment_ids())
                  : std::array<Matrix, Dim>{};

    const size_t num_points = mesh.number_of_grid_points();
    return std::make_tuple(
        ::Initialization::merge_into_databox<
            InitializeElement,
            db::AddSimpleTags<
                Tags::ParentElementId<Dim>, Tags::ChildElementIds<Dim>,
                Tags::ParentMesh<Dim>,
                observers::Tags::ObservationKeySuffix<Tags::MultigridLevel>,
                observers::Tags::ObservationKeySuffix<Tags::IsFinestLevel>,
                Tags::RestrictionOperator<Dim, OptionsGroup>,
                Tags::ProlongationOperator<Dim, OptionsGroup>,
                db::add_tag_prefix<Tags::PreSmoothingInitial, fields_tag>,
                db::add_tag_prefix<Tags::PreSmoothingSource, fields_tag>,
                db::add_tag_prefix<Tags::PreSmoothingResult, fields_tag>,
                db::add_tag_prefix<Tags::PreSmoothingResidual, fields_tag>,
                db::add_tag_prefix<Tags::PostSmoothingInitial, fields_tag>,
                db::add_tag_prefix<Tags::PostSmoothingSource, fields_tag>,
                db::add_tag_prefix<Tags::PostSmoothingResult, fields_tag>,
                db::add_tag_prefix<Tags::PostSmoothingResidual, fields_tag>>,
            db::AddComputeTags<
                domain::Tags::JacobianCompute<
                    domain::Tags::ElementMap<Dim>,
                    domain::Tags::Coordinates<Dim, Frame::Logical>>,
                domain::Tags::DetJacobianCompute<Dim, Frame::Logical,
                                                 Frame::Inertial>>>(
            std::move(box), std::move(parent_id), std::move(child_ids),
            std::move(parent_mesh), std::move(observation_key_suffix),
            get<Tags::MultigridLevel>(box) == 0
                ? observation_key_suffix
                : std::optional<std::string>{std::nullopt},
            std::move(restriction_operator), std::move(prolongation_operator),
            // The `PrepareSolve` action populates these
            // tags with initial values.
            db::item_type<
                db::add_tag_prefix<Tags::PreSmoothingInitial, fields_tag>>{
                num_points, 0.},
            db::item_type<
                db::add_tag_prefix<Tags::PreSmoothingSource, fields_tag>>{
                num_points, 0.},
            db::item_type<
                db::add_tag_prefix<Tags::PreSmoothingResult, fields_tag>>{
                num_points, 0.},
            db::item_type<
                db::add_tag_prefix<Tags::PreSmoothingResidual, fields_tag>>{
                num_points, 0.},
            db::item_type<
                db::add_tag_prefix<Tags::PostSmoothingInitial, fields_tag>>{
                num_points, 0.},
            db::item_type<
                db::add_tag_prefix<Tags::PostSmoothingSource, fields_tag>>{
                num_points, 0.},
            db::item_type<
                db::add_tag_prefix<Tags::PostSmoothingResult, fields_tag>>{
                num_points, 0.},
            db::item_type<
                db::add_tag_prefix<Tags::PostSmoothingResidual, fields_tag>>{
                num_points, 0.}));
  }
};

template <size_t Dim, typename SourceTag>
struct DataFromChildrenInboxTag
    : public Parallel::InboxInserters::Map<
          DataFromChildrenInboxTag<Dim, SourceTag>> {
  using temporal_id = size_t;
  using type =
      std::map<temporal_id, FixedHashMap<two_to_the(Dim), ElementId<Dim>,
                                         db::item_type<SourceTag>,
                                         boost::hash<ElementId<Dim>>>>;
};

template <typename FieldsTag, typename OptionsGroup, typename SourceTag>
struct PrepareSolve {
 private:
  using fields_tag = FieldsTag;

 public:
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            size_t Dim, typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&, bool> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& inboxes,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& /*element_id*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    // Terminate the algorithm on coarser multigrid levels. The algorithm will
    // be re-started on these levels whenever they receive data. This allows
    // the finest level to control when to terminate the algorithm. We make sure
    // to not terminate the algorithm if data has already been received from
    // children, because that means the algorithm should certainly continue.
    const bool terminate =
        (not db::get<LinearSolver::multigrid::Tags::IsFinestLevel>(box)) and
        (not(tuples::get<DataFromChildrenInboxTag<Dim, SourceTag>>(inboxes)
                 .count(0) > 0));
    return {std::move(box), terminate};
  }
};

template <typename FieldsTag, typename OptionsGroup, typename SourceTag>
struct SendResidualToParent {
 private:
  using fields_tag = FieldsTag;
  using source_tag = SourceTag;
  using residual_tag =
      db::add_tag_prefix<LinearSolver::Tags::Residual, FieldsTag>;

 public:
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            size_t Dim, typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    const auto& parent_id = get<Tags::ParentElementId<Dim>>(box);
    // We should always have a `parent_id` at this point because we skip this
    // part of the algorithm on the coarsest grid with the
    // `SkipPostsmoothingAtBottom` action
    ASSERT(parent_id,
           "Trying to send data to parent but no parent is set on element "
               << element_id << ".");

    db::mutate<db::add_tag_prefix<Tags::PreSmoothingResult, fields_tag>>(
        make_not_null(&box),
        [](const auto x, const auto& y) noexcept { *x = y; },
        db::get<fields_tag>(box));
    db::mutate<db::add_tag_prefix<Tags::PreSmoothingResidual, fields_tag>>(
        make_not_null(&box),
        [](const auto x, const auto& y) noexcept { *x = y; },
        db::get<residual_tag>(box));

    // TODO: Move jacobian ratio into restriction operator
    // TODO: Make sure the jacobians are handled correctly on curved meshes
    // TODO: Resolve the jacobian by numerically integrating the mass matrix
    // including the jacobian.
    auto residual = db::get<residual_tag>(box);
    if constexpr (not Metavariables::massive_operator) {
      residual *= get(
          db::get<domain::Tags::DetJacobian<Frame::Logical, Frame::Inertial>>(
              box));
    }

    // Restrict the residual to the coarser (parent) grid and treat as source.
    // We restrict before sending the data so the restriction operation is
    // parellelized. The parent only needs to sum up all child contributions.
    // TODO: Do nothing when parent is the the same element
    auto restricted_residual = apply_matrices(
        db::get<Tags::RestrictionOperator<Dim, OptionsGroup>>(box),
        db::item_type<source_tag>(std::move(residual)),
        // TODO: make sure apply_matrices works for non-square matrices
        db::get<domain::Tags::Mesh<Dim>>(box).extents());

    auto& receiver_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);
    const auto& temporal_id =
        db::get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
    Parallel::receive_data<DataFromChildrenInboxTag<Dim, SourceTag>>(
        receiver_proxy[*parent_id], temporal_id,
        std::make_pair(element_id, std::move(restricted_residual)),
        // We re-start the algorithm on coarser levels when they `receive_data`,
        // since it is terminated in `PrepareSolve`.
        true);
    return {std::move(box)};
  }
};

template <size_t Dim, typename FieldsTag, typename OptionsGroup,
          typename SourceTag>
struct RestrictResidualFromChildren {
 private:
  using fields_tag = FieldsTag;
  using operator_applied_to_fields_tag =
      db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, fields_tag>;
  using source_tag = SourceTag;

 public:
  using inbox_tags = tmpl::list<DataFromChildrenInboxTag<Dim, SourceTag>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& /*element_id*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    const auto& child_ids = get<Tags::ChildElementIds<Dim>>(box);

    if (child_ids.empty()) {
      // On finest grid the source should be set by the problem to solve
      db::mutate<db::add_tag_prefix<Tags::PreSmoothingInitial, fields_tag>>(
          make_not_null(&box),
          [](const auto x, const auto& y) noexcept { *x = y; },
          db::get<fields_tag>(box));
      db::mutate<db::add_tag_prefix<Tags::PreSmoothingSource, fields_tag>>(
          make_not_null(&box),
          [](const auto x, const auto& y) noexcept { *x = y; },
          db::get<source_tag>(box));
      return {std::move(box)};
    }
    auto& inbox =
        tuples::get<DataFromChildrenInboxTag<Dim, SourceTag>>(inboxes);
    const auto& temporal_id =
        db::get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
    const auto temporal_received = inbox.find(temporal_id);
    auto& children_data = temporal_received->second;

    // Assemble restricted data from children
    // TODO: Specialize for when single child is the same element
    db::mutate<source_tag>(
        make_not_null(&box),
        [&children_data](const gsl::not_null<db::item_type<source_tag>*> source,
                         const Mesh<Dim>& mesh,
                         const Scalar<DataVector>& det_jacobian) noexcept {
          *source = db::item_type<source_tag>{mesh.number_of_grid_points(), 0.};
          for (auto& child_id_and_data : children_data) {
            *source += child_id_and_data.second;
          }
          if constexpr (not Metavariables::massive_operator) {
            *source /= get(det_jacobian);
          }
        },
        db::get<domain::Tags::Mesh<Dim>>(box),
        db::get<domain::Tags::DetJacobian<Frame::Logical, Frame::Inertial>>(
            box));

    // On coarser grids we solve for a correction, so we set the initial guess
    // to zero
    db::mutate<fields_tag>(
        make_not_null(&box),
        [](const gsl::not_null<db::item_type<fields_tag>*> fields,
           const db::const_item_type<source_tag>& source) noexcept {
          *fields = make_with_value<db::item_type<fields_tag>>(source, 0.);
        },
        db::get<source_tag>(box));

    db::mutate<db::add_tag_prefix<Tags::PreSmoothingInitial, fields_tag>>(
        make_not_null(&box),
        [](const auto x, const auto& y) noexcept { *x = y; },
        db::get<fields_tag>(box));
    db::mutate<db::add_tag_prefix<Tags::PreSmoothingSource, fields_tag>>(
        make_not_null(&box),
        [](const auto x, const auto& y) noexcept { *x = y; },
        db::get<source_tag>(box));

    inbox.erase(temporal_received);
    return {std::move(box)};
  }

  template <typename DbTags, typename... InboxTags, typename Metavariables>
  static bool is_ready(const db::DataBox<DbTags>& box,
                       const tuples::TaggedTuple<InboxTags...>& inboxes,
                       const Parallel::GlobalCache<Metavariables>& /*cache*/,
                       const ElementId<Dim>& /*element_id*/) noexcept {
    const auto& child_ids = get<Tags::ChildElementIds<Dim>>(box);
    if (child_ids.empty()) {
      return true;
    }
    const auto& inbox =
        tuples::get<DataFromChildrenInboxTag<Dim, SourceTag>>(inboxes);
    const auto& temporal_id =
        db::get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
    const auto temporal_received = inbox.find(temporal_id);
    if (temporal_received == inbox.end()) {
      return false;
    }
    const auto& received_children_data = temporal_received->second;
    for (const auto& child_id : child_ids) {
      if (received_children_data.find(child_id) ==
          received_children_data.end()) {
        return false;
      }
    }
    return true;
  }
};

template <typename FieldsTag>
struct DataFromParentInboxTag : public Parallel::InboxInserters::Value<
                                    DataFromParentInboxTag<FieldsTag>> {
  using temporal_id = size_t;
  using type = std::map<temporal_id, db::item_type<FieldsTag>>;
};

template <typename OptionsGroup>
struct PostSmoothingLogFormatter {
  static std::string apply(const size_t iteration_id,
                           const size_t mg_level) noexcept {
    return "'" + Options::name<OptionsGroup>() + "' iteration " +
           get_output(iteration_id + 1) +
           " post-smoothing done on multigrid-level " + get_output(mg_level) +
           ".";
  }
};

template <typename FieldsTag, typename OptionsGroup, typename SourceTag>
struct SendCorrectionToChildren {
 private:
  using fields_tag = FieldsTag;

 public:
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            size_t Dim, typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    const bool is_coarsest_level =
        not get<Tags::ParentElementId<Dim>>(box).has_value();

    if (not is_coarsest_level and
        UNLIKELY(static_cast<int>(
                     get<LinearSolver::Tags::Verbosity<OptionsGroup>>(cache)) >=
                 static_cast<int>(::Verbosity::Verbose))) {
      Parallel::contribute_to_reduction<
          logging::Actions::Log<PostSmoothingLogFormatter<OptionsGroup>>,
          ParallelComponent, Tags::MultigridLevel>(
          Parallel::ReductionData<
              Parallel::ReductionDatum<size_t, funcl::AssertEqual<>>,
              Parallel::ReductionDatum<size_t, funcl::AssertEqual<>>>{
              get<LinearSolver::Tags::IterationId<OptionsGroup>>(box),
              get<Tags::MultigridLevel>(box)},
          Parallel::get_parallel_component<ParallelComponent>(
              cache)[element_id],
          Parallel::get_parallel_component<logging::Logger<Metavariables>>(
              cache),
          *get<Parallel::Tags::SectionBase<Tags::MultigridLevel>>(box),
          get<Tags::MultigridLevel>(box));
    }

    const auto& child_ids = get<Tags::ChildElementIds<Dim>>(box);
    db::mutate<db::add_tag_prefix<Tags::PostSmoothingResult, fields_tag>>(
        make_not_null(&box),
        [](const auto x, const auto& y) noexcept { *x = y; },
        db::get<fields_tag>(box));
    db::mutate<db::add_tag_prefix<Tags::PostSmoothingResidual, fields_tag>>(
        make_not_null(&box),
        [](const auto x, const auto& y) noexcept { *x = y; },
        db::get<db::add_tag_prefix<LinearSolver::Tags::Residual, fields_tag>>(
            box));
    if (child_ids.empty()) {
      return {std::move(box)};
    }
    auto& receiver_proxy =
        Parallel::get_parallel_component<ParallelComponent>(cache);
    const auto& temporal_id =
        db::get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
    for (const auto& child_id : child_ids) {
      Parallel::receive_data<DataFromParentInboxTag<FieldsTag>>(
          receiver_proxy[child_id], temporal_id, get<fields_tag>(box));
    }
    return {std::move(box)};
  }
};

template <size_t Dim, typename FieldsTag, typename OptionsGroup,
          typename SourceTag>
struct ProlongateCorrectionFromParent {
 private:
  using fields_tag = FieldsTag;

 public:
  using inbox_tags = tmpl::list<DataFromParentInboxTag<FieldsTag>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box, tuples::TaggedTuple<InboxTags...>& inboxes,
      const Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ElementId<Dim>& element_id, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    const auto& parent_id = get<Tags::ParentElementId<Dim>>(box);
    // We should always have a `parent_id` at this point because we skip this
    // part of the algorithm on the coarsest grid with the
    // `SkipPostsmoothingAtBottom` action
    ASSERT(parent_id,
           "Trying to receive data from parent but no parent is set on element "
               << element_id << ".");

    auto& inbox = tuples::get<DataFromParentInboxTag<FieldsTag>>(inboxes);
    const auto& temporal_id =
        db::get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
    const auto temporal_received = inbox.find(temporal_id);
    const auto& parent_correction = temporal_received->second;

    // Apply prolongation operator
    // TODO: Do nothing when parent is the same element
    db::mutate<fields_tag>(
        make_not_null(&box),
        [&](const gsl::not_null<db::item_type<fields_tag>*> fields,
            const std::array<Matrix, Dim>& prolongation_operator,
            const Mesh<Dim>& parent_mesh) noexcept {
          *fields += apply_matrices(prolongation_operator, parent_correction,
                                    parent_mesh.extents());
        },
        db::get<Tags::ProlongationOperator<Dim, OptionsGroup>>(box),
        db::get<Tags::ParentMesh<Dim>>(box));

    db::mutate<db::add_tag_prefix<Tags::PostSmoothingInitial, fields_tag>>(
        make_not_null(&box),
        [](const auto x, const auto& y) noexcept { *x = y; },
        db::get<fields_tag>(box));
    db::mutate<db::add_tag_prefix<Tags::PostSmoothingSource, fields_tag>>(
        make_not_null(&box),
        [](const auto x, const auto& y) noexcept { *x = y; },
        db::get<SourceTag>(box));

    inbox.erase(temporal_received);
    return {std::move(box)};
  }

  template <typename DbTags, typename... InboxTags, typename Metavariables>
  static bool is_ready(const db::DataBox<DbTags>& box,
                       const tuples::TaggedTuple<InboxTags...>& inboxes,
                       const Parallel::GlobalCache<Metavariables>& /*cache*/,
                       const ElementId<Dim>& /*element_id*/) noexcept {
    const auto& parent_id = get<Tags::ParentElementId<Dim>>(box);
    if (not parent_id) {
      return true;
    }
    const auto& inbox = tuples::get<DataFromParentInboxTag<FieldsTag>>(inboxes);
    const auto& temporal_id =
        db::get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);
    const auto temporal_received = inbox.find(temporal_id);
    if (temporal_received == inbox.end()) {
      return false;
    }
    return true;
  }
};

template <typename OptionsGroup>
struct PreSmoothingLogFormatter {
  static std::string apply(const size_t iteration_id, const size_t mg_level,
                           const bool is_coarsest_level) noexcept {
    return "'" + Options::name<OptionsGroup>() + "' iteration " +
           get_output(iteration_id + 1) +
           (is_coarsest_level ? " bottom-smoothing " : " pre-smoothing ") +
           "done on multigrid-level " + get_output(mg_level) + ".";
  }
};

template <typename FieldsTag, typename OptionsGroup, typename SourceTag>
struct SkipPostsmoothingAtBottom {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            size_t Dim, typename ActionList, typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&, bool, size_t> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& cache,
      const ElementId<Dim>& element_id, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    const bool is_coarsest_level =
        not get<Tags::ParentElementId<Dim>>(box).has_value();

    if (UNLIKELY(static_cast<int>(
                     get<LinearSolver::Tags::Verbosity<OptionsGroup>>(cache)) >=
                 static_cast<int>(::Verbosity::Verbose))) {
      Parallel::contribute_to_reduction<
          logging::Actions::Log<PreSmoothingLogFormatter<OptionsGroup>>,
          ParallelComponent, Tags::MultigridLevel>(
          Parallel::ReductionData<
              Parallel::ReductionDatum<size_t, funcl::AssertEqual<>>,
              Parallel::ReductionDatum<size_t, funcl::AssertEqual<>>,
              Parallel::ReductionDatum<bool, funcl::AssertEqual<>>>{
              get<LinearSolver::Tags::IterationId<OptionsGroup>>(box),
              get<Tags::MultigridLevel>(box), is_coarsest_level},
          Parallel::get_parallel_component<ParallelComponent>(
              cache)[element_id],
          Parallel::get_parallel_component<logging::Logger<Metavariables>>(
              cache),
          *get<Parallel::Tags::SectionBase<Tags::MultigridLevel>>(box),
          get<Tags::MultigridLevel>(box));
    }

    // On the coarsest grid skip the second smoothing step
    const size_t index_of_next_action =
        is_coarsest_level
            ? tmpl::index_of<ActionList,
                             SendCorrectionToChildren<FieldsTag, OptionsGroup,
                                                      SourceTag>>::value
            : tmpl::index_of<ActionList, SkipPostsmoothingAtBottom>::value + 1;
    return {std::move(box), false, index_of_next_action};
  }
};

}  // namespace LinearSolver::multigrid::detail