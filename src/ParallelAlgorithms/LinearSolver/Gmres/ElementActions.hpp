// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <tuple>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DenseVector.hpp"
#include "Parallel/ConstGlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Invoke.hpp"
#include "Parallel/Reduction.hpp"
#include "ParallelAlgorithms/LinearSolver/Gmres/ResidualMonitorActions.hpp"
#include "ParallelAlgorithms/LinearSolver/InnerProduct.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Requires.hpp"

#include "Parallel/Printf.hpp"

/// \cond
namespace tuples {
template <typename...>
class TaggedTuple;
}  // namespace tuples
namespace LinearSolver {
namespace gmres_detail {
template <typename Metavariables, typename FieldsTa, typename OptionsGroup>
struct ResidualMonitor;
}  // namespace gmres_detail
}  // namespace LinearSolver
/// \endcond

namespace LinearSolver {
namespace gmres_detail {

template <typename FieldsTag, typename OptionsGroup>
struct PrepareSolve {
 private:
  using fields_tag = FieldsTag;
  using initial_fields_tag =
      db::add_tag_prefix<LinearSolver::Tags::Initial, fields_tag>;
  using source_tag = db::add_tag_prefix<::Tags::FixedSource, fields_tag>;
  using operator_applied_to_fields_tag =
      db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, fields_tag>;
  using operand_tag =
      db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
  using preconditioned_operand_tag =
      db::add_tag_prefix<LinearSolver::Tags::Preconditioned, operand_tag>;
  using basis_history_tag =
      LinearSolver::Tags::KrylovSubspaceBasis<operand_tag>;
  using preconditioned_basis_history_tag =
      LinearSolver::Tags::KrylovSubspaceBasis<preconditioned_operand_tag>;

 public:
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&, bool> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::ConstGlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    // Parallel::printf("PrepareSolve\n");
    db::mutate<LinearSolver::Tags::IterationId<OptionsGroup>, operand_tag,
               initial_fields_tag, basis_history_tag,
               preconditioned_basis_history_tag>(
        make_not_null(&box),
        [](const gsl::not_null<size_t*> iteration_id,
           const gsl::not_null<db::item_type<operand_tag>*> operand,
           const gsl::not_null<db::item_type<initial_fields_tag>*>
               initial_fields,
           const gsl::not_null<db::item_type<basis_history_tag>*> basis_history,
           const gsl::not_null<db::item_type<preconditioned_basis_history_tag>*>
               preconditioned_basis_history,
           const db::item_type<source_tag>& source,
           const db::item_type<operator_applied_to_fields_tag>&
               operator_applied_to_fields,
           const db::item_type<fields_tag>& fields) noexcept {
          // We have not started iterating yet, so we initialize the current
          // iteration ID such that the _next_ iteration ID is zero.
          *iteration_id = std::numeric_limits<size_t>::max();
          *operand = source - operator_applied_to_fields;
          *initial_fields = fields;
          *basis_history = db::item_type<basis_history_tag>{};
          *preconditioned_basis_history =
              db::item_type<preconditioned_basis_history_tag>{};
        },
        get<source_tag>(box), get<operator_applied_to_fields_tag>(box),
        get<fields_tag>(box));

    Parallel::contribute_to_reduction<gmres_detail::InitializeResidualMagnitude<
        FieldsTag, OptionsGroup, ParallelComponent>>(
        Parallel::ReductionData<
            Parallel::ReductionDatum<double, funcl::Plus<>, funcl::Sqrt<>>>{
            inner_product(get<operand_tag>(box), get<operand_tag>(box))},
        Parallel::get_parallel_component<ParallelComponent>(cache)[array_index],
        Parallel::get_parallel_component<
            ResidualMonitor<Metavariables, FieldsTag, OptionsGroup>>(cache));

    return {
        std::move(box),
        // Terminate algorithm for now. The `ResidualMonitor` will receive the
        // reduction that is performed above and then broadcast to the following
        // action, which is responsible for restarting the algorithm.
        true};
  }
};

template <typename FieldsTag, typename OptionsGroup>
struct NormalizeInitialOperand {
 private:
  using fields_tag = FieldsTag;
  using operand_tag =
      db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
  using basis_history_tag =
      LinearSolver::Tags::KrylovSubspaceBasis<operand_tag>;

 public:
  template <
      typename ParallelComponent, typename DbTagsList, typename Metavariables,
      typename ArrayIndex, typename DataBox = db::DataBox<DbTagsList>,
      Requires<db::tag_is_retrievable_v<operand_tag, DataBox> and
               db::tag_is_retrievable_v<basis_history_tag, DataBox> and
               db::tag_is_retrievable_v<
                   LinearSolver::Tags::HasConverged<OptionsGroup>, DataBox>> =
          nullptr>
  static void apply(db::DataBox<DbTagsList>& box,
                    Parallel::ConstGlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index,
                    const double residual_magnitude,
                    const Convergence::HasConverged& has_converged) noexcept {
    // Parallel::printf("NormalizeInitialOperand\n");
    db::mutate<operand_tag, basis_history_tag,
               LinearSolver::Tags::HasConverged<OptionsGroup>>(
        make_not_null(&box),
        [residual_magnitude, &has_converged](
            const gsl::not_null<db::item_type<operand_tag>*> operand,
            const gsl::not_null<db::item_type<basis_history_tag>*>
                basis_history,
            const gsl::not_null<Convergence::HasConverged*>
                local_has_converged) noexcept {
          *operand /= residual_magnitude;
          basis_history->push_back(*operand);
          *local_has_converged = has_converged;
        });

    // Proceed with algorithm
    Parallel::get_parallel_component<ParallelComponent>(cache)[array_index]
        .perform_algorithm(true);
  }
};

template <typename FieldsTag, typename OptionsGroup>
struct PrepareStep {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      const Parallel::ConstGlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const ActionList /*meta*/,
      const ParallelComponent* const /*meta*/) noexcept {
    // Parallel::printf("PrepareStep\n");
    using fields_tag = FieldsTag;
    using operand_tag =
        db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
    using preconditioned_operand_tag =
        db::add_tag_prefix<LinearSolver::Tags::Preconditioned, operand_tag>;
    using orthogonalization_iteration_id_tag =
        db::add_tag_prefix<LinearSolver::Tags::Orthogonalization,
                           LinearSolver::Tags::IterationId<OptionsGroup>>;

    db::mutate<LinearSolver::Tags::IterationId<OptionsGroup>,
               orthogonalization_iteration_id_tag, preconditioned_operand_tag>(
        make_not_null(&box),
        [](const gsl::not_null<size_t*> iteration_id,
           const gsl::not_null<
               db::item_type<orthogonalization_iteration_id_tag>*>
               orthogonalization_iteration_id,
           const gsl::not_null<db::item_type<preconditioned_operand_tag>*>
               preconditioned_operand,
           const size_t& next_iteration_id,
           const db::const_item_type<operand_tag>& operand) noexcept {
          *iteration_id = next_iteration_id;
          *orthogonalization_iteration_id = 0;
          // Start the preconditioner with an initial guess equal to its source,
          // so not running the preconditioner at all means it is the identity.
          // Another possibility would be to start the preconditioner at zero
          // because we have no reason to expect the remaining residual to have
          // a particular form.
          *preconditioned_operand =
              db::item_type<preconditioned_operand_tag>(operand);
        },
        get<::Tags::Next<LinearSolver::Tags::IterationId<OptionsGroup>>>(box),
        get<operand_tag>(box));
    return {std::move(box)};
  }
};

template <typename FieldsTag, typename OptionsGroup>
struct PerformStep {
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static std::tuple<db::DataBox<DbTagsList>&&, bool> apply(
      db::DataBox<DbTagsList>& box,
      const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
      Parallel::ConstGlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const ActionList /*meta*/,
      // NOLINTNEXTLINE(readability-avoid-const-params-in-decls)
      const ParallelComponent* const /*meta*/) noexcept {
    // Parallel::printf("PerformStep\n");
    using fields_tag = FieldsTag;
    using operand_tag =
        db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
    using preconditioned_operand_tag =
        db::add_tag_prefix<LinearSolver::Tags::Preconditioned, operand_tag>;
    using operator_applied_to_preconditioned_operand_tag =
        db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo,
                           preconditioned_operand_tag>;
    using basis_history_tag =
        LinearSolver::Tags::KrylovSubspaceBasis<operand_tag>;
    using preconditioned_basis_history_tag =
        LinearSolver::Tags::KrylovSubspaceBasis<preconditioned_operand_tag>;

    ASSERT(get<LinearSolver::Tags::IterationId<OptionsGroup>>(box) !=
               std::numeric_limits<size_t>::max(),
           "'" + option_name<OptionsGroup>() +
               "' iteration ID is at initial state. Did you forget to "
               "invoke 'PrepareStep'?");

    db::mutate<preconditioned_basis_history_tag, operand_tag>(
        make_not_null(&box),
        [](const gsl::not_null<db::item_type<preconditioned_basis_history_tag>*>
               preconditioned_basis_history,
           const gsl::not_null<db::item_type<operand_tag>*> operand,
           const db::const_item_type<preconditioned_operand_tag>&
               preconditioned_operand,
           const db::const_item_type<
               operator_applied_to_preconditioned_operand_tag>&
               operator_applied_to_preconditioned_operand) noexcept {
          preconditioned_basis_history->push_back(preconditioned_operand);
          *operand = db::item_type<operand_tag>(
              operator_applied_to_preconditioned_operand);
        },
        get<preconditioned_operand_tag>(box),
        get<operator_applied_to_preconditioned_operand_tag>(box));

    Parallel::contribute_to_reduction<
        StoreOrthogonalization<FieldsTag, OptionsGroup, ParallelComponent>>(
        Parallel::ReductionData<
            Parallel::ReductionDatum<double, funcl::Plus<>>>{inner_product(
            get<basis_history_tag>(box)[0], get<operand_tag>(box))},
        Parallel::get_parallel_component<ParallelComponent>(cache)[array_index],
        Parallel::get_parallel_component<
            ResidualMonitor<Metavariables, FieldsTag, OptionsGroup>>(cache));

    // Terminate algorithm for now. The `ResidualMonitor` will receive the
    // reduction that is performed above and then broadcast to the following
    // action, which is responsible for restarting the algorithm.
    return {std::move(box), true};
  }
};

template <typename FieldsTag, typename OptionsGroup>
struct OrthogonalizeOperand {
 private:
  using fields_tag = FieldsTag;
  using operand_tag =
      db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
  using orthogonalization_iteration_id_tag =
      db::add_tag_prefix<LinearSolver::Tags::Orthogonalization,
                         LinearSolver::Tags::IterationId<OptionsGroup>>;
  using basis_history_tag =
      LinearSolver::Tags::KrylovSubspaceBasis<operand_tag>;

 public:
  template <
      typename ParallelComponent, typename DbTagsList, typename Metavariables,
      typename ArrayIndex, typename DataBox = db::DataBox<DbTagsList>,
      Requires<db::tag_is_retrievable_v<fields_tag, DataBox> and
               db::tag_is_retrievable_v<operand_tag, DataBox> and
               db::tag_is_retrievable_v<orthogonalization_iteration_id_tag,
                                        DataBox> and
               db::tag_is_retrievable_v<basis_history_tag, DataBox> and
               db::tag_is_retrievable_v<
                   LinearSolver::Tags::IterationId<OptionsGroup>, DataBox>> =
          nullptr>
  static void apply(db::DataBox<DbTagsList>& box,
                    Parallel::ConstGlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index,
                    const double orthogonalization) noexcept {
    // Parallel::printf("OrthogonalizeOperand\n");
    db::mutate<operand_tag, orthogonalization_iteration_id_tag>(
        make_not_null(&box),
        [orthogonalization](
            const gsl::not_null<db::item_type<operand_tag>*> operand,
            const gsl::not_null<
                db::item_type<orthogonalization_iteration_id_tag>*>
                orthogonalization_iteration_id,
            const db::const_item_type<basis_history_tag>&
                basis_history) noexcept {
          *operand -= orthogonalization *
                      gsl::at(basis_history, *orthogonalization_iteration_id);
          (*orthogonalization_iteration_id)++;
        },
        get<basis_history_tag>(box));

    const auto& next_orthogonalization_iteration_id =
        get<orthogonalization_iteration_id_tag>(box);
    const auto& iteration_id =
        get<LinearSolver::Tags::IterationId<OptionsGroup>>(box);

    if (next_orthogonalization_iteration_id <= iteration_id) {
      Parallel::contribute_to_reduction<
          StoreOrthogonalization<FieldsTag, OptionsGroup, ParallelComponent>>(
          Parallel::ReductionData<
              Parallel::ReductionDatum<double, funcl::Plus<>>>{
              inner_product(gsl::at(get<basis_history_tag>(box),
                                    next_orthogonalization_iteration_id),
                            get<operand_tag>(box))},
          Parallel::get_parallel_component<ParallelComponent>(
              cache)[array_index],
          Parallel::get_parallel_component<
              ResidualMonitor<Metavariables, FieldsTag, OptionsGroup>>(cache));
    } else {
      Parallel::contribute_to_reduction<StoreFinalOrthogonalization<
          FieldsTag, OptionsGroup, ParallelComponent>>(
          Parallel::ReductionData<
              Parallel::ReductionDatum<double, funcl::Plus<>>>{
              inner_product(get<operand_tag>(box), get<operand_tag>(box))},
          Parallel::get_parallel_component<ParallelComponent>(
              cache)[array_index],
          Parallel::get_parallel_component<
              ResidualMonitor<Metavariables, FieldsTag, OptionsGroup>>(cache));
    }
  }
};

template <typename FieldsTag, typename OptionsGroup>
struct NormalizeOperandAndUpdateField {
 private:
  using fields_tag = FieldsTag;
  using initial_fields_tag =
      db::add_tag_prefix<LinearSolver::Tags::Initial, fields_tag>;
  using operand_tag =
      db::add_tag_prefix<LinearSolver::Tags::Operand, fields_tag>;
  using preconditioned_operand_tag =
      db::add_tag_prefix<LinearSolver::Tags::Preconditioned, operand_tag>;
  using basis_history_tag =
      LinearSolver::Tags::KrylovSubspaceBasis<operand_tag>;
  using preconditioned_basis_history_tag =
      LinearSolver::Tags::KrylovSubspaceBasis<preconditioned_operand_tag>;

 public:
  template <
      typename ParallelComponent, typename DbTagsList, typename Metavariables,
      typename ArrayIndex, typename DataBox = db::DataBox<DbTagsList>,
      Requires<db::tag_is_retrievable_v<fields_tag, DataBox> and
               db::tag_is_retrievable_v<initial_fields_tag, DataBox> and
               db::tag_is_retrievable_v<operand_tag, DataBox> and
               db::tag_is_retrievable_v<basis_history_tag, DataBox> and
               db::tag_is_retrievable_v<
                   LinearSolver::Tags::HasConverged<OptionsGroup>, DataBox>> =
          nullptr>
  static void apply(db::DataBox<DbTagsList>& box,
                    Parallel::ConstGlobalCache<Metavariables>& cache,
                    const ArrayIndex& array_index, const double normalization,
                    const DenseVector<double>& minres,
                    const Convergence::HasConverged& has_converged) noexcept {
    // Parallel::printf("NormalizeOperandAndUpdateField\n");
    db::mutate<operand_tag, basis_history_tag, fields_tag,
               LinearSolver::Tags::HasConverged<OptionsGroup>>(
        make_not_null(&box),
        [normalization, &minres, &has_converged](
            const gsl::not_null<db::item_type<operand_tag>*> operand,
            const gsl::not_null<db::item_type<basis_history_tag>*>
                basis_history,
            const gsl::not_null<db::item_type<fields_tag>*> field,
            const gsl::not_null<Convergence::HasConverged*> local_has_converged,
            const db::const_item_type<initial_fields_tag>& initial_field,
            const db::item_type<preconditioned_basis_history_tag>&
                preconditioned_basis_history) noexcept {
          *operand /= normalization;
          basis_history->push_back(*operand);
          *field = initial_field;
          for (size_t i = 0; i < minres.size(); i++) {
            *field += minres[i] * gsl::at(preconditioned_basis_history, i);
          }
          *local_has_converged = has_converged;
        },
        get<initial_fields_tag>(box),
        get<preconditioned_basis_history_tag>(box));

    // Proceed with algorithm
    Parallel::get_parallel_component<ParallelComponent>(cache)[array_index]
        .perform_algorithm(true);
  }
};

}  // namespace gmres_detail
}  // namespace LinearSolver
