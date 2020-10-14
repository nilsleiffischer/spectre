// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/EagerMath/Norms.hpp"
#include "DataStructures/Variables.hpp"
#include "ErrorHandling/Assert.hpp"
#include "Evolution/Initialization/DiscontinuousGalerkin.hpp"
#include "Evolution/Initialization/Evolution.hpp"
#include "Evolution/Systems/ScalarWave/Constraints.hpp"
#include "Evolution/Systems/ScalarWave/System.hpp"
#include "Evolution/Systems/ScalarWave/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "ParallelAlgorithms/Initialization/MergeIntoDataBox.hpp"

namespace ScalarWave {
namespace Actions {
/// \ingroup InitializationGroup
/// \brief Initialize items related to constraints of the ScalarWave system
///
/// We add both constraints and the constraint damping parameter to the
/// evolution databox.
///
/// DataBox changes:
/// - Adds:
///   * `ScalarWave::Tags::ConstraintGamma2`
///   * `ScalarWave::Tags::OneIndexConstraint<Dim>`
///   * `ScalarWave::Tags::TwoIndexConstraint<Dim>`
///   * `::Tags::PointwiseL2Norm<ScalarWave::Tags::OneIndexConstraint<Dim>>`
///   * `::Tags::PointwiseL2Norm<ScalarWave::Tags::TwoIndexConstraint<Dim>>`
/// - Removes: nothing
/// - Modifies: nothing
///
/// \note This action relies on the `SetupDataBox` aggregated initialization
/// mechanism, so `Actions::SetupDataBox` must be present in the
/// `Initialization` phase action list prior to this action.
template <size_t Dim>
struct InitializeConstraints {
  using simple_tags = tmpl::list<ScalarWave::Tags::ConstraintGamma2>;

  using compute_tags = tmpl::list<
      ScalarWave::Tags::OneIndexConstraintCompute<Dim>,
      ScalarWave::Tags::TwoIndexConstraintCompute<Dim>,
      ::Tags::PointwiseL2NormCompute<ScalarWave::Tags::OneIndexConstraint<Dim>>,
      ::Tags::PointwiseL2NormCompute<
          ScalarWave::Tags::TwoIndexConstraint<Dim>>>;

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ActionList,
            typename ParallelComponent>
  static auto apply(db::DataBox<DbTagsList>& box,
                    const tuples::TaggedTuple<InboxTags...>& /*inboxes*/,
                    const Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const ActionList /*meta*/,
                    const ParallelComponent* const /*meta*/) noexcept {
    const auto& mesh = db::get<domain::Tags::Mesh<Dim>>(box);
    Scalar<DataVector> gamma_2{mesh.number_of_grid_points(), 0.};

    db::mutate_assign(make_not_null(&box), simple_tags{}, std::move(gamma_2));

    return std::make_tuple(std::move(box));
  }
};
}  // namespace Actions
}  // namespace ScalarWave
