// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/Actions/ComputeOperatorAction.hpp"
#include "Elliptic/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Elliptic/DiscontinuousGalerkin/InitializeElement.hpp"
#include "Elliptic/Systems/Poisson/Actions/Observe.hpp"
#include "Elliptic/Systems/Poisson/FirstOrderSystem.hpp"
#include "ErrorHandling/FloatingPointExceptions.hpp"
#include "IO/Observer/Actions.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Actions/ComputeNormalDotFluxes.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Actions/FluxCommunication2.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/BoundarySchemes/StrongFirstOrder.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/PopulateBoundaryMortars.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags.hpp"
#include "NumericalAlgorithms/LinearSolver/Actions/TerminateIfConverged.hpp"
#include "NumericalAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "NumericalAlgorithms/LinearSolver/Tags.hpp"
#include "Options/Options.hpp"
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/ConstGlobalCache.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Poisson/ProductOfSinusoids.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
template <size_t Dim>
struct Metavariables {
  static constexpr OptionString help{
      "Find the solution to a Poisson problem in Dim spatial dimensions.\n"
      "Analytic solution: ProductOfSinusoids\n"
      "Linear solver: GMRES\n"
      "Numerical flux: FirstOrderInternalPenaltyFlux"};

  // The system provides all equations specific to the problem.
  using system = Poisson::FirstOrderSystem<Dim>;

  // The analytic solution and corresponding source to solve the Poisson
  // equation for
  using analytic_solution_tag =
      OptionTags::AnalyticSolution<Poisson::Solutions::ProductOfSinusoids<Dim>>;

  // The linear solver algorithm. We must use GMRES since the operator is
  // not positive-definite for the first-order system.
  using linear_solver = LinearSolver::Gmres<Metavariables>;
  using temporal_id = LinearSolver::Tags::IterationId;

  // Parse numerical flux parameters from the input file to store in the cache.
  using normal_dot_numerical_flux =
      OptionTags::NumericalFlux<Poisson::FirstOrderInternalPenaltyFlux<Dim>>;
  using boundary_scheme = dg::BoundarySchemes::StrongFirstOrder<
      Dim, typename system::variables_tag, typename system::normal_dot_fluxes,
      ::Tags::NormalDotNumericalFluxComputer<
          typename normal_dot_numerical_flux::type>,
      LinearSolver::Tags::IterationId>;

  // Set up the domain creator from the input file.
  using domain_creator_tag = OptionTags::DomainCreator<Dim, Frame::Inertial>;

  // Collect all items to store in the cache.
  using const_global_cache_tag_list =
      tmpl::list<normal_dot_numerical_flux, analytic_solution_tag>;

  // Collect all reduction tags for observers
  using observed_reduction_data_tags = observers::collect_reduction_data_tags<
      tmpl::list<Poisson::Actions::Observe, linear_solver>>;

  // Specify all global synchronization points.
  enum class Phase { Initialization, RegisterWithObserver, Solve, Exit };

  // Specify all parallel components that will execute actions at some point.
  using component_list = tmpl::append<
      tmpl::list<Elliptic::DgElementArray<
          Metavariables,
          tmpl::list<
              Parallel::PhaseActions<
                  Phase, Phase::Initialization,
                  tmpl::list<Elliptic::dg::Actions::InitializeElement<Dim>>>,

              Parallel::PhaseActions<
                  Phase, Phase::RegisterWithObserver,
                  tmpl::list<observers::Actions::RegisterWithObservers<
                                 Poisson::Actions::Observe>,
                             Parallel::Actions::TerminatePhase>>,

              Parallel::PhaseActions<
                  Phase, Phase::Solve,
                  tmpl::list<
                      Poisson::Actions::Observe,
                      LinearSolver::Actions::TerminateIfConverged,
                      dg::Actions::ComputeNormalDotFluxes<
                          Tags::InternalDirections<Dim>,
                          typename system::variables_tag,
                          typename system::normal_dot_fluxes>,
                      dg::Actions::SendDataForFluxes<boundary_scheme>,
                      Elliptic::Actions::ComputeOperatorAction,
                      dg::Actions::ComputeNormalDotFluxes<
                          Tags::BoundaryDirectionsInterior<Dim>,
                          typename system::variables_tag,
                          typename system::normal_dot_fluxes>,
                      ::Actions::MutateApply<
                          ::dg::PopulateBoundaryMortars<boundary_scheme>>,
                      dg::Actions::ReceiveDataForFluxes<boundary_scheme>,
                      ::Actions::MutateApply<boundary_scheme>,
                      typename linear_solver::perform_step>>>,
          typename Elliptic::dg::Actions::InitializeElement<
              Dim>::AddOptionsToDataBox>>,
      typename linear_solver::component_list,
      tmpl::list<observers::Observer<Metavariables>,
                 observers::ObserverWriter<Metavariables>>>;

  // Specify the transitions between phases.
  static Phase determine_next_phase(
      const Phase& current_phase,
      const Parallel::CProxy_ConstGlobalCache<
          Metavariables>& /*cache_proxy*/) noexcept {
    switch (current_phase) {
      case Phase::Initialization:
        return Phase::RegisterWithObserver;
      case Phase::RegisterWithObserver:
        return Phase::Solve;
      case Phase::Solve:
        return Phase::Exit;
      case Phase::Exit:
        ERROR(
            "Should never call determine_next_phase with the current phase "
            "being 'Exit'");
      default:
        ERROR(
            "Unknown type of phase. Did you static_cast<Phase> an integral "
            "value?");
    }
  }
};

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling, &domain::creators::register_derived_with_charm};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};
/// \endcond
