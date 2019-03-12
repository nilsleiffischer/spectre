// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/Actions/ComputeOperatorAction.hpp"
#include "Elliptic/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Elliptic/DiscontinuousGalerkin/ImposeBoundaryConditions.hpp"
#include "Elliptic/Systems/Elasticity/Actions/Observe.hpp"
#include "Elliptic/Systems/Elasticity/FirstOrderSystem.hpp"
#include "ErrorHandling/FloatingPointExceptions.hpp"
#include "IO/Observer/Actions.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Actions/ApplyFluxes.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Actions/ComputeNonconservativeBoundaryFluxes.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Actions/FluxCommunication.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags.hpp"
#include "NumericalAlgorithms/LinearSolver/Actions/TerminateIfConverged.hpp"
#include "NumericalAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "NumericalAlgorithms/LinearSolver/IterationId.hpp"
#include "NumericalAlgorithms/LinearSolver/Tags.hpp"
#include "Options/Options.hpp"
#include "Parallel/ConstGlobalCache.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Elasticity/BentBeam.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/TMPL.hpp"

template <size_t Dim>
struct Metavariables {
  static constexpr OptionString help{
      "Find the solution to an elasticity problem in Dim spatial dimensions.\n"
      "Analytic solution: BentBeam\n"
      "Linear solver: GMRES\n"
      "Numerical flux: FirstOrderInternalPenaltyFlux"};

  // The system provides all equations specific to the problem.
  using system = Elasticity::FirstOrderSystem<Dim>;

  // The analytic solution and corresponding source to solve the elasticity
  // equations for
  using analytic_solution_tag =
      OptionTags::AnalyticSolution<Elasticity::Solutions::BentBeam>;

  // The linear solver algorithm. We must use GMRES since the operator is
  // not positive-definite for the first-order system.
  using linear_solver = LinearSolver::Gmres<Metavariables>;
  using temporal_id = LinearSolver::Tags::IterationId;

  // Parse numerical flux parameters from the input file to store in the cache.
  using normal_dot_numerical_flux = OptionTags::NumericalFluxParams<
      Elasticity::FirstOrderInternalPenaltyFlux<Dim>>;

  // Set up the domain creator from the input file.
  using domain_creator_tag = OptionTags::DomainCreator<Dim, Frame::Inertial>;

  // Collect all items to store in the cache.
  using const_global_cache_tag_list = tmpl::list<
      analytic_solution_tag,
      Elasticity::Tags::ConstitutiveRelation<
          typename analytic_solution_tag::type::constitutive_relation_type>>;

  using observed_reduction_data_tags = observers::collect_reduction_data_tags<
      tmpl::list<Elasticity::Actions::Observe, linear_solver>>;
  struct element_observation_type {};

  // Specify all parallel components that will execute actions at some point.
  using component_list = tmpl::append<
      tmpl::list<Elliptic::DgElementArray<
          Metavariables,
          tmpl::list<
              Elasticity::Actions::Observe,
              LinearSolver::Actions::TerminateIfConverged,
              dg::Actions::ComputeNonconservativeBoundaryFluxes<
                  Tags::InternalDirections<Dim>>,
              dg::Actions::SendDataForFluxes<Metavariables>,
              Elliptic::Actions::ComputeOperatorAction,
              dg::Actions::ComputeNonconservativeBoundaryFluxes<
                  Tags::BoundaryDirectionsInterior<Dim>>,
              Elliptic::dg::Actions::
                  ImposeHomogeneousDirichletBoundaryConditions<Metavariables>,
              dg::Actions::ReceiveDataForFluxes<Metavariables>,
              dg::Actions::ApplyFluxes, typename linear_solver::perform_step>>>,
      typename linear_solver::component_list,
      tmpl::list<observers::Observer<Metavariables>,
                 observers::ObserverWriter<Metavariables>>>;

  // Specify all global synchronization points.
  enum class Phase { Initialization, RegisterWithObserver, Solve, Exit };

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
