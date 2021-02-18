// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/Actions/ApplyLinearOperatorToInitialFields.hpp"
#include "Elliptic/Actions/InitializeAnalyticSolution.hpp"
#include "Elliptic/Actions/InitializeBackgroundFields.hpp"
#include "Elliptic/Actions/InitializeFields.hpp"
#include "Elliptic/Actions/InitializeFixedSources.hpp"
#include "Elliptic/DiscontinuousGalerkin/Actions/ApplyOperator.hpp"
#include "Elliptic/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Elliptic/DiscontinuousGalerkin/SubdomainOperator/InitializeSubdomain.hpp"
#include "Elliptic/DiscontinuousGalerkin/SubdomainOperator/SubdomainOperator.hpp"
#include "Elliptic/Systems/Xcts/FirstOrderSystem.hpp"
#include "Elliptic/Tags.hpp"
#include "Elliptic/Triggers/EveryNIterations.hpp"
#include "IO/Observer/Actions/RegisterEvents.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "NumericalAlgorithms/Convergence/Tags.hpp"
#include "Options/Options.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "ParallelAlgorithms/DiscontinuousGalerkin/InitializeDomain.hpp"
#include "ParallelAlgorithms/Events/ObserveErrorNorms.hpp"
#include "ParallelAlgorithms/Events/ObserveFields.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Actions/RunEventsAndTriggers.hpp"
#include "ParallelAlgorithms/Initialization/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "ParallelAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/Actions/RestrictFields.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/ElementsAllocator.hpp"
#include "ParallelAlgorithms/LinearSolver/Multigrid/Multigrid.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Actions/CommunicateOverlapFields.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/Schwarz.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"
#include "ParallelAlgorithms/NonlinearSolver/NewtonRaphson/NewtonRaphson.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/Flatness.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/Schwarzschild.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ErrorHandling/FloatingPointExceptions.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/TMPL.hpp"

namespace SolveXcts::OptionTags {
struct NonlinearSolverGroup {
  static std::string name() noexcept { return "NonlinearSolver"; }
  static constexpr Options::String help = "The iterative nonlinear solver";
};
struct NewtonRaphsonGroup {
  static std::string name() noexcept { return "NewtonRaphson"; }
  static constexpr Options::String help =
      "Options for the Newton-Raphson nonlinear solver";
  using group = NonlinearSolverGroup;
};
struct LinearSolverGroup {
  static std::string name() noexcept { return "LinearSolver"; }
  static constexpr Options::String help =
      "The iterative Krylov-subspace linear solver";
};
struct GmresGroup {
  static std::string name() noexcept { return "Gmres"; }
  static constexpr Options::String help = "Options for the GMRES linear solver";
  using group = LinearSolverGroup;
};
struct SchwarzSmootherGroup {
  static std::string name() noexcept { return "SchwarzSmoother"; }
  static constexpr Options::String help = "Options for the Schwarz smoother";
  using group = LinearSolverGroup;
};
struct MultigridGroup {
  static std::string name() noexcept { return "Multigrid"; }
  static constexpr Options::String help = "Options for the multigrid";
  using group = LinearSolverGroup;
};
}  // namespace SolveXcts::OptionTags

/// \cond
struct Metavariables {
  static constexpr size_t volume_dim = 3;
  using system =
      Xcts::FirstOrderSystem<Xcts::Equations::HamiltonianLapseAndShift,
                             Xcts::Geometry::FlatCartesian>;

  // List the possible backgrounds, i.e. the variable-independent part of the
  // equations that define the problem to solve (along with the boundary
  // conditions)
  using analytic_solution_registrars =
      tmpl::list<Xcts::Solutions::Registrars::Schwarzschild>;
  using background_tag = elliptic::Tags::Background<
      ::AnalyticData<3, analytic_solution_registrars>>;

  // List the possible initial guesses
  using initial_guess_registrars =
      tmpl::append<tmpl::list<Xcts::Solutions::Registrars::Flatness>,
                   analytic_solution_registrars>;
  using initial_guess_tag = elliptic::Tags::InitialGuess<
      ::AnalyticData<volume_dim, initial_guess_registrars>>;

  static constexpr Options::String help{
      "Find the solution to an XCTS problem."};

  // These are the fields we solve for
  using fields_tag = ::Tags::Variables<typename system::primal_fields>;
  // These are the fluxes corresponding to the fields, i.e. essentially their
  // first derivatives. These are background fields for the linearized sources.
  using fluxes_tag = ::Tags::Variables<typename system::primal_fluxes>;
  // These are the fixed sources, i.e. the RHS of the equations
  using fixed_sources_tag = db::add_tag_prefix<::Tags::FixedSource, fields_tag>;
  using operator_applied_to_fields_tag =
      db::add_tag_prefix<NonlinearSolver::Tags::OperatorAppliedTo, fields_tag>;

  using nonlinear_solver = NonlinearSolver::newton_raphson::NewtonRaphson<
      Metavariables, fields_tag, SolveXcts::OptionTags::NewtonRaphsonGroup,
      fixed_sources_tag, LinearSolver::multigrid::Tags::IsFinestLevel>;
  using nonlinear_solver_iteration_id =
      Convergence::Tags::IterationId<typename nonlinear_solver::options_group>;

  // The linear solver algorithm. We must use GMRES since the operator is
  // not guaranteed to be symmetric. It can be made symmetric by multiplying by
  // the DG mass matrix.
  using linear_solver = LinearSolver::gmres::Gmres<
      Metavariables, typename nonlinear_solver::linear_solver_fields_tag,
      SolveXcts::OptionTags::GmresGroup, true,
      typename nonlinear_solver::linear_solver_source_tag,
      LinearSolver::multigrid::Tags::IsFinestLevel>;
  using linear_solver_iteration_id =
      Convergence::Tags::IterationId<typename linear_solver::options_group>;
  // Precondition each linear solver iteration with a multigrid V-cycle
  using multigrid = LinearSolver::multigrid::Multigrid<
      volume_dim, typename linear_solver::operand_tag,
      SolveXcts::OptionTags::MultigridGroup, elliptic::dg::Tags::Massive,
      typename linear_solver::preconditioner_source_tag>;
  // Smooth each multigrid level with a number of Schwarz smoothing steps
  using subdomain_operator =
      elliptic::dg::subdomain_operator::SubdomainOperator<
          system, SolveXcts::OptionTags::SchwarzSmootherGroup>;
  // This data needs to be communicated on subdomain overlap regions
  using communicated_overlap_tags = tmpl::list<
      // For linearized sources
      fields_tag, fluxes_tag/*,
      // For AH boundary conditions
      domain::Tags::Interface<
          domain::Tags::BoundaryDirectionsInterior<volume_dim>,
          Xcts::Tags::ConformalFactor<DataVector>>,
      domain::Tags::Interface<
          domain::Tags::BoundaryDirectionsInterior<volume_dim>,
          Xcts::Tags::LapseTimesConformalFactor<DataVector>>,
      domain::Tags::Interface<
          domain::Tags::BoundaryDirectionsInterior<volume_dim>,
          ::Tags::NormalDotFlux<Xcts::Tags::ShiftExcess<DataVector, volume_dim,
                                                        Frame::Inertial>>>*/>;
  using schwarz_smoother = LinearSolver::Schwarz::Schwarz<
      typename multigrid::smooth_fields_tag,
      SolveXcts::OptionTags::SchwarzSmootherGroup, subdomain_operator,
      typename multigrid::smooth_source_tag,
      LinearSolver::multigrid::Tags::MultigridLevel>;
  // For the GMRES linear solver we need to apply the DG operator to its
  // internal "operand" in every iteration of the algorithm.
  using correction_vars_tag = typename linear_solver::operand_tag;
  using operator_applied_to_correction_vars_tag =
      db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo,
                         correction_vars_tag>;
  // The correction fluxes can be stored in an arbitrary tag. We don't need to
  // access them anywhere, they're just a memory buffer for the linearized
  // operator.
  using correction_fluxes_tag =
      db::add_tag_prefix<NonlinearSolver::Tags::Correction, fluxes_tag>;

  // Collect events and triggers
  // (public for use by the Charm++ registration code)
  using analytic_solution_fields = typename system::primal_fields;
  using observe_fields = analytic_solution_fields;
  using events =
      tmpl::list<dg::Events::Registrars::ObserveFields<
                     volume_dim, nonlinear_solver_iteration_id, observe_fields,
                     analytic_solution_fields,
                     LinearSolver::multigrid::Tags::MultigridLevel>,
                 dg::Events::Registrars::ObserveErrorNorms<
                     nonlinear_solver_iteration_id, analytic_solution_fields,
                     LinearSolver::multigrid::Tags::IsFinestLevel>>;
  using triggers = tmpl::list<elliptic::Triggers::Registrars::EveryNIterations<
      nonlinear_solver_iteration_id>>;

  // Collect all items to store in the cache.
  using const_global_cache_tags =
      tmpl::list<background_tag, initial_guess_tag,
                 Tags::EventsAndTriggers<events, triggers>>;

  // Collect all reduction tags for observers
  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::flatten<tmpl::list<
          typename Event<events>::creatable_classes, nonlinear_solver,
          linear_solver, multigrid, schwarz_smoother>>>;

  // Specify all global synchronization points.
  enum class Phase { Initialization, RegisterWithObserver, Solve, Exit };

  using initialization_actions = tmpl::list<
      Actions::SetupDataBox, dg::Actions::InitializeDomain<volume_dim>,
      typename nonlinear_solver::initialize_element,
      typename linear_solver::initialize_element,
      typename multigrid::initialize_element,
      typename schwarz_smoother::initialize_element,
      elliptic::Actions::InitializeFields<system, initial_guess_tag>,
      elliptic::Actions::InitializeFixedSources<system, background_tag>,
      elliptic::Actions::InitializeBackgroundFields<system, background_tag>,
      elliptic::Actions::InitializeOptionalAnalyticSolution<
          background_tag,
          tmpl::append<typename system::primal_fields,
                       typename system::primal_fluxes>,
          Xcts::Solutions::AnalyticSolution<analytic_solution_registrars>>,
      elliptic::dg::Actions::initialize_operator<
          system, nonlinear_solver_iteration_id, fields_tag,
          operator_applied_to_fields_tag, fluxes_tag>,
      elliptic::dg::Actions::initialize_operator<
          system, linear_solver_iteration_id, correction_vars_tag,
          operator_applied_to_correction_vars_tag, correction_fluxes_tag>,
      elliptic::dg::Actions::InitializeSubdomain<
          system, background_tag, typename schwarz_smoother::options_group>,
      Initialization::Actions::RemoveOptionsAndTerminatePhase>;

  template <bool Linearized>
  using build_operator_actions = elliptic::dg::Actions::apply_operator<
      system, Linearized,
      tmpl::conditional_t<Linearized, linear_solver_iteration_id,
                          nonlinear_solver_iteration_id>,
      tmpl::conditional_t<Linearized, correction_vars_tag, fields_tag>,
      tmpl::conditional_t<Linearized, operator_applied_to_correction_vars_tag,
                          operator_applied_to_fields_tag>,
      tmpl::conditional_t<Linearized, correction_fluxes_tag, fluxes_tag>>;

  using register_actions =
      tmpl::list<observers::Actions::RegisterEventsWithObservers,
                 typename nonlinear_solver::register_element,
                 typename schwarz_smoother::register_element,
                 typename multigrid::register_element,
                 Parallel::Actions::TerminatePhase>;

  template <typename Label>
  using smooth_actions = tmpl::list<build_operator_actions<true>,
                                    typename schwarz_smoother::template solve<
                                        build_operator_actions<true>, Label>>;

  using solve_actions = tmpl::list<
      // TODO: Only build nonlinear operator on finest grid
      build_operator_actions<false>, Actions::RunEventsAndTriggers,
      typename nonlinear_solver::template solve<
          build_operator_actions<false>,
          tmpl::list<
              LinearSolver::multigrid::Actions::SendFieldsToCoarserGrid<
                  fields_tag, typename multigrid::options_group, void>,
              LinearSolver::multigrid::Actions::SendFieldsToCoarserGrid<
                  fluxes_tag, typename multigrid::options_group, void>,
              LinearSolver::multigrid::Actions::ReceiveFieldsFromFinerGrid<
                  volume_dim, fields_tag, typename multigrid::options_group>,
              LinearSolver::multigrid::Actions::ReceiveFieldsFromFinerGrid<
                  volume_dim, fluxes_tag, typename multigrid::options_group>,
              LinearSolver::Schwarz::Actions::SendOverlapFields<
                  communicated_overlap_tags,
                  typename schwarz_smoother::options_group, false>,
              LinearSolver::Schwarz::Actions::ReceiveOverlapFields<
                  volume_dim, communicated_overlap_tags,
                  typename schwarz_smoother::options_group>,
              typename linear_solver::template solve<
                  typename multigrid::template solve<
                      smooth_actions<LinearSolver::multigrid::VcycleDownLabel>,
                      smooth_actions<LinearSolver::multigrid::VcycleUpLabel>>>>,
          tmpl::list<Actions::RunEventsAndTriggers>>,
      Parallel::Actions::TerminatePhase>;

  using dg_element_array = elliptic::DgElementArray<
      Metavariables,
      tmpl::list<Parallel::PhaseActions<Phase, Phase::Initialization,
                                        initialization_actions>,
                 Parallel::PhaseActions<Phase, Phase::RegisterWithObserver,
                                        register_actions>,
                 Parallel::PhaseActions<Phase, Phase::Solve, solve_actions>>,
      LinearSolver::multigrid::ElementsAllocator<
          volume_dim, typename multigrid::options_group>>;

  // Specify all parallel components that will execute actions at some point.
  using component_list = tmpl::flatten<
      tmpl::list<dg_element_array, typename nonlinear_solver::component_list,
                 typename linear_solver::component_list,
                 typename multigrid::component_list,
                 typename schwarz_smoother::component_list,
                 observers::Observer<Metavariables>,
                 observers::ObserverWriter<Metavariables>>>;

  // Specify the transitions between phases.
  static Phase determine_next_phase(
      const Phase& current_phase,
      const Parallel::CProxy_GlobalCache<
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
    &setup_error_handling,
    &disable_openblas_multithreading,
    &domain::creators::register_derived_with_charm,
    &Parallel::register_derived_classes_with_charm<
        metavariables::background_tag::type::element_type>,
    &Parallel::register_derived_classes_with_charm<
        metavariables::initial_guess_tag::type::element_type>,
    &Parallel::register_derived_classes_with_charm<
        metavariables::system::boundary_conditions_base>,
    &Parallel::register_derived_classes_with_charm<
        metavariables::schwarz_smoother::subdomain_solver>,
    &Parallel::register_derived_classes_with_charm<
        Event<metavariables::events>>,
    &Parallel::register_derived_classes_with_charm<
        Trigger<metavariables::triggers>>};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};
/// \endcond
