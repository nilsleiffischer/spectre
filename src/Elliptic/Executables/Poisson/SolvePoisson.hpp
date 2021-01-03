// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "Domain/Creators/RegisterDerivedWithCharm.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/Actions/InitializeAnalyticSolution.hpp"
#include "Elliptic/Actions/InitializeSystem.hpp"
#include "Elliptic/DiscontinuousGalerkin/DgElementArray.hpp"
#include "Elliptic/DiscontinuousGalerkin/ImposeBoundaryConditions.hpp"
// #include
// "Elliptic/DiscontinuousGalerkin/ImposeInhomogeneousBoundaryConditionsOnSource.hpp"
#include "Elliptic/BoundaryConditions/AnalyticSolution.hpp"
#include "Elliptic/BoundaryConditions/BoundaryCondition.hpp"
#include "Elliptic/BoundaryConditions/Zero.hpp"
#include "Elliptic/DiscontinuousGalerkin/InitializeFirstOrderOperator.hpp"
#include "Elliptic/DiscontinuousGalerkin/NumericalFluxes/FirstOrderInternalPenalty.hpp"
#include "Elliptic/FirstOrderOperator.hpp"
#include "Elliptic/Systems/Poisson/FirstOrderSystem.hpp"
#include "Elliptic/Tags.hpp"
#include "Elliptic/Triggers/EveryNIterations.hpp"
#include "ErrorHandling/FloatingPointExceptions.hpp"
#include "IO/Observer/Actions/RegisterEvents.hpp"
#include "IO/Observer/Helpers.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "NumericalAlgorithms/Convergence/Tags.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/BoundarySchemes/FirstOrder/FirstOrderScheme.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/Tags.hpp"
#include "Options/Options.hpp"
#include "Parallel/Actions/SetupDataBox.hpp"
#include "Parallel/Actions/TerminatePhase.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/RegisterDerivedClassesWithCharm.hpp"
#include "ParallelAlgorithms/Actions/MutateApply.hpp"
#include "ParallelAlgorithms/DiscontinuousGalerkin/CollectDataForFluxes.hpp"
#include "ParallelAlgorithms/DiscontinuousGalerkin/FluxCommunication.hpp"
#include "ParallelAlgorithms/DiscontinuousGalerkin/InitializeDomain.hpp"
#include "ParallelAlgorithms/DiscontinuousGalerkin/InitializeInterfaces.hpp"
#include "ParallelAlgorithms/DiscontinuousGalerkin/InitializeMortars.hpp"
#include "ParallelAlgorithms/Events/ObserveErrorNorms.hpp"
#include "ParallelAlgorithms/Events/ObserveFields.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Actions/RunEventsAndTriggers.hpp"
#include "ParallelAlgorithms/Initialization/Actions/AddComputeTags.hpp"
#include "ParallelAlgorithms/Initialization/Actions/RemoveOptionsAndTerminatePhase.hpp"
#include "ParallelAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"
#include "PointwiseFunctions/AnalyticData/AnalyticData.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Poisson/AnalyticSolution.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Poisson/Lorentzian.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Poisson/Moustache.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Poisson/ProductOfSinusoids.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Poisson/Zero.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/TMPL.hpp"

namespace SolvePoisson::OptionTags {
struct LinearSolverGroup {
  static std::string name() noexcept { return "LinearSolver"; }
  static constexpr Options::String help =
      "The iterative Krylov-subspace linear solver";
};
struct GmresGroup {
  static std::string name() noexcept { return "GMRES"; }
  static constexpr Options::String help = "Options for the GMRES linear solver";
  using group = LinearSolverGroup;
};
}  // namespace SolvePoisson::OptionTags

/// \cond
template <size_t Dim>
struct Metavariables {
  static constexpr size_t volume_dim = Dim;

  using system = Poisson::FirstOrderSystem<Dim, Poisson::Geometry::Euclidean>;
  // We'll probably always have an analytic solution for Poisson problems, so
  // we don't bother supporting non-solution backgrounds
  using analytic_solution_registrars = tmpl::flatten<tmpl::list<
      Poisson::Solutions::Registrars::ProductOfSinusoids<Dim>,
      tmpl::conditional_t<Dim == 1 or Dim == 2,
                          Poisson::Solutions::Registrars::Moustache<Dim>,
                          tmpl::list<>>,
      tmpl::conditional_t<Dim == 3,
                          Poisson::Solutions::Registrars::Lorentzian<Dim>,
                          tmpl::list<>>>>;
  // We use the analytic solution to impose boundary conditions. Also allow Zero
  // boundary conditions to (and perhaps other kinds in the future).
  using boundary_condition_registrars =
      tmpl::list<elliptic::BoundaryConditions::Registrars::AnalyticSolution<
                     volume_dim, typename system::primal_fields>,
                 elliptic::BoundaryConditions::Registrars::Zero<
                     volume_dim, typename system::primal_fields>>;
  // We currently only support the trivial "zero" initial guess. This will be
  // generalized ASAP.
  using initial_guess_registrars =
      tmpl::list<Poisson::Solutions::Registrars::Zero<Dim>>;

  static constexpr Options::String help{
      "Find the solution to a Poisson problem."};

  using fluxes_computer_tag =
      elliptic::Tags::FluxesComputer<typename system::fluxes>;
  using analytic_solution_tag = elliptic::Tags::Background<
      Poisson::Solutions::AnalyticSolution<Dim, analytic_solution_registrars>>;
  using boundary_conditions_tag = elliptic::Tags::BoundaryConditions<
      elliptic::BoundaryConditions::BoundaryCondition<
          Dim, boundary_condition_registrars>>;
  using linearized_boundary_conditions_tag =
      elliptic::Tags::LinearizedBoundaryConditions<
          typename boundary_conditions_tag::type::element_type>;
  using initial_guess_tag = elliptic::Tags::InitialGuess<
      ::AnalyticData<Dim, initial_guess_registrars>>;

  // The linear solver algorithm. We must use GMRES since the operator is
  // not positive-definite for the first-order system.
  using linear_solver =
      LinearSolver::gmres::Gmres<Metavariables, typename system::fields_tag,
                                 SolvePoisson::OptionTags::LinearSolverGroup,
                                 false>;
  using linear_solver_iteration_id =
      Convergence::Tags::IterationId<typename linear_solver::options_group>;
  // For the GMRES linear solver we need to apply the DG operator to its
  // internal "operand" in every iteration of the algorithm.
  using linear_operand_tag = db::add_tag_prefix<LinearSolver::Tags::Operand,
                                                typename system::fields_tag>;
  using primal_variables = db::wrap_tags_in<LinearSolver::Tags::Operand,
                                            typename system::primal_fields>;
  using auxiliary_variables =
      db::wrap_tags_in<LinearSolver::Tags::Operand,
                       typename system::auxiliary_fields>;

  // Parse numerical flux parameters from the input file to store in the cache.
  using normal_dot_numerical_flux = Tags::NumericalFlux<
      elliptic::dg::NumericalFluxes::FirstOrderInternalPenalty<
          volume_dim, fluxes_computer_tag, primal_variables,
          auxiliary_variables>>;
  // Specify the DG boundary scheme. We use the strong first-order scheme here
  // that only requires us to compute normals dotted into the first-order
  // fluxes.
  using boundary_scheme = dg::FirstOrderScheme::FirstOrderScheme<
      volume_dim, linear_operand_tag,
      db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo,
                         linear_operand_tag>,
      normal_dot_numerical_flux, linear_solver_iteration_id>;

  // Collect events and triggers
  // (public for use by the Charm++ registration code)
  using observe_fields = typename system::fields_tag::tags_list;
  using analytic_solution_fields = observe_fields;
  using events =
      tmpl::list<dg::Events::Registrars::ObserveFields<
                     volume_dim, linear_solver_iteration_id, observe_fields,
                     analytic_solution_fields>,
                 dg::Events::Registrars::ObserveErrorNorms<
                     linear_solver_iteration_id, analytic_solution_fields>>;
  using triggers = tmpl::list<elliptic::Triggers::Registrars::EveryNIterations<
      linear_solver_iteration_id>>;

  // Collect all items to store in the cache.
  using const_global_cache_tags =
      tmpl::list<analytic_solution_tag, boundary_conditions_tag,
                 linearized_boundary_conditions_tag, initial_guess_tag,
                 fluxes_computer_tag, normal_dot_numerical_flux,
                 Tags::EventsAndTriggers<events, triggers>>;

  // Collect all reduction tags for observers
  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::flatten<tmpl::list<
          typename Event<events>::creatable_classes, linear_solver>>>;

  // Specify all global synchronization points.
  enum class Phase { Initialization, RegisterWithObserver, Solve, Exit };

  using initialization_actions = tmpl::list<
      Actions::SetupDataBox, dg::Actions::InitializeDomain<volume_dim>,
      dg::Actions::InitializeInterfaces<
          system, dg::Initialization::slice_tags_to_face<>,
          dg::Initialization::slice_tags_to_exterior<>,
          dg::Initialization::face_compute_tags<>,
          dg::Initialization::exterior_compute_tags<>, false, false>,
      typename linear_solver::initialize_element,
      elliptic::Actions::InitializeSystem<system, analytic_solution_tag>,
      elliptic::Actions::InitializeAnalyticSolution<analytic_solution_tag,
                                                    analytic_solution_fields>,
      //   elliptic::dg::Actions::ImposeInhomogeneousBoundaryConditionsOnSource<
      //       Metavariables>,
      dg::Actions::InitializeMortars<boundary_scheme>,
      elliptic::dg::Actions::InitializeFirstOrderOperator<
          volume_dim, typename system::fluxes, typename system::sources,
          linear_operand_tag, primal_variables, auxiliary_variables>,
      Initialization::Actions::RemoveOptionsAndTerminatePhase>;

  using build_linear_operator_actions = tmpl::list<
      dg::Actions::CollectDataForFluxes<
          boundary_scheme, domain::Tags::InternalDirections<volume_dim>>,
      dg::Actions::SendDataForFluxes<boundary_scheme>,
      Actions::MutateApply<elliptic::FirstOrderOperator<
          volume_dim, LinearSolver::Tags::OperatorAppliedTo,
          linear_operand_tag>>,
      elliptic::dg::Actions::ImposeBoundaryConditions<
          linearized_boundary_conditions_tag, linear_operand_tag,
          primal_variables, auxiliary_variables, fluxes_computer_tag>,
      dg::Actions::CollectDataForFluxes<
          boundary_scheme,
          domain::Tags::BoundaryDirectionsInterior<volume_dim>>,
      dg::Actions::ReceiveDataForFluxes<boundary_scheme>,
      Actions::MutateApply<boundary_scheme>>;

  using register_actions =
      tmpl::list<observers::Actions::RegisterEventsWithObservers,
                 Parallel::Actions::TerminatePhase>;

  using solve_actions = tmpl::list<
      typename linear_solver::template solve<tmpl::list<
          Actions::RunEventsAndTriggers, build_linear_operator_actions>>,
      Actions::RunEventsAndTriggers, Parallel::Actions::TerminatePhase>;

  using dg_element_array = elliptic::DgElementArray<
      Metavariables,
      tmpl::list<Parallel::PhaseActions<Phase, Phase::Initialization,
                                        initialization_actions>,
                 Parallel::PhaseActions<Phase, Phase::RegisterWithObserver,
                                        register_actions>,
                 Parallel::PhaseActions<Phase, Phase::Solve, solve_actions>>>;

  // Specify all parallel components that will execute actions at some point.
  using component_list = tmpl::flatten<
      tmpl::list<dg_element_array, typename linear_solver::component_list,
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
        metavariables::analytic_solution_tag::type::element_type>,
    &Parallel::register_derived_classes_with_charm<
        metavariables::initial_guess_tag::type::element_type>,
    &Parallel::register_derived_classes_with_charm<
        metavariables::boundary_conditions_tag::type::element_type>,
    &Parallel::register_derived_classes_with_charm<
        metavariables::linearized_boundary_conditions_tag::type::element_type>,
    &Parallel::register_derived_classes_with_charm<
        Event<metavariables::events>>,
    &Parallel::register_derived_classes_with_charm<
        Trigger<metavariables::triggers>>};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};
/// \endcond
