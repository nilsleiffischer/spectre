// Distributed under the MIT License.
// See LICENSE.txt for details.

#define CATCH_CONFIG_RUNNER

#include <vector>

#include "ErrorHandling/FloatingPointExceptions.hpp"
#include "IO/Observer/Helpers.hpp"            // IWYU pragma: keep
#include "IO/Observer/ObserverComponent.hpp"  // IWYU pragma: keep
#include "ParallelAlgorithms/LinearSolver/ConjugateGradient/ConjugateGradient.hpp"
// #include "NumericalAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "Helpers/ParallelAlgorithms/NonlinearSolver/NonlinearSolverAlgorithmTestHelpers.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Main.hpp"
#include "ParallelAlgorithms/NonlinearSolver/Globalization/LineSearch/LineSearch.hpp"
#include "ParallelAlgorithms/NonlinearSolver/NewtonRaphson/NewtonRaphson.hpp"
#include "Utilities/TMPL.hpp"

namespace helpers = NonlinearSolverAlgorithmTestHelpers;

namespace {

struct LinearSolverGroup {
  static std::string name() noexcept { return "LinearSolver"; }
  static constexpr OptionString help = "Options for the linear solver";
};

struct Metavariables {
  using nonlinear_solver = NonlinearSolver::NewtonRaphson<
      Metavariables, helpers::fields_tag,
      NonlinearSolver::Globalization::LineSearch>;
  using linear_solver = LinearSolver::cg::ConjugateGradient<
      Metavariables, helpers::correction_tag, LinearSolverGroup>;

  using component_list =
      tmpl::append<tmpl::list<helpers::ElementArray<Metavariables>,
                              observers::ObserverWriter<Metavariables>,
                              helpers::OutputCleaner<Metavariables>>,
                   typename nonlinear_solver::component_list,
                   typename linear_solver::component_list>;

  using observed_reduction_data_tags =
      observers::collect_reduction_data_tags<tmpl::list<linear_solver>>;

  static constexpr const char* const help{
      "Test the Newton-Raphson nonlinear solver algorithm"};
  static constexpr bool ignore_unrecognized_command_line_options = false;

  enum class Phase {
    Initialization,
    RegisterWithObserver,
    Solve,
    TestResult,
    CleanOutput,
    Exit
  };

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
        return Phase::TestResult;
      case Phase::TestResult:
        return Phase::CleanOutput;
      default:
        return Phase::Exit;
    }
  }
};

}  // namespace

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};

using charmxx_main_component = Parallel::Main<Metavariables>;

#include "Parallel/CharmMain.tpp"  // IWYU pragma: keep
