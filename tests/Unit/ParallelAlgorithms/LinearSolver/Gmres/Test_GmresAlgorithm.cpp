// Distributed under the MIT License.
// See LICENSE.txt for details.

#define CATCH_CONFIG_RUNNER

#include <vector>

#include "ErrorHandling/FloatingPointExceptions.hpp"
#include "Helpers/ParallelAlgorithms/LinearSolver/LinearSolverAlgorithmTestHelpers.hpp"
#include "Parallel/InitializationFunctions.hpp"
#include "Parallel/Main.hpp"
#include "ParallelAlgorithms/LinearSolver/Gmres/Gmres.hpp"
#include "Utilities/TMPL.hpp"

namespace helpers = LinearSolverAlgorithmTestHelpers;

namespace {

struct SerialGmres {
  static constexpr OptionString help =
      "Options for the iterative linear solver";
};

struct Metavariables {
  static constexpr const char* const help{
      "Test the GMRES linear solver algorithm"};

  using linear_solver = LinearSolver::Gmres<Metavariables, helpers::fields_tag,
                                            SerialGmres, false>;
  using preconditioner = void;

  using component_list = helpers::component_list<Metavariables>;
  using element_observation_type = helpers::element_observation_type;
  using observed_reduction_data_tags =
      helpers::observed_reduction_data_tags<Metavariables>;
  static constexpr bool ignore_unrecognized_command_line_options = false;
  using Phase = helpers::Phase;
  static constexpr auto determine_next_phase =
      helpers::determine_next_phase<Metavariables>;
};

}  // namespace

static const std::vector<void (*)()> charm_init_node_funcs{
    &setup_error_handling};
static const std::vector<void (*)()> charm_init_proc_funcs{
    &enable_floating_point_exceptions};

using charmxx_main_component = Parallel::Main<Metavariables>;

#include "Parallel/CharmMain.tpp"  // IWYU pragma: keep
