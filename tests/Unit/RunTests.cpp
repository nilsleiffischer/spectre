// Distributed under the MIT License.
// See LICENSE.txt for details.

#define CATCH_CONFIG_RUNNER

#include "RunTests.hpp"

#include "Framework/TestingFramework.hpp"

#include <charm++.h>
#include <cstddef>
#include <exception>
#include <limits>
#include <memory>
#include <string>

#include "ErrorHandling/Abort.hpp"
#include "ErrorHandling/FloatingPointExceptions.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Informer/InfoFromBuild.hpp"
#include "Parallel/Exit.hpp"
#include "Parallel/Printf.hpp"
#include "tests/Unit/RunTestsRegister.hpp"

RunTests::RunTests(CkArgMsg* msg) {
  std::set_terminate([]() { abort("Called terminate. Aborting..."); });
  register_run_tests_libs();
  Parallel::printf("%s", info_from_build().c_str());
  enable_floating_point_exceptions();
  Catch::StringMaker<double>::precision =
      std::numeric_limits<double>::max_digits10;
  Catch::StringMaker<float>::precision =
      std::numeric_limits<float>::max_digits10;
  const int result = Catch::Session().run(msg->argc, msg->argv);
  // In the case where we run all the non-failure tests at once we must ensure
  // that we only initialize and finalize the python env once. Initialization is
  // done in the constructor of SetupLocalPythonEnvironment, while finalization
  // is done in the constructor of RunTests.
  pypp::SetupLocalPythonEnvironment::finalize_env();
  if (0 == result) {
    Parallel::exit();
  }
  abort("A catch test has failed.");
}

#include "tests/Unit/RunTests.def.h"  /// IWYU pragma: keep

// Needed for tests that use the GlobalCache since it registers itself with
// Charm++. However, since Parallel/CharmMain.tpp isn't included in the RunTests
// executable, no actual registration is done, the GlobalCache is only
// queued for registration.
namespace Parallel::charmxx {
class RegistrationHelper;
/// \cond
std::unique_ptr<RegistrationHelper>* charm_register_list = nullptr;
size_t charm_register_list_capacity = 0;
size_t charm_register_list_size = 0;
/// \endcond
}  // namespace Parallel::charmxx
