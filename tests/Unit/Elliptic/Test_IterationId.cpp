// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "tests/Unit/TestingFramework.hpp"

#include <functional>
#include <string>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Elliptic/IterationId.hpp"
#include "Elliptic/Tags.hpp"
#include "NumericalAlgorithms/LinearSolver/Tags.hpp"
#include "NumericalAlgorithms/NonlinearSolver/Tags.hpp"
#include "Utilities/GetOutput.hpp"
#include "tests/Unit/TestHelpers.hpp"

SPECTRE_TEST_CASE("Unit.Elliptic.IterationId", "[Unit][Elliptic]") {
  using iteration_id = Elliptic::IterationId<NonlinearSolver::Tags::IterationId,
                                             LinearSolver::Tags::IterationId>;
  {
    INFO("Comparison logic");
    check_cmp(iteration_id{1_st, 1_st}, iteration_id{1_st, 2_st});
    check_cmp(iteration_id{1_st, 1_st}, iteration_id{2_st, 1_st});
    check_cmp(iteration_id{1_st, 1_st}, iteration_id{2_st, 0_st});
  }
  {
    INFO("Hash");
    using Hash = std::hash<iteration_id>;
    CHECK(Hash{}(iteration_id{1_st, 2_st}) == Hash{}(iteration_id{1_st, 2_st}));
    CHECK(Hash{}(iteration_id{1_st, 2_st}) != Hash{}(iteration_id{2_st, 2_st}));
    CHECK(Hash{}(iteration_id{1_st, 2_st}) != Hash{}(iteration_id{1_st, 1_st}));
  }
}
