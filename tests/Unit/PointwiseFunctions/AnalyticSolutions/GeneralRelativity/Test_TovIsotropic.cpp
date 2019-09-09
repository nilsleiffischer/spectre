// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "tests/Unit/TestingFramework.hpp"

#include <cmath>
#include <cstddef>
#include <memory>
#include <pup.h>

#include "DataStructures/Tensor/Tensor.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/TovIsotropic.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "tests/Unit/TestHelpers.hpp"

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.AnalyticSolutions.Gr.TovIsotropic",
                  "[Unit][PointwiseFunctions]") {
  EquationsOfState::PolytropicFluid<true> equation_of_state{497503.72443934507,
                                                            3.0};

  const gr::Solutions::TovIsotropic solution(equation_of_state,
                                             0.0005123881685413372, 0.);

  Approx numerical_approx = Approx::custom().epsilon(1.e-5).scale(1.);
  CHECK(solution.outer_radius() == numerical_approx(9.749166076324654));
}
