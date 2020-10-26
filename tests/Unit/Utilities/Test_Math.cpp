// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cmath>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "Utilities/Math.hpp"
#include "Utilities/TypeTraits.hpp"

namespace {
template <size_t N>
void test_smoothstep() {
  CAPTURE(N);
  CHECK(smoothstep<N>(0., 1., 0.) == approx(0.));
  CHECK(smoothstep<N>(0., 1., 0.5) == approx(0.5));
  CHECK(smoothstep<N>(0., 1., 1.) == approx(1.));
  CHECK_ITERABLE_APPROX(
      smoothstep<N>(0., 1., DataVector({-1., 0., 0.5, 1., 2.})),
      DataVector({0., 0., 0.5, 1., 1.}));
  CHECK_ITERABLE_APPROX(
      smoothstep<N>(-1., 1., DataVector({-2., -1., 0., 1., 2.})),
      DataVector({0., 0., 0.5, 1., 1.}));
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Utilities.Math", "[Unit][Utilities]") {
  {
    INFO("Test number_of_digits");
    CHECK(2 == number_of_digits(10));
    CHECK(1 == number_of_digits(0));
    CHECK(1 == number_of_digits(-1));
    CHECK(1 == number_of_digits(9));
    CHECK(2 == number_of_digits(-99));
  }

  {
    INFO("Test evaluate_polynomial");
    const std::vector<double> poly_coeffs{1., 2.5, 0.3, 1.5};
    CHECK_ITERABLE_APPROX(evaluate_polynomial(poly_coeffs, 0.5), 2.5125);
    CHECK_ITERABLE_APPROX(
        evaluate_polynomial(std::array<double, 4>{1., 2.5, 0.3, 1.5}, 0.5),
        2.5125);
    CHECK_ITERABLE_APPROX(
        evaluate_polynomial(poly_coeffs,
                            DataVector({-0.5, -0.1, 0., 0.8, 1., 12.})),
        DataVector({-0.3625, 0.7515, 1., 3.96, 5.3, 2666.2}));
    const std::vector<DataVector> poly_variable_coeffs{DataVector{1., 0., 2.},
                                                       DataVector{0., 2., 1.}};
    CHECK_ITERABLE_APPROX(
        evaluate_polynomial(poly_variable_coeffs, DataVector({0., 0.5, 1.})),
        DataVector({1., 1., 3.}));
  }

  {
    INFO("Test smoothstep");
    test_smoothstep<0>();
    test_smoothstep<1>();
    test_smoothstep<2>();
    test_smoothstep<3>();
    // Test a case that failed in Release builds when a static_cast was missing
    CHECK_ITERABLE_APPROX(
        smoothstep<2>(0., 2.,
                      DataVector({-2., -1.5, -1., -0.5, -0.25, 0., 0.25, 0.5,
                                  1., 1.5, 2.})),
        DataVector({0., 0., 0., 0., 0., 0., 1.605224609375e-02, 1.03515625e-01,
                    0.5, 8.96484375e-01, 1.}));
  }

  {
    INFO("Test inverse roots and step_function");
    CHECK(step_function(1.0) == 1.0);
    CHECK(step_function(0.5) == 1.0);
    CHECK(step_function(-10) == 0);
    CHECK(step_function(0.0) == 1.0);
    CHECK(step_function(0) == 1);
    CHECK(invsqrt(4.0) == 0.5);
    CHECK(invsqrt(10) == 1.0 / sqrt(10));
    CHECK(approx(invcbrt(27.0)) == (1 / 3.0));
    CHECK(approx(invcbrt(1.0 / 64.0)) == 4.0);
  }

  {
    INFO("Test sign function");
    CHECK(sgn(2) == 1);
    CHECK(sgn(0) == 0);
    CHECK(sgn(-2) == -1);
    static_assert(std::is_same_v<decltype(sgn(-2)), int>,
                  "Failed testing type of sgn");

    CHECK(sgn(2.14) == 1.0);
    CHECK(sgn(0.0) == 0.0);
    CHECK(sgn(-3.87) == -1.0);
    static_assert(std::is_same_v<decltype(sgn(2.14)), double>,
                  "Failed testing type of sgn");

    CHECK(sgn(static_cast<unsigned>(2)) == 1);
    CHECK(sgn(static_cast<unsigned>(0)) == 0);
    CHECK(sgn(static_cast<unsigned>(-1)) == 1);
    static_assert(
        std::is_same_v<decltype(sgn(static_cast<unsigned>(2))), unsigned>,
        "Failed testing type of sgn");
  }
}
