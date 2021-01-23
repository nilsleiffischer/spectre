// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <functional>
#include <utility>

#include "DataStructures/ApplyMatrices.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/DenseMatrix.hpp"
#include "DataStructures/DenseVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/NumericalAlgorithms/LinearSolver/TestHelpers.hpp"
#include "NumericalAlgorithms/LinearSolver/ExplicitInverse.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/ElementCenteredSubdomainData.hpp"
#include "ParallelAlgorithms/LinearSolver/Schwarz/OverlapHelpers.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"

namespace helpers = TestHelpers::LinearSolver;

namespace {
struct ScalarFieldTag {
  using type = Scalar<DataVector>;
};
}  // namespace

namespace LinearSolver::Serial {

SPECTRE_TEST_CASE("Unit.LinearSolver.Serial.ExplicitInverse",
                  "[Unit][NumericalAlgorithms][LinearSolver]") {
  {
    INFO("Solve a simple matrix");
    const DenseMatrix<double> matrix{{4., 1.}, {3., 1.}};
    const helpers::ApplyMatrix linear_operator{matrix};
    const DenseVector<double> source{1., 2.};
    const DenseVector<double> expected_solution{-1., 5.};
    DenseVector<double> solution(2);
    const ExplicitInverse<> solver{true};
    const auto has_converged =
        solver.solve(make_not_null(&solution), linear_operator, source);
    REQUIRE(has_converged);
    DenseMatrix<double> matrix_repr{solver.matrix_representation()};
    CHECK_MATRIX_APPROX(matrix_repr, blaze::inv(matrix));
    CHECK_ITERABLE_APPROX(solution, expected_solution);
    {
      INFO("Resetting");
      ExplicitInverse<> resetting_solver{true};
      resetting_solver.solve(make_not_null(&solution), linear_operator, source);
      // Solving a different operator after resetting should work
      resetting_solver.reset();
      const DenseMatrix<double> matrix2{{4., 1.}, {1., 3.}};
      const helpers::ApplyMatrix linear_operator2{matrix2};
      const DenseVector<double> expected_solution2{0.0909090909090909,
                                                   0.6363636363636364};
      resetting_solver.solve(make_not_null(&solution), linear_operator2,
                             source);
      matrix_repr = resetting_solver.matrix_representation();
      CHECK_MATRIX_APPROX(matrix_repr, blaze::inv(matrix2));
      CHECK_ITERABLE_APPROX(solution, expected_solution2);
      // When resetting is disabled, the solver should keep applying the cached
      // inverse even when solving a different operator
      ExplicitInverse<> non_resetting_solver{false};
      non_resetting_solver.solve(make_not_null(&solution), linear_operator,
                                 source);
      non_resetting_solver.reset();
      non_resetting_solver.solve(make_not_null(&solution), linear_operator2,
                                 source);
      // Still the inverse of the operator we solved first
      matrix_repr = non_resetting_solver.matrix_representation();
      CHECK_MATRIX_APPROX(matrix_repr, blaze::inv(matrix));
      CHECK_ITERABLE_APPROX(solution, blaze::inv(matrix) * source);
    }
  }
  {
    INFO("Solve a heterogeneous data structure");
    using SubdomainData = ::LinearSolver::Schwarz::ElementCenteredSubdomainData<
        1, tmpl::list<ScalarFieldTag>>;

    const Matrix matrix_element{{4., 1., 1.}, {1., 1., 3.}, {0., 2., 0.}};
    const Matrix matrix_overlap{{4., 1.}, {3., 1.}};
    const ::LinearSolver::Schwarz::OverlapId<1> overlap_id{
        Direction<1>::lower_xi(), ElementId<1>{0}};
    const std::array<std::reference_wrapper<const Matrix>, 1> matrices_element{
        matrix_element};
    const std::array<std::reference_wrapper<const Matrix>, 1> matrices_overlap{
        matrix_overlap};
    const auto linear_operator = [&matrices_element, &matrices_overlap,
                                  &overlap_id](
                                     const gsl::not_null<SubdomainData*> result,
                                     const SubdomainData& operand) noexcept {
      apply_matrices(make_not_null(&result->element_data), matrices_element,
                     operand.element_data, Index<1>{3});
      apply_matrices(make_not_null(&result->overlap_data.at(overlap_id)),
                     matrices_overlap, operand.overlap_data.at(overlap_id),
                     Index<1>{2});
    };

    SubdomainData source{3};
    get(get<ScalarFieldTag>(source.element_data)) = DataVector{1., 2., 1.};
    source.overlap_data.emplace(overlap_id,
                                typename SubdomainData::OverlapData{2});
    get(get<ScalarFieldTag>(source.overlap_data.at(overlap_id))) =
        DataVector{1., 2.};
    auto expected_solution = make_with_value<SubdomainData>(source, 0.);
    get(get<ScalarFieldTag>(expected_solution.element_data)) =
        DataVector{0., 0.5, 0.5};
    get(get<ScalarFieldTag>(expected_solution.overlap_data.at(overlap_id))) =
        DataVector{-1., 5.};

    const ExplicitInverse<> solver{true};
    auto solution = make_with_value<SubdomainData>(source, 0.);
    solver.solve(make_not_null(&solution), linear_operator, source);
    CHECK(solver.size() == 5);
    DenseMatrix<double> expected_matrix(5, 5, 0.);
    blaze::submatrix(expected_matrix, 0, 0, 3, 3) = matrix_element;
    blaze::submatrix(expected_matrix, 3, 3, 2, 2) = matrix_overlap;
    blaze::invert(expected_matrix);
    const DenseMatrix<double> matrix_repr{solver.matrix_representation()};
    CHECK_MATRIX_APPROX(matrix_repr, expected_matrix);
    CHECK_VARIABLES_APPROX(solution.element_data,
                           expected_solution.element_data);
    CHECK_VARIABLES_APPROX(solution.overlap_data.at(overlap_id),
                           expected_solution.overlap_data.at(overlap_id));
  }
}

}  // namespace LinearSolver::Serial
