// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "tests/Unit/TestingFramework.hpp"

#include <cstddef>
#include <tuple>

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"  // IWYU pragma: keep
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Domain/CoordinateMaps/Affine.hpp"
#include "Domain/CoordinateMaps/CoordinateMap.hpp"
#include "Domain/CoordinateMaps/ProductMaps.hpp"
#include "Domain/Mesh.hpp"
#include "Elliptic/Systems/Poisson/FirstOrderSystem.hpp"
#include "Elliptic/Systems/Poisson/Tags.hpp"  // IWYU pragma: keep
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Poisson/Moustache.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"
#include "tests/Unit/PointwiseFunctions/AnalyticSolutions/FirstOrderEllipticSolutionsTestHelpers.hpp"
#include "tests/Unit/Pypp/CheckWithRandomValues.hpp"
#include "tests/Unit/Pypp/SetupLocalPythonEnvironment.hpp"
#include "tests/Unit/TestCreation.hpp"
#include "tests/Unit/TestHelpers.hpp"

namespace {

template <size_t Dim>
struct MoustacheProxy : Poisson::Solutions::Moustache<Dim> {
  using Poisson::Solutions::Moustache<Dim>::Moustache;

  using field_tags = tmpl::list<
      Poisson::Tags::Field,
      ::Tags::deriv<Poisson::Tags::Field, tmpl::size_t<Dim>, Frame::Inertial>>;
  using source_tags = tmpl::list<Tags::FixedSource<Poisson::Tags::Field>>;

  tuples::tagged_tuple_from_typelist<field_tags> field_variables(
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x) const noexcept {
    return Poisson::Solutions::Moustache<Dim>::variables(x, field_tags{});
  }

  tuples::tagged_tuple_from_typelist<source_tags> source_variables(
      const tnsr::I<DataVector, Dim, Frame::Inertial>& x) const noexcept {
    return Poisson::Solutions::Moustache<Dim>::variables(x, source_tags{});
  }
};

template <size_t Dim, typename... Maps>
void test_solution(const Mesh<Dim>& mesh,
                   const domain::CoordinateMap<Frame::Logical, Frame::Inertial,
                                               Maps...>& coord_map,
                   const double tolerance) {
  const MoustacheProxy<Dim> solution{};
  pypp::check_with_random_values<
      1, tmpl::list<Poisson::Tags::Field,
                    ::Tags::deriv<Poisson::Tags::Field, tmpl::size_t<Dim>,
                                  Frame::Inertial>>>(
      &MoustacheProxy<Dim>::field_variables, solution, "Moustache",
      {"field", "field_gradient"}, {{{0., 1.}}}, std::make_tuple(),
      DataVector(5));
  pypp::check_with_random_values<
      1, tmpl::list<Tags::FixedSource<Poisson::Tags::Field>>>(
      &MoustacheProxy<Dim>::source_variables, solution, "Moustache", {"source"},
      {{{0., 1.}}}, std::make_tuple(), DataVector(5));

  Poisson::Solutions::Moustache<Dim> created_solution =
      test_creation<Poisson::Solutions::Moustache<Dim>>("  ");
  CHECK(created_solution == solution);
  test_serialization(solution);

  FirstOrderEllipticSolutionsTestHelpers::verify_solution<
      Poisson::FirstOrderSystem<Dim>>(
      created_solution, typename Poisson::FirstOrderSystem<Dim>::fluxes{}, mesh,
      coord_map, tolerance);
}

}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.AnalyticSolutions.Poisson.Moustache",
                  "[PointwiseFunctions][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/AnalyticSolutions/Poisson"};

  using AffineMap = domain::CoordinateMaps::Affine;
  test_solution(
      Mesh<1>{10, Spectral::Basis::Legendre,
              Spectral::Quadrature::GaussLobatto},
      domain::CoordinateMap<Frame::Logical, Frame::Inertial, AffineMap>{
          {-1., 1., 0., 1.}},
      1.e-1);

  using AffineMap2D =
      domain::CoordinateMaps::ProductOf2Maps<AffineMap, AffineMap>;
  test_solution(
      Mesh<2>{10, Spectral::Basis::Legendre,
              Spectral::Quadrature::GaussLobatto},
      domain::CoordinateMap<Frame::Logical, Frame::Inertial, AffineMap2D>{
          {{-1., 1., 0., 1.}, {-1., 1., 0., 1.}}},
      1.e-1);
}
