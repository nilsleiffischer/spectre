// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <string>
#include <tuple>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/Schwarzschild.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace {

using field_tags =
    tmpl::list<Xcts::Tags::ConformalFactor<DataVector>,
               Xcts::Tags::LapseTimesConformalFactor<DataVector>,
               Xcts::Tags::ShiftExcess<DataVector, 3, Frame::Inertial>>;
using auxiliary_field_tags =
    tmpl::list<::Tags::deriv<Xcts::Tags::ConformalFactor<DataVector>,
                             tmpl::size_t<3>, Frame::Inertial>,
               ::Tags::deriv<Xcts::Tags::LapseTimesConformalFactor<DataVector>,
                             tmpl::size_t<3>, Frame::Inertial>,
               Xcts::Tags::ShiftStrain<DataVector, 3, Frame::Inertial>>;
using background_tags =
    tmpl::list<Xcts::Tags::ConformalMetric<DataVector, 3, Frame::Inertial>,
               gr::Tags::TraceExtrinsicCurvature<DataVector>,
               ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataVector>,
                             tmpl::size_t<3>, Frame::Inertial>,
               Xcts::Tags::ShiftBackground<DataVector, 3, Frame::Inertial>>;
using matter_source_tags =
    tmpl::list<gr::Tags::EnergyDensity<DataVector>,
               gr::Tags::StressTrace<DataVector>,
               gr::Tags::MomentumDensity<3, Frame::Inertial, DataVector>>;
using fixed_source_tags = db::wrap_tags_in<Tags::FixedSource, field_tags>;

template <Xcts::Solutions::SchwarzschildCoordinates Coords>
struct SchwarzschildProxy : Xcts::Solutions::Schwarzschild<Coords> {
  using Xcts::Solutions::Schwarzschild<Coords>::Schwarzschild;
  tuples::tagged_tuple_from_typelist<
      tmpl::append<field_tags, auxiliary_field_tags>>
  field_variables(const tnsr::I<DataVector, 3, Frame::Inertial>& x) const
      noexcept {
    return Xcts::Solutions::Schwarzschild<Coords>::variables(
        x, tmpl::append<field_tags, auxiliary_field_tags>{});
  }
  tuples::tagged_tuple_from_typelist<background_tags> background_variables(
      const tnsr::I<DataVector, 3, Frame::Inertial>& x) const noexcept {
    return Xcts::Solutions::Schwarzschild<Coords>::variables(x,
                                                             background_tags{});
  }
  tuples::tagged_tuple_from_typelist<matter_source_tags>
  matter_source_variables(
      const tnsr::I<DataVector, 3, Frame::Inertial>& x) const noexcept {
    return Xcts::Solutions::Schwarzschild<Coords>::variables(
        x, matter_source_tags{});
  }
  tuples::tagged_tuple_from_typelist<fixed_source_tags> fixed_source_variables(
      const tnsr::I<DataVector, 3, Frame::Inertial>& x) const noexcept {
    return Xcts::Solutions::Schwarzschild<Coords>::variables(
        x, fixed_source_tags{});
  }
};

template <Xcts::Solutions::SchwarzschildCoordinates Coords>
void test_solution(const double radius_at_horizon,
                   const std::string& py_functions_suffix) {
  CAPTURE(Coords);
  const SchwarzschildProxy<Coords> solution{};
  CHECK(solution.radius_at_horizon() == approx(radius_at_horizon));
  const double inner_radius = 0.5 * radius_at_horizon;
  const double outer_radius = 2. * radius_at_horizon;
  pypp::check_with_random_values<
      1, tmpl::append<field_tags, auxiliary_field_tags>>(
      &SchwarzschildProxy<Coords>::field_variables, solution, "Schwarzschild",
      {"conformal_factor_" + py_functions_suffix,
       "lapse_times_conformal_factor_" + py_functions_suffix,
       "shift_" + py_functions_suffix,
       "conformal_factor_gradient_" + py_functions_suffix,
       "lapse_times_conformal_factor_gradient_" + py_functions_suffix,
       "shift_strain_" + py_functions_suffix},
      {{{inner_radius, outer_radius}}}, std::tuple<>{}, DataVector(5));
  pypp::check_with_random_values<1, background_tags>(
      &SchwarzschildProxy<Coords>::background_variables, solution,
      "Schwarzschild",
      {"conformal_spatial_metric_" + py_functions_suffix,
       "extrinsic_curvature_trace_" + py_functions_suffix,
       "extrinsic_curvature_trace_gradient_" + py_functions_suffix,
       "shift_background"},
      {{{inner_radius, outer_radius}}}, std::tuple<>{}, DataVector(5));
  pypp::check_with_random_values<1, matter_source_tags>(
      &SchwarzschildProxy<Coords>::matter_source_variables, solution,
      "Schwarzschild", {"energy_density", "stress_trace", "momentum_density"},
      {{{inner_radius, outer_radius}}}, std::tuple<>{}, DataVector(5));
  pypp::check_with_random_values<1, fixed_source_tags>(
      &SchwarzschildProxy<Coords>::fixed_source_variables, solution,
      "Schwarzschild",
      {"conformal_factor_fixed_source",
       "lapse_times_conformal_factor_fixed_source", "shift_fixed_source"},
      {{{inner_radius, outer_radius}}}, std::tuple<>{}, DataVector(5));

  const auto created_solution =
      TestHelpers::test_creation<Xcts::Solutions::Schwarzschild<Coords>>("");
  CHECK(created_solution == solution);
  test_serialization(solution);
}

}  // namespace

SPECTRE_TEST_CASE(
    "Unit.PointwiseFunctions.AnalyticSolutions.Xcts.Schwarzschild",
    "[PointwiseFunctions][Unit]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/AnalyticSolutions/Xcts"};
  test_solution<Xcts::Solutions::SchwarzschildCoordinates::Isotropic>(
      0.5, "isotropic");
  test_solution<Xcts::Solutions::SchwarzschildCoordinates::PainleveGullstrand>(
      2., "painleve_gullstrand");
  test_solution<Xcts::Solutions::SchwarzschildCoordinates::KerrSchildIsotropic>(
      1.2727410334221052, "kerr_schild_isotropic");
}
