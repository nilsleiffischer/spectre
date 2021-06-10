// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>
#include <utility>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/FaceNormal.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/BoundaryConditions/ApplyBoundaryCondition.hpp"
#include "Elliptic/BoundaryConditions/BoundaryCondition.hpp"
#include "Elliptic/Systems/Xcts/BoundaryConditions/ApparentHorizon.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/NormalDotFlux.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/Kerr.hpp"
#include "Utilities/TMPL.hpp"

namespace Xcts::BoundaryConditions {

namespace {
// Make a metric approximately Riemannian
void make_metric_riemannian(
    const gsl::not_null<tnsr::II<DataVector, 3>*> inv_conformal_metric) {
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < i; ++j) {
      inv_conformal_metric->get(i, j) *= 1.e-2;
    }
    inv_conformal_metric->get(i, i) = abs(inv_conformal_metric->get(i, i));
  }
}

// Generate a face normal pretending the surface is a coordinate sphere
std::tuple<tnsr::i<DataVector, 3>, tnsr::ij<DataVector, 3>, Scalar<DataVector>>
make_spherical_face_normal(
    tnsr::I<DataVector, 3> x, const std::array<double, 3>& center,
    const tnsr::II<DataVector, 3>& inv_conformal_metric) {
  for (size_t d = 0; d < 3; ++d) {
    x.get(d) -= gsl::at(center, d);
  }
  Scalar<DataVector> proper_radius{x.begin()->size(), 0.};
  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      get(proper_radius) +=
          inv_conformal_metric.get(i, j) * x.get(i) * x.get(j);
    }
  }
  get(proper_radius) = sqrt(get(proper_radius));
  tnsr::i<DataVector, 3> face_normal{x.begin()->size()};
  get<0>(face_normal) = -get<0>(x) / get(proper_radius);
  get<1>(face_normal) = -get<1>(x) / get(proper_radius);
  get<2>(face_normal) = -get<2>(x) / get(proper_radius);
  tnsr::ij<DataVector, 3> deriv_unnormalized_face_normal{x.begin()->size(), 0.};
  get<0, 0>(deriv_unnormalized_face_normal) = -1.;
  get<1, 1>(deriv_unnormalized_face_normal) = -1.;
  get<2, 2>(deriv_unnormalized_face_normal) = -1.;
  return {std::move(face_normal), std::move(deriv_unnormalized_face_normal),
          std::move(proper_radius)};
}
std::tuple<tnsr::i<DataVector, 3>, tnsr::ij<DataVector, 3>, Scalar<DataVector>>
make_spherical_face_normal_flat_cartesian(tnsr::I<DataVector, 3> x,
                                          const std::array<double, 3>& center) {
  for (size_t d = 0; d < 3; ++d) {
    x.get(d) -= gsl::at(center, d);
  }
  Scalar<DataVector> euclidean_radius = magnitude(x);
  tnsr::i<DataVector, 3> face_normal{x.begin()->size()};
  get<0>(face_normal) = -get<0>(x) / get(euclidean_radius);
  get<1>(face_normal) = -get<1>(x) / get(euclidean_radius);
  get<2>(face_normal) = -get<2>(x) / get(euclidean_radius);
  tnsr::ij<DataVector, 3> deriv_unnormalized_face_normal{x.begin()->size(), 0.};
  get<0, 0>(deriv_unnormalized_face_normal) = -1.;
  get<1, 1>(deriv_unnormalized_face_normal) = -1.;
  get<2, 2>(deriv_unnormalized_face_normal) = -1.;
  return {std::move(face_normal), std::move(deriv_unnormalized_face_normal),
          std::move(euclidean_radius)};
}

template <Xcts::Geometry ConformalGeometry, bool Linearized, typename... Args>
void apply_boundary_condition_impl(
    const gsl::not_null<Scalar<DataVector>*> n_dot_conformal_factor_gradient,
    const gsl::not_null<Scalar<DataVector>*>
        n_dot_lapse_times_conformal_factor_gradient,
    const gsl::not_null<tnsr::I<DataVector, 3>*> shift_excess,
    Scalar<DataVector> conformal_factor,
    Scalar<DataVector> lapse_times_conformal_factor,
    tnsr::I<DataVector, 3> n_dot_longitudinal_shift_excess,
    const std::array<double, 3>& center, const std::array<double, 3>& spin,
    Args&&... args) {
  const ApparentHorizon<ConformalGeometry> boundary_condition{center, spin,
                                                              std::nullopt};
  const auto direction = Direction<3>::lower_xi();
  const auto box = db::create<make_interface_tags<
      tmpl::conditional_t<
          Linearized,
          typename ApparentHorizon<ConformalGeometry>::argument_tags_linearized,
          typename ApparentHorizon<ConformalGeometry>::argument_tags>,
      domain::Tags::BoundaryDirectionsInterior<3>,
      tmpl::conditional_t<
          Linearized,
          typename ApparentHorizon<ConformalGeometry>::volume_tags_linearized,
          typename ApparentHorizon<ConformalGeometry>::volume_tags>>>(
      std::unordered_map<Direction<3>, std::decay_t<decltype(args)>>{
          {direction, std::move(args)}}...);
  elliptic::apply_boundary_condition<Linearized>(
      boundary_condition, box, direction, make_not_null(&conformal_factor),
      make_not_null(&lapse_times_conformal_factor), shift_excess,
      n_dot_conformal_factor_gradient,
      n_dot_lapse_times_conformal_factor_gradient,
      make_not_null(&n_dot_longitudinal_shift_excess));
}

void apply_boundary_condition(
    const gsl::not_null<Scalar<DataVector>*> n_dot_conformal_factor_gradient,
    const gsl::not_null<Scalar<DataVector>*>
        n_dot_lapse_times_conformal_factor_gradient,
    const gsl::not_null<tnsr::I<DataVector, 3>*> shift_excess,
    const Scalar<DataVector>& conformal_factor,
    const Scalar<DataVector>& lapse_times_conformal_factor,
    const tnsr::I<DataVector, 3>& n_dot_longitudinal_shift_excess,
    const std::array<double, 3>& center, const std::array<double, 3>& spin,
    const tnsr::I<DataVector, 3>& x,
    const Scalar<DataVector>& extrinsic_curvature_trace,
    const tnsr::I<DataVector, 3>& shift_background,
    const tnsr::II<DataVector, 3>& longitudinal_shift_background,
    tnsr::II<DataVector, 3> inv_conformal_metric,
    const tnsr::Ijj<DataVector, 3>& conformal_christoffel_second_kind) {
  make_metric_riemannian(make_not_null(&inv_conformal_metric));
  auto [face_normal, deriv_unnormalized_face_normal, face_normal_magnitude] =
      make_spherical_face_normal(x, center, inv_conformal_metric);
  apply_boundary_condition_impl<Xcts::Geometry::Curved, false>(
      n_dot_conformal_factor_gradient,
      n_dot_lapse_times_conformal_factor_gradient, shift_excess,
      conformal_factor, lapse_times_conformal_factor,
      n_dot_longitudinal_shift_excess, center, spin, std::move(face_normal),
      std::move(deriv_unnormalized_face_normal),
      std::move(face_normal_magnitude), x, extrinsic_curvature_trace,
      shift_background, longitudinal_shift_background,
      std::move(inv_conformal_metric), conformal_christoffel_second_kind);
}

void apply_boundary_condition_flat_cartesian(
    const gsl::not_null<Scalar<DataVector>*> n_dot_conformal_factor_gradient,
    const gsl::not_null<Scalar<DataVector>*>
        n_dot_lapse_times_conformal_factor_gradient,
    const gsl::not_null<tnsr::I<DataVector, 3>*> shift_excess,
    const Scalar<DataVector>& conformal_factor,
    const Scalar<DataVector>& lapse_times_conformal_factor,
    const tnsr::I<DataVector, 3>& n_dot_longitudinal_shift_excess,
    const std::array<double, 3>& center, const std::array<double, 3>& spin,
    const tnsr::I<DataVector, 3>& x,
    const Scalar<DataVector>& extrinsic_curvature_trace,
    const tnsr::I<DataVector, 3>& shift_background,
    const tnsr::II<DataVector, 3>& longitudinal_shift_background) {
  auto [face_normal, deriv_unnormalized_face_normal, face_normal_magnitude] =
      make_spherical_face_normal_flat_cartesian(x, center);
  apply_boundary_condition_impl<Xcts::Geometry::FlatCartesian, false>(
      n_dot_conformal_factor_gradient,
      n_dot_lapse_times_conformal_factor_gradient, shift_excess,
      conformal_factor, lapse_times_conformal_factor,
      n_dot_longitudinal_shift_excess, center, spin, std::move(face_normal),
      std::move(deriv_unnormalized_face_normal),
      std::move(face_normal_magnitude), x, extrinsic_curvature_trace,
      shift_background, longitudinal_shift_background);
}

void apply_boundary_condition_linearized(
    const gsl::not_null<Scalar<DataVector>*>
        n_dot_conformal_factor_gradient_correction,
    const gsl::not_null<Scalar<DataVector>*>
        n_dot_lapse_times_conformal_factor_gradient_correction,
    const gsl::not_null<tnsr::I<DataVector, 3>*> shift_excess_correction,
    const Scalar<DataVector>& conformal_factor_correction,
    const Scalar<DataVector>& lapse_times_conformal_factor_correction,
    const tnsr::I<DataVector, 3>& n_dot_longitudinal_shift_excess_correction,
    const std::array<double, 3>& center, const std::array<double, 3>& spin,
    const tnsr::I<DataVector, 3>& x,
    const Scalar<DataVector>& extrinsic_curvature_trace,
    const tnsr::II<DataVector, 3>& longitudinal_shift_background,
    const Scalar<DataVector>& conformal_factor,
    const Scalar<DataVector>& lapse_times_conformal_factor,
    const tnsr::I<DataVector, 3>& n_dot_longitudinal_shift_excess,
    tnsr::II<DataVector, 3> inv_conformal_metric,
    const tnsr::Ijj<DataVector, 3>& conformal_christoffel_second_kind) {
  make_metric_riemannian(make_not_null(&inv_conformal_metric));
  auto [face_normal, deriv_unnormalized_face_normal, face_normal_magnitude] =
      make_spherical_face_normal(x, center, inv_conformal_metric);
  apply_boundary_condition_impl<Xcts::Geometry::Curved, true>(
      n_dot_conformal_factor_gradient_correction,
      n_dot_lapse_times_conformal_factor_gradient_correction,
      shift_excess_correction, conformal_factor_correction,
      lapse_times_conformal_factor_correction,
      n_dot_longitudinal_shift_excess_correction, center, spin,
      std::move(face_normal), std::move(deriv_unnormalized_face_normal),
      std::move(face_normal_magnitude), extrinsic_curvature_trace,
      longitudinal_shift_background, conformal_factor,
      lapse_times_conformal_factor, n_dot_longitudinal_shift_excess,
      std::move(inv_conformal_metric), conformal_christoffel_second_kind);
}

void apply_boundary_condition_linearized_flat_cartesian(
    const gsl::not_null<Scalar<DataVector>*>
        n_dot_conformal_factor_gradient_correction,
    const gsl::not_null<Scalar<DataVector>*>
        n_dot_lapse_times_conformal_factor_gradient_correction,
    const gsl::not_null<tnsr::I<DataVector, 3>*> shift_excess_correction,
    const Scalar<DataVector>& conformal_factor_correction,
    const Scalar<DataVector>& lapse_times_conformal_factor_correction,
    const tnsr::I<DataVector, 3>& n_dot_longitudinal_shift_excess_correction,
    const std::array<double, 3>& center, const std::array<double, 3>& spin,
    const tnsr::I<DataVector, 3>& x,
    const Scalar<DataVector>& extrinsic_curvature_trace,
    const tnsr::II<DataVector, 3>& longitudinal_shift_background,
    const Scalar<DataVector>& conformal_factor,
    const Scalar<DataVector>& lapse_times_conformal_factor,
    const tnsr::I<DataVector, 3>& n_dot_longitudinal_shift_excess) {
  auto [face_normal, deriv_unnormalized_face_normal, face_normal_magnitude] =
      make_spherical_face_normal_flat_cartesian(x, center);
  apply_boundary_condition_impl<Xcts::Geometry::FlatCartesian, true>(
      n_dot_conformal_factor_gradient_correction,
      n_dot_lapse_times_conformal_factor_gradient_correction,
      shift_excess_correction, conformal_factor_correction,
      lapse_times_conformal_factor_correction,
      n_dot_longitudinal_shift_excess_correction, center, spin,
      std::move(face_normal), std::move(deriv_unnormalized_face_normal),
      std::move(face_normal_magnitude), extrinsic_curvature_trace,
      longitudinal_shift_background, conformal_factor,
      lapse_times_conformal_factor, n_dot_longitudinal_shift_excess);
}

}  // namespace

SPECTRE_TEST_CASE("Unit.Xcts.BoundaryConditions.ApparentHorizon",
                  "[Unit][Elliptic]") {
  // Test factory-creation
  const auto created = TestHelpers::test_creation<
      std::unique_ptr<elliptic::BoundaryConditions::BoundaryCondition<
          3, tmpl::list<Registrars::ApparentHorizon<Xcts::Geometry::Curved>>>>>(
      "ApparentHorizon:\n"
      "  Center: [1., 2., 3.]\n"
      "  Spin: [0.1, 0.2, 0.3]\n"
      "  Mass: Auto\n");
  REQUIRE(dynamic_cast<const ApparentHorizon<Xcts::Geometry::Curved>*>(
              created.get()) != nullptr);
  const auto& boundary_condition =
      dynamic_cast<const ApparentHorizon<Xcts::Geometry::Curved>&>(*created);
  {
    INFO("Properties");
    CHECK(boundary_condition.center() == std::array<double, 3>{{1., 2., 3.}});
    CHECK(boundary_condition.spin() == std::array<double, 3>{{0.1, 0.2, 0.3}});
    CHECK(boundary_condition.mass() == std::nullopt);
  }
  {
    INFO("Semantics");
    test_serialization(boundary_condition);
    test_copy_semantics(boundary_condition);
    auto move_boundary_condition = boundary_condition;
    test_move_semantics(std::move(move_boundary_condition), boundary_condition);
  }
  {
    INFO("Random-value tests");
    pypp::SetupLocalPythonEnvironment local_python_env(
        "Elliptic/Systems/Xcts/BoundaryConditions/");
    pypp::check_with_random_values<11>(
        &apply_boundary_condition, "ApparentHorizon",
        {"normal_dot_conformal_factor_gradient",
         "normal_dot_lapse_times_conformal_factor_gradient", "shift_excess"},
        {{{0.5, 2.},
          {0.5, 2.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.}}},
        DataVector{3});
    pypp::check_with_random_values<9>(
        &apply_boundary_condition_flat_cartesian, "ApparentHorizon",
        {"normal_dot_conformal_factor_gradient_flat_cartesian",
         "normal_dot_lapse_times_conformal_factor_gradient",
         "shift_excess_flat_cartesian"},
        {{{0.5, 2.},
          {0.5, 2.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.}}},
        DataVector{3});
    pypp::check_with_random_values<13>(
        &apply_boundary_condition_linearized, "ApparentHorizon",
        {"normal_dot_conformal_factor_gradient_correction",
         "normal_dot_lapse_times_conformal_factor_gradient",
         "shift_excess_correction"},
        {{{0.5, 2.},
          {0.5, 2.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {0.5, 2.},
          {0.5, 2.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.}}},
        DataVector{3});
    pypp::check_with_random_values<11>(
        &apply_boundary_condition_linearized_flat_cartesian, "ApparentHorizon",
        {"normal_dot_conformal_factor_gradient_correction_flat_cartesian",
         "normal_dot_lapse_times_conformal_factor_gradient",
         "shift_excess_correction_flat_cartesian"},
        {{{0.5, 2.},
          {0.5, 2.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {-1., 1.},
          {0.5, 2.},
          {0.5, 2.},
          {-1., 1.}}},
        DataVector{3});
  }
  {
    INFO("Consistency with Kerr solution");
    const double mass = 1.2;
    const std::array<double, 3> center{{0., 0., 0.}};
    // TODO: Add a spin
    const std::array<double, 3> spin{{0., 0., 0.}};
    const Solutions::Kerr<> solution{mass, spin, center};
    const ApparentHorizon<Xcts::Geometry::Curved> kerr_horizon{center, spin,
                                                               mass};
    const DataVector used_for_size{3};
    // Choose an arbitrary set of points on the horizon surface
    MAKE_GENERATOR(generator, 3636894815);
    std::uniform_real_distribution<> dist_phi(0., 2. * M_PI);
    std::uniform_real_distribution<> dist_theta(0., M_PI);
    auto phi = make_with_random_values<DataVector>(
        make_not_null(&generator), make_not_null(&dist_phi), used_for_size);
    CAPTURE(phi);
    auto theta = make_with_random_values<DataVector>(
        make_not_null(&generator), make_not_null(&dist_theta), used_for_size);
    CAPTURE(theta);
    const double horizon_coord_radius = 2. * mass;
    tnsr::I<DataVector, 3> x{used_for_size.size()};
    get<0>(x) = horizon_coord_radius * cos(phi) * sin(theta);
    get<1>(x) = horizon_coord_radius * sin(phi) * sin(theta);
    get<2>(x) = horizon_coord_radius * cos(theta);
    // Get background fields from the solution
    const auto background_fields = solution.variables(
        x,
        tmpl::list<Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>,
                   Tags::ConformalChristoffelSecondKind<DataVector, 3,
                                                        Frame::Inertial>,
                   gr::Tags::TraceExtrinsicCurvature<DataVector>,
                   Tags::ShiftBackground<DataVector, 3, Frame::Inertial>,
                   Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
                       DataVector, 3, Frame::Inertial>>{});
    const auto& inv_conformal_metric =
        get<Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>>(
            background_fields);
    // Set up the face normal on the horizon surface. It points _into_ the
    // horizon because the computational domain fills the space outside of it
    // and the normal always points away from the computational domain. Since
    // the surface is a coordinate sphere the normalized face normal is just n_i
    // = -x^i / sqrt(gamma^ij x^i x^j)
    const auto [face_normal, deriv_unnormalized_face_normal,
                face_normal_magnitude] =
        make_spherical_face_normal(x, center, inv_conformal_metric);
    // Retrieve the expected surface vars and fluxes from the solution
    const auto surface_vars_expected =
        variables_from_tagged_tuple(solution.variables(
            x,
            tmpl::list<Tags::ConformalFactor<DataVector>,
                       Tags::LapseTimesConformalFactor<DataVector>,
                       Tags::ShiftExcess<DataVector, 3, Frame::Inertial>>{}));
    const auto surface_fluxes_expected =
        variables_from_tagged_tuple(solution.variables(
            x,
            tmpl::list<::Tags::Flux<Tags::ConformalFactor<DataVector>,
                                    tmpl::size_t<3>, Frame::Inertial>,
                       ::Tags::Flux<Tags::LapseTimesConformalFactor<DataVector>,
                                    tmpl::size_t<3>, Frame::Inertial>,
                       Tags::LongitudinalShiftExcess<DataVector, 3,
                                                     Frame::Inertial>>{}));
    Variables<tmpl::list<
        ::Tags::NormalDotFlux<Tags::ConformalFactor<DataVector>>,
        ::Tags::NormalDotFlux<Tags::LapseTimesConformalFactor<DataVector>>,
        ::Tags::NormalDotFlux<
            Tags::ShiftExcess<DataVector, 3, Frame::Inertial>>>>
        n_dot_surface_fluxes_expected{used_for_size.size()};
    normal_dot_flux(make_not_null(&n_dot_surface_fluxes_expected), face_normal,
                    surface_fluxes_expected);
    // Apply the boundary conditions, passing garbage for the data that the
    // boundary conditions are expected to fill
    auto surface_vars = surface_vars_expected;
    auto n_dot_surface_fluxes = n_dot_surface_fluxes_expected;
    get(get<::Tags::NormalDotFlux<Tags::ConformalFactor<DataVector>>>(
        n_dot_surface_fluxes)) = std::numeric_limits<double>::signaling_NaN();
    get<0>(get<Tags::ShiftExcess<DataVector, 3, Frame::Inertial>>(
        surface_vars)) = std::numeric_limits<double>::signaling_NaN();
    get<1>(get<Tags::ShiftExcess<DataVector, 3, Frame::Inertial>>(
        surface_vars)) = std::numeric_limits<double>::signaling_NaN();
    get<2>(get<Tags::ShiftExcess<DataVector, 3, Frame::Inertial>>(
        surface_vars)) = std::numeric_limits<double>::signaling_NaN();
    kerr_horizon.apply(
        make_not_null(&get<Tags::ConformalFactor<DataVector>>(surface_vars)),
        make_not_null(
            &get<Tags::LapseTimesConformalFactor<DataVector>>(surface_vars)),
        make_not_null(&get<Tags::ShiftExcess<DataVector, 3, Frame::Inertial>>(
            surface_vars)),
        make_not_null(
            &get<::Tags::NormalDotFlux<Tags::ConformalFactor<DataVector>>>(
                n_dot_surface_fluxes)),
        make_not_null(&get<::Tags::NormalDotFlux<
                          Tags::LapseTimesConformalFactor<DataVector>>>(
            n_dot_surface_fluxes)),
        make_not_null(&get<::Tags::NormalDotFlux<
                          Tags::ShiftExcess<DataVector, 3, Frame::Inertial>>>(
            n_dot_surface_fluxes)),
        face_normal, deriv_unnormalized_face_normal, face_normal_magnitude, x,
        get<gr::Tags::TraceExtrinsicCurvature<DataVector>>(background_fields),
        get<Tags::ShiftBackground<DataVector, 3, Frame::Inertial>>(
            background_fields),
        get<Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
            DataVector, 3, Frame::Inertial>>(background_fields),
        get<Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>>(
            background_fields),
        get<Tags::ConformalChristoffelSecondKind<DataVector, 3,
                                                 Frame::Inertial>>(
            background_fields));
    // Check the result.
    CHECK_VARIABLES_APPROX(surface_vars, surface_vars_expected);
    CHECK_VARIABLES_APPROX(n_dot_surface_fluxes, n_dot_surface_fluxes_expected);
  }
}

}  // namespace Xcts::BoundaryConditions
