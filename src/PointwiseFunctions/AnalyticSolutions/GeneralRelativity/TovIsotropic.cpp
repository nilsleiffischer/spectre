// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/TovIsotropic.hpp"

// Need Boost MultiArray because it is used internally by ODEINT
#include "DataStructures/BoostMultiArray.hpp"  // IWYU pragma: keep

#include <algorithm>
#include <array>
#include <boost/numeric/odeint.hpp>  // IWYU pragma: keep
#include <cmath>
#include <cstddef>
#include <functional>
#include <ostream>
#include <pup.h>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "ErrorHandling/Assert.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Tov.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"

#include "Parallel/Printf.hpp"

namespace gr {
namespace Solutions {

namespace {

template <typename F>
void dr_conformal_factor_potential(
    gsl::not_null<double*> dr_conformal_factor_potential,
    const double areal_radius, const F& mass_over_areal_radius_at_areal_radius,
    const double central_mass_density) noexcept {
  // See class documentation for choice of cutoff
  const double radial_expansion_cutoff =
      cbrt(std::numeric_limits<double>::epsilon());
  if (UNLIKELY(areal_radius < radial_expansion_cutoff)) {
    *dr_conformal_factor_potential =
        -4. / 3. * M_PI * central_mass_density * areal_radius;
  } else {
    const double sqrt_one_minus_two_m_over_r =
        sqrt(1. - 2. * mass_over_areal_radius_at_areal_radius(areal_radius));
    *dr_conformal_factor_potential =
        (sqrt_one_minus_two_m_over_r - 1.) /
        (areal_radius * sqrt_one_minus_two_m_over_r);
  }
}

class IntegralObserver {
 public:
  void operator()(const double current_conformal_factor_potential,
                  const double current_areal_radius) noexcept {
    areal_radius.push_back(current_areal_radius);
    conformal_factor_potential.push_back(current_conformal_factor_potential);
  }
  std::vector<double> areal_radius;
  std::vector<double> conformal_factor_potential;
};

template <typename DataType>
TovIsotropic::RadialVariables<DataType> interior_solution(
    const EquationsOfState::EquationOfState<true, 1>& equation_of_state,
    const DataType& radius, const DataType& conformal_factor,
    const DataType& log_specific_enthalpy, const DataType& mass_within_radius,
    const double log_lapse_at_outer_radius) noexcept {
  TovIsotropic::RadialVariables<DataType> result{radius};
  result.conformal_factor = Scalar<DataType>{conformal_factor};
  DataType areal_radius = square(conformal_factor) * radius;
  DataType m_over_areal_radius = mass_within_radius / areal_radius;
  DataType sqrt_one_minus_two_m_over_areal_radius =
      sqrt(1. - 2. * m_over_areal_radius);
  result.dr_conformal_factor =
      Scalar<DataType>{0.5 * conformal_factor / radius *
                       (sqrt_one_minus_two_m_over_areal_radius - 1.)};
  result.specific_enthalpy = Scalar<DataType>{exp(log_specific_enthalpy)};
  result.rest_mass_density = equation_of_state.rest_mass_density_from_enthalpy(
      result.specific_enthalpy);
  result.pressure =
      equation_of_state.pressure_from_density(result.rest_mass_density);
  result.specific_internal_energy =
      equation_of_state.specific_internal_energy_from_density(
          result.rest_mass_density);
  result.lapse =
      Scalar<DataType>{exp(log_lapse_at_outer_radius - log_specific_enthalpy)};
  result.dr_lapse =
      Scalar<DataType>{get(result.lapse) / radius *
                       (m_over_areal_radius + 4 * M_PI * get(result.pressure) *
                                                  square(areal_radius)) /
                       sqrt_one_minus_two_m_over_areal_radius};
  return result;
}

template <typename DataType>
TovIsotropic::RadialVariables<DataType> vacuum_solution(
    const DataType& radius, const double total_mass) noexcept {
  TovIsotropic::RadialVariables<DataType> result{radius};
  const DataType m_over_two_r = 0.5 * total_mass / radius;
  const DataType one_plus_m_over_two_r = 1. + m_over_two_r;
  result.conformal_factor = Scalar<DataType>{one_plus_m_over_two_r};
  result.dr_conformal_factor =
      Scalar<DataType>{-0.5 * total_mass / square(radius)};
  result.specific_enthalpy = make_with_value<Scalar<DataType>>(radius, 1.);
  result.rest_mass_density = make_with_value<Scalar<DataType>>(radius, 0.);
  result.pressure = make_with_value<Scalar<DataType>>(radius, 0.);
  result.specific_internal_energy =
      make_with_value<Scalar<DataType>>(radius, 0.);
  result.lapse = Scalar<DataType>{(1. - m_over_two_r) / one_plus_m_over_two_r};
  result.dr_lapse = Scalar<DataType>(
      m_over_two_r / radius / one_plus_m_over_two_r * (1. + get(result.lapse)));
  return result;
}

}  // namespace

TovIsotropic::TovIsotropic(
    const EquationsOfState::EquationOfState<true, 1>& equation_of_state,
    const double central_mass_density,
    const double log_enthalpy_at_outer_radius, const double absolute_tolerance,
    const double relative_tolerance) noexcept
    : areal_solution_{equation_of_state, central_mass_density,
                      log_enthalpy_at_outer_radius, absolute_tolerance,
                      relative_tolerance} {
  const double areal_outer_radius = areal_solution_.outer_radius();
  total_mass_ = areal_solution_.mass(areal_outer_radius);
  log_lapse_at_outer_radius_ =
      0.5 * log(1. - 2. * total_mass_ / areal_outer_radius);
  const double outer_radius =
      0.5 * (areal_outer_radius - total_mass_ +
             sqrt(square(areal_outer_radius) -
                  2. * total_mass_ * areal_outer_radius));
  // Is it possible to just reference the member function?
  const auto mass_over_areal_radius_at_areal_radius = [this](
      const double areal_radius) noexcept {
    return this->areal_solution_.mass_over_radius(areal_radius);
  };

  double current_areal_radius = 0.;
  double current_conformal_factor_potential = 0.;  // ?
  double current_dr_conformal_factor_potential{};
  dr_conformal_factor_potential(
      &current_dr_conformal_factor_potential, current_areal_radius,
      mass_over_areal_radius_at_areal_radius, central_mass_density);

  const double initial_step = areal_outer_radius / 100.;  // ?
  using StateDopri5 = boost::numeric::odeint::runge_kutta_dopri5<double>;
  boost::numeric::odeint::dense_output_runge_kutta<
      boost::numeric::odeint::controlled_runge_kutta<StateDopri5>>
      dopri5 = make_dense_output(absolute_tolerance, relative_tolerance,
                                 StateDopri5{});
  IntegralObserver integral_observer{};
  boost::numeric::odeint::integrate_adaptive(
      dopri5,
      [&mass_over_areal_radius_at_areal_radius, &central_mass_density ](
          const double /*local_conformal_factor_potential*/,
          double& local_dr_conformal_factor_potential,
          const double local_areal_radius) noexcept {
        return dr_conformal_factor_potential(
            &local_dr_conformal_factor_potential, local_areal_radius,
            mass_over_areal_radius_at_areal_radius, central_mass_density);
      },
      current_conformal_factor_potential, current_areal_radius,
      areal_outer_radius, initial_step, std::ref(integral_observer));

  const double matching_constant =
      outer_radius / areal_outer_radius *
      exp(integral_observer.conformal_factor_potential.back());

  const size_t num_points = integral_observer.areal_radius.size();
  std::vector<double> isotropic_radius(num_points);
  std::vector<double> conformal_factor(num_points);
  for (size_t i = 0; i < num_points; i++) {
    const double conformal_factor_square =
        exp(integral_observer.conformal_factor_potential[i]) /
        matching_constant;
    conformal_factor[i] = sqrt(conformal_factor_square);
    isotropic_radius[i] =
        integral_observer.areal_radius[i] / conformal_factor_square;
  }
  conformal_factor_interpolant_ =
      intrp::BarycentricRational(isotropic_radius, conformal_factor, 5);
  // This should be equal to the outer_radius calculated before, but differs
  // by a numerical error. We need it to be exact for checking bounds of the
  // interpolation.
  outer_radius_ = isotropic_radius.back();
}

double TovIsotropic::outer_radius() const noexcept { return outer_radius_; }

template <>
TovIsotropic::RadialVariables<double> TovIsotropic::radial_variables(
    const EquationsOfState::EquationOfState<true, 1>& equation_of_state,
    const tnsr::I<double, 3>& x) const noexcept {
  // add small number to avoid FPEs at origin
  const double radius = get(magnitude(x)) + 1.e-30 * outer_radius_;
  if (radius >= outer_radius_) {
    return vacuum_solution(radius, total_mass_);
  }
  return interior_solution(equation_of_state, radius, conformal_factor(radius),
                           log_specific_enthalpy(radius), mass(radius),
                           log_lapse_at_outer_radius_);
}

template <>
TovIsotropic::RadialVariables<DataVector> TovIsotropic::radial_variables(
    const EquationsOfState::EquationOfState<true, 1>& equation_of_state,
    const tnsr::I<DataVector, 3>& x) const noexcept {
  // add small number to avoid FPEs at origin
  const DataVector radius = get(magnitude(x)) + 1.e-30 * outer_radius_;
  if (min(radius) >= outer_radius_) {
    return vacuum_solution(radius, total_mass_);
  }
  const size_t num_points = radius.size();
  if (max(radius) <= outer_radius_) {
    DataVector conformal_factor_at_radii{num_points};
    DataVector log_of_specific_enthalpy{num_points};
    DataVector mass_within_radius{num_points};
    for (size_t i = 0; i < num_points; i++) {
      const double r = radius[i];
      conformal_factor_at_radii[i] = conformal_factor(r);
      log_of_specific_enthalpy[i] = log_specific_enthalpy(r);
      mass_within_radius[i] = mass(r);
    }
    return interior_solution(equation_of_state, radius,
                             conformal_factor_at_radii,
                             log_of_specific_enthalpy, mass_within_radius,
                             log_lapse_at_outer_radius_);
  }
  RadialVariables<DataVector> result{radius};
  for (size_t i = 0; i < num_points; i++) {
    const double r = radius[i];
    auto radial_vars_at_r =
        (r <= outer_radius_
             ? interior_solution(equation_of_state, r, conformal_factor(r),
                                 log_specific_enthalpy(r), mass(r),
                                 log_lapse_at_outer_radius_)
             : vacuum_solution(r, total_mass_));
    get(result.conformal_factor)[i] = get(radial_vars_at_r.conformal_factor);
    get(result.dr_conformal_factor)[i] =
        get(radial_vars_at_r.dr_conformal_factor);
    get(result.specific_enthalpy)[i] = get(radial_vars_at_r.specific_enthalpy);
    get(result.rest_mass_density)[i] = get(radial_vars_at_r.rest_mass_density);
    get(result.pressure)[i] = get(radial_vars_at_r.pressure);
    get(result.specific_internal_energy)[i] =
        get(radial_vars_at_r.specific_internal_energy);
    get(result.lapse)[i] = get(radial_vars_at_r.lapse);
    get(result.dr_lapse)[i] = get(radial_vars_at_r.dr_lapse);
  }
  return result;
}

void TovIsotropic::pup(PUP::er& p) noexcept {  // NOLINT
  p | areal_solution_;
  p | outer_radius_;
  p | total_mass_;
  p | log_lapse_at_outer_radius_;
  p | conformal_factor_interpolant_;
}

}  // namespace Solutions
}  // namespace gr
