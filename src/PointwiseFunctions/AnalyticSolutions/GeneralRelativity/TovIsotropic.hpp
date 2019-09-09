// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>

#include "DataStructures/Tensor/Tensor.hpp"
#include "NumericalAlgorithms/Interpolation/BarycentricRational.hpp"
#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/Tov.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace gr {
namespace Solutions {

class TovIsotropic {
 public:
  TovIsotropic(
      const EquationsOfState::EquationOfState<true, 1>& equation_of_state,
      double central_mass_density, double log_enthalpy_at_outer_radius = 0.0,
      double absolute_tolerance = 1.0e-14,
      double relative_tolerance = 1.0e-14) noexcept;

  TovIsotropic() = default;
  TovIsotropic(const TovIsotropic& /*rhs*/) = delete;
  TovIsotropic& operator=(const TovIsotropic& /*rhs*/) = delete;
  TovIsotropic(TovIsotropic&& /*rhs*/) noexcept = default;
  TovIsotropic& operator=(TovIsotropic&& /*rhs*/) noexcept = default;
  ~TovIsotropic() = default;

  /// \brief The outer isotropic radius of the solution.
  ///
  /// \note This is the radius at which `log_specific_enthalpy` is equal
  /// to the value of `log_enthalpy_at_outer_radius` that was given when
  /// constructing this TovSolution
  double outer_radius() const noexcept;

  double conformal_factor(const double r) const noexcept {
    ASSERT(
        r >= 0.0 and r <= outer_radius_,
        "Invalid radius: " << r << " not in [0.0, " << outer_radius_ << "]\n");
    return conformal_factor_interpolant_(r);
  }

  double log_specific_enthalpy(const double r) const noexcept {
    ASSERT(
        r >= 0.0 and r <= outer_radius_,
        "Invalid radius: " << r << " not in [0.0, " << outer_radius_ << "]\n");
    return areal_solution_.log_specific_enthalpy(square(conformal_factor(r)) *
                                                 r);
  }

  double mass(const double r) const noexcept {
    ASSERT(
        r >= 0.0 and r <= outer_radius_,
        "Invalid radius: " << r << " not in [0.0, " << outer_radius_ << "]\n");
    return areal_solution_.mass(square(conformal_factor(r)) * r);
  }

  template <typename DataType>
  struct RadialVariables {
    explicit RadialVariables(DataType isotropic_radius_in) noexcept
        : isotropic_radius(std::move(isotropic_radius_in)),
          specific_enthalpy(
              make_with_value<Scalar<DataType>>(isotropic_radius, 0.)),
          rest_mass_density(
              make_with_value<Scalar<DataType>>(isotropic_radius, 0.)),
          pressure(make_with_value<Scalar<DataType>>(isotropic_radius, 0.)),
          specific_internal_energy(
              make_with_value<Scalar<DataType>>(isotropic_radius, 0.)),
          lapse(make_with_value<Scalar<DataType>>(isotropic_radius, 0.)),
          dr_lapse(make_with_value<Scalar<DataType>>(isotropic_radius, 0.)),
          conformal_factor(
              make_with_value<Scalar<DataType>>(isotropic_radius, 0.)),
          dr_conformal_factor(
              make_with_value<Scalar<DataType>>(isotropic_radius, 0.)) {}
    DataType isotropic_radius{};
    Scalar<DataType> specific_enthalpy{};
    Scalar<DataType> rest_mass_density{};
    Scalar<DataType> pressure{};
    Scalar<DataType> specific_internal_energy{};
    Scalar<DataType> lapse{};
    Scalar<DataType> dr_lapse{};
    Scalar<DataType> conformal_factor{};
    Scalar<DataType> dr_conformal_factor{};
  };

  /// \brief The radial variables from which the hydrodynamic quantities and
  /// spacetime metric can be computed.
  ///
  /// For radii greater than the outer_radius, this returns the appropriate
  /// vacuum spacetime.
  ///
  /// \note This solution of the TOV equations is a function of areal radius.
  template <typename DataType>
  RadialVariables<DataType> radial_variables(
      const EquationsOfState::EquationOfState<true, 1>& equation_of_state,
      const tnsr::I<DataType, 3>& x) const noexcept;

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) noexcept;

 private:
  TovSolution areal_solution_{};
  double outer_radius_{std::numeric_limits<double>::signaling_NaN()};
  double total_mass_{std::numeric_limits<double>::signaling_NaN()};
  double log_lapse_at_outer_radius_{
      std::numeric_limits<double>::signaling_NaN()};
  intrp::BarycentricRational conformal_factor_interpolant_;
};

}  // namespace Solutions
}  // namespace gr
