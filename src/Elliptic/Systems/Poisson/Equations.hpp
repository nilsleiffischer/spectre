// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/TMPL.hpp"

#include "PointwiseFunctions/Elasticity/ConstitutiveRelations/Tags.hpp"

/// \cond
class DataVector;
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace Poisson {

/*!
 * \brief Compute the fluxes \f$F^i=\partial_i u(x)\f$ for the Poisson
 * equation on a flat spatial metric in Cartesian coordinates.
 */
template <size_t Dim>
void euclidean_fluxes(
    gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux_for_field,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& field_gradient) noexcept;

/*!
 * \brief Compute the fluxes \f$F^i=\sqrt{\gamma}\gamma^{ij}\partial_j u(x)\f$
 * for the curved-space Poisson equation on a spatial metric \f$\gamma_{ij}\f$.
 */
template <size_t Dim>
void noneuclidean_fluxes(
    gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*> flux_for_field,
    const tnsr::II<DataVector, Dim, Frame::Inertial>& inv_spatial_metric,
    const Scalar<DataVector>& det_spatial_metric,
    const tnsr::i<DataVector, Dim, Frame::Inertial>& field_gradient) noexcept;

/*!
 * \brief Compute the fluxes \f$F^i_j=\delta^i_j u(x)\f$ for the auxiliary
 * field in the first-order formulation of the Poisson equation.
 *
 * \see Poisson::FirstOrderSystem
 */
template <size_t Dim>
void auxiliary_fluxes(gsl::not_null<tnsr::Ij<DataVector, Dim, Frame::Inertial>*>
                          flux_for_gradient,
                      const Scalar<DataVector>& field) noexcept;

/*!
 * \brief Compute the fluxes \f$F^i_A\f$ for the Poisson equation on a flat
 * metric in Cartesian coordinates.
 *
 * \see Poisson::FirstOrderSystem
 */
template <size_t Dim>
struct EuclideanFluxes {
  using argument_tags = tmpl::list<Elasticity::Tags::ConstitutiveRelation<Dim>>;
  // Doesn't work yet since elliptic::Tags::FirstOrderFluxesCompute doesn't
  // support `volume_tags`.
  using volume_tags = tmpl::list<Elasticity::Tags::ConstitutiveRelation<Dim>>;

  static void apply(
      const gsl::not_null<tnsr::I<DataVector, Dim, Frame::Inertial>*>
          flux_for_field,
      const Elasticity::ConstitutiveRelations::ConstitutiveRelation<
          Dim>& /*meta*/,
      const tnsr::i<DataVector, Dim, Frame::Inertial>&
          field_gradient) noexcept {
    euclidean_fluxes(flux_for_field, field_gradient);
  }
  static void apply(
      const gsl::not_null<tnsr::Ij<DataVector, Dim, Frame::Inertial>*>
          flux_for_gradient,
      const Elasticity::ConstitutiveRelations::ConstitutiveRelation<
          Dim>& /*meta*/,
      const Scalar<DataVector>& field) noexcept {
    auxiliary_fluxes(flux_for_gradient, field);
  }
  // clang-tidy: no runtime references
  void pup(PUP::er& /*p*/) noexcept {}  // NOLINT
};

/*!
 * \brief Compute the sources \f$S_A\f$ for the Poisson equation.
 *
 * \see Poisson::FirstOrderSystem
 */
struct Sources {
  using argument_tags = tmpl::list<>;
  static void apply(const gsl::not_null<Scalar<DataVector>*> source_for_field,
                    const Scalar<DataVector>& /*field*/) noexcept {
    get(*source_for_field) = 0.;
  }
};

}  // namespace Poisson
