// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Elliptic/Systems/Elasticity/Equations.hpp"

#include <cstddef>

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Elasticity/FirstOrderSystem.hpp"
#include "Elliptic/Systems/Elasticity/Tags.hpp"
#include "PointwiseFunctions/Elasticity/ConstitutiveRelations/ConstitutiveRelation.hpp"
#include "PointwiseFunctions/GeneralRelativity/IndexManipulation.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/TMPL.hpp"

namespace Elasticity {

template <size_t Dim>
void primal_fluxes(
    const gsl::not_null<tnsr::IJ<DataVector, Dim>*> flux_for_displacement,
    const tnsr::ii<DataVector, Dim>& strain,
    const ConstitutiveRelations::ConstitutiveRelation<Dim>&
        constitutive_relation,
    const tnsr::I<DataVector, Dim>& coordinates) noexcept {
  const auto stress = constitutive_relation.stress(strain, coordinates);
  // To set the components of the flux each component of the symmetric stress
  // tensor is used twice. So the tensor can't be moved in its entirety.
  for (size_t d = 0; d < Dim; d++) {
    for (size_t e = 0; e < Dim; e++) {
      // Also, the stress has lower and the flux upper indices. The minus sign
      // originates in the definition of the stress \f$T^{ij} = -Y^{ijkl}
      // S_{kl}\f$.
      flux_for_displacement->get(d, e) = -stress.get(d, e);
    }
  }
}

template <size_t Dim>
void add_non_euclidean_sources(
    const gsl::not_null<tnsr::I<DataVector, Dim>*> source_for_displacement,
    const gsl::not_null<tnsr::ii<DataVector, Dim>*> source_for_strain,
    const tnsr::ijj<DataVector, Dim>& christoffel_first_kind,
    const tnsr::Ijj<DataVector, Dim>& christoffel_second_kind,
    const tnsr::i<DataVector, Dim>& christoffel_contracted,
    const tnsr::I<DataVector, Dim>& displacement,
    const tnsr::II<DataVector, Dim>& stress) noexcept {
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j < Dim; ++j) {
      source_for_displacement->get(j) -=
          christoffel_contracted.get(i) * stress.get(i, j);
      for (size_t k = 0; k < Dim; ++k) {
        source_for_displacement->get(j) -=
            christoffel_second_kind.get(j, i, k) * stress.get(i, k);
      }
    }
  }
  for (size_t i = 0; i < Dim; ++i) {
    for (size_t j = 0; j <= i; ++j) {
      for (size_t k = 0; k < Dim; ++k) {
        source_for_strain->get(i, j) +=
            christoffel_first_kind.get(k, i, j) * displacement.get(k);
      }
    }
  }
}

template <size_t Dim>
void auxiliary_fluxes(
    const gsl::not_null<tnsr::Ijj<DataVector, Dim>*> flux_for_strain,
    const tnsr::I<DataVector, Dim>& displacement) noexcept {
  std::fill(flux_for_strain->begin(), flux_for_strain->end(), 0.);
  // The off-diagonal elements are calculated by going over the upper triangular
  // matrix (the lower triangular matrix, excluding the diagonal elements, is
  // set by virtue of the tensor being symmetric in its last two indices) and
  // the symmetrisation is completed by going over the diagonal elements again.
  for (size_t d = 0; d < Dim; d++) {
    flux_for_strain->get(d, d, d) += 0.5 * displacement.get(d);
    for (size_t e = 0; e < Dim; e++) {
      flux_for_strain->get(d, e, d) += 0.5 * displacement.get(e);
    }
  }
}

template <size_t Dim>
void non_euclidean_auxiliary_fluxes(
    const gsl::not_null<tnsr::Ijj<DataVector, Dim>*> flux_for_strain,
    const tnsr::ii<DataVector, Dim>& metric,
    const tnsr::I<DataVector, Dim>& displacement) noexcept {
  const auto co_displacement = raise_or_lower_index(displacement, metric);
  std::fill(flux_for_strain->begin(), flux_for_strain->end(), 0.);
  for (size_t d = 0; d < Dim; d++) {
    flux_for_strain->get(d, d, d) += 0.5 * co_displacement.get(d);
    for (size_t e = 0; e < Dim; e++) {
      flux_for_strain->get(d, e, d) += 0.5 * co_displacement.get(e);
    }
  }
}

}  // namespace Elasticity

#define DIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATE(_, data)                                             \
  template void Elasticity::primal_fluxes<DIM(data)>(                    \
      gsl::not_null<tnsr::IJ<DataVector, DIM(data)>*>,                   \
      const tnsr::ii<DataVector, DIM(data)>&,                            \
      const Elasticity::ConstitutiveRelations::ConstitutiveRelation<DIM( \
          data)>&,                                                       \
      const tnsr::I<DataVector, DIM(data)>&) noexcept;                   \
  template void Elasticity::add_non_euclidean_sources<DIM(data)>(        \
      gsl::not_null<tnsr::I<DataVector, DIM(data)>*>,                    \
      gsl::not_null<tnsr::ii<DataVector, DIM(data)>*>,                   \
      const tnsr::ijj<DataVector, DIM(data)>&,                           \
      const tnsr::Ijj<DataVector, DIM(data)>&,                           \
      const tnsr::i<DataVector, DIM(data)>&,                             \
      const tnsr::I<DataVector, DIM(data)>&,                             \
      const tnsr::II<DataVector, DIM(data)>&) noexcept;                  \
  template void Elasticity::auxiliary_fluxes<DIM(data)>(                 \
      gsl::not_null<tnsr::Ijj<DataVector, DIM(data)>*>,                  \
      const tnsr::I<DataVector, DIM(data)>&) noexcept;                   \
  template void Elasticity::non_euclidean_auxiliary_fluxes<DIM(data)>(   \
      gsl::not_null<tnsr::Ijj<DataVector, DIM(data)>*>,                  \
      const tnsr::ii<DataVector, DIM(data)>&,                            \
      const tnsr::I<DataVector, DIM(data)>&) noexcept;

GENERATE_INSTANTIATIONS(INSTANTIATE, (2, 3))

#undef INSTANTIATE
#undef DIM
