// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Tags.hpp"
#include "Elliptic/Systems/Elasticity/Tags.hpp"

namespace Elasticity {

template <size_t Dim>
void strain(gsl::not_null<tnsr::ii<DataVector, Dim>*> strain,
            const tnsr::I<DataVector, Dim>& displacement, const Mesh<Dim>& mesh,
            const InverseJacobian<DataVector, Dim, Frame::Logical,
                                  Frame::Inertial>& inv_jacobian) noexcept;

namespace Tags {

template <size_t Dim>
struct StrainCompute : Elasticity::Tags::Strain<Dim>, db::ComputeTag {
  using base = Elasticity::Tags::Strain<Dim>;
  using return_type = tnsr::ii<DataVector, Dim>;
  using argument_tags = tmpl::list<
      Elasticity::Tags::Displacement<Dim>, domain::Tags::Mesh<Dim>,
      domain::Tags::InverseJacobian<Dim, Frame::Logical, Frame::Inertial>>;
  static constexpr auto function = &strain<Dim>;
};

}  // namespace Tags
}  // namespace Elasticity
