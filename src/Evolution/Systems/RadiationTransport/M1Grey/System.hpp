// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Variables.hpp"
#include "Evolution/Conservative/ConservativeDuDt.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Characteristics.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Fluxes.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Sources.hpp"
#include "Evolution/Systems/RadiationTransport/M1Grey/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace Tags {
template <class>
class Variables;
}  // namespace Tags

/// \ingroup EvolutionSystemsGroup
/// \brief Items related to general relativistic radiation transport
namespace RadiationTransport {
/// The M1 scheme for radiation transport
///
/// References:
/// - Post-merger evolution of a neutron star-black hole binary with
///   neutrino transport, \cite Foucart2015vpa
namespace M1Grey {

template <typename NeutrinoSpeciesList>
struct System;

template <typename... NeutrinoSpecies>
struct System<tmpl::list<NeutrinoSpecies...>> {
  static constexpr bool is_in_flux_conservative_form = true;
  static constexpr bool has_primitive_and_conservative_vars = false;
  static constexpr size_t volume_dim = 3;
  // If coupling to hydro, we'll want 3D equations of state
  // i.e. P(rho,T,Ye)... but this is not implemented yet.
  // For early tests of M1, we'll ignore coupling to the fluid
  // and provide analytical expressions for its 4-velocity / LorentzFactor
  //static constexpr size_t thermodynamic_dim = 3;

  using variables_tag = ::Tags::Variables<
      tmpl::list<Tags::TildeE<Frame::Physical, NeutrinoSpecies>...,
                 Tags::TildeS<Frame::Physical, NeutrinoSpecies>...>>;

  using spacetime_variables_tag = ::Tags::Variables<tmpl::list<
      gr::Tags::Lapse<DataVector>,
      gr::Tags::Shift<3, Frame::Physical, DataVector>,
      gr::Tags::SpatialMetric<3, Frame::Physical, DataVector>,
      gr::Tags::InverseSpatialMetric<3, Frame::Physical, DataVector>,
      ::Tags::deriv<gr::Tags::Lapse<DataVector>, tmpl::size_t<3>,
                    Frame::Physical>,
      ::Tags::deriv<gr::Tags::Shift<3, Frame::Physical, DataVector>,
                    tmpl::size_t<3>, Frame::Physical>,
      ::Tags::deriv<gr::Tags::SpatialMetric<3, Frame::Physical, DataVector>,
                    tmpl::size_t<3>, Frame::Physical>,
      gr::Tags::ExtrinsicCurvature<3, Frame::Physical, DataVector>>>;

  using hydro_variables_tag = ::Tags::Variables<
      tmpl::list<hydro::Tags::LorentzFactor<DataVector>,
                 hydro::Tags::SpatialVelocity<DataVector, 3, Frame::Physical>>>;

  using primitive_variables_tag = ::Tags::Variables<tmpl::list<
      Tags::ClosureFactor<NeutrinoSpecies>...,
      Tags::TildeP<Frame::Physical, NeutrinoSpecies>...,
      Tags::TildeJ<NeutrinoSpecies>..., Tags::TildeHNormal<NeutrinoSpecies>...,
      Tags::TildeHSpatial<Frame::Physical, NeutrinoSpecies>...>>;

  template <typename Tag>
  using magnitude_tag = ::Tags::NonEuclideanMagnitude<
      Tag, gr::Tags::InverseSpatialMetric<3, Frame::Physical, DataVector>>;

  using char_speeds_tag = Tags::CharacteristicSpeedsCompute;

  using volume_fluxes = ComputeFluxes<NeutrinoSpecies...>;

  using volume_sources = ComputeSources<NeutrinoSpecies...>;

  using compute_time_derivative = ConservativeDuDt<System>;

  using sourced_variables =
      tmpl::list<Tags::TildeE<Frame::Physical, NeutrinoSpecies>...,
                 Tags::TildeS<Frame::Physical, NeutrinoSpecies>...>;
};
}  // namespace M1Grey
}  // namespace RadiationTransport
