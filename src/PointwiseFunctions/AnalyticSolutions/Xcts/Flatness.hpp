// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <limits>
#include <ostream>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "Options/Options.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/AnalyticSolution.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace Xcts::Solutions {

/// \cond
template <typename Registrars>
struct Flatness;

namespace Registrars {
struct Flatness {
  template <typename Registrars>
  using f = Solutions::Flatness<Registrars>;
};
}  // namespace Registrars
/// \endcond

/// Flat spacetime in general relativity. Useful as initial guess.
template <typename Registrars = tmpl::list<Solutions::Registrars::Flatness>>
class Flatness : public AnalyticSolution<Registrars> {
 private:
  using Base = AnalyticSolution<Registrars>;

 public:
  using options = tmpl::list<>;
  static constexpr Options::String help{
      "Flat spacetime, useful as initial guess."};

  Flatness() = default;
  Flatness(const Flatness&) noexcept = default;
  Flatness& operator=(const Flatness&) noexcept = default;
  Flatness(Flatness&&) noexcept = default;
  Flatness& operator=(Flatness&&) noexcept = default;
  ~Flatness() noexcept = default;

  /// \cond
  explicit Flatness(CkMigrateMessage* m) noexcept : Base(m) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(Flatness);
  /// \endcond

  template <typename DataType, typename... RequestedTags>
  tuples::TaggedTuple<RequestedTags...> variables(
      const tnsr::I<DataType, 3, Frame::Inertial>& x,
      tmpl::list<RequestedTags...> /*meta*/) const noexcept {
    using supported_tags_zero = tmpl::list<
        gr::Tags::TraceExtrinsicCurvature<DataType>,
        ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataType>,
                      tmpl::size_t<3>, Frame::Inertial>,
        ::Tags::deriv<Tags::ConformalFactor<DataType>, tmpl::size_t<3>,
                      Frame::Inertial>,
        ::Tags::deriv<Tags::LapseTimesConformalFactor<DataType>,
                      tmpl::size_t<3>, Frame::Inertial>,
        Tags::ShiftBackground<DataType, 3, Frame::Inertial>,
        Tags::ShiftExcess<DataType, 3, Frame::Inertial>,
        Tags::ShiftStrain<DataType, 3, Frame::Inertial>,
        gr::Tags::EnergyDensity<DataType>, gr::Tags::StressTrace<DataType>,
        gr::Tags::MomentumDensity<3, Frame::Inertial, DataType>,
        ::Tags::FixedSource<Tags::ConformalFactor<DataType>>,
        ::Tags::FixedSource<Tags::LapseTimesConformalFactor<DataType>>,
        ::Tags::FixedSource<Tags::ShiftExcess<DataType, 3, Frame::Inertial>>>;
    using supported_tags_one =
        tmpl::list<Tags::ConformalFactor<DataType>,
                   Tags::LapseTimesConformalFactor<DataType>>;
    using supported_tags = tmpl::append<
        supported_tags_zero, supported_tags_one,
        tmpl::list<Tags::ConformalMetric<DataType, 3, Frame::Inertial>>>;
    static_assert(tmpl::size<tmpl::list_difference<tmpl::list<RequestedTags...>,
                                                   supported_tags>>::value == 0,
                  "The requested tag is not supported");
    const auto make_value = [&x](auto tag_v) noexcept {
      using tag = std::decay_t<decltype(tag_v)>;
      if constexpr (tmpl::list_contains_v<supported_tags_zero, tag>) {
        return make_with_value<typename tag::type>(x, 0.);
      } else if constexpr (tmpl::list_contains_v<supported_tags_one, tag>) {
        return make_with_value<typename tag::type>(x, 1.);
      } else if constexpr (std::is_same_v<
                               tag, Tags::ConformalMetric<DataType, 3,
                                                          Frame::Inertial>>) {
        auto flat_metric = make_with_value<typename tag::type>(x, 0.);
        get<0, 0>(flat_metric) = 1.;
        get<1, 1>(flat_metric) = 1.;
        get<2, 2>(flat_metric) = 1.;
        return flat_metric;
      }
    };
    return {make_value(RequestedTags{})...};
  }
};

template <typename Registrars>
bool operator==(const Flatness<Registrars>& /*lhs*/,
                const Flatness<Registrars>& /*rhs*/) noexcept {
  return true;
}

template <typename Registrars>
bool operator!=(const Flatness<Registrars>& lhs,
                const Flatness<Registrars>& rhs) noexcept {
  return not(lhs == rhs);
}

/// \cond
template <typename Registrars>
PUP::able::PUP_ID Flatness<Registrars>::my_PUP_ID = 0;  // NOLINT
/// \endcond

}  // namespace Xcts::Solutions
