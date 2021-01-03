// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Parallel/CharmPupable.hpp"
#include "PointwiseFunctions/AnalyticData/Elasticity/AnalyticData.hpp"
#include "Utilities/FakeVirtual.hpp"
#include "Utilities/Registration.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Elasticity {
/// Analytic solutions of the linear Elasticity equations
namespace Solutions {

/// Base class for analytic solutions of the linear Elasticity equations
template <size_t Dim, typename Registrars>
class AnalyticSolution
    : public Elasticity::AnalyticData::AnalyticData<Dim, Registrars> {
 private:
  using Base = Elasticity::AnalyticData::AnalyticData<Dim, Registrars>;

 protected:
  /// \cond
  AnalyticSolution() = default;
  AnalyticSolution(const AnalyticSolution&) = default;
  AnalyticSolution(AnalyticSolution&&) = default;
  AnalyticSolution& operator=(const AnalyticSolution&) = default;
  AnalyticSolution& operator=(AnalyticSolution&&) = default;
  /// \endcond

 public:
  ~AnalyticSolution() override = default;

  using registrars = Registrars;
  using creatable_classes = Registration::registrants<registrars>;

  /// \cond
  explicit AnalyticSolution(CkMigrateMessage* m) noexcept : Base(m) {}
  WRAPPED_PUPable_abstract(AnalyticSolution);  // NOLINT
  /// \endcond

  template <typename DataType, typename... RequestedTags>
  tuples::TaggedTuple<RequestedTags...> variables(
      const tnsr::I<DataType, Dim, Frame::Inertial>& x,
      tmpl::list<RequestedTags...> /*meta*/) const noexcept {
    return call_with_dynamic_type<
        tuples::TaggedTuple<RequestedTags...>,
        tmpl::filter<
            creatable_classes,
            std::is_base_of<AnalyticSolution<Dim, Registrars>, tmpl::_1>>>(
        this, [&x](auto* const derived) noexcept {
          return derived->variables(x, tmpl::list<RequestedTags...>{});
        });
  }
};

}  // namespace Solutions
}  // namespace Elasticity
