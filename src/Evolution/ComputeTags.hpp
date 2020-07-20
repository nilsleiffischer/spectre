// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/Tags.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Tags.hpp"
#include "Time/Tags.hpp"
#include "Utilities/TMPL.hpp"

namespace evolution {
namespace Tags {
/*!
 * \brief Use the `AnalyticSolutionTag` to compute the analytic solution of the
 * tags in `AnalyticFieldsTagList`.
 */
template <size_t Dim, typename AnalyticSolutionTag,
          typename AnalyticFieldsTagList>
struct AnalyticCompute
    : db::add_tag_prefix<::Tags::Analytic,
                         ::Tags::Variables<AnalyticFieldsTagList>>,
      db::ComputeTag {
  using base = db::add_tag_prefix<::Tags::Analytic,
                                  ::Tags::Variables<AnalyticFieldsTagList>>;
  using return_type = typename base::type;
  using argument_tags =
      tmpl::list<AnalyticSolutionTag,
                 domain::Tags::Coordinates<Dim, Frame::Inertial>, ::Tags::Time>;
  static void function(
      const gsl::not_null<return_type*> analytic_solution,
      const db::const_item_type<AnalyticSolutionTag>&
          analytic_solution_computer,
      const tnsr::I<DataVector, Dim, Frame::Inertial>& inertial_coords,
      const double time) noexcept {
    *analytic_solution =
        variables_from_tagged_tuple(analytic_solution_computer.variables(
            inertial_coords, time, AnalyticFieldsTagList{}));
  }
};
}  // namespace Tags
}  // namespace evolution
