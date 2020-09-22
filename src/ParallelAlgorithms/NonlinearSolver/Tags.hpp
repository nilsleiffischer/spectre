// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines DataBox tags for the linear solver

#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "Informer/Verbosity.hpp"
#include "NumericalAlgorithms/Convergence/Criteria.hpp"
#include "NumericalAlgorithms/Convergence/HasConverged.hpp"
#include "ParallelAlgorithms/LinearSolver/Tags.hpp"
#include "Utilities/Functional.hpp"
#include "Utilities/Numeric.hpp"
#include "Utilities/Requires.hpp"
#include "Utilities/TypeTraits.hpp"

namespace NonlinearSolver::Tags {

/*
 * \brief The correction \f$\delta x\f$ to improve a solution \f$x_0\f$
 *
 * \details A linear problem \f$Ax=b\f$ can be equivalently formulated as the
 * problem \f$A\delta x=b-A x_0\f$ for the correction \f$\delta x\f$ to an
 * initial guess \f$x_0\f$. More importantly, we can use a correction scheme
 * to solve a nonlinear problem \f$A_\mathrm{nonlinear}(x)=b\f$ by repeatedly
 * solving a linearization of it. For instance, a Newton-Raphson scheme
 * iteratively refines an initial guess \f$x_0\f$ by repeatedly solving the
 * linearized problem
 * \f[
 * \frac{\delta A_\mathrm{nonlinear}}{\delta x}(x_k)\delta x_k =
 * b-A_\mathrm{nonlinear}(x_k)
 * \f]
 * for the correction \f$\delta x_k\f$ and then updating the solution as
 * \f$x_{k+1}=x_k + \delta x_k\f$.
 */
template <typename Tag>
struct Correction : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    return "Correction(" + Tag::name() + ")";
  }
  using type = typename Tag::type;
  using tag = Tag;
};

/*!
 * \brief The nonlinear operator \f$A_\mathrm{nonlinear}\f$ applied to the data
 * in `Tag`
 */
template <typename Tag>
struct OperatorAppliedTo : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Nonlinear" prefix to abbreviate the namespace for uniqueness
    return "NonlinearOperatorAppliedTo(" + db::tag_name<Tag>() + ")";
  }
  using type = typename Tag::type;
  using tag = Tag;
};

/*!
 * \brief The nonlinear residual
 * \f$r_\mathrm{nonlinear} = b - A_\mathrm{nonlinear}(\delta x)\f$
 */
template <typename Tag>
struct Residual : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Nonlinear" prefix to abbreviate the namespace for uniqueness
    return "NonlinearResidual(" + db::tag_name<Tag>() + ")";
  }
  using type = typename Tag::type;
  using tag = Tag;
};

/// Compute the residual \f$r=b - Ax\f$ from the `SourceTag` \f$b\f$ and the
/// `db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo, FieldsTag>`
/// \f$Ax\f$.
template <typename FieldsTag, typename SourceTag>
struct ResidualCompute : db::add_tag_prefix<Residual, FieldsTag>,
                         db::ComputeTag {
  using base = db::add_tag_prefix<Residual, FieldsTag>;
  using argument_tags =
      tmpl::list<SourceTag, db::add_tag_prefix<OperatorAppliedTo, FieldsTag>>;
  using return_type = typename base::type;
  static void function(
      const gsl::not_null<return_type*> residual,
      const typename SourceTag::type& source,
      const typename db::add_tag_prefix<OperatorAppliedTo, FieldsTag>::type&
          operator_applied_to_fields) noexcept {
    *residual = source - operator_applied_to_fields;
  }
};

template <typename OptionsGroup>
struct StepLength : db::SimpleTag {
  using type = double;
  static std::string name() noexcept {
    return "StepLength(" + Options::name<OptionsGroup>() + ")";
  }
};

template <typename OptionsGroup>
struct GlobalizationIterationId : db::SimpleTag {
  using type = size_t;
  static std::string name() noexcept {
    return "GlobalizationIterationId(" + Options::name<OptionsGroup>() + ")";
  }
};

}  // namespace NonlinearSolver::Tags
