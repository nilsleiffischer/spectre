// Distributed under the MIT License.
// See LICENSE.txt for details.

/// \file
/// Defines DataBox tags for the linear solver

#pragma once

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataBox/TagName.hpp"
#include "DataStructures/DenseMatrix.hpp"
#include "Informer/Verbosity.hpp"
#include "Utilities/Gsl.hpp"

/*!
 * \ingroup LinearSolverGroup
 * \brief Functionality for solving linear systems of equations
 */
namespace LinearSolver {

/*!
 * \ingroup LinearSolverGroup
 * \brief The \ref DataBoxGroup tags associated with the linear solver
 */
namespace Tags {

/*!
 * \brief The operand that the local linear operator \f$A\f$ is applied to
 *
 * \details The result of the operation should be wrapped in
 * `LinearSolver::Tags::OperatorAppliedTo`.
 */
template <typename Tag>
struct Operand : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Linear" prefix to abbreviate the namespace for uniqueness
    return "LinearOperand(" + db::tag_name<Tag>() + ")";
  }
  using type = typename Tag::type;
  using tag = Tag;
};

/*!
 * \brief The linear operator \f$A\f$ applied to the data in `Tag`
 */
template <typename Tag>
struct OperatorAppliedTo : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Linear" prefix to abbreviate the namespace for uniqueness
    return "LinearOperatorAppliedTo(" + db::tag_name<Tag>() + ")";
  }
  using type = typename Tag::type;
  using tag = Tag;
};

/*!
 * \brief The residual \f$r=b - Ax\f$
 */
template <typename Tag>
struct Residual : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Linear" prefix to abbreviate the namespace for uniqueness
    return "LinearResidual(" + db::tag_name<Tag>() + ")";
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

template <typename Tag>
struct Initial : db::PrefixTag, db::SimpleTag {
  using type = typename Tag::type;
  using tag = Tag;
};

/*!
 * \brief The magnitude square \f$\langle \cdot,\cdot\rangle\f$ w.r.t.
 * the `LinearSolver::inner_product`
 */
template <typename Tag>
struct MagnitudeSquare : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Linear" prefix to abbreviate the namespace for uniqueness
    return "LinearMagnitudeSquare(" + db::tag_name<Tag>() + ")";
  }
  using type = double;
  using tag = Tag;
};

/*!
 * \brief The magnitude \f$\sqrt{\langle \cdot,\cdot\rangle}\f$ w.r.t.
 * the `LinearSolver::inner_product`
 */
template <typename Tag>
struct Magnitude : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Linear" prefix to abbreviate the namespace for uniqueness
    return "LinearMagnitude(" + db::tag_name<Tag>() + ")";
  }
  using type = double;
  using tag = Tag;
};

/*!
 * \brief The prefix for tags related to an orthogonalization procedure
 */
template <typename Tag>
struct Orthogonalization : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Linear" prefix to abbreviate the namespace for uniqueness
    return "LinearOrthogonalization(" + db::tag_name<Tag>() + ")";
  }
  using type = typename Tag::type;
  using tag = Tag;
};

/*!
 * \brief A Hessenberg matrix built up during an orthogonalization procedure
 */
template <typename Tag>
struct OrthogonalizationHistory : db::PrefixTag, db::SimpleTag {
  static std::string name() noexcept {
    // Add "Linear" prefix to abbreviate the namespace for uniqueness
    return "LinearOrthogonalizationHistory(" + db::tag_name<Tag>() + ")";
  }
  using type = DenseMatrix<double>;
  using tag = Tag;
};

/*!
 * \brief A set of \f$n\f$ vectors that form a basis of the \f$n\f$-th Krylov
 * subspace \f$K_n(A,b)\f$
 *
 * \details The Krylov subspace \f$K_n(A,b)\f$ spanned by this basis is the one
 * generated by the linear operator \f$A\f$ and source \f$b\f$ that are
 * represented by the tags
 * `db::add_tag_prefix<LinearSolver::Tags::OperatorAppliedTo,
 * db::add_tag_prefix<LinearSolver::Tags::Operand, Tag>>` and
 * `db::add_tag_prefix<::Tags::FixedSource, Tag>`, respectively. Therefore, each
 * basis vector is of the type db::add_tag_prefix<Operand, Tag>::type.
 */
template <typename Tag>
struct KrylovSubspaceBasis : db::PrefixTag, db::SimpleTag {
  // Use automatically generated name because a Krylov subspace always refers to
  // a linear operator
  using type = std::vector<typename Tag::type>;
  using tag = Tag;
};

/// Indicates the `Tag` is related to preconditioning of the linear solve
template <typename Tag>
struct Preconditioned : db::PrefixTag, db::SimpleTag {
  using type = typename Tag::type;
  using tag = Tag;
};

}  // namespace Tags
}  // namespace LinearSolver
