// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Matrix.hpp"  // IWYU pragma: keep

/// \cond
template <size_t Dim>
class Mesh;
/// \endcond

namespace Spectral {

/// The portion of an element covered by a child element.
enum class ChildSize { Full, UpperHalf, LowerHalf };

/// The portion of an element covered by a mortar.
using MortarSize = ChildSize;

std::ostream& operator<<(std::ostream& os, ChildSize mortar_size) noexcept;

/*!
 * \brief The projection matrix from a child element to its parent.
 *
 * The projection matrices returned by this function (and by
 * projection_matrix_parent_to_child()) define orthogonal projection operators
 * between the spaces of functions on a parent element and its children. These
 * projections are usually the correct way to transfer data between elements in
 * a mesh-refinement hierarchy, as well as between an element face and its
 * adjacent mortars.
 *
 * These functions assume that the `child_mesh` is at least as fine as the
 * `parent_mesh`, i.e. functions on the `parent_mesh` can be represented exactly
 * on the `child_mesh`. In practice this means that functions can be projected
 * to a mortar (the `child_mesh`) from both adjacent element faces (the
 * `parent_mesh`) without losing accuracy. Similarly, functions in a
 * mesh-refinement hierarchy don't lose accuracy when an element is split
 * (h-refined). For this reason, the `projection_matrix_child_to_parent` is
 * sometimes referred to as a "restriction operator" and the
 * `projection_matrix_parent_to_child` as a "prolongation operator".
 *
 * The projections implemented here are \f$L^2\f$ projections, or "Galerkin"
 * projections. They are defined such that the projection from a parent to one
 * of its children ("prolongation") is an (exact) interpolation, and the
 * projection from a child to its parent ("restriction") is its inverse.
 * Specifically, for every function \f$u_f\f$ on the "fine" child-mesh \f$e_f\f$
 * and every function \f$u_c\f$ on the "coarse" parent-mesh \$e_c\f$ the
 * restriction operator \f$R\f$ satisfies the condition:
 *
 * \f{equation}
 * \int_{e_c} R(u_f) u_c \mathrm{d}x = \int_{e_f} u_f I(u_c) \mathrm{d}x
 * \f}
 *
 * where \f$I\f$ denotes the interpolation operator from the coarse to the fine
 * mesh. When we choose a set of basis functions \f$\phi_k(x)\f$ on the coarse
 * mesh and another on the fine mesh, we can express the restriction operator in
 * terms of mass matrices \f$M_{ij}=\int_e \phi_i(x) \phi_k(x) \mathrm{d}x\f$:
 *
 * \f{equation}
 * \quad R = M_c^{-1} I^T M_f
 * \text{.}
 * \f}
 *
 * The restriction operation is an \f$L^2\f$ projection in the sense that
 * \f$R(u_f)=u_c\f$ minimizes the quantity \f$||u_c-u_f||^2\f$. It also has an
 * interpretation of constructing a coarse-mesh representation of the fine-mesh
 * function \f$u_f\f$ by "averaging" over its modes in a mass-conservative way.
 *
 * The half-interval projections are based on an equation derived by
 * Saul.  This shows that the projection from the spectral basis for
 * the entire interval to the spectral basis for the upper half
 * interval is
 * \f{equation*}
 * T_{jk} = \frac{2 j + 1}{2} 2^j \sum_{n=0}^{j-k} \binom{j}{k+n}
 * \binom{(j + k + n - 1)/2}{j} \frac{(k + n)!^2}{(2 k + n + 1)! n!}
 * \f}
 */
const Matrix& projection_matrix_child_to_parent(const Mesh<1>& child_mesh,
                                                const Mesh<1>& parent_mesh,
                                                ChildSize size) noexcept;

/// The projection matrix from a parent element to one of its children.
///
/// \see projection_matrix_child_to_parent()
const Matrix& projection_matrix_parent_to_child(const Mesh<1>& parent_mesh,
                                                const Mesh<1>& child_mesh,
                                                ChildSize size) noexcept;

}  // namespace Spectral
