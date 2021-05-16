// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <cstddef>
#include <pybind11/pybind11.h>
#include <string>

#include "Elliptic/DiscontinuousGalerkin/Python/ApplyOperator.hpp"
#include "Elliptic/Systems/Poisson/FirstOrderSystem.hpp"
#include "Elliptic/Systems/Poisson/Geometry.hpp"
#include "Utilities/GetOutput.hpp"

namespace py = pybind11;

namespace elliptic::dg::py_bindings {

namespace detail {
template <size_t Dim>
void bind_apply_poisson_operator_impl(py::module& m) {  // NOLINT
  using system =
      Poisson::FirstOrderSystem<Dim, Poisson::Geometry::FlatCartesian>;
  bind_apply_operator<system>(
      m, "apply_poisson_operator_flat_cartesian_" + get_output(Dim) + "d");
}
}  // namespace detail

void bind_apply_poisson_operator(py::module& m) {  // NOLINT
  detail::bind_apply_poisson_operator_impl<1>(m);
  //   detail::bind_apply_poisson_operator_impl<2>(m);
  //   detail::bind_apply_poisson_operator_impl<3>(m);
}

}  // namespace elliptic::dg::py_bindings
