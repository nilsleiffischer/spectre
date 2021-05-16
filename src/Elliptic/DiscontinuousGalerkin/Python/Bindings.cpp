// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace elliptic::dg {

namespace py_bindings {
void bind_apply_poisson_operator(py::module& m);  // NOLINT
}  // namespace py_bindings

PYBIND11_MODULE(_PyEllipticDg, m) {  // NOLINT
  py_bindings::bind_apply_poisson_operator(m);
}

}  // namespace domain
