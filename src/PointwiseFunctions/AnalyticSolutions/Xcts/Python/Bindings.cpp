// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace Xcts {
namespace Solutions {
namespace py_bindings {
void bind_constant_density_star(py::module& m);  // NOLINT
void bind_schwarzschild(py::module& m);          // NOLINT
}  // namespace py_bindings

PYBIND11_MODULE(_PyXctsSolutions, m) {  // NOLINT
  py_bindings::bind_constant_density_star(m);
  py_bindings::bind_schwarzschild(m);
}

}  // namespace Solutions
}  // namespace Poisson