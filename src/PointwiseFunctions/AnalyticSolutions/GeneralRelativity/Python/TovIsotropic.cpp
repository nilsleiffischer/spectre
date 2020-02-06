// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <string>

#include "PointwiseFunctions/AnalyticSolutions/GeneralRelativity/TovIsotropic.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"

namespace py = pybind11;

namespace gr {
namespace Solutions {
namespace py_bindings {

void bind_tov_isotropic(py::module& m) {  // NOLINT
  py::class_<TovIsotropic>(m, "TovIsotropic")
      .def(py::init<const EquationsOfState::EquationOfState<true, 1>&, double,
                    double, double, double>(),
           py::arg("equation_of_state"), py::arg("central_mass_density"),
           py::arg("log_enthalpy_at_outer_radius") = 0.,
           py::arg("absolute_tolerance") = 1.e-14,
           py::arg("relative_tolerance") = 1.e-14)
      .def("outer_radius", &TovIsotropic::outer_radius)
      .def("conformal_factor", py::vectorize(&TovIsotropic::conformal_factor))
      .def("mass", py::vectorize(&TovIsotropic::mass))
      .def("log_specific_enthalpy",
           py::vectorize(&TovIsotropic::log_specific_enthalpy));
}

}  // namespace py_bindings
}  // namespace Solutions
}  // namespace gr
