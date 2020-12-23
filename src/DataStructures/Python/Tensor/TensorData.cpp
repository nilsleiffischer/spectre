// Distributed under the MIT License.
// See LICENSE.txt for details.

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TensorData.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/StdHelpers.hpp"

namespace py = pybind11;

namespace py_bindings {
void bind_tensordata(py::module& m) {  // NOLINT
  // Wrapper for TensorComponent
  py::class_<TensorComponent>(m, "TensorComponent")
      .def(py::init<std::string, DataVector>(), py::arg("name"),
           py::arg("data"))
      .def_readwrite("name", &TensorComponent::name)
      .def_readwrite("data", &TensorComponent::data)
      .def("__str__", get_output<TensorComponent>)
      .def("__repr__", get_output<TensorComponent>)
      // NOLINTNEXTLINE(misc-redundant-expression)
      .def(py::self == py::self)
      // NOLINTNEXTLINE(misc-redundant-expression)
      .def(py::self != py::self);

  /*
  py::enum_<Spectral::Basis>(m, "SpectralBasis")
      .value("Legendre", Spectral::Basis::Legendre)
      .value("Chebyshev", Spectral::Basis::Chebyshev)
      .value("FiniteDifference", Spectral::Basis::FiniteDifference)
      .export_values();
  py::enum_<Spectral::Quadrature>(m, "SpectralQuadrature")
      .value("Gauss", Spectral::Quadrature::Gauss)
      .value("GaussLobatto", Spectral::Quadrature::GaussLobatto)
      .value("CellCentered", Spectral::Quadrature::CellCentered)
      .value("FaceCentered", Spectral::Quadrature::FaceCentered)
      .export_values();
      */

  // Wrapper for ExtentsAndTensorVolumeData
  py::class_<ExtentsAndTensorVolumeData>(m, "ExtentsAndTensorVolumeData")
      .def(py::init<std::vector<size_t>, std::vector<TensorComponent>>(),
           py::arg("extents"), py::arg("components"))
      .def_readwrite("extents", &ExtentsAndTensorVolumeData::extents)
      .def_readwrite("tensor_components",
                     &ExtentsAndTensorVolumeData::tensor_components)
      .def("__str__",
           [](const ExtentsAndTensorVolumeData& extents_and_data) {
             return "(" + get_output(extents_and_data.extents) + "," +
                    get_output(extents_and_data.tensor_components) + ")";
           })
      .def("__repr__", [](const ExtentsAndTensorVolumeData& extents_and_data) {
        return "(" + get_output(extents_and_data.extents) + "," +
               get_output(extents_and_data.tensor_components) + ")";
      });
  py::class_<ElementVolumeData, ExtentsAndTensorVolumeData>(m,
                                                            "ElementVolumeData")
      .def(py::init<std::vector<size_t>, std::vector<TensorComponent>,
                    std::vector<Spectral::Basis>,
                    std::vector<Spectral::Quadrature>>(),
           py::arg("extents"), py::arg("components"), py::arg("basis"),
           py::arg("quadrature"))
      .def_readwrite("basis", &ElementVolumeData::basis)
      .def_readwrite("quadrature", &ElementVolumeData::quadrature);
}
}  // namespace py_bindings
