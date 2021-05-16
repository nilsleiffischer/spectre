# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "PyEllipticDg")

spectre_python_add_module(
  dg
  MODULE_PATH elliptic
  LIBRARY_NAME ${LIBRARY}
  SOURCES
  ApplyPoissonOperator.cpp
  Bindings.cpp
  )

spectre_python_link_libraries(
  ${LIBRARY}
  PRIVATE
  DomainStructure
  Elliptic
  EllipticDg
  Poisson
  PoissonBoundaryConditions
  pybind11::module
  Utilities
  )
