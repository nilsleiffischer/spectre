# Distributed under the MIT License.
# See LICENSE.txt for details.

option(BUILD_PYTHON_BINDINGS "Build the python bindings for SpECTRE" OFF)

if(BUILD_PYTHON_BINDINGS)
  if(NOT ${MEMORY_ALLOCATOR} STREQUAL SYSTEM)
    message(FATAL_ERROR "Python bindings currently require using the system "
      "memory allocator. Set MEMORY_ALLOCATOR=SYSTEM.")
  endif()

  # Make sure to find the Python interpreter first, so it is consistent with
  # the one that pybind11 uses
  find_package(PythonInterp)
  find_package(PythonLibs)
  # Uses `Findpybind11.cmake` to find the headers. Since we can't rely on the
  # corresponding cmake files to be installed as well we bundle them in
  # `external/pybind11`.
  find_package(pybind11 REQUIRED)

  # Load the CMake files from `external/pybind11`
  list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/external/pybind11/tools")
  include(pybind11Tools)

  set_property(
    GLOBAL APPEND PROPERTY SPECTRE_THIRD_PARTY_LIBS
    pybind11
    )

  message(STATUS "Pybind11 include: ${pybind11_INCLUDE_DIR}")
endif()
