# Distributed under the MIT License.
# See LICENSE.txt for details.

find_package(Git)

# Get the current working branch and commit hash
if(EXISTS ${CMAKE_SOURCE_DIR}/.git AND Git_FOUND)
  execute_process(
    COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_BRANCH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe --abbrev=0 --always --tags
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
else()
  set(GIT_BRANCH "NOT_IN_GIT_REPO")
  set(GIT_COMMIT_HASH "NOT_IN_GIT_REPO")
endif()

message(STATUS "\nUseful Information:")
message(STATUS "Git Branch: " ${GIT_BRANCH})
message(STATUS "Git Hash: " ${GIT_COMMIT_HASH})
message(STATUS "Build Directory: " ${CMAKE_BINARY_DIR})
message(STATUS "Source Directory: " ${CMAKE_SOURCE_DIR})
message(STATUS "Bin Directory: " ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
message(STATUS "CMake Modules Path: " ${CMAKE_MODULE_PATH})
message(STATUS "CMAKE_CXX_FLAGS: " ${CMAKE_CXX_FLAGS})
message(STATUS "CMAKE_CXX_LINK_FLAGS: " ${CMAKE_CXX_LINK_FLAGS})
message(STATUS "CMAKE_CXX_FLAGS_DEBUG: " ${CMAKE_CXX_FLAGS_DEBUG})
message(STATUS "CMAKE_CXX_FLAGS_RELEASE: " ${CMAKE_CXX_FLAGS_RELEASE})
message(STATUS "CMAKE_EXE_LINKER_FLAGS: " ${CMAKE_EXE_LINKER_FLAGS})
message(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})
message(STATUS "CMAKE_CXX_COMPILER: " ${CMAKE_CXX_COMPILER})
message(STATUS "BUILD_SHARED_LIBS: " ${BUILD_SHARED_LIBS})
message(STATUS "USE_PCH: ${USE_PCH}")
message(STATUS "ASAN: ${ASAN}")
message(STATUS "UBSAN_UNDEFINED: ${UBSAN_UNDEFINED}")
message(STATUS "UBSAN_INTEGER: ${UBSAN_INTEGER}")
message(STATUS "USE_SYSTEM_INCLUDE: ${USE_SYSTEM_INCLUDE}")
if (PYTHONINTERP_FOUND)
  message(STATUS "Python: " ${PYTHON_EXECUTABLE})
  message(STATUS "Python Version: ${PYTHON_VERSION_STRING}")
else()
  message(STATUS "Python: Not found")
endif()
message(STATUS "BUILD_PYTHON_BINDINGS: " ${BUILD_PYTHON_BINDINGS})
if(CLANG_TIDY_BIN)
  message(STATUS "Found clang-tidy: ${CLANG_TIDY_BIN}")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(
      STATUS
      "Could not find clang-tidy even though LLVM clang is installed"
  )
endif()

if (CODE_COVERAGE)
  message(STATUS "Code coverage enabled. All prerequisites found:")
  message(STATUS "  gcov: ${GCOV}")
  message(STATUS "  lcov: ${LCOV}")
  message(STATUS "  genhtml: ${GENHTML}")
  message(STATUS "  sed: ${SED}")
endif()

if(DOXYGEN_FOUND)
  message(STATUS "Doxygen: " ${DOXYGEN_EXECUTABLE})
else()
  message(STATUS "Doxygen: Not found, documentation cannot be built.")
endif()

if(BUILD_PYTHON_BINDINGS AND "${JEMALLOC_LIB_TYPE}" STREQUAL SHARED)
  message(STATUS
    "When using the python bindings you must run python as:\n"
    "   LD_PRELOAD=${JEMALLOC_LIBRARIES} python ..."
    "Alternatively, use the system allocator by setting "
    "-D MEMORY_ALLOCATOR=SYSTEM")
endif()
