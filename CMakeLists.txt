# CMake build system for the Enzyme Fortran benchmarks.
#
# This project does not build the benchmarks themselves; the CMake project
# exists to configure LLVM's lit site configuration and to drive LLVM's
# ``lit`` test runner against the Fortran benchmarks, mimicking the
# build/test infrastructure used by Enzyme itself (enzyme/test/CMakeLists.txt)
# and by enzyme-benchmarks (the C++ benchmark suite).

cmake_minimum_required(VERSION 3.20)
# NOTE: C needed for LLVMConfig.cmake's internal checks
project(enzyme-fortran-benchmarks LANGUAGES C Fortran)

# Locate the same LLVM installation that Enzyme was built against. Required
# so that we can find LLVM's CMake helpers (e.g. ``configure_lit_site_cfg``,
# ``add_lit_testsuite``) and the LLVM tools (``opt``, ``FileCheck``, ...).
find_package(LLVM REQUIRED CONFIG)
message(STATUS "Found LLVM ${LLVM_PACKAGE_VERSION}: ${LLVM_DIR}")

list(APPEND CMAKE_MODULE_PATH "${LLVM_CMAKE_DIR}")
include(AddLLVM)

# Locate the Enzyme build tree so that the benchmarks can reference the
# Enzyme LLVM/flang plugins and the Fortran enzyme module directory. Either
# point this at the build directory directly (the directory containing
# ``Enzyme/LLVMEnzyme-<ver>.so`` etc.) or at a standard Enzyme install
# prefix.
set(Enzyme_BINARY_DIR "" CACHE PATH "Path to the Enzyme build directory")
if(NOT Enzyme_BINARY_DIR)
  message(FATAL_ERROR "Enzyme_BINARY_DIR must be set to the Enzyme build directory")
endif()

# Locate the ``lit`` test runner. AddLLVM's helpers expect this in the form
# of the LLVM_EXTERNAL_LIT cache FILEPATH; we always override it from the
# environment, picking up whatever 'lit' is on PATH (typically the one in
# the active Python virtualenv).
find_program(ENZYME_BENCH_LIT_EXECUTABLE NAMES llvm-lit lit REQUIRED)
set(LLVM_EXTERNAL_LIT "${ENZYME_BENCH_LIT_EXECUTABLE}"
    CACHE FILEPATH "Path to the lit test runner" FORCE)
message(STATUS "Found lit: ${LLVM_EXTERNAL_LIT}")

enable_testing()

# The SUNDIALS benchmarks link against SUNDIALS (including its Fortran module
# interface). Build the pinned SUNDIALS submodule once and expose its Fortran
# module include dir and link flags to the lit configuration as the
# %sundialsFortran and %sundialsLibs substitutions. Set SUNDIALS_PREFIX to an
# existing SUNDIALS install (with the Fortran module interface enabled) to skip
# building the submodule.
set(SUNDIALS_PREFIX "" CACHE PATH "Path to an existing SUNDIALS install with the Fortran module interface; skips building dependencies/sundials")
set(SUNDIALS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/sundials")
set(SUNDIALS_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/dependencies/sundials")
set(SUNDIALS_INSTALL_DIR "${CMAKE_CURRENT_BINARY_DIR}/dependencies/sundials-install")
set(BENCH_SUNDIALS_LIB_DIR "")
if(SUNDIALS_PREFIX)
  set(BENCH_SUNDIALS_FORTRAN_INCLUDE "${SUNDIALS_PREFIX}/fortran")
  set(BENCH_SUNDIALS_LIB_DIR "${SUNDIALS_PREFIX}/lib")
  set(BENCH_SUNDIALS_DEP_TARGET "")
else()
  if(NOT EXISTS "${SUNDIALS_SOURCE_DIR}/CMakeLists.txt")
    message(FATAL_ERROR
      "SUNDIALS source not found at ${SUNDIALS_SOURCE_DIR}; "
      "init the dependencies/sundials submodule or set SUNDIALS_PREFIX")
  endif()
  include(ExternalProject)
  ExternalProject_Add(bench-sundials
    SOURCE_DIR "${SUNDIALS_SOURCE_DIR}"
    BINARY_DIR "${SUNDIALS_BINARY_DIR}"
    INSTALL_DIR "${SUNDIALS_INSTALL_DIR}"
    CMAKE_ARGS
      -G ${CMAKE_GENERATOR}
      -DCMAKE_BUILD_TYPE=Release
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DBUILD_FORTRAN_MODULE_INTERFACE=ON
      -DBUILD_STATIC_LIBS=OFF
      -DBUILD_ARKODE=OFF -DBUILD_CVODES=OFF -DBUILD_IDA=OFF -DBUILD_IDAS=OFF
      -DBUILD_KINSOL=OFF
      -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_CXX=OFF
      -DEXAMPLES_ENABLE_F2003=OFF -DEXAMPLES_INSTALL=OFF
      -DBUILD_TESTING=OFF
      -DCMAKE_INSTALL_PREFIX=${SUNDIALS_INSTALL_DIR}
    BUILD_ALWAYS ON)
  set(BENCH_SUNDIALS_FORTRAN_INCLUDE "${SUNDIALS_INSTALL_DIR}/fortran")
  set(BENCH_SUNDIALS_DEP_TARGET bench-sundials)
endif()
# SUNDIALS installs its libs under lib or lib64 depending on the platform;
# after building the submodule, pick whichever exists. (With SUNDIALS_PREFIX
# the install is expected to use lib, as SUNDIALS' own installs do.)
if(NOT BENCH_SUNDIALS_LIB_DIR)
  set(BENCH_SUNDIALS_LIB_DIR "${SUNDIALS_INSTALL_DIR}/lib")
endif()
if(NOT SUNDIALS_PREFIX)
  if(NOT EXISTS "${SUNDIALS_INSTALL_DIR}/lib" AND EXISTS "${SUNDIALS_INSTALL_DIR}/lib64")
    set(BENCH_SUNDIALS_LIB_DIR "${SUNDIALS_INSTALL_DIR}/lib64")
  endif()
endif()

# Libraries the Robertson CVODE benchmarks link against, in link order.
set(BENCH_SUNDIALS_LIBS
    fcvode_mod cvode fnvecserial_mod nvecserial fsunlinsolspgmr_mod sunlinsolspgmr
    fsunlinsoldense_mod sunlinsoldense fsunmatrixdense_mod sunmatrixdense fcore_mod core)
set(BENCH_SUNDIALS_LIB_FLAGS
    "-L${BENCH_SUNDIALS_LIB_DIR} -Wl,-rpath,${BENCH_SUNDIALS_LIB_DIR}")
foreach(lib ${BENCH_SUNDIALS_LIBS})
  string(APPEND BENCH_SUNDIALS_LIB_FLAGS " -lsundials_${lib}")
endforeach()
add_custom_target(bench-sundials-deps
  COMMAND ${CMAKE_COMMAND} -E echo
          "SUNDIALS available in ${BENCH_SUNDIALS_LIB_DIR}, Fortran modules in ${BENCH_SUNDIALS_FORTRAN_INCLUDE}"
  DEPENDS ${BENCH_SUNDIALS_DEP_TARGET})

# The OpenAD benchmarks build examples shipped with the OpenAD/F tool
# (benchmarks/openad). Download the last OpenAD snapshot once at configure time
# and unpack those examples into the build tree; they are exposed to lit as
# %openad. Set BENCH_OPENAD=OFF to skip them, e.g. without network access.
option(BENCH_OPENAD "Download the OpenAD/F examples for the OpenAD benchmarks" ON)
set(BENCH_OPENAD_EXAMPLES "")
if(BENCH_OPENAD)
  set(OPENAD_TARBALL "${CMAKE_CURRENT_BINARY_DIR}/dependencies/OpenAD_2014-03-15.tgz")
  file(DOWNLOAD
    "http://www.mcs.anl.gov/~utke/OpenAD_tars/493/OpenAD_2014-03-15.tgz"
    "${OPENAD_TARBALL}"
    EXPECTED_HASH SHA256=1e41bc4388e3c7f99d64dde11e22e801d6705a0c6049b30c04b48f89d8b28800
    SHOW_PROGRESS)
  set(OPENAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/dependencies/openad")
  file(ARCHIVE_EXTRACT INPUT "${OPENAD_TARBALL}" DESTINATION "${OPENAD_DIR}"
    PATTERNS OpenAD/COPYRIGHT OpenAD/Examples/BoxModel OpenAD/Examples/FDC
             OpenAD/Examples/ShallowWater)
  set(BENCH_OPENAD_EXAMPLES "${OPENAD_DIR}/OpenAD/Examples")
endif()

# The OpenAD shallow water example reads its input with netCDF-Fortran's
# Fortran 77 interface (netcdf.inc).
find_path(BENCH_NETCDF_INCLUDE_DIR netcdf.inc)
find_library(BENCH_NETCDFF_LIBRARY netcdff)
find_library(BENCH_NETCDF_LIBRARY netcdf)
set(BENCH_NETCDF_INCLUDE "")
set(BENCH_NETCDF_LIBS "")
if(BENCH_NETCDF_INCLUDE_DIR AND BENCH_NETCDFF_LIBRARY AND BENCH_NETCDF_LIBRARY)
  set(BENCH_NETCDF_INCLUDE "-I${BENCH_NETCDF_INCLUDE_DIR}")
  set(BENCH_NETCDF_LIBS "${BENCH_NETCDFF_LIBRARY} ${BENCH_NETCDF_LIBRARY}")
endif()

# These variables are referenced by lit.site.cfg.py.in.
set(BENCH_Fortran_COMPILER "${CMAKE_Fortran_COMPILER}")
set(BENCH_ENZYME_BINARY_DIR "${Enzyme_BINARY_DIR}")
set(BENCH_LLVM_TOOLS_DIR "${LLVM_TOOLS_BINARY_DIR}")
set(BENCH_LLVM_VERSION_MAJOR "${LLVM_VERSION_MAJOR}")
set(BENCH_LLVM_SHLIBEXT "${LLVM_PLUGIN_EXT}")

configure_lit_site_cfg(
  ${CMAKE_CURRENT_SOURCE_DIR}/lit.site.cfg.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/lit.site.cfg.py MAIN_CONFIG
  ${CMAKE_CURRENT_SOURCE_DIR}/lit.cfg.py)

# Unlike enzyme-benchmarks, which is built from within Enzyme's own build
# tree, this project is standalone, so the Enzyme plugin targets do not exist
# here. Instead, check that the plugin files the lit substitutions reference
# are present and expose the check as a target the lit suite can depend on.
set(ENZYME_BENCH_PLUGIN
    "${Enzyme_BINARY_DIR}/Enzyme/LLVMEnzyme-${LLVM_VERSION_MAJOR}${LLVM_PLUGIN_EXT}")
set(ENZYME_BENCH_MODULES_DIR "${Enzyme_BINARY_DIR}/modules")
add_custom_target(bench-enzyme-deps
  COMMAND ${CMAKE_COMMAND} -E echo
          "Enzyme plugins found in ${Enzyme_BINARY_DIR}"
  DEPENDS "${ENZYME_BENCH_PLUGIN}" bench-sundials-deps)
if(NOT EXISTS "${ENZYME_BENCH_PLUGIN}")
  message(FATAL_ERROR
    "LLVMEnzyme plugin not found at ${ENZYME_BENCH_PLUGIN}; "
    "point Enzyme_BINARY_DIR at an Enzyme build with ENZYME_FORTRAN=ON")
endif()
if(NOT EXISTS "${ENZYME_BENCH_MODULES_DIR}/enzyme.mod")
  message(FATAL_ERROR
    "Fortran enzyme module not found in ${ENZYME_BENCH_MODULES_DIR}; "
    "point Enzyme_BINARY_DIR at an Enzyme build with ENZYME_FORTRAN=ON")
endif()

# Run the benchmark suite
add_lit_testsuite(bench-enzyme-fortran "Running Enzyme Fortran benchmarks"
    ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS bench-enzyme-deps
    ARGS -v -j 1
)

set_target_properties(bench-enzyme-fortran PROPERTIES FOLDER "bench Tests")

# Also register with CTest so ``ctest`` picks the suite up from the build
# directory. (add_lit_testsuite only creates a build target.)
add_test(
  NAME bench-enzyme-fortran-lit
  COMMAND ${LLVM_EXTERNAL_LIT} -sv ${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(benchmarks)
