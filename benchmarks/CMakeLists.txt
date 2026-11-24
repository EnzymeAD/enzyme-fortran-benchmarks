# Fortran benchmarks for Enzyme.
#
# Each subdirectory holds one benchmark: ``.f90`` files carrying lit
# ``! RUN:`` / ``! CHECK:`` directives, registered with the lit testsuite
# declared in the parent CMakeLists.txt.

add_subdirectory(square)
add_subdirectory(sundials)
