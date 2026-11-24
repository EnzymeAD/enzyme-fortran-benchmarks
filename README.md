# Enzyme Fortran benchmarks

This repository collects a variety of benchmarks for Enzyme's Fortran support.

## Building

The benchmarks are driven by LLVM's [lit](https://llvm.org/docs/CommandGuide/lit.html)
test runner, configured via CMake. Point the build at an Enzyme build tree
(with `-DENZYME_FORTRAN=ON`) and at the LLVM that Enzyme was built against:

```sh
cmake -G Ninja -S . -B build \
  -DCMAKE_Fortran_COMPILER=flang \
  -DLLVM_DIR="$(spack find --format '{prefix}' llvm)/lib/cmake/llvm" \
  -DEnzyme_BINARY_DIR=../enzyme/build/flang
```

## Running

```sh
cmake --build build --target bench-enzyme-fortran
```

The suite is also registered with CTest, so `ctest --test-dir build` works too.

Each benchmark is a `.f90` file under `benchmarks/` carrying `! RUN:` /
`! CHECK:` lit directives. The lit configuration provides the same
substitutions as Enzyme's own Fortran tests (`%fc`, `%loadFortran`,
`%loadEnzyme`, `%enzyme`, `%loadFlangEnzyme`, `%opt`, `%FileCheck`, ...), so
benchmarks can pipe flang-emitted IR through `opt` with Enzyme loaded, or let
flang run Enzyme directly via the FlangEnzyme plugin (`REQUIRES: flangenzyme`).

## Benchmarks

- [x] Square
- [ ] SUNDIALS
- [ ] Fiats
- [ ] MITgcm
