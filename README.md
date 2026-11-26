# Enzyme Fortran benchmarks

This repository collects a variety of benchmarks for Enzyme's Fortran support.

## Dependencies

- **Enzyme**: a build configured with `-DENZYME_FORTRAN=ON`, which provides
  the `LLVMEnzyme` plugin and the Fortran `enzyme` module. The benchmarks that
  let flang run Enzyme also need `-DENZYME_FLANG=ON` (the `FlangEnzyme` plugin
  and its config file), and the OpenAD benchmarks, which link with full LTO,
  additionally need the `LLDEnzyme` plugin and `ld.lld`.
- **LLVM**: the LLVM that Enzyme was built against, with flang, `FileCheck`
  and `lit`.
- **Downloaded or bundled by the build**:
  - [SUNDIALS](https://github.com/LLNL/sundials), the pinned
    `dependencies/sundials` submodule, built by CMake (or an existing install
    via `SUNDIALS_PREFIX`).
  - [OpenAD/F](https://www.mcs.anl.gov/OpenAD/) examples, from the
    `OpenAD_2014-03-15.tgz` snapshot downloaded at configure time
    (`-DBENCH_OPENAD=OFF` skips it).
- **System libraries**: [netCDF-Fortran](https://github.com/Unidata/netcdf-fortran)
  (with netCDF-C) for the OpenAD shallow water benchmark, found by CMake; the
  benchmark is unsupported without it.

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

### SUNDIALS benchmarks

The benchmarks under `benchmarks/sundials/` link against SUNDIALS' CVODE and
its Fortran module interface. The pinned `dependencies/sundials` submodule is
built automatically as part of the `bench-enzyme-fortran` target (CVODE only,
Fortran module interface on) and installed into the build tree; the lit
configuration exposes it through the `%sundialsFortran` (Fortran `.mod`
include path) and `%sundialsLibs` (link flags, with rpath) substitutions. To
use an existing SUNDIALS install instead, configure with
`-DSUNDIALS_PREFIX=/path/to/sundials/install` (built with
`-DBUILD_FORTRAN_MODULE_INTERFACE=ON`).

### Timing benchmarks

The Fortran program itself times repeated calls with `system_clock` (repetition
count passed as the first command-line argument, default 9) and prints
average/minimum wall times and the adjoint/primal ratio, e.g.:

```
repetitions: 9
primal average wall time (ms):     0.000038
adjoint average wall time (ms):     0.000027
primal minimum wall time (ms):     0.000026
adjoint minimum wall time (ms):     0.000026
average adjoint / primal time ratio:   0.7185
```

The lit `RUN:` line compiles the benchmark with Enzyme loaded, runs it, and
`FileCheck`s the printed timings.

The helper script `bench.sh` runs the compiled benchmarks directly from their
lit `Output` directories. It discovers every benchmark automatically and
accepts name filters (`./bench.sh bratu`) and a repetition count
(`-r/--repetitions N`, default 9). See `./bench.sh --help`.

### OpenAD benchmarks

The benchmarks under `benchmarks/openad/` differentiate examples that ship
with the OpenAD/F tool, unmodified: a three-box ocean model (BoxModel), the
MINPACK-2 flow in a driven cavity (FDC), and wind-driven channel flow over a
sill (ShallowWater, the SWM case of Vlasenko, Köhl and Stammer, Comput. Phys.
Commun. 199, 2016). CMake downloads the last OpenAD snapshot
(`OpenAD_2014-03-15.tgz`, checked by SHA-256) into the build tree and lit
exposes its examples as `%openad`; configure with `-DBENCH_OPENAD=OFF` to skip
the download. Each benchmark checks its gradient against finite differences,
and the shallow water one also against the gradient OpenAD computed, before
timing it. They build with full LTO through FlangEnzyme's config file
(`REQUIRES: flangenzyme_lld`), and the shallow water benchmark needs
netCDF-Fortran.

## Benchmarks

- [x] Square
- [x] OpenAD (box model, driven cavity, shallow water)
- [x] SUNDIALS
- [ ] Fiats
- [ ] MITgcm
