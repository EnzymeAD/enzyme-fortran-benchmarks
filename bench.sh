#!/bin/bash
# Run all benchmarks

cd build/benchmarks/benchmarks
for BENCHMARK in square; do
  cd ${BENCHMARK}/Output
  for CASE in 1 2; do
    echo "${BENCHMARK}, case ${CASE}"
    echo
    "./${BENCHMARK}.f90.tmp${CASE}"
    echo
  done
  cd ../..
done
