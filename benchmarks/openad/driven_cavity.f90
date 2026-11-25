! REQUIRES: fortran, flangenzyme_lld, openad
! RUN: mkdir -p %t.mod && %fc %flangEnzymeConfig -flto -O2 %loadFortran -module-dir %t.mod %openad/FDC/head.f %openad/FDC/init.f %S/Inputs/openad_harness.f90 %s -o %t2 && %t2 9 | FileCheck %s

! Timing benchmark for OpenAD/F's flow in a driven cavity example
! (Examples/FDC), the MINPACK-2 dfdcfj problem, on a 10 x 10 grid with
! Reynolds number 10 at the problem's standard starting point. The gradient of
! the residual, weighted by a fixed vector, is taken with respect to the
! stream function.

module driven_cavity_wrapper
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: driven_cavity, driven_cavity_gradient, nx, ny

  integer, parameter :: nx = 10, ny = 10
  real(real64), parameter :: reynolds = 10

contains

  subroutine driven_cavity(n, x, m, y)
    integer, intent(in) :: n, m
    real(real64), intent(in) :: x(n)
    real(real64), intent(out) :: y(m)
    real(real64) :: r
    external :: head

    r = reynolds
    call head(nx, ny, x, y, r)
  end subroutine driven_cavity

  subroutine driven_cavity_gradient(n, x, dx, m, w)
    use enzyme, only: enzyme_autodiff, enzyme_const, enzyme_dup
    integer, intent(in) :: n, m
    real(real64), intent(in) :: x(n), w(m)
    real(real64), intent(out) :: dx(n)
    real(real64) :: y(m), dy(m)

    dx = 0
    dy = w
    call enzyme_autodiff(driven_cavity, enzyme_const, n, enzyme_dup, x, dx, &
                         enzyme_const, m, enzyme_dup, y, dy)
  end subroutine driven_cavity_gradient

end module driven_cavity_wrapper

program main
  use, intrinsic :: iso_fortran_env, only: real64
  use driven_cavity_wrapper, only: driven_cavity, driven_cavity_gradient, nx, ny
  use openad_harness, only: read_repetitions, fill_deterministic, check_and_time
  implicit none
  integer, parameter :: n = nx * ny, m = nx * ny
  real(real64) :: x(n), dx(n), w(m), v(n)
  external :: init

  call init(nx, ny, x)
  call fill_deterministic(w, 1)
  call fill_deterministic(v, 2)
  call check_and_time(driven_cavity, driven_cavity_gradient, n, x, m, w, v, &
                      1.0e-6_real64, 1.0e-6_real64, read_repetitions(), dx)
end program main

! CHECK: gradient check: PASS
! CHECK: repetitions: 9
! CHECK: primal average wall time (s):
! CHECK: adjoint average wall time (s):
! CHECK: primal minimum wall time (s):
! CHECK: adjoint minimum wall time (s):
! CHECK: average adjoint / primal time ratio:
