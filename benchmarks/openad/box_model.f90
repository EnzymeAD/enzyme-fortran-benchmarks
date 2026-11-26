! REQUIRES: fortran, flangenzyme_lld, openad
! RUN: mkdir -p %t.mod && %fc %flangEnzymeConfig -flto -O2 %loadFortran -module-dir %t.mod -Wl,-mllvm,-enzyme-global-activity %openad/BoxModel/all_globals_mod.f %openad/BoxModel/head.f %S/Inputs/openad_harness.f90 %s -o %t2 && %t2 9 | FileCheck %s

! Timing benchmark for OpenAD/F's box model (Examples/BoxModel): a three-box
! model of the thermohaline circulation, integrated for 50 years in 3650 time
! steps. The gradient of the final temperatures and salinities, weighted by a
! fixed vector, is taken with respect to the perturbations xx of the initial
! state. The model keeps its whole state in module all_globals_mod.

module box_model_wrapper
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: box_model, box_model_gradient

contains

  subroutine box_model(n, x, m, y)
    use all_globals_mod, only: xx, tnew, snew, ndim
    integer, intent(in) :: n, m
    real(real64), intent(in) :: x(n)
    real(real64), intent(out) :: y(m)
    external :: box_model_body

    xx = x
    call box_model_body
    y(1:ndim) = tnew
    y(ndim + 1:2 * ndim) = snew
  end subroutine box_model

  subroutine box_model_gradient(n, x, dx, m, w)
    use enzyme, only: enzyme_autodiff, enzyme_const, enzyme_dup
    integer, intent(in) :: n, m
    real(real64), intent(in) :: x(n), w(m)
    real(real64), intent(out) :: dx(n)
    real(real64) :: y(m), dy(m)

    dx = 0
    dy = w
    call enzyme_autodiff(box_model, enzyme_const, n, enzyme_dup, x, dx, &
                         enzyme_const, m, enzyme_dup, y, dy)
  end subroutine box_model_gradient

end module box_model_wrapper

program main
  use, intrinsic :: iso_fortran_env, only: real64
  use box_model_wrapper, only: box_model, box_model_gradient
  use openad_harness, only: read_repetitions, fill_deterministic, check_and_time
  implicit none
  integer, parameter :: n = 6, m = 6
  ! Finite difference step along v, and the largest relative difference
  ! between the adjoint and the finite difference directional derivatives.
  real(real64), parameter :: h = 1.0e-6_real64, tol = 1.0e-5_real64
  real(real64) :: x(n), dx(n), w(m), v(n)
  external :: box_ini_params

  call box_ini_params
  x = 0
  call fill_deterministic(w, 1)
  call fill_deterministic(v, 2)
  call check_and_time(box_model, box_model_gradient, n, x, m, w, v, &
                      h, tol, read_repetitions(), dx)
end program main

! CHECK: gradient check: PASS
! CHECK: repetitions: 9
! CHECK: primal average wall time (s):
! CHECK: adjoint average wall time (s):
! CHECK: primal minimum wall time (s):
! CHECK: adjoint minimum wall time (s):
! CHECK: average adjoint / primal time ratio:
