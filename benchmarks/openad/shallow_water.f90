! REQUIRES: fortran, flangenzyme_lld, openad, netcdf
! RUN: rm -rf %t.run && mkdir -p %t.mod %t.run && cp %openad/ShallowWater/Data/* %t.run/
! RUN: sed -i "/^# fullio suppressio/{n;s/.*/.false. .true./}" %t.run/data
! RUN: sed -n 's/^grad(\(.*\))=/\1/p' %openad/ShallowWater/shallow_water_ad.out > %t.run/grad_reference.txt
! RUN: %fc %flangEnzymeConfig -flto -O2 -fdefault-real-8 -fdefault-double-8 %loadFortran -module-dir %t.mod \
! RUN:   -I%openad/ShallowWater/Inc %netcdfInclude -Wl,-mllvm,-enzyme-global-activity \
! RUN:   %shallow_water_sources %S/Inputs/openad_harness.f90 %s -o %t2 %netcdfLibs
! RUN: cd %t.run && %t2 9 | FileCheck %s

! Timing benchmark for OpenAD/F's shallow water example
! (Examples/ShallowWater): wind-driven channel flow over a meridional sill on a
! 20 x 20 grid (Losch and Wunsch, J. Atmos. Oceanic Technol. 20, 2003,
! https://doi.org/10.1175/1520-0426(2003)020<1685:BTAACV>2.0.CO;2), the SWM
! case of Vlasenko, Koehl and Stammer, Comput. Phys. Commun. 199 (2016),
! https://doi.org/10.1016/j.cpc.2015.10.008. The
! gradient of the zonal volume transport after 40000 time steps is taken with
! respect to the 160 bathymetry controls, and compared both with finite
! differences and with the gradient OpenAD computed (shallow_water_ad.out).
!
! Compiled with -fdefault-real-8 as the example's own build does, so the
! default real is real64 here too.

module shallow_water_wrapper
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: shallow_water, shallow_water_gradient

contains

  subroutine shallow_water(n, x, m, y)
    integer, intent(in) :: n, m
    real(real64), intent(in) :: x(n)
    real(real64), intent(out) :: y(m)
    real(real64) :: cost
    external :: forward_model

    ! The zonal volume transport is the only output, so m is 1.
    call forward_model(n, x, cost)
    y(1) = cost
  end subroutine shallow_water

  subroutine shallow_water_gradient(n, x, dx, m, w)
    use enzyme, only: enzyme_autodiff, enzyme_const, enzyme_dup
    integer, intent(in) :: n, m
    real(real64), intent(in) :: x(n), w(m)
    real(real64), intent(out) :: dx(n)
    real(real64) :: y(m), dy(m)

    dx = 0
    dy = w
    call enzyme_autodiff(shallow_water, enzyme_const, n, enzyme_dup, x, dx, &
                         enzyme_const, m, enzyme_dup, y, dy)
  end subroutine shallow_water_gradient

end module shallow_water_wrapper

program main
  use, intrinsic :: iso_fortran_env, only: real64
  use adsize, only: ncmax
  use shallow_water_wrapper, only: shallow_water, shallow_water_gradient
  use openad_harness, only: read_repetitions, fill_deterministic, check_and_time
  implicit none
  integer, parameter :: m = 1
  ! Finite difference step along v, and the largest relative difference
  ! between the adjoint and the finite difference directional derivatives.
  real(real64), parameter :: h = 1.0e-6_real64, tol = 1.0e-5_real64
  real(real64), parameter :: reference_tolerance = 1.0e-8_real64
  real(real64) :: xc(ncmax), w(m), t0, max_err
  real(real64), allocatable :: v(:), dxc(:), reference(:)
  integer :: nc, i, index, unit
  external :: inifields, readparms, make_weights, read_eta_data
  external :: length_of_control_vector, map_to_control_vector

  call inifields
  call readparms
  call make_weights
  t0 = 0
  call read_eta_data(t0)
  call length_of_control_vector(nc)
  call map_to_control_vector(nc, xc)

  allocate (v(nc), dxc(nc), reference(nc))
  w = 1
  call fill_deterministic(v, 2)
  call check_and_time(shallow_water, shallow_water_gradient, nc, xc(1:nc), m, &
                      w, v, h, tol, read_repetitions(), dxc)

  open (newunit=unit, file="grad_reference.txt", status="old", action="read")
  do i = 1, nc
    read (unit, *) index, reference(i)
  end do
  close (unit)
  max_err = maxval(abs(dxc - reference)) / maxval(abs(reference))
  write (*, "(a40,es11.4)") "largest difference to OpenAD gradient: ", max_err
  if (.not. (max_err <= reference_tolerance)) error stop "OpenAD gradient mismatch"
  write (*, "(a)") "OpenAD gradient check: PASS"
end program main

! CHECK: gradient check: PASS
! CHECK: repetitions: 9
! CHECK: primal average wall time (s):
! CHECK: adjoint average wall time (s):
! CHECK: primal minimum wall time (s):
! CHECK: adjoint minimum wall time (s):
! CHECK: average adjoint / primal time ratio:
! CHECK: OpenAD gradient check: PASS
