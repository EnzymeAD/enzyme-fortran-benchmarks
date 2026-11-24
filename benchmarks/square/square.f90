! REQUIRES: fortran
! RUN: %fc -flto -O2 -c %loadFortran %s -o /dev/stdout | %opt %loadEnzyme %enzyme -o %t.ll && %fc -flto -O2 %t.ll -o %t1 && %t1 9 | FileCheck %s
! RUN: %if flangenzyme %{ %fc -O2 %loadFortran %loadFlangEnzyme %s -o %t2 && %t2 9 | FileCheck %s %}

! A minimal timing benchmark, timing repeated calls to the primal
! square(x) = x**2 and its reference and Enzyme-based derivatives. Average and
! minimum wall times are printed for each case, as well as various ratios.

program main
  use, intrinsic :: iso_fortran_env, only: int64
  use enzyme, only: enzyme_dup, enzyme_autodiff
  implicit none
  integer :: args, repetition, repetitions, ios
  integer(int64) :: start_count, end_count, clock_rate
  real :: x, ref_dx, ad_dx, y, primal, primal_acc, ref_acc, ad_acc
  real, allocatable :: primal_times(:), ref_times(:), ad_times(:)
  character(len=32) :: repetitions_arg

  repetitions = 9
  args = command_argument_count()
  if (args > 1) error stop "usage: timing [repetitions]"
  if (args == 1) then
    call get_command_argument(1, repetitions_arg)
    read (repetitions_arg, *, iostat=ios) repetitions
    if (ios /= 0 .or. repetitions < 1) error stop "repetitions must be a positive integer"
  end if
  allocate(primal_times(repetitions))
  allocate(ref_times(repetitions))
  allocate(ad_times(repetitions))

  x = 3.0
  ref_dx = 0.0
  ad_dx = 0.0
  primal_acc = 0.0
  ref_acc = 0.0
  ad_acc = 0.0
  call system_clock(count_rate=clock_rate)

  ! Time the primal calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    primal = square(x)
    call system_clock(end_count)
    primal_times(repetition) = real(end_count - start_count, kind(1.0)) / real(clock_rate, kind(1.0))
    primal_acc = primal_acc + primal
  end do

  ! Time the reference calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    ref_dx = dsquare(x)
    call system_clock(end_count)
    ref_times(repetition) = real(end_count - start_count, kind(1.0)) / real(clock_rate, kind(1.0))
    ref_acc = ref_acc + ref_dx
    ref_dx = 0.0
  end do

  ! Time the Enzyme-differentiated calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    call enzyme_autodiff(square, enzyme_dup, x, ad_dx)
    call system_clock(end_count)
    ad_times(repetition) = real(end_count - start_count, kind(1.0)) / real(clock_rate, kind(1.0))
    ad_acc = ad_acc + ad_dx
    ad_dx = 0.0
  end do

  ! Guard against the optimizer removing the calls
  if (primal_acc /= 9.0 * 9.0) error stop "unexpected primal result"
  if (ref_acc /= 9.0 * 6.0) error stop "unexpected reference result"
  if (ad_acc /= 9.0 * 6.0) error stop "unexpected adjoint result"

  write (*, "(a,i0)") "repetitions: ", repetitions
  write (*, "(a40,es11.4)") "primal average wall time (s): ", sum(primal_times) / repetitions
  write (*, "(a40,es11.4)") "reference average wall time (s): ", sum(ref_times) / repetitions
  write (*, "(a40,es11.4)") "adjoint average wall time (s): ", sum(ad_times) / repetitions
  write (*, "(a40,es11.4)") "primal minimum wall time (s): ", minval(primal_times)
  write (*, "(a40,es11.4)") "reference minimum wall time (s): ", minval(ref_times)
  write (*, "(a40,es11.4)") "adjoint minimum wall time (s): ", minval(ad_times)
  write (*, "(a40,es11.4)") "average adjoint / reference time ratio: ", sum(ad_times) / sum(ref_times)
  write (*, "(a40,es11.4)") "average adjoint / primal time ratio: ", sum(ad_times) / sum(primal_times)

contains

  ! Function to differentiate
  real function square(x)
    real, intent(in) :: x
    square = x**2
  end function square

  ! Reference implementation for derivative
  real function dsquare(x)
    real, intent(in) :: x
    dsquare = 2 * x
  end function dsquare

end program main

! CHECK: repetitions: 9
! CHECK: primal average wall time (s):
! CHECK: adjoint average wall time (s):
! CHECK: primal minimum wall time (s):
! CHECK: adjoint minimum wall time (s):
! CHECK: average adjoint / primal time ratio:
