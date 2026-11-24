! REQUIRES: fortran
! RUN: %fc -flto -O2 -c %loadFortran %s -o /dev/stdout | %opt %loadEnzyme %enzyme -o %t.ll && %fc -flto -O2 %t.ll -o %t1 && %t1 9 | FileCheck %s
! RUN: %if flangenzyme %{ %fc -O2 %loadFortran %loadFlangEnzyme %s -o %t2 && %t2 9 | FileCheck %s %}

! Timing benchmark for the time-dependent 1D Bratu problem
!
!   u_t = u_xx + lambda * exp(u) on (0, 1),  u(0, t) = u(1, t) = 0,  u(x, 0) = 0,
!
! discretized in space with central differences (method of lines) and
! integrated with backward Euler. Each implicit step solves
!
!   G(u) = u - u_old - dt * f(u) = 0
!
! with Newton's method using an assembled tridiagonal Jacobian, built either
! column by column with Enzyme forward mode on the residual routine (adjoint)
! or by a reference hand-derived formula, followed by a direct tridiagonal
! solve. This mirrors what a stiff ODE solver such as DVODE with MF = 21
! asks of a user jac routine: Enzyme differentiates the model, not the
! integrator. The primal run times the residual evaluation itself. For
! lambda = 1 the solution relaxes to the lower steady state of the Bratu
! equation, whose closed form is the reference solution.

module bratu_implicit_euler
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: dp, residual, jacobian_adjoint, jacobian_ref, integrate_tridiagonal
  integer, parameter :: dp = real64
contains

  subroutine residual(n, dx, lambda, u, f)
    integer, intent(in) :: n
    real(dp), intent(in) :: dx, lambda
    real(dp), intent(in) :: u(n)
    real(dp), intent(out) :: f(n)
    real(dp) :: ul, ur
    integer :: i
    do i = 1, n
      ul = 0.0_dp
      ur = 0.0_dp
      if (i > 1) ul = u(i - 1)
      if (i < n) ur = u(i + 1)
      f(i) = (ul - 2.0_dp * u(i) + ur) / dx**2 + lambda * exp(u(i))
    end do
  end subroutine residual

  subroutine jacobian_vector_product_adjoint(n, dx, lambda, u, v, jv)
    use enzyme, only: enzyme_const, enzyme_dup, enzyme_fwddiff
    integer, intent(in) :: n
    real(dp), intent(in) :: dx, lambda
    real(dp), intent(in) :: u(n), v(n)
    real(dp), intent(out) :: jv(n)
    real(dp) :: f(n)
    f = 0.0_dp
    jv = 0.0_dp
    call enzyme_fwddiff(residual, enzyme_const, n, enzyme_const, dx, &
                        enzyme_const, lambda, enzyme_dup, u, v, &
                        enzyme_dup, f, jv)
  end subroutine jacobian_vector_product_adjoint

  subroutine jacobian_adjoint(n, dx, lambda, u, jac)
    integer, intent(in) :: n
    real(dp), intent(in) :: dx, lambda
    real(dp), intent(in) :: u(n)
    real(dp), intent(out) :: jac(n, n)
    real(dp) :: v(n), column(n)
    integer :: j

    do j = 1, n
      v = 0.0_dp
      v(j) = 1.0_dp
      call jacobian_vector_product_adjoint(n, dx, lambda, u, v, column)
      jac(:, j) = column
    end do
  end subroutine jacobian_adjoint

  ! Reference implementation for the Jacobian
  subroutine jacobian_ref(n, dx, lambda, u, jac)
    integer, intent(in) :: n
    real(dp), intent(in) :: dx, lambda
    real(dp), intent(in) :: u(n)
    real(dp), intent(out) :: jac(n, n)
    integer :: i

    jac = 0.0_dp
    do i = 1, n
      jac(i, i) = -2.0_dp / dx**2 + lambda * exp(u(i))
      if (i > 1) jac(i, i - 1) = 1.0_dp / dx**2
      if (i < n) jac(i, i + 1) = 1.0_dp / dx**2
    end do
  end subroutine jacobian_ref

  subroutine integrate_tridiagonal(n, dx, lambda, dt, nsteps, use_adjoint, u, newton_steps)
    integer, intent(in) :: n, nsteps
    real(dp), intent(in) :: dx, lambda, dt
    logical, intent(in) :: use_adjoint
    real(dp), intent(out) :: u(n)
    integer, intent(out) :: newton_steps
    real(dp) :: u_old(n), f(n), g(n), jac(n, n)
    real(dp) :: diagonal(n), lower(n - 1), upper(n - 1), d(n)
    real(dp) :: norm0, norm, tol
    integer :: step, newton, i

    u = 0.0_dp
    newton_steps = 0
    do step = 1, nsteps
      u_old = u
      call residual(n, dx, lambda, u, f)
      g = u - u_old - dt * f
      norm0 = sqrt(dot_product(g, g))
      norm = norm0
      tol = 1.0e-10_dp * norm0 + 1.0e-14_dp
      newton = 0
      do while (norm > tol .and. newton < 20)
        if (use_adjoint) then
          call jacobian_adjoint(n, dx, lambda, u, jac)
        else
          call jacobian_ref(n, dx, lambda, u, jac)
        end if
        diagonal = 1.0_dp - dt * [(jac(i, i), i = 1, n)]
        lower = -dt * [(jac(i + 1, i), i = 1, n - 1)]
        upper = -dt * [(jac(i, i + 1), i = 1, n - 1)]
        call solve_tridiagonal(diagonal, lower, upper, -1.0_dp * g, d)
        u = u + d
        call residual(n, dx, lambda, u, f)
        g = u - u_old - dt * f
        norm = sqrt(dot_product(g, g))
        newton = newton + 1
        if (norm <= 1.0e-12_dp) exit
      end do
      if (norm > tol .and. norm > 1.0e-12_dp) then
        error stop "implicit Euler Newton step did not converge"
      end if
      newton_steps = newton_steps + newton
    end do
  end subroutine integrate_tridiagonal

  subroutine solve_tridiagonal(diagonal, lower, upper, rhs, solution)
    real(dp), intent(in) :: diagonal(:), lower(:), upper(:), rhs(:)
    real(dp), intent(out) :: solution(:)
    real(dp) :: modified_upper(size(upper)), denominator
    integer :: i, n
    n = size(diagonal)
    modified_upper(1) = upper(1) / diagonal(1)
    solution(1) = rhs(1) / diagonal(1)
    do i = 2, n
      denominator = diagonal(i) - lower(i - 1) * modified_upper(i - 1)
      if (i < n) modified_upper(i) = upper(i) / denominator
      solution(i) = (rhs(i) - lower(i - 1) * solution(i - 1)) / denominator
    end do
    do i = n - 1, 1, -1
      solution(i) = solution(i) - modified_upper(i) * solution(i + 1)
    end do
  end subroutine solve_tridiagonal
end module bratu_implicit_euler

program bratu_implicit_euler_jacobian_timing
  use, intrinsic :: iso_fortran_env, only: int64
  use bratu_implicit_euler, only: dp, integrate_tridiagonal, residual, &
    jacobian_adjoint, jacobian_ref
  implicit none
  integer, parameter :: n = 200
  integer, parameter :: nsteps = 50
  real(dp), parameter :: lambda = 1.0_dp
  real(dp), parameter :: dt = 0.05_dp
  real(dp), parameter :: theta = 1.5171646_dp
  integer :: args, i, repetition, repetitions, ios
  integer :: newton_steps(2)
  integer(int64) :: start_count, end_count, clock_rate
  real(dp) :: dx, x(n), f(n), solution_ref(n), u_sample(n)
  real(dp) :: jac_ad(n, n), jac_ref(n, n)
  real(dp) :: u_ad(n), u_ref(n)
  real(dp) :: primal_acc
  real(dp), allocatable :: primal_times(:), ref_times(:), ad_times(:)
  character(len=32) :: repetitions_arg

  repetitions = 9
  args = command_argument_count()
  if (args > 1) error stop "usage: bratu_implicit_euler_jacobian [repetitions]"
  if (args == 1) then
    call get_command_argument(1, repetitions_arg)
    read (repetitions_arg, *, iostat=ios) repetitions
    if (ios /= 0 .or. repetitions < 1) error stop "repetitions must be a positive integer"
  end if
  allocate (primal_times(repetitions), ref_times(repetitions), ad_times(repetitions))

  dx = 1.0_dp / (n + 1)
  x = [(real(i, dp) * dx, i = 1, n)]
  solution_ref = -2.0_dp * log(cosh(theta * (x - 0.5_dp) / 2.0_dp) / cosh(theta / 4.0_dp))
  u_sample = 0.5_dp * sin(3.141592653589793_dp * x)
  call jacobian_adjoint(n, dx, lambda, u_sample, jac_ad)
  call jacobian_ref(n, dx, lambda, u_sample, jac_ref)
  if (maxval(abs(jac_ad - jac_ref)) * dx**2 > 1.0e-12_dp) then
    error stop "reference Jacobian is incorrect"
  end if
  primal_acc = 0.0_dp
  call system_clock(count_rate=clock_rate)

  ! Time the primal calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    call residual(n, dx, lambda, x, f)
    call system_clock(end_count)
    primal_times(repetition) = real(end_count - start_count, dp) / real(clock_rate, dp)
    primal_acc = primal_acc + sum(f)
  end do

  ! Time the reference calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    call integrate_tridiagonal(n, dx, lambda, dt, nsteps, .false., u_ref, newton_steps(1))
    call system_clock(end_count)
    ref_times(repetition) = real(end_count - start_count, dp) / real(clock_rate, dp)
  end do

  ! Time the Enzyme-differentiated (adjoint) calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    call integrate_tridiagonal(n, dx, lambda, dt, nsteps, .true., u_ad, newton_steps(2))
    call system_clock(end_count)
    ad_times(repetition) = real(end_count - start_count, dp) / real(clock_rate, dp)
  end do

  ! Guard against the optimizer removing the calls
  if (primal_acc == 0.0_dp) error stop "unexpected primal result"
  if (maxval(abs(u_ad - solution_ref)) >= 1.0e-5_dp) error stop "unexpected adjoint result"
  if (maxval(abs(u_ref - solution_ref)) >= 1.0e-5_dp) error stop "unexpected reference result"
  if (maxval(abs(u_ad - u_ref)) >= 1.0e-10_dp) error stop "solutions differ"

  write (*, "(a,f4.1)") "both setups reached t = ", nsteps * dt
  write (*, "(a,l1)") "solutions agree: ", maxval(abs(u_ad - u_ref)) < 1.0e-10_dp
  write (*, "(a,i0,1x,i0)") "Newton iterations (adjoint, reference): ", newton_steps
  write (*, "(a,i0)") "repetitions: ", repetitions
  write (*, "(a,f12.3)") "primal average wall time (ms): ", 1.0e3_dp * sum(primal_times) / repetitions
  write (*, "(a,f12.3)") "reference average wall time (ms): ", 1.0e3_dp * sum(ref_times) / repetitions
  write (*, "(a,f12.3)") "adjoint average wall time (ms): ", 1.0e3_dp * sum(ad_times) / repetitions
  write (*, "(a,f12.3)") "primal minimum wall time (ms): ", 1.0e3_dp * minval(primal_times)
  write (*, "(a,f12.3)") "reference minimum wall time (ms): ", 1.0e3_dp * minval(ref_times)
  write (*, "(a,f12.3)") "adjoint minimum wall time (ms): ", 1.0e3_dp * minval(ad_times)
  write (*, "(a,f8.4)") "average adjoint / reference time ratio: ", sum(ad_times) / sum(ref_times)
  write (*, "(a,f10.4)") "average adjoint / primal time ratio: ", sum(ad_times) / sum(primal_times)
end program bratu_implicit_euler_jacobian_timing

! CHECK: both setups reached t =
! CHECK: solutions agree: T
! CHECK: Newton iterations (adjoint, reference):
! CHECK: repetitions: 9
! CHECK: primal average wall time (ms):
! CHECK: reference average wall time (ms):
! CHECK: adjoint average wall time (ms):
! CHECK: primal minimum wall time (ms):
! CHECK: reference minimum wall time (ms):
! CHECK: adjoint minimum wall time (ms):
! CHECK: average adjoint / reference time ratio:
! CHECK: average adjoint / primal time ratio:
