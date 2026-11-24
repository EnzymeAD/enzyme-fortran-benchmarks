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
! with Newton's method, and each Newton system (I - dt f'(u)) d = -G(u) with
! conjugate gradients. The only derivative information needed is the product
! f'(u) v, obtained either from Enzyme forward mode on the residual routine
! (adjoint) or from a reference hand-derived JVP; the solver itself is ordinary
! Fortran and is not differentiated. This is the same division of labour as
! handing a stiff ODE solver such as DVODE a user Jacobian: Enzyme
! differentiates the model, not the integrator. The primal run times the
! residual evaluation itself. For lambda = 1 the solution relaxes to the
! lower steady state of the Bratu equation, whose closed form is the
! reference solution.

module bratu_implicit_euler
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: dp, residual, jacobian_vector_product_ad, jacobian_vector_product_ref, &
    integrate_krylov
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

  subroutine jacobian_vector_product_ad(n, dx, lambda, u, v, jv)
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
  end subroutine jacobian_vector_product_ad

  ! Reference implementation for the Jacobian-vector product
  subroutine jacobian_vector_product_ref(n, dx, lambda, u, v, jv)
    integer, intent(in) :: n
    real(dp), intent(in) :: dx, lambda
    real(dp), intent(in) :: u(n), v(n)
    real(dp), intent(out) :: jv(n)
    real(dp) :: vl, vr
    integer :: i
    do i = 1, n
      vl = 0.0_dp
      vr = 0.0_dp
      if (i > 1) vl = v(i - 1)
      if (i < n) vr = v(i + 1)
      jv(i) = (vl - 2.0_dp * v(i) + vr) / dx**2 + lambda * exp(u(i)) * v(i)
    end do
  end subroutine jacobian_vector_product_ref

  subroutine integrate_krylov(n, dx, lambda, dt, nsteps, use_ad, u, newton_steps, krylov_steps)
    integer, intent(in) :: n, nsteps
    real(dp), intent(in) :: dx, lambda, dt
    logical, intent(in) :: use_ad
    real(dp), intent(out) :: u(n)
    integer, intent(out) :: newton_steps, krylov_steps
    real(dp) :: u_old(n)
    integer :: step

    u = 0.0_dp
    newton_steps = 0
    krylov_steps = 0
    do step = 1, nsteps
      u_old = u
      call implicit_euler_step(n, dx, lambda, dt, u_old, u, use_ad, newton_steps, krylov_steps)
    end do
  end subroutine integrate_krylov

  subroutine implicit_euler_step(n, dx, lambda, dt, u_old, u, use_ad, newton_steps, &
                                 krylov_steps)
    integer, intent(in) :: n
    real(dp), intent(in) :: dx, lambda, dt
    real(dp), intent(in) :: u_old(n)
    real(dp), intent(inout) :: u(n)
    logical, intent(in) :: use_ad
    integer, intent(inout) :: newton_steps, krylov_steps
    real(dp), parameter :: eta = 1.0e-3_dp
    real(dp) :: f(n), g(n), d(n), r(n), p(n), jp(n), ap(n)
    real(dp) :: norm0, norm, tol, alpha, beta, rr, rr_new
    integer :: newton, it

    u = u_old
    call residual(n, dx, lambda, u, f)
    g = u - u_old - dt * f
    norm0 = sqrt(dot_product(g, g))
    norm = norm0
    tol = 1.0e-10_dp * norm0 + 1.0e-14_dp
    newton = 0
    do while (norm > tol .and. newton < 20)
      d = 0.0_dp
      r = -g
      p = r
      rr = dot_product(r, r)
      do it = 1, 10 * n
        if (use_ad) then
          call jacobian_vector_product_ad(n, dx, lambda, u, p, jp)
        else
          call jacobian_vector_product_ref(n, dx, lambda, u, p, jp)
        end if
        krylov_steps = krylov_steps + 1
        ap = p - dt * jp
        alpha = rr / dot_product(p, ap)
        d = d + alpha * p
        r = r - alpha * ap
        rr_new = dot_product(r, r)
        if (sqrt(rr_new) <= eta * norm) exit
        beta = rr_new / rr
        p = r + beta * p
        rr = rr_new
      end do
      u = u + d
      call residual(n, dx, lambda, u, f)
      g = u - u_old - dt * f
      norm = sqrt(dot_product(g, g))
      newton = newton + 1
    end do
    newton_steps = newton_steps + newton
  end subroutine implicit_euler_step
end module bratu_implicit_euler

program bratu_implicit_euler_jvp_timing
  use, intrinsic :: iso_fortran_env, only: int64
  use bratu_implicit_euler, only: dp, integrate_krylov, residual, &
    jacobian_vector_product_ad, jacobian_vector_product_ref
  implicit none
  integer, parameter :: n = 200
  integer, parameter :: nsteps = 50
  real(dp), parameter :: lambda = 1.0_dp
  real(dp), parameter :: dt = 0.05_dp
  real(dp), parameter :: theta = 1.5171646_dp
  integer :: args, i, repetition, repetitions, ios
  integer :: newton_steps(2), krylov_steps(2)
  integer(int64) :: start_count, end_count, clock_rate
  real(dp) :: dx, x(n), f(n), solution_ref(n), u_sample(n), v(n)
  real(dp) :: jv_ad(n), jv_ref(n)
  real(dp) :: u_ad(n), u_ref(n)
  real(dp) :: primal_acc
  real(dp), allocatable :: primal_times(:), ref_times(:), ad_times(:)
  character(len=32) :: repetitions_arg

  repetitions = 9
  args = command_argument_count()
  if (args > 1) error stop "usage: bratu_implicit_euler_jvp [repetitions]"
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
  v = sin(2.0_dp * 3.141592653589793_dp * x)
  call jacobian_vector_product_ad(n, dx, lambda, u_sample, v, jv_ad)
  call jacobian_vector_product_ref(n, dx, lambda, u_sample, v, jv_ref)
  if (maxval(abs(jv_ad - jv_ref)) * dx**2 > 1.0e-12_dp) then
    error stop "reference JVP is incorrect"
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
    call integrate_krylov(n, dx, lambda, dt, nsteps, .false., u_ref, newton_steps(1), &
                          krylov_steps(1))
    call system_clock(end_count)
    ref_times(repetition) = real(end_count - start_count, dp) / real(clock_rate, dp)
  end do

  ! Time the Enzyme-differentiated (adjoint) calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    call integrate_krylov(n, dx, lambda, dt, nsteps, .true., u_ad, newton_steps(2), &
                          krylov_steps(2))
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
  write (*, "(a,i0,1x,i0)") "Krylov iterations (adjoint, reference): ", krylov_steps
  write (*, "(a,i0)") "repetitions: ", repetitions
  write (*, "(a40,es11.4)") "primal average wall time (s): ", sum(primal_times) / repetitions
  write (*, "(a40,es11.4)") "reference average wall time (s): ", sum(ref_times) / repetitions
  write (*, "(a40,es11.4)") "adjoint average wall time (s): ", sum(ad_times) / repetitions
  write (*, "(a40,es11.4)") "primal minimum wall time (s): ", minval(primal_times)
  write (*, "(a40,es11.4)") "reference minimum wall time (s): ", minval(ref_times)
  write (*, "(a40,es11.4)") "adjoint minimum wall time (s): ", minval(ad_times)
  write (*, "(a40,es11.4)") "average adjoint / reference time ratio: ", sum(ad_times) / sum(ref_times)
  write (*, "(a40,es11.4)") "average adjoint / primal time ratio: ", sum(ad_times) / sum(primal_times)
end program bratu_implicit_euler_jvp_timing

! CHECK: both setups reached t =
! CHECK: solutions agree: T
! CHECK: Newton iterations (adjoint, reference):
! CHECK: Krylov iterations (adjoint, reference):
! CHECK: repetitions: 9
! CHECK: primal average wall time (s):
! CHECK: reference average wall time (s):
! CHECK: adjoint average wall time (s):
! CHECK: primal minimum wall time (s):
! CHECK: reference minimum wall time (s):
! CHECK: adjoint minimum wall time (s):
! CHECK: average adjoint / reference time ratio:
! CHECK: average adjoint / primal time ratio:
