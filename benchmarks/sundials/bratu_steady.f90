! REQUIRES: fortran
! RUN: %if flangenzyme %{ %fc -O2 %loadFortran %loadFlangEnzyme %s -o %t2 && %t2 9 | FileCheck %s %}

! Timing benchmark for the steady 1D Bratu problem
!
!   u'' + lambda * exp(u) = 0 on (0, 1),  u(0) = u(1) = 0,
!
! discretized with second-order central differences and solved by a
! Jacobian-free Newton-Krylov method: Newton's method on the residual F(u),
! each Newton system F'(u) d = -F(u) solved by conjugate gradients, and every
! product F'(u) v obtained either from Enzyme forward mode on the residual
! routine (adjoint) or from a reference hand-derived JVP. The timer encloses
! the complete solve; only the JVP callback differs between the two timed runs.
! The primal run times the residual evaluation itself. With lambda = 1 the
! discrete Jacobian is well conditioned and the iteration converges to the
! lower branch, whose closed form is the reference solution.

module bratu_model
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: dp, residual, jacobian_vector_product_ad, jacobian_vector_product_ref
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
end module bratu_model

module newton_krylov
  use bratu_model, only: dp, residual
  implicit none
  private
  public :: solve_bratu
contains

  subroutine solve_bratu(n, dx, lambda, u, converged, newton_steps, krylov_steps, jvp)
    integer, intent(in) :: n
    real(dp), intent(in) :: dx, lambda
    real(dp), intent(out) :: u(n)
    logical, intent(out) :: converged
    integer, intent(out) :: newton_steps, krylov_steps
    interface
      subroutine jvp(n, dx, lambda, u, v, jv)
        import :: dp
        integer, intent(in) :: n
        real(dp), intent(in) :: dx, lambda
        real(dp), intent(in) :: u(n), v(n)
        real(dp), intent(out) :: jv(n)
      end subroutine jvp
    end interface
    real(dp), parameter :: pi = 3.141592653589793_dp
    real(dp), parameter :: eta = 1.0e-3_dp
    real(dp) :: f(n), d(n), r(n), p(n), ap(n)
    real(dp) :: norm0, norm, tol, alpha, beta, rr, rr_new
    integer :: i, it

    u = sin(pi * ([(real(i, dp), i = 1, n)] / (n + 1)))
    call residual(n, dx, lambda, u, f)
    norm0 = sqrt(dot_product(f, f))
    norm = norm0
    tol = 1.0e-8_dp * norm0
    newton_steps = 0
    krylov_steps = 0

    do while (norm > tol .and. newton_steps < 20)
      d = 0.0_dp
      r = -f
      p = r
      rr = dot_product(r, r)
      do it = 1, 10 * n
        call jvp(n, dx, lambda, u, p, ap)
        krylov_steps = krylov_steps + 1
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
      norm = sqrt(dot_product(f, f))
      newton_steps = newton_steps + 1
    end do
    converged = norm <= tol
  end subroutine solve_bratu
end module newton_krylov

program bratu_steady
  use, intrinsic :: iso_fortran_env, only: int64
  use bratu_model, only: dp, residual, jacobian_vector_product_ad, jacobian_vector_product_ref
  use newton_krylov, only: solve_bratu
  implicit none
  integer, parameter :: n = 200
  real(dp), parameter :: lambda = 1.0_dp
  real(dp), parameter :: theta = 1.5171646_dp
  integer :: args, repetition, repetitions, ios
  integer(int64) :: start_count, end_count, clock_rate
  integer :: newton_steps(2), krylov_steps(2)
  real(dp) :: dx, x(n), f(n), solution_ref(n), u_ad(n), u_ref(n)
  real(dp) :: primal_acc
  real(dp), allocatable :: primal_times(:), ref_times(:), ad_times(:)
  logical :: converged(2)
  character(len=32) :: repetitions_arg

  repetitions = 9
  args = command_argument_count()
  if (args > 1) error stop "usage: bratu_steady [repetitions]"
  if (args == 1) then
    call get_command_argument(1, repetitions_arg)
    read (repetitions_arg, *, iostat=ios) repetitions
    if (ios /= 0 .or. repetitions < 1) error stop "repetitions must be a positive integer"
  end if
  allocate (primal_times(repetitions), ref_times(repetitions), ad_times(repetitions))

  dx = 1.0_dp / (n + 1)
  x = [(real(repetition, dp) * dx, repetition = 1, n)]
  solution_ref = -2.0_dp * log(cosh(theta * (x - 0.5_dp) / 2.0_dp) / cosh(theta / 4.0_dp))
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
    call solve_bratu(n, dx, lambda, u_ref, converged(1), newton_steps(1), &
                     krylov_steps(1), jacobian_vector_product_ref)
    call system_clock(end_count)
    ref_times(repetition) = real(end_count - start_count, dp) / real(clock_rate, dp)
  end do

  ! Time the Enzyme-differentiated (adjoint) calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    call solve_bratu(n, dx, lambda, u_ad, converged(2), newton_steps(2), &
                     krylov_steps(2), jacobian_vector_product_ad)
    call system_clock(end_count)
    ad_times(repetition) = real(end_count - start_count, dp) / real(clock_rate, dp)
  end do

  ! Guard against the optimizer removing the calls
  if (primal_acc == 0.0_dp) error stop "unexpected primal result"
  if (.not. all(converged)) error stop "a Newton-Krylov solve did not converge"
  if (any(newton_steps /= newton_steps(1))) error stop "Newton iteration counts differ"
  if (any(krylov_steps /= krylov_steps(1))) error stop "Krylov iteration counts differ"
  if (maxval(abs(u_ad - solution_ref)) >= 1.0e-5_dp) error stop "unexpected adjoint result"
  if (maxval(abs(u_ref - solution_ref)) >= 1.0e-5_dp) error stop "unexpected reference result"
  if (maxval(abs(u_ad - u_ref)) >= 1.0e-12_dp) error stop "JVP results differ"

  write(*, "(a,l1)") "both setups converged: ", all(converged)
  write(*, "(a,l1)") "solutions agree: ", maxval(abs(u_ad - u_ref)) < 1.0e-12_dp
  write(*, "(a,i0)") "Newton iterations per solve: ", newton_steps(1)
  write(*, "(a,i0)") "Krylov iterations per solve: ", krylov_steps(1)
  write(*, "(a,i0)") "repetitions: ", repetitions
  write(*, "(a40,es11.4)") "primal average wall time (s): ", sum(primal_times) / repetitions
  write(*, "(a40,es11.4)") "reference average wall time (s): ", sum(ref_times) / repetitions
  write(*, "(a40,es11.4)") "adjoint average wall time (s): ", sum(ad_times) / repetitions
  write(*, "(a40,es11.4)") "primal minimum wall time (s): ", minval(primal_times)
  write(*, "(a40,es11.4)") "reference minimum wall time (s): ", minval(ref_times)
  write(*, "(a40,es11.4)") "adjoint minimum wall time (s): ", 1 * minval(ad_times)
  write(*, "(a40,es11.4)") "average adjoint / reference time ratio: ", sum(ad_times) / sum(ref_times)
  write(*, "(a40,es11.4)") "average adjoint / primal time ratio: ", sum(ad_times) / sum(primal_times)
  write(*,*)
  write(*,*) "Table row:"
  write(*,"(a41,2(1x,'&',1x,i6),3(1x,'&',es11.4),1x,'&',f6.2,1x,'\%\\')") "1D Bratu (Jacobian-free)", &
    krylov_steps, sum(primal_times) / repetitions, sum(ref_times) / repetitions, &
    sum(ad_times) / repetitions, (sum(ad_times) / sum(ref_times) - 1.0) * 100
end program bratu_steady

! CHECK: both setups converged: T
! CHECK: solutions agree: T
! CHECK: Newton iterations per solve:
! CHECK: Krylov iterations per solve:
! CHECK: repetitions: 9
! CHECK: primal average wall time (s):
! CHECK: reference average wall time (s):
! CHECK: adjoint average wall time (s):
! CHECK: primal minimum wall time (s):
! CHECK: reference minimum wall time (s):
! CHECK: adjoint minimum wall time (s):
! CHECK: average adjoint / reference time ratio:
! CHECK: average adjoint / primal time ratio:
