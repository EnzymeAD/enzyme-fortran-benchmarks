! Shared driver for the OpenAD/F example benchmarks. Each benchmark wraps its
! example as y = f(x) and its Enzyme gradient as dx = J^T w; this module
! checks the gradient against central finite differences (a dot product test)
! and times repeated primal and gradient calls.
module openad_harness
  use, intrinsic :: iso_fortran_env, only: int64, real64
  implicit none
  private
  public :: read_repetitions, fill_deterministic, check_and_time

  abstract interface
    subroutine primal(n, x, m, y)
      import :: real64
      implicit none
      integer, intent(in) :: n, m
      real(real64), intent(in) :: x(n)
      real(real64), intent(out) :: y(m)
    end subroutine primal

    subroutine gradient(n, x, dx, m, w)
      import :: real64
      implicit none
      integer, intent(in) :: n, m
      real(real64), intent(in) :: x(n), w(m)
      real(real64), intent(out) :: dx(n)
    end subroutine gradient
  end interface

contains

  ! Repetition count from the first command-line argument, 9 by default.
  integer function read_repetitions()
    integer :: ios
    character(len=32) :: argument

    read_repetitions = 9
    if (command_argument_count() > 1) error stop "usage: benchmark [repetitions]"
    if (command_argument_count() == 1) then
      call get_command_argument(1, argument)
      read (argument, *, iostat=ios) read_repetitions
      if (ios /= 0 .or. read_repetitions < 1) &
        error stop "repetitions must be a positive integer"
    end if
  end function read_repetitions

  ! Fills v with values in [0.5, 1.5) that do not depend on the compiler's
  ! random number generator.
  subroutine fill_deterministic(v, seed)
    real(real64), intent(out) :: v(:)
    integer, intent(in) :: seed
    integer :: i

    do i = 1, size(v)
      v(i) = 0.5_real64 + modulo(0.6180339887498949_real64 * (i + seed), 1.0_real64)
    end do
  end subroutine fill_deterministic

  ! Compares <J^T w, v> with (<w, f(x + h v)> - <w, f(x - h v)>) / (2 h), stops
  ! when the relative error exceeds tol, and then times the primal and the
  ! gradient. dx returns the gradient.
  subroutine check_and_time(f, df, n, x, m, w, v, h, tol, repetitions, dx)
    procedure(primal) :: f
    procedure(gradient) :: df
    integer, intent(in) :: n, m, repetitions
    real(real64), intent(in) :: x(n), w(m), v(n), h, tol
    real(real64), intent(out) :: dx(n)
    real(real64) :: y(m), yp(m), ym(m), fd, ad, err
    real(real64) :: primal_times(repetitions), adjoint_times(repetitions)
    integer(int64) :: start_count, end_count, clock_rate
    integer :: repetition

    call df(n, x, dx, m, w)
    call f(n, x + h * v, m, yp)
    call f(n, x - h * v, m, ym)
    fd = (dot_product(w, yp) - dot_product(w, ym)) / (2 * h)
    ad = dot_product(dx, v)
    err = abs(ad - fd) / max(abs(fd), tiny(fd))
    write (*, "(a40,es11.4)") "adjoint directional derivative: ", ad
    write (*, "(a40,es11.4)") "finite difference: ", fd
    write (*, "(a40,es11.4)") "relative error: ", err
    if (.not. (err <= tol)) error stop "gradient check failed"
    write (*, "(a)") "gradient check: PASS"

    call system_clock(count_rate=clock_rate)
    do repetition = 1, repetitions
      call system_clock(start_count)
      call f(n, x, m, y)
      call system_clock(end_count)
      primal_times(repetition) = real(end_count - start_count, real64) / clock_rate
    end do
    do repetition = 1, repetitions
      call system_clock(start_count)
      call df(n, x, dx, m, w)
      call system_clock(end_count)
      adjoint_times(repetition) = real(end_count - start_count, real64) / clock_rate
    end do

    write (*, "(a,i0)") "repetitions: ", repetitions
    write (*, "(a40,es11.4)") "primal average wall time (s): ", sum(primal_times) / repetitions
    write (*, "(a40,es11.4)") "adjoint average wall time (s): ", sum(adjoint_times) / repetitions
    write (*, "(a40,es11.4)") "primal minimum wall time (s): ", minval(primal_times)
    write (*, "(a40,es11.4)") "adjoint minimum wall time (s): ", minval(adjoint_times)
    write (*, "(a40,es11.4)") "average adjoint / primal time ratio: ", &
      sum(adjoint_times) / sum(primal_times)
  end subroutine check_and_time

end module openad_harness
