! REQUIRES: fortran
! RUN: %fc -flto -O2 -c %loadFortran %sundialsFortran %s -o /dev/stdout | %opt %loadEnzyme %enzyme -o %t.ll && %fc -flto -O2 %t.ll %sundialsLibs -o %t1 && %t1 9 | FileCheck %s
! RUN: %if flangenzyme %{ %fc -O2 %loadFortran %loadFlangEnzyme %sundialsFortran %s %sundialsLibs -o %t2 && %t2 9 | FileCheck %s %}

! Timing benchmark for the Robertson chemical kinetics problem
!
!   y1' = -k1 y1 + k2 y2 y3
!   y2' =  k1 y1 - k2 y2 y3 - k3 y2^2
!   y3' =  k3 y2^2
!
! integrated with SUNDIALS' CVODE (BDF) using a Jacobian-free Newton-Krylov
! linear solver (SPGMR). The Jacobian-vector products CVODE asks for through
! its JacTimes callback are obtained either from Enzyme forward mode on the
! right-hand side routine (adjoint) or from a reference hand-derived JVP. The
! complete CVODE solve is timed; only the derivative callback differs between
! the two timed runs. The primal run times the right-hand side evaluation
! itself.

module robertson_model
  use, intrinsic :: iso_c_binding, only: c_int64_t, c_double
  implicit none
  private
  public :: neq, rates, jacobian_vector_product_ad, jacobian_vector_product_ref
  integer(c_int64_t), parameter :: neq = 3
  real(c_double), parameter :: k1 = 0.04d0, k2 = 1.0d4, k3 = 3.0d7
contains

  subroutine rates(n, y, ydot)
    integer, intent(in) :: n
    real(c_double), intent(in) :: y(n)
    real(c_double), intent(out) :: ydot(n)

    ydot(1) = -k1 * y(1) + k2 * y(2) * y(3)
    ydot(3) = k3 * y(2)**2
    ydot(2) = -ydot(1) - ydot(3)
  end subroutine rates

  subroutine jacobian_vector_product_ad(n, y, v, jv)
    use enzyme, only: enzyme_const, enzyme_dup, enzyme_fwddiff
    integer, intent(in) :: n
    real(c_double), intent(in) :: y(n), v(n)
    real(c_double), intent(out) :: jv(n)
    real(c_double) :: ydot(n)

    ydot = 0.0d0
    jv = 0.0d0
    call enzyme_fwddiff(rates, enzyme_const, n, enzyme_dup, y, v, &
                         enzyme_dup, ydot, jv)
  end subroutine jacobian_vector_product_ad

  ! Reference implementation for the Jacobian-vector product
  subroutine jacobian_vector_product_ref(n, y, v, jv)
    integer, intent(in) :: n
    real(c_double), intent(in) :: y(n), v(n)
    real(c_double), intent(out) :: jv(n)

    jv(1) = -k1 * v(1) + k2 * (y(3) * v(2) + y(2) * v(3))
    jv(3) = 2.0d0 * k3 * y(2) * v(2)
    jv(2) = -jv(1) - jv(3)
  end subroutine jacobian_vector_product_ref
end module robertson_model

module robertson_cvode
  use, intrinsic :: iso_c_binding, only: c_int, c_double, c_long, c_ptr, c_funptr, &
    c_associated, c_funloc, c_null_funptr
  use fsundials_core_mod, only: N_Vector, SUNMatrix, SUNLinearSolver, FSUNContext_Create, &
    FSUNContext_Free, FSUNLinSolFree, FN_VDestroy, FN_VGetArrayPointer, &
    SUN_COMM_NULL, SUN_PREC_NONE
  use fcvode_mod, only: FCVodeCreate, FCVodeInit, FCVodeSVtolerances, FCVodeSetLinearSolver, &
    FCVodeSetJacTimes, FCVode, FCVodeGetNumSteps, FCVodeGetNumJtimesEvals, &
    FCVodeFree, CV_BDF, CV_NORMAL, CV_SUCCESS
  use fnvector_serial_mod, only: FN_VMake_Serial
  use fsunlinsol_spgmr_mod, only: FSUNLinSol_SPGMR
  use robertson_model, only: neq, rates, jacobian_vector_product_ad, jacobian_vector_product_ref
  implicit none
  private
  public :: solve_robertson_jtimes
  integer(c_long) :: counts(1)

! The complete CVODE solve is timed; only its JacTimes callback changes.
contains

  integer(c_int) function rhs(t, sunvec_y, sunvec_f, user_data) result(ierr) bind(C)
    real(c_double), value :: t
    type(N_Vector), intent(inout) :: sunvec_y, sunvec_f
    type(c_ptr), value :: user_data
    real(c_double), pointer, contiguous :: y(:), f(:)

    y => FN_VGetArrayPointer(sunvec_y)
    f => FN_VGetArrayPointer(sunvec_f)
    call rates(int(neq), y, f)
    ierr = 0
  end function rhs

  integer(c_int) function jtimes_ad(sunvec_v, sunvec_jv, t, sunvec_y, sunvec_fy, &
                                   user_data, sunvec_tmp) result(ierr) bind(C)
    type(N_Vector), intent(inout) :: sunvec_v, sunvec_jv, sunvec_y, sunvec_fy, sunvec_tmp
    real(c_double), value :: t
    type(c_ptr), value :: user_data
    real(c_double), pointer, contiguous :: v(:), jv(:), y(:)

    v => FN_VGetArrayPointer(sunvec_v)
    jv => FN_VGetArrayPointer(sunvec_jv)
    y => FN_VGetArrayPointer(sunvec_y)
    call jacobian_vector_product_ad(int(neq), y, v, jv)
    ierr = 0
  end function jtimes_ad

  integer(c_int) function jtimes_ref(sunvec_v, sunvec_jv, t, sunvec_y, sunvec_fy, &
                                       user_data, sunvec_tmp) result(ierr) bind(C)
    type(N_Vector), intent(inout) :: sunvec_v, sunvec_jv, sunvec_y, sunvec_fy, sunvec_tmp
    real(c_double), value :: t
    type(c_ptr), value :: user_data
    real(c_double), pointer, contiguous :: v(:), jv(:), y(:)

    v => FN_VGetArrayPointer(sunvec_v)
    jv => FN_VGetArrayPointer(sunvec_jv)
    y => FN_VGetArrayPointer(sunvec_y)
    call jacobian_vector_product_ref(int(neq), y, v, jv)
    ierr = 0
  end function jtimes_ref

  subroutine solve_robertson_jtimes(use_ad, yval, nst, jvp_evals)
    logical, intent(in) :: use_ad
    real(c_double), intent(out) :: yval(neq)
    integer(c_long), intent(out) :: nst, jvp_evals
    integer, parameter :: nout = 12
    real(c_double) :: rtol, tout, tret(1), avtol(neq)
    integer(c_int) :: retval, iout
    type(N_Vector), pointer :: sunvec_y, sunvec_av
    type(SUNMatrix), pointer :: sunmat
    type(SUNLinearSolver), pointer :: sunls
    type(c_ptr) :: cvode_mem, sunctx
    type(c_funptr) :: jtimes_callback

    retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)
    yval = [1.0d0, 0.0d0, 0.0d0]
    avtol = [1.0d-12, 1.0d-18, 1.0d-12]
    rtol = 1.0d-8
    tout = 0.4d0
    if (.not. c_associated(sunctx)) error stop "Cannot create SUNContext"
    sunvec_y => FN_VMake_Serial(neq, yval, sunctx)
    sunvec_av => FN_VMake_Serial(neq, avtol, sunctx)
    if (.not. (associated(sunvec_y) .and. associated(sunvec_av))) stop 1
    cvode_mem = FCVodeCreate(CV_BDF, sunctx)
    if (.not. c_associated(cvode_mem)) stop 1
    retval = FCVodeInit(cvode_mem, c_funloc(rhs), 0.0d0, sunvec_y)
    if (retval /= 0) stop 1
    retval = FCVodeSVtolerances(cvode_mem, rtol, sunvec_av)
    if (retval /= 0) stop 1
    sunmat => null()
    sunls => FSUNLinSol_SPGMR(sunvec_y, SUN_PREC_NONE, 5, sunctx)
    if (.not. associated(sunls)) stop 1
    retval = FCVodeSetLinearSolver(cvode_mem, sunls, sunmat)
    if (retval /= 0) stop 1
    if (use_ad) then
      jtimes_callback = c_funloc(jtimes_ad)
    else
      jtimes_callback = c_funloc(jtimes_ref)
    end if
    retval = FCVodeSetJacTimes(cvode_mem, c_null_funptr, jtimes_callback)
    if (retval /= 0) stop 1

    do iout = 1, nout
      retval = FCVode(cvode_mem, tout, sunvec_y, tret(1), CV_NORMAL)
      if (retval /= CV_SUCCESS) error stop "FCVode failed"
      tout = tout * 10.0d0
    end do
    retval = FCVodeGetNumSteps(cvode_mem, counts)
    nst = counts(1)
    retval = FCVodeGetNumJtimesEvals(cvode_mem, counts)
    jvp_evals = counts(1)
    if (retval /= 0) error stop "Cannot read CVODE statistics"
    call FCVodeFree(cvode_mem)
    retval = FSUNLinSolFree(sunls)
    call FN_VDestroy(sunvec_y)
    call FN_VDestroy(sunvec_av)
    retval = FSUNContext_Free(sunctx)
  end subroutine solve_robertson_jtimes
end module robertson_cvode

program robertson_cvode_jtimes
  use, intrinsic :: iso_c_binding, only: c_double, c_long
  use, intrinsic :: iso_fortran_env, only: int64
  use robertson_model, only: neq, rates, jacobian_vector_product_ad, jacobian_vector_product_ref
  use robertson_cvode, only: solve_robertson_jtimes
  implicit none
  integer :: args, repetition, repetitions, ios
  integer(c_long) :: nst(2), jvp_evals(2)
  integer(int64) :: start_count, end_count, clock_rate
  real(c_double) :: y(neq), v(neq), ydot(neq), jv_ad(neq), jv_ref(neq)
  real(c_double) :: final_y(neq, 2)
  real(c_double) :: primal_acc
  real(c_double), allocatable :: primal_times(:), ref_times(:), ad_times(:)
  character(len=32) :: repetitions_arg

  repetitions = 9
  args = command_argument_count()
  if (args > 1) error stop "usage: robertson_cvode_jtimes [repetitions]"
  if (args == 1) then
    call get_command_argument(1, repetitions_arg)
    read (repetitions_arg, *, iostat=ios) repetitions
    if (ios /= 0 .or. repetitions < 1) error stop "repetitions must be a positive integer"
  end if
  allocate (primal_times(repetitions), ref_times(repetitions), ad_times(repetitions))

  y = [0.9d0, 1.0d-5, 0.1d0]
  v = [0.1d0, 0.2d0, 0.3d0]
  call jacobian_vector_product_ad(int(neq), y, v, jv_ad)
  call jacobian_vector_product_ref(int(neq), y, v, jv_ref)
  if (maxval(abs(jv_ad - jv_ref)) > 1.0d-12 * maxval(abs(jv_ref))) then
    error stop "reference JVP is incorrect"
  end if
  primal_acc = 0.0d0
  call system_clock(count_rate=clock_rate)

  ! Time the primal (right-hand side) calls
  do repetition = 1, repetitions
    call system_clock(start_count)
    call rates(int(neq), y, ydot)
    call system_clock(end_count)
    primal_times(repetition) = real(end_count - start_count, c_double) / real(clock_rate, c_double)
    primal_acc = primal_acc + sum(ydot)
  end do

  ! Time the reference (hand-coded JVP) solves
  do repetition = 1, repetitions
    call system_clock(start_count)
    call solve_robertson_jtimes(.false., final_y(:, 1), nst(1), jvp_evals(1))
    call system_clock(end_count)
    ref_times(repetition) = real(end_count - start_count, c_double) / real(clock_rate, c_double)
  end do

  ! Time the Enzyme-differentiated (adjoint) solves
  do repetition = 1, repetitions
    call system_clock(start_count)
    call solve_robertson_jtimes(.true., final_y(:, 2), nst(2), jvp_evals(2))
    call system_clock(end_count)
    ad_times(repetition) = real(end_count - start_count, c_double) / real(clock_rate, c_double)
  end do

  ! Guard against the optimizer removing the calls
  if (primal_acc == 0.0d0) error stop "unexpected primal result"
  if (any(nst <= 0)) error stop "a CVODE solve did not run"
  if (any(jvp_evals <= 0)) error stop "a solver did not use the supplied JVP"
  if (any(abs(sum(final_y, dim=1) - 1.0d0) > 1.0d-8)) error stop "mass is not conserved"
  if (maxval(abs(final_y(:, 1) - final_y(:, 2))) > 1.0d-8) error stop "solutions differ"

  write (*, "(a,l1)") "both CVODE setups completed: ", all(nst > 0)
  write (*, "(a,l1)") "final solutions agree: ", &
    maxval(abs(final_y(:, 1) - final_y(:, 2))) <= 1.0d-8
  write (*, "(a,i0,1x,i0)") "CVODE steps (reference, adjoint): ", nst
  write (*, "(a,i0,1x,i0)") "JVP evaluations (reference, adjoint): ", jvp_evals
  write (*, "(a,i0)") "repetitions: ", repetitions
  write (*, "(a40,es11.4)") "primal average wall time (s): ", sum(primal_times) / repetitions
  write (*, "(a40,es11.4)") "reference average wall time (s): ", sum(ref_times) / repetitions
  write (*, "(a40,es11.4)") "adjoint average wall time (s): ", sum(ad_times) / repetitions
  write (*, "(a40,es11.4)") "primal minimum wall time (s): ", minval(primal_times)
  write (*, "(a40,es11.4)") "reference minimum wall time (s): ", minval(ref_times)
  write (*, "(a40,es11.4)") "adjoint minimum wall time (s): ", minval(ad_times)
  write (*, "(a40,es11.4)") "average adjoint / reference time ratio: ", sum(ad_times) / sum(ref_times)
  write (*, "(a40,es11.4)") "average adjoint / primal time ratio: ", sum(ad_times) / sum(primal_times)
end program robertson_cvode_jtimes

! CHECK: both CVODE setups completed: T
! CHECK: final solutions agree: T
! CHECK: CVODE steps (reference, adjoint):
! CHECK: JVP evaluations (reference, adjoint):
! CHECK: repetitions: 9
! CHECK: primal average wall time (s):
! CHECK: reference average wall time (s):
! CHECK: adjoint average wall time (s):
! CHECK: primal minimum wall time (s):
! CHECK: reference minimum wall time (s):
! CHECK: adjoint minimum wall time (s):
! CHECK: average adjoint / reference time ratio:
! CHECK: average adjoint / primal time ratio:
