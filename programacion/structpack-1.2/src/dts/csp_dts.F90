!> \brief Interface for using the sp_dts module routines from C code.
!>
!> Appropriate C function prototypes for using these routines are provided
!> in the csp_dts.h header file.
!> \author StructPack Team
!> \date April, 2011
module csp_dts
  use iso_c_binding
  use sp_dts
  implicit none

  public :: csp_dts_sv, csp_dts_nrm1, csp_dts_gemv

  private
  contains

  !> \brief C interface to the sp_dts_sv subroutine of the sp_dts module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         t       (double*) Pointer to the first element of the first column (row) of the Toeplitz matrix.
  !> \param[out]        x       (double*) Pointer to the first element of solution vector.
  !> \param[in,out]     b       (double*) Pointer to the first element of the right hand side vector. Changed on output.
  !> \param[in]         nb      (int) Block size.
  !> \param[in]         piv     (int) Uses local diagonal pivoting if nonzero.
  !> \param[in]         affi    (int) Uses local core affinity.
  subroutine csp_dts_sv(n, t, x, b, nb, piv, affi ) bind(c)
    integer(kind=c_int), value, intent(in) :: n
    real(c_double), dimension(n), intent(in) :: t
    real(c_double), dimension(n), intent(out) :: x
    real(c_double), dimension(n), intent(inout) :: b
    integer(kind=c_int), value, intent(in) :: nb
    integer(kind=c_int), value, intent(in) :: piv
    integer(kind=c_int), value, intent(in) :: affi

    logical :: af  
    af = affi.ne.0

    if (piv .ne. 0) then
      call sp_dts_sv(t, x, b, nb, .true., af )
    else
      call sp_dts_sv(t, x, b, nb, .false., af )
    end if
  end subroutine csp_dts_sv


  !> \brief C interface to the sp_dts_nrm1 function of the sp_dts module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         t       (double*) Pointer to the first element of the first column (row) of the Toeplitz matrix.
  !> \return            (double) Returned value of ctpsynrm1 (1-norm of the symmetric Toeplitz matrix).
  function csp_dts_nrm1(n, t) bind(c)
    real(kind=c_double) :: csp_dts_nrm1
    integer(kind=c_int), value, intent(in) :: n
    real(kind=c_double), dimension(n), intent(in) :: t

    csp_dts_nrm1 = sp_dts_nrm1(t)
  end function csp_dts_nrm1


  !> \brief C interface to the sp_dts_gemv subroutine of the sp_dts module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         alpha   (double) Value alpha.
  !> \param[in]         t       (double*) Pointer to the first element of the first column (row) of the Toeplitz matrix.
  !> \param[out]        x       (double*) Pointer to the first element of vector x.
  !> \param[in,out]     b       (double*) Pointer to the first element of vector b.
  subroutine csp_dts_gemv(n, alpha, t, x, b) bind(c)
    integer(kind=c_int), value, intent(in) :: n
    real(kind=c_double), value, intent(in) :: alpha
    real(kind=c_double), dimension(n), intent(in) :: t, x
    real(kind=c_double), dimension(n), intent(inout) :: b

    call sp_dts_gemv( alpha, t, x, b )
  end subroutine csp_dts_gemv


end module csp_dts
