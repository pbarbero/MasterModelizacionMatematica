!> \brief Interface for using the sp_dpis module routines from C code.
!>
!> Appropriate C function prototypes for using these routines are provided
!> in the csp_dpis.h header file.
!> \author StructPack Team
!> \date April, 2011
module csp_dpis
  use iso_c_binding
  use sp_dpis
  implicit none

  public :: csp_dpis_rojosv, csp_dpis_dstsv, &
    csp_dpis_ldltsv, csp_dpis_maxsvd, &
    csp_dpis_minsvd, csp_dpis_gemv

  private
  contains

  !> \brief C interface to the sp_dpis_rojosv subroutine of the sp_dpis module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         t0      (double) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (double) Sub/Super--diagonal of the Toeplitz matrix.
  !> \param[out]        x       (double*) Pointer to the first element of solution vector.
  subroutine csp_dpis_rojosv(n, t0, t1, x) bind(c)
    integer(kind=c_int), value, intent(in) :: n
    real(kind=c_double), value, intent(in) :: t0, t1
    real(c_double), dimension(*), intent(inout) :: x

    call sp_dpis_rojosv(n, t0, t1, x)
  end subroutine csp_dpis_rojosv


  !> \brief C interface to the sp_dpis_dstsv subroutine of the sp_dpis module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         t0      (double) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (double) Sub/Super--diagonal of the Toeplitz matrix.
  !> \param[out]        x       (double*) Pointer to the first element of solution vector.
  subroutine csp_dpis_dstsv(n, t0, t1, x) bind(c)
    integer(kind=c_int), value, intent(in) :: n
    real(kind=c_double), value, intent(in) :: t0, t1
    real(c_double), dimension(*), intent(inout) :: x

    call sp_dpis_dstsv(n, t0, t1, x)
  end subroutine csp_dpis_dstsv


  !> \brief C interface to the sp_dpis_ldltsv subroutine of the sp_dpis module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         t0      (double) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (double) Sub/Super--diagonal of the Toeplitz matrix.
  !> \param[out]        x       (double*) Pointer to the first element of solution vector.
  !> \return                    0 if everything went fine, NON 0 if an error took place
  function csp_dpis_ldltsv(n, t0, t1, x) bind(c)
    integer(kind=c_int) :: csp_dpis_ldltsv
    integer(kind=c_int), value, intent(in) :: n
    real(kind=c_double), value, intent(in) :: t0, t1
    real(c_double), dimension(*), intent(inout) :: x

    integer(kind=c_int) :: status

    call sp_dpis_ldltsv(n, t0, t1, x, status)
    csp_dpis_ldltsv = status

  end function csp_dpis_ldltsv


  !> \brief C interface to the sp_dpis_maxsvd function of the sp_dpis module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         t0      (double) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (double) Sub/Super--diagonal of the Toeplitz matrix.
  !> \return                    (double) Maximum SVD of the symmetric tridiagonal Toeplitz matrix.
  function csp_dpis_maxsvd(n, t0, t1) bind(c)
    real(kind=c_double) :: csp_dpis_maxsvd
    integer(kind=c_int), value, intent(in) :: n
    real(kind=c_double), value, intent(in) :: t0, t1

    csp_dpis_maxsvd = sp_dpis_maxsvd(n, t0, t1)
  end function csp_dpis_maxsvd


  !> \brief C interface to the sp_dpis_minsvd function of the sp_dpis module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         t0      (double) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (double) Sub/Super--diagonal of the Toeplitz matrix.
  !> \return                    (double) Minimum SVD of the symmetric tridiagonal Toeplitz matrix.
  function csp_dpis_minsvd(n, t0, t1) bind(c)
    real(kind=c_double) :: csp_dpis_minsvd
    integer(kind=c_int), value, intent(in) :: n
    real(kind=c_double), value, intent(in) :: t0, t1

    csp_dpis_minsvd = sp_dpis_minsvd(n, t0, t1)
  end function csp_dpis_minsvd


  !> \brief C interface to the sp_dpis_gemv subroutine of the sp_dpis module.
  !> \param[in]         n       (int) Size of matrix.
  !> \param[in]         alpha   (double) Value alpha.
  !> \param[in]         t0      (double) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (double) Sub/Super--diagonal of the Toeplitz matrix.
  !> \param[out]        x       (double*) Pointer to the first element of vector x.
  !> \param[in,out]     b       (double*) Pointer to the first element of vector b.
  subroutine csp_dpis_gemv(n, alpha, t0, t1, x, b) bind(c)
    integer(kind=c_int), value, intent(in) :: n
    real(kind=c_double), value, intent(in) :: alpha, t0, t1
    real(kind=c_double), dimension(*), intent(in) :: x
    real(kind=c_double), dimension(*), intent(inout) :: b

    call sp_dpis_gemv( n, alpha, t0, t1, x, b )
  end subroutine csp_dpis_gemv


end module csp_dpis
