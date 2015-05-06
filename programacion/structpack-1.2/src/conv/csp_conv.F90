!> \brief Interface for using the sp_conv module routines from C code.
!>
!> Appropriate C function prototypes for using these routines are provided
!> in the csp_conv.h header file.
!> \author StructPack Team
!> \date April, 2012
module csp_conv
    use iso_c_binding
    use sp_conv
    implicit none
  
    public :: csp_convn, csp_convs

    private 
    contains
    
    !> \brief C interface to the sp_convn subroutine of the sp_conv module which computes y = beta*y+alpha*T*x with T non-symmetric.
    !> \param[in]         n       (int) Size of matrix.
    !> \param[in]         alpha   (double)  Scalar alpha. 
    !> \param[in]         t       (double*) Pointer to the first element of the first column of the Toeplitz matrix.
    !> \param[in]         u       (double*) Pointer to the first element of the first row of the Toeplitz matrix.
    !> \param[out]        x       (double*) Pointer to the first element of vector x.
    !> \param[in]         beta    (double)  Scalar beta. 
    !> \param[out]        y       (double*) Pointer to the first element of vector y. 
    subroutine csp_convn( n, alpha, t, u, x, beta, y ) bind(c)
        integer(kind=c_int), value,        intent(in)  :: n
        real(kind=c_double), dimension(n), intent(in)  :: t, u, x
        real(kind=c_double), dimension(n), intent(out) :: y
        real(kind=c_double), value       , intent(in)  :: alpha, beta
        
        call conv( alpha, t, x, beta, y, u )
    end subroutine csp_convn
    
    !> \brief C interface to the sp_convn subroutine of the sp_conv module which computes y = beta*y+alpha*T*x with T symmetric.
    !> \param[in]         n       (int) Size of matrix.
    !> \param[in]         alpha   (double)  Scalar alpha. 
    !> \param[in]         t       (double*) Pointer to the first element of the first column (row) of the Toeplitz matrix.
    !> \param[out]        x       (double*) Pointer to the first element of vector x.
    !> \param[in]         beta    (double)  Scalar beta. 
    !> \param[out]        y       (double*) Pointer to the first element of vector y. 
    subroutine csp_convs( n, alpha, t, x, beta, y ) bind(c)
        integer(kind=c_int), value,        intent(in)  :: n
        real(kind=c_double), dimension(n), intent(in)  :: t, x
        real(kind=c_double), dimension(n), intent(out) :: y
        real(kind=c_double), value       , intent(in)  :: alpha, beta
        
        call conv( alpha, t, x, beta, y )
    end subroutine csp_convs
    
end module csp_conv
