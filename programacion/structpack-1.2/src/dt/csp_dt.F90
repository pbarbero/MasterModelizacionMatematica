!> \brief Interface for using the sp_dt module routines from C code.
!>
!> Appropriate C function prototypes for using these routines are provided
!> in the csp_dt.h header file.
!> \author Daniel Arguelles Martino
!> \date February, 2012
module csp_dt
    use iso_c_binding
    use sp_dt
    implicit none

    public :: csp_dt_sv, csp_dt_nrm1, csp_dt_gemv

    private
    contains
    
    !> \brief C interface to the sp_dt_sv subroutine of the sp_dt module.
    !> \param[in]         n       (int) Size of matrix.
    !> \param[in]         u       (double*) Pointer to the first element of the first column of the Toeplitz matrix.
    !> \param[in]         v       (double*) Pointer to the first element of the first row of the Toeplitz matrix.
    !> \param[out]        x       (double*) Pointer to the first element of solution vector.
    !> \param[in,out]     b       (double*) Pointer to the first element of the right hand side vector. Changed on output.
    !> \param[in]         nb      (int) Block size.
    !> \param[in]         piv     (int) Uses local diagonal pivoting if nonzero.
    !> \param[in]         ref     (int) Number of iterative refinement steps.   
    !> \param[out]        INFO    (int*) = 0: successful exit
    !>                              = 1: first element of u and first element of v are different
    !>                              = 2: u and v have different length
    !>                              = 3: u and b have different length
    subroutine csp_dt_sv(n, u, v, x, b, nb, piv, ref, INFO) bind(c)
        integer(kind=c_int), value,   intent(in   ) :: n
        real(c_double), dimension(n), intent(in   ) :: u, v
        real(c_double), dimension(n), intent(out  ) :: x
        real(c_double), dimension(n), intent(inout) :: b
        integer(kind=c_int), value,   intent(in   ) :: nb
        integer(kind=c_int), value,   intent(in   ) :: piv
        integer(kind=c_int), value,   intent(in   ) :: ref
        integer(kind=c_int),          intent(inout) :: INFO
        
        if (piv .ne. 0) then
            call sp_dt_sv( u, v, x, b, nb, .true., ref, INFO )
        else
            call sp_dt_sv( u, v, x, b, nb, .false., ref, INFO )
        end if
    end subroutine csp_dt_sv


    !> \brief C interface to the sp_dt_nrm1 function of the sp_dt module.
    !> \param[in]         n       (int) Size of matrix.
    !> \param[in]         u       (double*) Pointer to the first element of the first column of the Toeplitz matrix.
    !> \param[in]         v       (double*) Pointer to the first element of the first row of the Toeplitz matrix.
    !> \return            (double) Returned value of 1-norm of the general Toeplitz matrix.
    function csp_dt_nrm1(n, u, v) bind(c)
        real   ( kind=c_double )                           :: csp_dt_nrm1
        integer( kind=c_int ),    value,        intent(in) :: n
        real   ( kind=c_double ), dimension(n), intent(in) :: u, v

        csp_dt_nrm1 = sp_dt_nrm1( u, v )
    end function csp_dt_nrm1

    
    !> \brief C interface to the sp_dt_gemv subroutine of the sp_dt module.
    !> \param[in]         n       (int) Size of matrix.
    !> \param[in]         alpha   (double) Value alpha.
    !> \param[in]         u       (double*) Pointer to the first element of the first column of the Toeplitz matrix.
    !> \param[in]         v       (double*) Pointer to the first element of the first row of the Toeplitz matrix.
    !> \param[out]        x       (double*) Pointer to the first element of vector x.
    !> \param[in,out]     b       (double*) Pointer to the first element of vector b.
    subroutine csp_dt_gemv(n, alpha, u, v, x, beta, b) bind(c)
        integer(kind=c_int), value,        intent(in   ) :: n
        real(kind=c_double), value,        intent(in   ) :: alpha, beta
        real(kind=c_double), dimension(n), intent(in   ) :: u, v, x
        real(kind=c_double), dimension(n), intent(inout) :: b

        call sp_dt_gemv( alpha, u, v, x, beta, b )
    end subroutine csp_dt_gemv

end module csp_dt
