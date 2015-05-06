module csp_dtspg
    use iso_c_binding
    use sp_dtspg
    
    implicit none
    public ::csp_dtspg_sv, csp_dtspg_rg,csp_dtspg_gemv,csp_dtspg_nrm1
    
    private
    contains
    
    !> M integer number of rows of op( T )
    subroutine csp_dtspg_sv( M, N, T, b, x, K )  bind(c)
        integer(kind=c_int), value,               intent(in)    :: M , N, K
        real(c_double), dimension(N,M), intent(in)    :: T
        real(c_double), dimension(max(M,N),K), intent(in)    :: b
        real(c_double), dimension(max(M,N),K), intent(inout) :: x
        
        CALL sp_dtspg_sv( T, b, x )
        
    end subroutine
    
    subroutine csp_dtspg_rg( nBlocks, bSize, seed, T ) bind(c)
        integer(kind=c_int), value, intent(in) :: nBlocks,bSize,seed
        real(c_double), dimension(nBlocks*bSize,nBlocks*bSize), intent(inout)    :: T
        
        CALL sp_dtspg_rg( nBlocks, bSize, seed, T )
       ! T(:,1:bSize) = TRANSPOse(T(1:bSize,:))
    end subroutine
    
    real(c_double) function csp_dtspg_nrm1( nBlocks, bSize, T ) bind(c)
        integer(kind=c_int), value, intent(in) :: nBlocks, bSize
        real(c_double), dimension(bSize,bSize*nBlocks), intent(in) :: T

        csp_dtspg_nrm1 = sp_dtspg_nrm1( T )
        
    end function
    
    subroutine csp_dtspg_gemv( nBlocks, bSize, alfa, T, x, beta, b ) bind(c)
        integer(kind=c_int), value, intent(in) :: nBlocks, bSize
        real(c_double), dimension(bSize,bSize*nBlocks), intent(in) :: T
        real(c_double), dimension(bSize*nBlocks), intent(in) :: x
        real(c_double), dimension(bSize*nBlocks), intent(inout) :: b
        real(c_double), value, intent(in) :: alfa, beta
        
       
       CALL sp_dtspg_gemv( alfa, T, x, beta, b )
        
    end subroutine    
        
        
end module
