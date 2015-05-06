module csp_toegen
    use iso_c_binding
    use sp_toegen
    implicit none
  
    public
    contains
    
    
    subroutine csp_dkms( a, eps, n ) bind(c)
        integer(c_int),                   intent(inout) :: n
        real(c_double), dimension(2*n-1), intent(inout) :: a
        real(c_double),                   intent(inout) :: eps
        
        call sp_dkms( a, eps, n )
    end subroutine csp_dkms
    
    subroutine csp_dsweet1( a, eps, n ) bind(c)
        integer(c_int),                   intent(inout) :: n
        real(c_double), dimension(2*n-1), intent(inout) :: a
        real(c_double),                   intent(inout) :: eps
        
        call sp_dsweet1( a, eps, n )
    end subroutine csp_dsweet1
    
    subroutine csp_dsweet2( a, eps, n ) bind(c)
        integer(c_int),                   intent(inout) :: n
        real(c_double), dimension(2*n-1), intent(inout) :: a
        real(c_double),                   intent(inout) :: eps
        
        call sp_dsweet2( a, eps, n )
    end subroutine csp_dsweet2
    
    subroutine csp_dsweet3( a, n ) bind(c)
        integer(c_int),                   intent(inout) :: n
        real(c_double), dimension(2*n-1), intent(inout) :: a
        
        call sp_dsweet3( a, n )
    end subroutine csp_dsweet3
    
    subroutine csp_dfreund0( a, n ) bind(c)
        integer(c_int),                   intent(inout) :: n
        real(c_double), dimension(2*n-1), intent(inout) :: a
        
        call sp_dfreund0( a, n )
    end subroutine csp_dfreund0
    
    subroutine csp_dfreund1( a, n ) bind(c)
        integer(c_int),                   intent(inout) :: n
        real(c_double), dimension(2*n-1), intent(inout) :: a
        
        call sp_dfreund1( a, n )
    end subroutine csp_dfreund1
    
    subroutine csp_dfreund2( a, n ) bind(c)
        integer(c_int),                   intent(inout) :: n
        real(c_double), dimension(2*n-1), intent(inout) :: a
        
        call sp_dfreund2( a, n )
    end subroutine csp_dfreund2

    
end module
