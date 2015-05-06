module sp_toegen
    use sp_utils
    implicit none
    
    REAL( kind=r_kind ), PARAMETER, PRIVATE ::  START_VALUE = 0.5d+0, &
                                                FACTOR      = 2.0d+0
                                                

    public    
    contains 
    
    subroutine sp_dkms( a, eps, n )
        REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: a
        REAL( kind=r_kind ),               INTENT( inout ) :: eps
        INTEGER,                           INTENT( inout ) :: n
        
        INTEGER :: i
        
        a(n) = eps
        if( 1 .le. n ) then
            a(n-1) = START_VALUE            
            a(n+1) = a(n-1)
            do i=2, n
                a(n-i) = a(n+i-1)/FACTOR
                a(n+i) = a(n-i)
            end do
        end if
    end subroutine sp_dkms
    
    subroutine sp_dsweet1( a, eps, n )
        REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: a
        REAL( kind=r_kind ),               INTENT( inout ) :: eps
        INTEGER,                           INTENT( inout ) :: n

        eps = 5d-8
        n=6
        a(1) = 3.0d+0
        a(2) = 2.0d+0
        a(3) = 6.0d+0
        a(4) = 1.0d+0
        a(5) = 8.0d+0
        a(6) = 4.0d+0
        a(7) = 6.0d+0
        a(8) = 71.0d+0/15.0d+0 + eps
        a(9) = 5.0d+0
        a(10)= 3.0d+0
        a(11)= 1.0d+0
    
    end subroutine sp_dsweet1
    
    subroutine sp_dsweet2( a, eps, n )
        REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: a
        REAL( kind=r_kind ),               INTENT( inout ) :: eps
        INTEGER,                           INTENT( inout ) :: n

        eps = 5d-13
        n=6
        a(1) = 3.0d+0
        a(2) = 2.0d+0
        a(3) = 6.0d+0
        a(4) = 1.0d+0
        a(5) = 4.0d+0
        a(6) = 8.0d+0
        a(7) = 4.0d+0
        a(8) = 34.0d+0 + eps
        a(9) = 5.0d+0
        a(10)= 3.0d+0
        a(11)= 1.0d+0
    end subroutine sp_dsweet2
    
    subroutine sp_dsweet3( a, n )
        REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: a
        INTEGER,                           INTENT( inout ) :: n

        n = 13
        a(1) = -15.0d+0
        a(2) = 10.0d+0
        a(3) = 1.0d+0
        a(4) = 7.0d+0
        a(5) = 2.0d+0
        a(6) = 5.0d+0
        a(7) = 3.0d+0
        a(8) = 5.850d+0
        a(9) = 5.697d+0
        a(10)= 2.0d+0
        a(11)= 6.0d+0
        a(12)= -1.0d+0
        a(13)= 5.0d+0
        a(14)= 1.0d+0
        a(15)= 3.0d+0
        a(16)= 12.755d+0
        a(17)= 19.656d+0
        a(18)= 28.361d+0
        a(19)= 7.0d+0
        a(20)= 1.0d+0
        a(21)= 2.0d+0
        a(22)= 1.0d+0
        a(23)= -6.0d+0
        a(24)= 1.0d+0
        a(25)= 0.5d+0
    end subroutine sp_dsweet3
    
    subroutine sp_dfreund0( a, n )
    REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: a
        INTEGER,                           INTENT( inout ) :: n
        n=7
        a(1)=1.0d+0
        a(2)=0.0d+0
        a(3)=1.0d+0
        a(4)=1.0d+0
        a(5)=0.0d+0
        a(6)=1.0d+0
        a(7)=0.0d+0
        a(8)=1.0d+0
        a(9)=0.0d+0
        a(10)=1.0d+0
        a(11)=1.0d+0
        a(12)=0.0d+0
        a(13)=1.0d+0
    end subroutine sp_dfreund0

    subroutine sp_dfreund1( a, n )
REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: a
        INTEGER,                           INTENT( inout ) :: n
        n=5
        a(1)=1.21785304238395d+0
        a(2)=1.06413364195684d+0
        a(3)=-1.62115749363923d+0
        a(4)=1.27324683138786d+0
        a(5)=-1.00000000000001d+0
        a(6)=0.78539366864947d+0
        a(7)=3.41046741401696d+0
        a(8)=-17.92422495778239d+0
        a(9)=38.20692196916536d+0
    end subroutine sp_dfreund1

    subroutine sp_dfreund2( a, n )
        REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: a
        INTEGER,                           INTENT( inout ) :: n
        n=6
        a(1)=4.51853189291597d+0
        a(2)=-2.22889818997626d+0
        a(3)=1.16717755769466d+0
        a(4)=-1.10855680502906d+0
        a(5)=1.05288024249153d+0
        a(6)=-0.99999999999998d+0
        a(7)=0.94977563415339d+0
        a(8)=3.85673107101965d+0
        a(9)=-13.61721591570147d+0
        a(10)=3.81850412563076d+0
        a(11)=73.05176317918625d+0
    end subroutine sp_dfreund2
    
end module
