! dtspg_example.c
!  A simple C program that solves a general symmetric block Toeplitz system using
!  the csp_dtspg_sv function.
! 
!  Author: Daniel Arguelles Martino <daniel.arguelles.martino@gmail.com>
!  Date: 1 August 2012
! 
!
! Solves the system
! 
!   Tx = b
! 
!  where
! 
!   T =
!    (  71, 101,  23,  30 )
!    ( 101, 194,  70,  96 )
!    (  23,  70,  71, 101 )
!    (  30,  96, 101, 194 )
! 
!   b =
!       (1  1  1  1)

PROGRAM dtsbg_example
  USE structpack
  IMPLICIT NONE
    INTEGER, PARAMETER :: n = 4, nb = 2
    REAL(KIND=r_kind), DIMENSION(n,1) :: x, b
    REAL(KIND=r_kind), DIMENSION(n,n) :: T
    
    T(1,1) = 71.0;  T(1,2) = 101.0; T(1,3) = 23.0;  T(1,4) = 30.0;
    T(2,1) = 101.0; T(2,2) = 194.0; T(2,3) = 70.0;  T(2,4) = 96.0;
    T(3,1) = 23.0;  T(3,2) = 70.0;  T(3,3) = 71.0;  T(3,4) = 101.0;
    T(4,1) = 30.0;  T(4,2) = 96.0;  T(4,3) = 101.0; T(4,4) = 194.0;
    
    x = 0.0
    b(1:4,1) = 1.0

    print *, "T:"
    print *, T(1,:)
    print *, T(2,:)
    print *, T(3,:)
    print *, T(4,:)
    
    PRINT*, "b:", b 

    CALL sp_dtspg_sv(T(1:2,:), b, x)

    PRINT*, "x:", x
    PRINT*, "b:", b
END PROGRAM dtsbg_example
 