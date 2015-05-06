! dts_example.F90
! A simple Fortran program that solves a general non symmetric Toeplitz system
! using the sp_dt_sv subroutine.
!
! Author: Daniel Arguelles Martino <daniel.arguelles.martino@gmail.com>
! Date: 1 August 2012

! Solves the system
!
!  Tx = b
!
! where
!
!  T =
!      (1  2  3  4)
!      (5  1  2  3)
!      (6  2  1  2)
!      (7  3  2  1)
!
!  b =
!      (1  1  1  1)
!
 
PROGRAM dt_example
  USE structpack
  IMPLICIT NONE
    INTEGER, PARAMETER :: n = 4, nb = 2
    INTEGER :: i, INFO
    REAL(KIND=r_kind), DIMENSION(n) :: u,v, x, b
    u = (/(i,i=1,4)/)
    v = (/(i,i=4,7)/)
    v(1) = 1
    x = 0.0
    b(1:4) = 1.0

    PRINT*, "u:", u
    PRINT*, "v:", v
    PRINT*, "b:", b 

    CALL sp_dt_sv(u,v, x, b, nb, .TRUE., 1, INFO)

    PRINT*, "x:", x
    PRINT*, "b:", b
END PROGRAM dt_example
