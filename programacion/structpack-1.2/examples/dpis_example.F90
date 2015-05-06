! tpystrid_example.F90
! A simple Fortran program that solves a tridiagonal symmetric Toeplitz system using
! the sp_dpis_rojosv subroutine.
!
! Author: Pablo Martinez Naredo <pmnaredo@gmail.com>
! Date: 27 june 2011

! Solves the system
!
!  Tx = b
!
! where
!
!  T =
!      (1  2  0  0)
!      (2  1  2  0)
!      (0  2  1  2)
!      (0  0  2  1)
!
!  b =
!      (1  1  1  1)
!
 
PROGRAM tpystrid_example
  USE structpack
  USE sp_dpis
  IMPLICIT NONE
    INTEGER, PARAMETER :: n = 4
    REAL(KIND=r_kind) :: t0 = 1, t1 = 2
    REAL(KIND=r_kind), DIMENSION(n) :: x
    x(1:4) = 1_r_kind

    PRINT*, "t0:", t0
    PRINT*, "t1:", t1
    PRINT*, "b:", x 

    CALL sp_dpis_rojosv(n, t0, t1, x)

    PRINT*, "x:", x
END PROGRAM tpystrid_example
