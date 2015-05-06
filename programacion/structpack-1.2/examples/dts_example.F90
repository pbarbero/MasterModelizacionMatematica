! dts_example.F90
! A simple Fortran program that solves a general symmetric Toeplitz system
! using the sp_dts_sv subroutine.
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
!      (1  2  3  4)
!      (2  1  2  3)
!      (3  2  1  2)
!      (4  3  2  1)
!
!  b =
!      (1  1  1  1)
!
 
PROGRAM dts_example
  USE structpack
  IMPLICIT NONE
    INTEGER, PARAMETER :: n = 4, nb = 20
    REAL(KIND=r_kind), DIMENSION(n) :: t, x, b
    t(1) = 1_r_kind
    t(2) = 2_r_kind
    t(3) = 3_r_kind
    t(4) = 4_r_kind
    x = 0_r_kind
    b(1:4) = 1_r_kind

    PRINT*, "t:", t
    PRINT*, "b:", b 

    CALL sp_dts_sv(t, x, b, nb, .TRUE.)

    PRINT*, "x:", x
    PRINT*, "b:", b
END PROGRAM dts_example
