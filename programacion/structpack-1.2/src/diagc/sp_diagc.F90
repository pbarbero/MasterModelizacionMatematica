!!> \file sp_diagc.F90
!!>
!!> This module file contains a subroutine that computes the diagonal of a
!!> Cauchy-like matrix obtained by means of certain product
!!>
!!> \author P. Alonso
!!>
!!> \date Dec-2010

module sp_diagc
  use sp_utils
  use sp_dst
  
  implicit none
  
  public :: diagc
  private
  
contains
 
    !> \brief Computes the diagonal of a Cauchy-like matrix
    !> 
    !>   Computes the diagonal of a Cauchy-like matrix obtained
    !>   by means of the product 
    !>
    !>         C = S * T * S
    !>
    !>   where S is the DST ans T is an order n Toeplitz matrix defined 
    !>   by Column and Row_. Row_ is an optional argument not referenced
    !>   if T is symmetric.
    !> \param[in] Column
    !> \param[out] Diagonal
    !> \param[in] Row_
  subroutine diagc( Column, Diagonal, Row_ )
    !
    !   diagc computes the diagonal of a Cauchy-like matrix obtained
    !   by means of the product 
    !
    !         C = S * T * S
    !
    !   where S is the DST ans T is an order n Toeplitz matrix defined 
    !   by Column and Row_. Row_ is an optional argument not referenced
    !   if T is symmetric.
    !
    real(kind=r_kind), dimension(:), intent(in) :: Column 
    real(kind=r_kind), dimension(:), intent(out) :: Diagonal 
    real(kind=r_kind), dimension(:), intent(in), optional :: Row_ 

    real(kind=r_kind), dimension(1) :: Row

    if( present(Row_) ) then 
      call diagcx( 'N', size(Column,1), Column, Row_, Diagonal )
    else
      call diagcx( 'Y', size(Column,1), Column, Row , Diagonal )
    end if

  end subroutine diagc

  subroutine diagcx( sim, n, U, V, D )
    !
    !     diagcx computes the diagonal of a Cauchy-like matrix obtained
    !     by means of the product 
    !
    !           C = S * T * S
    !
    !     where S is the DST ans T is an order n Toeplitz matrix.
    !     
    !     If sim='Y' or sim='y' is assumed that T is symmetric.
    !     If sim='N' or sim='n' is assumed that T is non-symmetric.
    !     U and V are the column and row of T, respectively.
    !     If the T is symmetric, V is not referenced.
    !     On output, D have the diagonal entries of C.
    !
    !     D MUST be of size n+1
    ! 
    !     Pedro Alonso
    !     September 2005
    ! 
    character, intent(in) :: sim
    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: U
    double precision, dimension(:), intent(in) :: V
    double precision, dimension(n), intent(out) :: D

    double precision, allocatable, dimension(:) :: W

    if( sim == 'N' .or. sim == 'n' ) then 
       !allocate( W(n+1) )
       allocate( W(n) )
    end if

    call ddiagcx( sim, n, U, V, D, W )

    if( sim == 'N' .or. sim == 'n' ) then 
       deallocate( W )
    end if

    return 

  end subroutine diagcx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

  SUBROUTINE DDIAGCX( SIM, N, U, V, D, W )

    !
    !  -- Pedro Alonso
    !     february 2005
    !     modified september 2005
    !
    !     .. Scalar Arguments ..
    CHARACTER, intent(in) :: SIM
    INTEGER, intent(in)   :: N
    !     ..
    !     .. Array Arguments ..
    real(kind=r_kind), intent(in)  :: U(:), V(:)
    real(kind=r_kind), intent(out) :: D(:), W(:)
    !
    !  Purpose
    !  =======
    !
    !  DDIAGCX is the driver routine for computing the diagonal of
    !  a Cauchy-like matrix C in which is converted a symmetric as non-symmetric
    !  Toeplitz matrix T with the DST S, C=S*T*S.
    !
    !  Arguments
    !  =========
    !
    !  SIM     (INPUT) CHARACTER*1
    !          SIM = 'Y' or SIM='N' indicates if the Toeplitz matrix is
    !          symmetric or not.
    !
    !  N       (input) INTEGER
    !          The number of elements of vectors T and C. N >= 4 because 
    !          DDIAGC can not be called with a lower dimension.
    !
    !  U       (input) DOUBLE PRECISION array, dimension N
    !          First column of the Toeplitz matrix.
    !
    !  V       (input) DOUBLE PRECISION array, dimension N
    !          First row of the Toeplitz matrix.
    !          If SIM = 'Y' or SIM = 'y' this array is not referenced.
    !
    !  D       (output) DOUBLE PRECISION array, dimension N
    !          Diagonal entries of a Cauchy-like matrix in which the 
    !          Toeplitz matrix has been transformed.
    !
    !  W       (workspace) DOUBLE PRECISION arrays, dimension N
    !          If SIM = 'Y' or SIM = 'y' this array is not referenced.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    real(kind=r_kind)   ONE         , TWO         
    PARAMETER        ( ONE = 1.0_r_kind , TWO = 2.0_r_kind  )
    !
    INTEGER            INFO
    real(kind=r_kind) :: AUX(n)
    !     ..
    !     .. External Functions ..
    LOGICAL            LSAME
    EXTERNAL           LSAME
    !     ..
    !     .. External Subroutines ..
    EXTERNAL           XERBLA, DSCAL, DAXPY
    !     
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF( .NOT.(LSAME(SIM,'Y') .OR. &
              LSAME(SIM,'y') .OR.& 
              LSAME(SIM,'N') .OR. &
              LSAME(SIM,'n')) ) THEN
       INFO = 1
    ELSE IF( N.LT.4 ) THEN
       INFO = 2
    END IF

    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DDIAGCX', -1 )
       RETURN
    END IF

    !     Return if possible
    IF( N.EQ.0 ) RETURN

    AUX(1)  = U(1) / TWO
    AUX(2:) = U(2:)
    !     Computing the diagonal of the lower triangular part of T
    CALL DDIAGC( N, AUX, D )

    IF( LSAME(SIM,'Y') .OR. LSAME(SIM,'y') ) THEN
       CALL DSCAL( N, TWO, D, 1 )
    ELSE
       !       Saving of the first entry of U
       AUX(1)  = V(1) / TWO
       AUX(2:) = V(2:)
       !       Computing the diagonal of the upper triangular part of T
       CALL DDIAGC( N, AUX, W )
       CALL DAXPY( N, ONE, W, 1, D, 1 )
     END IF

    RETURN

  END SUBROUTINE DDIAGCX


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  SUBROUTINE DDIAGC( N, T, D )
    !
    !  -- Pedro Alonso
    !     february 2005
    !
    !     .. Scalar Arguments ..
    INTEGER, intent(in) ::  N
    !     ..
    !     .. Array Arguments ..
    real(kind=r_kind), intent(in)  ::  T(:)
    real(kind=r_kind), intent(out) ::  D(:)
    !
    !  Purpose
    !  =======
    !
    !  DDIAGC computes the diagonal entries of a Cauchy-like matrix in which has
    !  been transformed matrix a lower triangular N by N Toeplitz 
    !  matrix by means of the discrete sine transform (DST).
    !
    !  Arguments
    !  =========
    !
    !  N       (input) INTEGER
    !          The number of elements of vectors T and C.  N >= 4 because
    !          the matrix G can not be formed with lower dimension in 
    !          this version.
    !
    !  T       (input) DOUBLE PRECISION array, dimension N
    !          First column of the Toeplitz matrix.
    !
    !  D       (output) DOUBLE PRECISION array, dimension N
    !          Diagonal entries of a Cauchy-like matrix in which 
    !          the Toeplitz matrix has been transformed.
    !
    !  =====================================================================
    !
    !     .. Parameters ..
    real(kind=r_kind)   ONE         , TWO         , FOUR
    PARAMETER        ( ONE = 1.0_r_kind , TWO = 2.0_r_kind , FOUR= 4.0_r_kind  )
    !
    !     .. Local Scalars ..
    INTEGER            INFO, I
    real(kind=r_kind)   AUX, BETAIP, BETAP, BETAII, BETAI, HP, HI
    !
    !     .. External Functions ..
    LOGICAL            LSAME
    EXTERNAL           LSAME
    !     ..
    !     .. External Subroutines ..
    EXTERNAL           XERBLA
    !     
    !     ..
    !     .. Intrinsic ..
    INTRINSIC          SIN, ATAN, MOD
    !     ..
    !     .. Executable Statements ..
    !
    !     Test the input parameters.
    !
    INFO = 0
    IF( N.LT.4 ) THEN
       INFO = 1
    END IF
    !
    IF( INFO.NE.0 ) THEN
       CALL XERBLA( 'DDIAGC', -1 )
       RETURN
    END IF
    !
    AUX = TWO*N+TWO
    !     % impares
    BETAIP  = N*T(1);
    BETAP   = BETAIP;
    HP      = -T(1);
    D(1) = 3*BETAP;
    D(3) = -BETAP;
    !     % impares
    BETAII  = (N-1)*T(2);
    BETAI   = BETAII;
    HI      = -T(2);
    D(2) = TWO*BETAI;
    D(4) = -BETAI;

    DO I=3,N-2,2
       !       % pares
       BETAIP    = (N-I+1)*T(I) + TWO*HP;
       BETAP     = BETAP + BETAIP;
       HP        = -T(I) + HP;
       D(I-2) = (D(I-2) - BETAP)/AUX
       D(I)   = D(I)   + TWO*BETAP;
       D(I+2) = -BETAP;
       !       % impares
       IF( (I+3).LE.N) THEN
          BETAII    = (N-I)*T(I+1) + TWO*HI;
          BETAI     = BETAI + BETAII;
          HI        = -T(I+1) + HI;
          D(I-1) = (D(I-1) - BETAI)/AUX
          D(I+1) = D(I+1) + TWO*BETAI;
          D(I+3) = -BETAI;
       END IF
    END DO

    IF( MOD(N-1,2).EQ.0 ) THEN
       BETAII    = TWO*T(N-1) + TWO*HI;
       BETAI     = BETAI + BETAII;
       D(N-3) = (D(N-3) - BETAI)/AUX
       D(N-1) = (D(N-1) + TWO*BETAI)/AUX
    ELSE
       BETAIP    = TWO*T(N-1) + TWO*HP;
       BETAP     = BETAP + BETAIP;
       D(N-3) = (D(N-3) - BETAP)/AUX
       D(N-1) = (D(N-1) + TWO*BETAP)/AUX
    END IF

    IF( MOD(N,2).EQ.0 ) THEN
       BETAII    = T(N) + TWO*HI;
       BETAI     = BETAI + BETAII;
       D(N-2) = (D(N-2) - BETAI)/AUX
       D(N)   = (D(N)   + 3*BETAI)/AUX
    ELSE
       BETAIP    = T(N) + TWO*HP;
       BETAP     = BETAP + BETAIP;
       D(N-2) = (D(N-2) - BETAP)/AUX
       D(N)   = (D(N)   + 3*BETAP)/AUX
    END IF

    CALL dst( D )

    AUX = FOUR*ATAN(ONE) / DBLE(N+1)
    DO I = 1, N
       D( I ) = D( I ) / ( TWO*SIN( I*AUX ) )
    END DO

    RETURN

  END SUBROUTINE DDIAGC
  
end module sp_diagc
