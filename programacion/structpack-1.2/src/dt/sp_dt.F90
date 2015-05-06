!> \date May, 2012
!> \author StructPack team
!> \brief Module for the solution of non-Symmetric Toeplitz Linear Systems

module sp_dt
    use omp_lib
    use sp_utils
    use sp_dst
    use sp_dctii
    use sp_conv
    
    IMPLICIT NONE

    ! precision
    !
    PUBLIC :: sp_dt_sv, sp_dt_nrm1, sp_dt_gemv

    PRIVATE
    
  
    type pivot_type
        INTEGER, DIMENSION(:), ALLOCATABLE :: v
    end type
  
    contains
    
    !> \brief Solves the linear system
    !>        T x = b,
    !> where T is a non-symmetric Toeplitz matrix.
    !> \param[in]         u         First colum of matrix T
    !> \param[in]         v         First row of matrix T
    !> \param[out]        x         Solution vector of Tx=b returned by this routine.
    !> \param[in,out]     b         Right hand side vector of Tx=b. Changed on output.
    !> \param[in]         nb        Block size.
    !> \param[in]         pivoting  Uses local diagonal pivoting if nonzero.
    !> \param[in]         ref       Number of iterative refinement steps.   
    !> \param[out]        INFO      = 0: successful exit
    !>                              = 1: first element of u and first element of v are different
    !>                              = 2: u and v have different length
    !>                              = 3: u and b have different length
    subroutine sp_dt_sv( u, v, x, b, nb,  pivoting, refs, INFO )
        INTEGER,                          INTENT(IN)    :: nb
        REAL(kind=r_kind),  DIMENSION(:), INTENT(IN)    :: u, v, b
        REAL(kind=r_kind),  DIMENSION(:), INTENT(INOUT) :: x
        LOGICAL,                          INTENT(IN)    :: pivoting
        INTEGER,                          INTENT(IN)    :: refs
        INTEGER,                          INTENT(INOUT) :: INFO

        REAL(kind=r_kind), DIMENSION(:),   ALLOCATABLE :: Omega, Lambda
        REAL(kind=r_kind), DIMENSION(:,:), ALLOCATABLE :: G, H, AUX
        type(block_type),  DIMENSION(:,:), ALLOCATABLE :: L, mU
        integer,           DIMENSION(:),   ALLOCATABLE :: P

        INTEGER          :: n, num_blocks, threads, ref

        n = size(u,1)

        IF ((size(u,1) /= n) .OR. (u(1)/=v(1)) .OR. (size(b,1) /= n)) THEN
            ! Error en los datos de entrada
            INFO = -1
            return
        END IF
        threads = omp_get_max_threads()
        call omp_set_num_threads(threads)
        !$OMP PARALLEL SECTIONS
        
        !$OMP SECTION
        num_blocks = n/nb
        IF (mod(n,nb).NE.0) num_blocks = num_blocks+1
        
        CALL ALLOCATE_BLOCKS(L,  n, num_blocks, nb)
        CALL ALLOCATE_BLOCKS(mU, n, num_blocks, nb)  
        allocate( P(n) )
        
        !$OMP SECTION
        x=b
        call DST(x)
        
        !$OMP SECTION
        ALLOCATE (G(n,4), H(n,4), Omega(n), Lambda(n), AUX(n,2))
        CALL gen(u, v, G, H, Omega, Lambda, threads)

        !$OMP END PARALLEL SECTIONS
       

        CALL luCB(nb, G, H, Omega, Lambda, pivoting, L, mU, P, AUX, threads)

       
        CALL DTRSV_Blocks( 'L', L,  x, num_blocks, nb, 'N', 'U', threads, P )
        CALL DTRSV_Blocks( 'L', mU, x, num_blocks, nb, 'T', 'N', threads )
        CALL IDCTII( x )
        
        !Iterative refinement
        do ref = 1, refs
          aux(:,1) = b
          call sp_dt_gemv( -1.0D+0, u, v, x, 1.0D+0, aux(:,1) )
          call DST( aux(:,1) )
          call DTRSV_Blocks( 'L', L,  aux(:,1), num_blocks, nb, 'N', 'U', threads, P )
          call DTRSV_Blocks( 'L', mU, aux(:,1), num_blocks, nb, 'T', 'N', threads )
          call IDCTII( aux(:,1) )
          x = x + aux(:,1)
        end do

        CALL DEALLOCATE_BLOCKS(L)
        CALL DEALLOCATE_BLOCKS(mU)  
        DEALLOCATE( P )
        DEALLOCATE(G)
        DEALLOCATE(H)
        DEALLOCATE(Omega)
        DEALLOCATE(Lambda)
        DEALLOCATE(AUX)

    END subroutine sp_dt_sv
        
            
    subroutine ALLOCATE_BLOCKS(LD, m, num_blocks, nb)
        type(block_type), dimension(:,:), allocatable, intent(inout) :: LD
        integer, intent(in):: num_blocks, m, nb
        integer :: i, j,modm

        allocate(LD(num_blocks,num_blocks))
        modm=mod(m,nb)
        do j = 1, num_blocks,1
            do i = j, num_blocks,1
                LD(i,j)%x= nb
                LD(i,j)%y= nb
                if (modm .ne. 0) then
                    if (i.eq.num_blocks) then
                        LD(i,j)%x= modm
                    end if
                    if (j.eq.num_blocks) then
                        LD(i,j)%y= modm
                    end if
                end if
            allocate(LD(i,j)%v(LD(i,j)%x,LD(i,j)%y))
            end do
        end do  
        
    end subroutine ALLOCATE_BLOCKS
    
    SUBROUTINE DEALLOCATE_BLOCKS( m )
        INTEGER :: I, J
        type(block_type),  DIMENSION(:,:),INTENT(INOUT)  :: m
        DO j=1, size(m,2)
            DO i=j, size(m,1)
                DEALLOCATE( m(i,j)%v )
            END DO
        END DO
    END SUBROUTINE DEALLOCATE_BLOCKS



SUBROUTINE luCB( nb, G, H, Omega, Lambda, pivoting , L, U, P, aux, threads )
    INTEGER,                           INTENT(IN) :: nb, threads
    LOGICAL,                           INTENT(IN) :: pivoting
    REAL(kind=r_kind), DIMENSION(:,:), INTENT(INOUT) :: G,H, aux
    REAL(kind=r_kind), DIMENSION(:),   INTENT(INOUT) :: Omega, Lambda
    type(block_type),  DIMENSION(:,:), INTENT(INOUT) :: L,U
    INTEGER,           DIMENSION(:),   INTENT(INOUT) :: P
    
    INTEGER :: n, i,k,o, q, diagonalSize, ib,kb

    n = size(G,1)
               
    IF ( (size(H,1) /= size(G,1)) .OR. (size(H,2) /= size(G,2)) ) THEN
        stop "luCB: Error en las dimensiones de entrada"
    END IF
    
    DO k=1, n, nb
        kb = 1+(k-1)/nb
        q = k+min( nb, n-k+1 )-1
        diagonalSize = q-k+1
        
        call luC( G(k:q,:),H(k:q,:),Omega(k:q),Lambda(k:q),pivoting, &
                  L( kb, kb )%v,U( kb, kb )%v, P( k:q ), aux, threads )
        

    !$OMP PARALLEL PRIVATE(ib, o)  num_threads(threads)
    
    !$OMP DO SCHEDULE(guided)
        DO i=k+nb, n, nb
            ib = 1+(i-1)/nb
            o = i+min( nb, n-i+1 )-1
            
            CALL updateL(U( kb, kb )%v, G(k:q,:),G(i:o,:),H(k:q,:), &
                 Omega(i:o),Lambda(k:q), L( ib, kb )%v, 1)
        END DO
    !$OMP END DO NOWAIT

    !$OMP DO SCHEDULE(guided)
        DO i=k+nb, n, nb
            ib = 1+(i-1)/nb
            o = i+min( nb, n-i+1 )-1
            
            CALL updateU(U( kb, kb )%v,H(k:q,:),H(i:o,:),G(k:q,:), &
                Lambda(i:o),Omega(k:q), U( ib, kb )%v, 1)
        END DO
    !$OMP END DO
    !$OMP END PARALLEL

     
     END DO

END SUBROUTINE luCB


SUBROUTINE luC( G,H,Omega,Lambda,pivoting,L,U,P, workspace, threads )

    LOGICAL,                           INTENT(IN)    :: pivoting
    REAL(kind=r_kind), DIMENSION(:),   INTENT(IN)    :: Lambda
    REAL(kind=r_kind), DIMENSION(:,:), INTENT(INOUT) :: L, U, G, H, workspace
    REAL(kind=r_kind), DIMENSION(:),   INTENT(INOUT) :: Omega
    INTEGER,           DIMENSION(:),   INTENT(INOUT) :: P
    INTEGER,                           INTENT(IN)    :: threads
    
    INTEGER :: n,r, j, k, columSize
    REAL(kind=r_kind) :: d, aux

    n = size(G,1)
    r = size(G,2)
    
    P = (/(j,j=1,n)/)
    
    DO j=1,n,1
        columSize = n-j+1;

        call DGEMV( 'N', n-j+1, r, 1.0d+0, G(j:n,:), columSize, H(J,:), 1, 0.0d+0, workspace(:,2), 1)
        workspace(1:columSize,2) = workspace(1:columSize,2) / ( Omega(j:n)-Lambda(j) )

        IF ( pivoting ) THEN
            CALL posmaxabs( workspace(1:columSize,2), k )
        ELSE
            k = 1
        END IF
        
        d = workspace(k,2)        
        
        IF ( k /= 1 ) THEN
            workspace(k,2) = workspace(1,2)
            workspace(1,2)= d
            k = k + j -1 

            aux = P( j )      !
            P( j ) = P( k )   ! swap( piv(j), piv(k) )
            P( k ) = aux      !

            IF( j > 1 ) THEN
                workspace(1:j-1,1) = L(j, 1:j-1)
                L( j,1:j-1 ) = L(k,1:j-1)
                L( k,1:j-1 ) = workspace(1:j-1,1)
            END IF

            aux = Omega(j)
            Omega(j) = Omega(k)
            Omega(k) = aux
            
            workspace(1:4,1) = G(j,:)
            G(j,:) = G(k,:)
            G(k,:) = workspace(1:4,1)

        END IF


        U(j,j) = d
        L(j,j) = 1.0d+0
        
        IF( j < n ) THEN

            L( j+1:n, j) = workspace( 2: columSize,2)  / d
            
            call DGEMV( 'N', n-j, r, 1.0d+0, H(j+1:n,:), n-j, G(j,:), 1, 0.0d+0, U( j+1:n, j ), 1)                  
            U( j+1:n, j ) = U( j+1:n, j ) / ( ( Omega(j)-Lambda(j+1:n)) )
            
            G(j+1:n,1 ) = G(j+1:n,1 ) - (L(j+1:n,j) * G(j,1))
            G(j+1:n,2 ) = G(j+1:n,2 ) - (L(j+1:n,j) * G(j,2))
            G(j+1:n,3 ) = G(j+1:n,3 ) - (L(j+1:n,j) * G(j,3))
            G(j+1:n,4 ) = G(j+1:n,4 ) - (L(j+1:n,j) * G(j,4))
            
            H( j+1:n, 1) = H(j+1:n,1) - ((U(j+1:n,j) * H(j,1)) / d)
            H( j+1:n, 2) = H(j+1:n,2) - ((U(j+1:n,j) * H(j,2)) / d)
            H( j+1:n, 3) = H(j+1:n,3) - ((U(j+1:n,j) * H(j,3)) / d)
            H( j+1:n, 4) = H(j+1:n,4) - ((U(j+1:n,j) * H(j,4)) / d)
        END IF

    END DO

END SUBROUTINE luC

    subroutine gen(u, v, G, H, Omega, Lambda, threads )
        REAL(kind=r_kind), DIMENSION(:),   INTENT(IN)    :: u, v
        REAL(kind=r_kind), DIMENSION(:,:), INTENT(INOUT) :: G, H
        REAL(kind=r_kind), DIMENSION(:),   INTENT(OUT)   :: Omega, Lambda
        INTEGER,                           INTENT(IN)    :: threads

        INTEGER                                          :: n, j
        n       = size(u,1)

                 G(:,1)  = (/-v(n),  v(n:2:-1)-v(n-1:1:-1)/)
                 G(:,2)  = (/1.0d+0, (0.0d+0,j=2,n)/)
                 G(:,3)  = (/u(2:n)-u(1:n-1), -u(n)/)
                 G(:,4)  = (/(0.0d+0,j=1,n-1), 1.0d+0/)
                 Omega   = 2.0d+0*DCOS(((/(j,j=1,n)/) * PI) / (n+1))
                 
                CALL DST( G(:,1) )
                CALL DST( G(:,2) )
                CALL DST( G(:,3) )
                CALL DST( G(:,4) )
                 H(:,1) = (/(0.0d+0,j=1,n-1), 0.5d+0/)
                 H(:,2) = (/-v(2:n)*0.5d+0, 0.0d+0/)
                 H(:,3) = (/0.5d+0, (0.0d+0,j=2,n)/)
                 H(:,4) = (/0.0d+0, -u(n:2:-1)*0.5d+0 /)
                 Lambda  = 2.0d+0*DCOS(((/(j,j=0,n-1)/)* PI) / n)
                 
                CALL DCTII( H(:,1))
                CALL DCTII( H(:,2))
                CALL DCTII( H(:,3))
                CALL DCTII( H(:,4))
        H(1,1) = H(1,1) - 1.0d+0
        H(1,2) = H(1,2) - sum( -v(2:n))
        H(1,3) = H(1,3) - 1.0d+0
        H(1,4) = H(1,4) - sum(-u(2:n))        
            
#ifdef EXTRAE
    call Extrae_user_function  (0)
#endif        

    END subroutine gen


SUBROUTINE updateL( d, dG, G, dH, Omega, Lambda, L, threads) 
    REAL(kind=r_kind), DIMENSION(:),   INTENT(IN)    :: Omega, Lambda
    REAL(kind=r_kind), DIMENSION(:,:), INTENT(IN)    :: d, dG, dH  !(:,4)
    INTEGER,                           INTENT(IN)    :: threads
    REAL(kind=r_kind), DIMENSION(:,:), INTENT(INOUT) :: G !(:,4)
    REAL(kind=r_kind), DIMENSION(:,:), INTENT(OUT)   :: L  
    INTEGER                                         :: n, m, j

    n = size(G,1)
    m = size(dH,1)
    L = 0.0d+0
    IF ( (size(dH,2) /= size(G,2) ) .OR. ( size(Omega,1) /= size(G,1)) &
        .OR. (size(Lambda,1) /= size(dH,1)) .OR. ( size(d,1) /= size(dG,1)) ) THEN
        stop "updateL: Error en las dimensiones de entrada"
    END IF
        
    DO j=1,m,1
        L(:,j) = MATMUL(G, dH(j,:)) / (d(j,j) * ( Omega - Lambda(j))) 
       
        G(:,1) = G(:,1) - L(:,j) * dG(j,1)
        G(:,2) = G(:,2) - L(:,j) * dG(j,2)
        G(:,3) = G(:,3) - L(:,j) * dG(j,3)
        G(:,4) = G(:,4) - L(:,j) * dG(j,4)

    END DO

END SUBROUTINE updateL


SUBROUTINE updateU( d, dG, G, dH, Omega, Lambda, U, threads) 
    REAL(kind=r_kind), DIMENSION(:),   INTENT(IN)    :: Omega, Lambda
    REAL(kind=r_kind), DIMENSION(:,:), INTENT(IN)    :: d, dG,dH  !(:,4)
    INTEGER,                           INTENT(IN)    :: threads
    
    REAL(kind=r_kind), DIMENSION(:,:), INTENT(INOUT) :: G
    REAL(kind=r_kind), DIMENSION(:,:), INTENT(OUT)   :: U
    
    INTEGER :: n, m, j

    n = size(G,  1 )
    m = size(dH, 1 )
    
    DO j=1, m, 1

        U(:,j) = MATMUL( G, dH(j,:) ) / (-( Omega - Lambda(j) ))
        G(:,1) = G(:,1) - U(:,j) * (dG(j,1) / d(j,j))
        G(:,2) = G(:,2) - U(:,j) * (dG(j,2) / d(j,j))
        G(:,3) = G(:,3) - U(:,j) * (dG(j,3) / d(j,j))
        G(:,4) = G(:,4) - U(:,j) * (dG(j,4) / d(j,j))

    END DO
    
END SUBROUTINE updateU

SUBROUTINE posmaxabs( v, p )

    REAL(kind=r_kind), DIMENSION(:), INTENT(in) :: v
    INTEGER,                         INTENT(out):: p
    
    INTEGER           :: i,n
    REAL(kind=r_kind) :: maxvalue,tmp
    maxvalue = DABS(v(1))
    p = 1
    n = size(v,1)
    DO i=2, n
        tmp = DABS(v(i))
        if( maxvalue .le. tmp ) then
            maxvalue = tmp
            p = i
        END IF
    END DO        
END SUBROUTINE posmaxabs


    !> \brief Computes the 1-norm of a non-symmetric Toeplitz matrix T.
    !> \param[in]         u         First colum of matrix T
    !> \param[in]         v         First row of matrix T
    !> \return                      1-norm of T.
    real(kind=r_kind) function sp_dt_nrm1( u, v )
        real(kind=r_kind), DIMENSION(:), INTENT(in) :: u,v
        integer :: n, i
        real(kind=r_kind) :: tmp, maximum 

        n = size(u)
        
        maximum = SUM( DABS(u) )
        tmp = maximum
        do i = n, 2, -1
            tmp = tmp - DABS(u(i)) + DABS(v( n-i+1 ))
            if ( tmp .gt. maximum ) maximum = tmp
        end do
      
        sp_dt_nrm1 = maximum

    end function sp_dt_nrm1


  !> \brief Performs the operation
  !>        b = b + alpha * T * x,
  !> where T is a non-symmetric Toeplitz matrix.
  !> \param[in]      alpha  Scalar alpha.
  !> \param[in]      u      First column of matrix T.
  !> \param[in]      v      First row of matrix T.
  !> \param[in]      x      Vector x.
  !> \param[in,out]  b      Vector b. Changed on output.
    subroutine sp_dt_gemv(alpha, u, v, x, beta, b)
        real(kind=r_kind), intent(in) :: alpha, beta
        real(kind=r_kind), dimension(:), intent(in) :: u, v, x
        real(kind=r_kind), dimension(:), intent(inout) :: b

        !integer :: n, i
       ! real(kind=r_kind), dimension(2*size(u)-1) :: a

        !n = size(u)
       ! a = (/v(n:2:-1),u/)
        
       ! do i = 1, n
       !     b = b + alpha * a(n-i+1:2*n-i) * x
       ! end do

        call conv( alpha, u, x, beta, b, v )
    end subroutine sp_dt_gemv


end module sp_dt



