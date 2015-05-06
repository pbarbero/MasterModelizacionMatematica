!> \date May, 2012
!> \author StructPack team
!> \brief Module for the solution of Symmetric block Toeplitz Linear Systems

module sp_dtspg
    use omp_lib
    use sp_utils
    
        PUBLIC :: sp_dtspg_sv, sp_dtspg_rg, sp_dtspg_nrm1,sp_dtspg_gemv

        PRIVATE

        contains
               
        !> \brief Build a random symmetric block toeplitz T
        !> \param[in]   nBlocks     Number of rows per dimension of matrix T
        !> \param[in]   bSize       Size of each block of matrix T
        !> \param[out]  T           Symmetric block toeplitz T
        subroutine sp_dtspg_rg( nBlocks, bSize, seed, T )
            INTEGER, INTENT( in ) :: nBlocks, bSize, seed
            REAL( kind=r_kind ), DIMENSION(:,:), INTENT( inout ) :: T

            INTEGER :: i,j, mseed, idist
            INTEGER,             DIMENSION(4)                :: iseed
            REAL( kind=r_kind ), DIMENSION(:,:), ALLOCATABLE :: H, BLOCK
            idist = 2
            mseed = mod( seed, 4095 )
            iseed(1) = mseed
            iseed(2) = mseed
            iseed(3) = mseed
            iseed(4) = mseed
            if ( mod(mseed,2) == 0) iseed(4) = mseed + 1

            allocate( H( (2*nBlocks-1)*bSize,nBlocks*bSize) )
            allocate( BLOCK( bSize, bSize ) )
            
            H = ZERO

            do i=1,nBlocks
                call DLARNV( idist, iseed, bSize*bSize, BLOCK )

                do j=1,nBlocks
                    H((i+j-2)*bSize+1:(i+j-2)*bSize+bSize, (j-1)*bSize+1:(j-1)*bSize+bSize) = BLOCK
                end do
            end do
    
                        
            CALL dgemm( 'T', 'N', bSize, nBlocks*bSize, (2*nBlocks-1)*bSize, &
                        ONE, H, (2*nBlocks-1)*bSize, H, (2*nBlocks-1)*bSize, ZERO, T, nBlocks*bSize )
            
            do j=1, nBlocks
                do i=2, nBlocks
                    if( i .gt. j ) then
                        T((i-1)*bSize+1:i*bSize, (j-1)*bSize+1:j*bSize) = &
                            TRANSPOSE(T(1:bSize, ((abs(i-j))*bSize)+1:(abs(i-j)+1)*bSize ))
                    else
                        T((i-1)*bSize+1:i*bSize, (j-1)*bSize+1:j*bSize) = &
                            T(1:bSize, ((abs(i-j))*bSize)+1:(abs(i-j)+1)*bSize )    
                    end if
                        
                end do
            end do
            
            DEALLOCATE( H, BLOCK )
        end subroutine  
        
    !> \brief Solves the linear system
    !>        T X = B,
    !> where T is a symmetric block Toeplitz matrix.
    !> \param[in]         T         First row of blocks of matrix T
    !> \param[in]         b         Right hand side Matrix of TX=B
    !> \param[out]        x         Solution matrix of TX=B returned by this routine.
        subroutine sp_dtspg_sv( T, b, x )
            REAL( kind=r_kind ), DIMENSION(:,:), INTENT(  in  ) :: T
            REAL( kind=r_kind ), DIMENSION(:,:), INTENT(  in  ) :: b
            REAL( kind=r_kind ), DIMENSION(:,:), INTENT( inout ):: x
            INTEGER :: n, bSize, nBlock, i,j,OMP_GET_MAX_THREADS,threads, nsystems
            type( block_type ), DIMENSION(:,:), ALLOCATABLE :: C
            REAL( kind=r_kind ), DIMENSION(:),   ALLOCATABLE :: y

            if( size(x,2) .ne. size(b,2) ) stop '2nd and 3rd params must have the same dimension lenght'

            n      = size( T, 2 ) ! Tamaño de M
            bSize  = size( T, 1 ) ! Tamaño de bloque
            nBlock = n/bSize      ! Numero de bloques
            nsystems = size(x,2)

            allocate( C(nBlock,nBlock), y(n) )
            do i = 1, nBlock
                do j=i, nBlock
                    allocate( C(i,j)%v(bSize, bSize) ) 
                end do
            end do
            
            call sp_dtspg_ch( T, C )
            threads = OMP_GET_MAX_THREADS()
            
            if( nsystems .eq. 1 ) then
                x(:,1) = b(:,1)
                CALL DTRSV_Blocks( 'U', C, x(:,1), nBlock, bSize, 'T', 'N', threads )
                CALL DTRSV_Blocks( 'U', C, x(:,1), nBlock, bSize, 'N', 'N', threads )
            else
                x = b
                CALL DTRSM_Blocks( C, x, nBlock, bSize, 'T', 'N' )
                CALL DTRSM_Blocks( C, x, nBlock, bSize, 'N', 'N' )
            end if


            do i = 1, nBlock
               do j=i,nBlock  
                deallocate( C(i,j)%v ) 
                end do
            end do
            deallocate( C, y )
            
        end subroutine sp_dtspg_sv
        

    


        
        subroutine sp_dtspg_ch( T, C )
            REAL( kind=r_kind ), DIMENSION(:,:), INTENT(  in  ) :: T
            type( block_type ),  DIMENSION(:,:),   INTENT( inout ):: C
                        
            INTEGER :: n, bSize, nBlock, INFO, i,j
            type( block_type ),  DIMENSION(:),   ALLOCATABLE :: G
            REAL( kind=r_kind ), DIMENSION(:,:), ALLOCATABLE :: H, M, Xin
            
            
            n      = size( T, 2 ) ! Tamaño de M
            bSize  = size( T, 1 ) ! Tamaño de bloque
            nBlock = n/bSize      ! Numero de bloques


            allocate ( H(2*bSize,2*bSize), M(2*bSize,n), Xin(2*bSize,n) )
            
            allocate( G(nBlock) )
            do i = 1, nBlock
                allocate( G(i)%v(2*bSize, bSize) ) 
            end do
           
            C( 1, 1 )%v = T(1:bSize,1:bSize)
            
            CALL DPOTRF( 'U', bSize, C( 1, 1 )%v, bSize, INFO )
            
            
            G(1)%v(1:bSize,1:bSize) = C( 1, 1 )%v
            
            CALL DTRTRI( 'U', 'N', bSize, G(1)%v(1:bSize,1:bSize), bSize, INFO )! inv

            !$OMP PARALLEL DO
            do i=2,nBlock
                C(1, i)%v = T(1:bSize, (i-1)*bSize+1:i*bSize)
                CALL DTRMM( 'L','U','T','N',bSize,bSize,ONE,G(1)%v(1:bSize,1:bSize),bSize,C(1,i)%v, bSize)
            end do
            !$OMP END PARALLEL DO
            
                        
            if( nBlock .eq. 1 ) return
            
            
            !$OMP PARALLEL DO
            do i=1,nBlock-1
                G(i)%v(1:bSize,         1:bSize ) = C( 1, i   )%v
                G(i)%v(bSize+1:2*bSize, 1:bSize ) = C( 1, i+1 )%v
            end do
            !$OMP END PARALLEL DO
            
            CALL normablock( G, n-bSize, G,M, H, Xin ) !%v(:,1:n-bSize)
            
            
            !$OMP PARALLEL DO
            do i=1,nBlock-1
                C( 2,i+1 )%v = G(i)%v(1:bSize,1:bSize)
            end do
            !$OMP END PARALLEL DO
            
            do i = 2,nBlock-1
            
            
                do j=1, nBlock-2
                    G(j)%v(bSize+1:2*bSize,:) = G(j+1)%v(bSize+1:2*bSize,:)
                end do
                
                CALL normablock(G, n-bSize*i, G,M, H, Xin)
                               
                !$OMP PARALLEL DO
                do j=i, nBlock-1
                
                    C(i+1,j+1)%v = G(j-i+1)%v(1:bSize,1:bSize)
                end do
                !$OMP END PARALLEL DO
            end do
        
            
            
            do i = 1, nBlock
                deallocate( G(i)%v ) 
            end do
            deallocate( G )
            
        end subroutine

            
        subroutine normablock( G, n, X, M, H, Xin)
            type( block_type ), DIMENSION(:), INTENT( inout ) :: G
            type( block_type ), DIMENSION(:), INTENT( inout ) :: X
            INTEGER, INTENT(IN) :: n
            
            INTEGER :: bSize, i, j, INFO
            REAL( kind=r_kind ), DIMENSION(:,:), INTENT( inout ) :: H, M, Xin
            bSize = size( G(1)%v(:,:), 2 )  ! Tamaño de bloque
            

            H = ZERO
            DO j=1, 2*bSize
                H(j,j) = ONE
            END DO

            ! ........  [Q,R] = qr(G(bSize+1:2*bSize,1:bSize)) ..........
            M(bSize+1:2*bSize, 1:bSize) = G(1)%v(bSize+1:2*bSize,:)!G(2)%v(:,1:bSize)
            CALL DGEQRF( bSize, bSize, M(bSize+1:2*bSize, 1:bSize), bSize, Xin(1:bSize,1), &
                        Xin(1:bSize,2), bsize, INFO )
            
            H( bSize+1:2*bSize, bSize+1:2*bSize ) = M(bSize+1:2*bSize, 1:bSize)
            CALL DORGQR( bSize, bSize, bSize, H( bSize+1:2*bSize, bSize+1:2*bSize ), bSize, &
                        Xin(1:bSize,1), Xin(1:bSize,2), bSize, INFO )
            !  ................  Fin QR ................................
            
            M(1:bSize, 1:bSize) = G(1)%v(1:bSize,:)!G(1)%v(:,1:bSize)
            
            DO j = 1,bSize
                DO i = 1,j 
                    CALL hyperblock( M(j,j:bSize), M(i+bSize,j:bSize), H(:,j), H(:,i+bSize), Xin(:,1) )
                END DO
            END DO

            !$OMP PARALLEL DO
            do j=2,n/bSize
                CALL DGEMM( 'T', 'N', 2*bSize, bSize,2*bSize, &
                            ONE, H, 2*bSize, G(j)%v, 2*bSize, &
                            ZERO, Xin(:,(j-1)*bSize+1:j*bSize), 2*bSize )
                X(j)%v = Xin(:,(j-1)*bSize+1:j*bSize)      
            end do
            !$OMP END PARALLEL DO
            X(1)%v = M(1:2*bSize,1:bSize)
            
        end subroutine normablock


        subroutine hyperblock( u, v, h1, h2, tmp1 )
        
            REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: u, v, h1, h2
            REAL( kind=r_kind ), DIMENSION(:), INTENT( inout ) :: tmp1
          
            INTEGER :: nb, n
            REAL( kind=r_kind ) :: x,y, modulo, u1, v1
            nb= size( u, 1 )  ! Tamaño de los vectores
            if( size( v, 1 ) .ne. nb )   stop 'hyperblock: Dimension de los vectores distinta'
            
            n = size( h1, 1 )  ! Tamaño de los vectores de H
            if( size( h2, 1 ) .ne. n ) stop 'hyperblock: Dimension de los vectores distinta'
              
            x = u(1)**2
            y = v(1)**2
            if( x < y ) stop 'No definida positiva'

            modulo = dsqrt( x - y );
            u1 = u(1) / modulo
            v1 = v(1) / modulo
            
            
            tmp1(1:size(u)) = ( u1*u  - v1*v )
            v    = ( u1*v  - v1*u ) 
            u    = tmp1(1:size(u))  
            
            tmp1 = ( u1*h1 - v1*h2) 
            h2   = ( u1*h2 - v1*h1) 
            h1   = tmp1
        end subroutine hyperblock

        
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!      Non optimal functions      !!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !> \brief Computes the 1-norm of a symmetric block Toeplitz matrix T.
        !> \param[in]         T         First row block of matrix T
        !> \return                      1-norm of T.
        real(kind=r_kind) function sp_dtspg_nrm1( T )
            REAL( kind=r_kind), DIMENSION(:,:), INTENT( in ) :: T
            REAL( kind=r_kind), DIMENSION(:,:), allocatable :: Tin
            integer n, nBlocks, bSize, i, j
            DOUBLE PRECISION ::DUMMY
            
            n = size(T,2)
            bSize = size(T,1)
            nBlocks = n/bSize
            
            allocate(Tin(n,n))
            Tin(1:bSize, :) = T
            
            do j=1, nBlocks
                do i=2, nBlocks
                    if( i .gt. j ) then
                        Tin((i-1)*bSize+1:i*bSize, (j-1)*bSize+1:j*bSize) = &
                            TRANSPOSE(Tin(1:bSize, ((abs(i-j))*bSize)+1:(abs(i-j)+1)*bSize ))
                    else
                        Tin((i-1)*bSize+1:i*bSize, (j-1)*bSize+1:j*bSize) = &
                            Tin(1:bSize, ((abs(i-j))*bSize)+1:(abs(i-j)+1)*bSize )    
                    end if
                        
                end do
            end do

            ! sp_dtspg_nrm1 = DLANGE( '1', n, n, Tin, n, DUMMY )
           
           sp_dtspg_nrm1 = sum(dabs(Tin(:,1)))
           do i=2, n
                sp_dtspg_nrm1 = max(sum(abs(Tin(:,i))),sp_dtspg_nrm1)
           end do
            
            deallocate( Tin )
            
        end function
        
        
        !> \brief Performs the operation
        !>        b = beta* b + alpha * T * x,
        !> where T is a symmetric block Toeplitz matrix.
        !> \param[in]      alpha  Scalar alpha.
        !> \param[in]      T      First row block of T
        !> \param[in]      x      Vector x.
        !> \param[in]      beta   Scalar beta
        !> \param[in,out]  b      Vector b. Changed on output.
        subroutine sp_dtspg_gemv( alfa, T, x, beta, b )
            REAL( kind=r_kind ), intent( in )                 :: alfa, beta
            REAL( kind=r_kind ), DIMENSION(:,:), intent( in ) :: T
            REAL( kind=r_kind ), DIMENSION(:), intent( in )   :: x
            REAL( kind=r_kind ), DIMENSION(:), intent( inout ):: b
            
            
            REAL( kind=r_kind), DIMENSION(:,:), allocatable :: Tin
            integer :: i,j, nBlocks, bSize, n
            
           
            n = size(T,2)
            bSize = size(T,1)
            nBlocks = n/bSize
            
            allocate(Tin(n,n))
            Tin(1:bSize, :) = T
            
            do j=1, nBlocks
                do i=2, nBlocks
                    if( i .gt. j ) then
                        Tin((i-1)*bSize+1:i*bSize, (j-1)*bSize+1:j*bSize) = &
                            TRANSPOSE(Tin(1:bSize, ((abs(i-j))*bSize)+1:(abs(i-j)+1)*bSize ))
                    else
                        Tin((i-1)*bSize+1:i*bSize, (j-1)*bSize+1:j*bSize) = &
                            Tin(1:bSize, ((abs(i-j))*bSize)+1:(abs(i-j)+1)*bSize )    
                    end if
                        
                end do
            end do

            b = alfa*MATMUL(Tin,x) + beta*b           
            
        end subroutine

end module
