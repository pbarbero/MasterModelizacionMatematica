!> \date April, 2013
!> \author StructPack team
!> \brief Module for the solution of Symmetric Toeplitz Linear Systems
module sp_dts
  use iso_c_binding
  use omp_lib
  use sp_utils
  use sp_dst
  use sp_diagc
  implicit none
  
  public :: sp_dts_sv, sp_dts_nrm1, sp_dts_gemv

  private
  
contains

  !> \brief Solves the linear system
  !>        T x = b,
  !> where T is a symmetric Toeplitz matrix.
  !> \param[in]      T    First column (row) of matrix T.
  !> \param[out]     x    Solution vector of Tx=b returned by this routine.
  !> \param[in,out]  b    Right hand side vector of Tx=b. Changed on output.
  !> \param[in]      nb   Block size.
  !> \param[in]      piv  Uses local diagonal pivoting if piv=.TRUE..
  !> \param[in]      affi Uses local core affinity if affi=.TRUE..
  subroutine sp_dts_sv(T, x, b, nb, piv, affi )
    integer, intent(in) :: nb
    real(kind=r_kind), dimension(:), intent(in) :: T
    real(kind=r_kind), dimension(:), intent(out):: x
    real(kind=r_kind), dimension(:), intent(inout):: b
    logical, intent(in):: piv, affi

    real(kind=r_kind), dimension(:), allocatable :: g11, g12, g21, g22, bsol, diag, lambda
    integer, dimension(:), allocatable :: piv1, piv2
    integer, dimension(:), allocatable :: piv1_t, piv2_t

    integer :: num_blocks1, num_blocks2
    type(block_type), dimension(:,:), allocatable :: LD1, LD2
    integer :: n, m1, m2, status, threads, max_threads

    n = size(T,1)
    m1 = size(T,1)/2+mod(size(T,1),2)
    m2 = size(T,1)/2

    num_blocks1 = m1 / nb
    if (mod(m1,nb) .ne. 0) num_blocks1 = num_blocks1+1
    num_blocks2 = m2 / nb
    if (mod(m2,nb) .ne. 0) num_blocks2 = num_blocks2+1

    CALL allocate_LD(LD1, num_blocks1, m1, nb)
    CALL allocate_LD(LD2, num_blocks2, m2, nb)

    allocate(g11(m1), g12(m1), g21(m2), g22(m2), bsol(n), diag(n), lambda(n), Stat=status)
    if ( piv ) then
      allocate( piv1(m1), piv2(m2), Stat=status)
      allocate( piv1_t(m1), piv2_t(m2), Stat=status)
    end if

    CALL cauchy_toep_generators(t, lambda(1:m1), lambda(m1+1:n), diag(1:m1), diag(m1+1:n), g11, g12, g21, g22)

    bsol = b
    CALL DST(bsol)
    bsol = (/ bsol(1:n:2), bsol(2:n:2) /)

    CALL OMP_SET_NESTED(.true.)

    max_threads = OMP_GET_MAX_THREADS()
    threads = 1
    if(max_threads.ge.2) threads = 2

    !$OMP PARALLEL num_threads(threads)
    if( affi .and. threads.ge.2 ) call affinity2( )
    !$OMP SECTIONS
    !$OMP SECTION
    CALL LDLt_Blocks(  g11, g12, lambda(1:m1), diag(1:m1), LD1, nb, piv1, affi )
    CALL DTRSV_Blocks( 'L', LD1, bsol(1:m1), num_blocks1, nb, 'N', 'U', threads/2, piv1 )
    CALL DTRSD_Blocks( LD1, bsol(1:m1), num_blocks1, nb )
    CALL DTRSV_Blocks( 'L', LD1, bsol(1:m1), num_blocks1, nb, 'T', 'U', threads/2, piv1 )
    !$OMP SECTION
    CALL LDLt_Blocks(  g21, g22, lambda(m1+1:n), diag(m1+1:n), LD2, nb, piv2, affi )
    CALL DTRSV_Blocks( 'L', LD2, bsol(m1+1:n), num_blocks2, nb, 'N', 'U', threads/2, piv2 )
    CALL DTRSD_Blocks( LD2, bsol(m1+1:n), num_blocks2, nb )
    CALL DTRSV_Blocks( 'L', LD2, bsol(m1+1:n), num_blocks2, nb, 'T', 'U', threads/2, piv2 )
    !$OMP END SECTIONS 
    !$OMP END PARALLEL

    x(1:n:2) = bsol(1:m1)
    x(2:n:2) = bsol(m1+1:n)

    CALL DST(x)

    x = x / (TWO*(n+1))

    deallocate(g11, g12, g21, g22, bsol, diag, lambda)
    if ( piv ) then
      deallocate(piv1,piv2)
      deallocate(piv1_t,piv2_t)
    end if
  
    CALL deallocate_LD(LD1, num_blocks1)
    CALL deallocate_LD(LD2, num_blocks2)

  end subroutine sp_dts_sv


  subroutine cauchy_toep_generators( t, lambda1, lambda2, diag1, diag2, g11, g12, g21, g22 )
    
    ! cauchy_toep_generators obtains the rank displacement arguments of the two Cauchy-like 
    ! linear systems yielding from a symmetric Toeplitz matrix
    !
    !   t: first column (row) of the toeplitz matrix
    !
    !   lambda1: diagonal of the displacement matrix of equation 1
    !
    !   lambda2: diagonal of the displacement matrix of equation 2
    !
    !   diag1: diagonal of the Cauchy-like matrix of equation 1
    !
    !   diag2: diagonal of the Cauchy-like matrix of equation 2
    !
    !   g11: First column of the generator of the displacement equation 1
    !
    !   g12: Second column of the generator of the displacement equation 1
    !
    !   g21: First column of the generator of the displacement equation 2
    !
    !   g22: Second column of the generator of the displacement equation 2
    !    
    real(kind=r_kind), dimension(:), intent(in)  :: t
    real(kind=r_kind), dimension(:), intent(out) :: lambda1, diag1
    real(kind=r_kind), dimension(:), intent(out) :: g11, g12
    real(kind=r_kind), dimension(:), intent(out) :: lambda2, diag2
    real(kind=r_kind), dimension(:), intent(out) :: g21, g22

    integer :: n, n1, n2, i
    real(kind=r_kind), dimension(size(t,1)) :: Poe, g1, g2, lambda, diag
    real(kind=r_kind) :: sqrtnp1, sqrtnp2
    real(kind=r_kind), parameter :: PI = FOUR * atan(ONE)

    n = size(t)
    n2 = n/2
    n1 = n2 + mod(n,2)

    sqrtnp1 = sqrt(TWO)
    sqrtnp2 = TWO * sqrtnp1

    !$OMP PARALLEL SECTIONS PRIVATE(i, Poe)

    !$OMP SECTION       
    ! First column of generator G of TS
    g1(1) = ZERO
    g1(2:n-1) = sqrtnp1 * T(3:n)
    g1(n) = ZERO
    call dst (g1)
    Poe = g1
    g1 = (/Poe(1:n:2),Poe(2:n:2)/)
    g11 = g1(1:n1)
    g21 = g1(n1+1:n)

    !$OMP SECTION       
    ! Second column of generator G of C directly
    g2(1:n/2) = (/ (sqrtnp2*sin( real(i)*PI / real(n+1)),i=1,n/2) /)
    g2(n:(n+1)/2+1:-1) = g2( 1:n/2 )
    if (mod(n,2) .eq. 1) g2( n/2+1 ) = sqrtnp2
    Poe = g2
    g2 = (/Poe(1:n:2),Poe(2:n:2)/)
    g12 = g2(1:n1)
    g22 = g2(n1+1:n)
    
    !$OMP SECTION       
    ! Computing vector lambda
    lambda(1:n/2) = (/ ( TWO * cos( real(i)*PI / real(n+1) ), i=1,n/2 ) /)
    lambda(n:(n+1)/2+1:-1) = -lambda( 1:n/2 )
    lambda = TWO * ( n + 1 ) * lambda
    if (mod(n,2) .eq. 1) lambda(n/2+1) = ZERO
    Poe = lambda
    lambda = (/Poe(1:n:2),Poe(2:n:2)/)
    lambda1 = lambda(1:n1)
    lambda2 = lambda(n1+1:n)

    !$OMP SECTION       
    ! Computing diagonal entries of C a priori
    call diagC( t, diag )
    Poe = diag
    diag = (/Poe(1:n:2),Poe(2:n:2)/)
    diag1 = diag(1:n1)
    diag2 = diag(n1+1:n)

    !$OMP END PARALLEL SECTIONS
  end subroutine cauchy_toep_generators


  subroutine LDLt_Blocks( g1, g2, lambda, diagC, L, nb, piv, affi ) 
    type(block_type), dimension(:,:), intent(inout) :: L
    real(kind=r_kind), dimension(:), intent(inout) :: diagC, g1, g2, lambda
    integer, intent(in) :: nb
    integer, dimension(:), allocatable, intent(inout) :: piv
    logical :: affi

    integer :: n, i, k, threads, bs, ik, fk, ii, fi
    logical :: pivoting
    n = size(L,1)
    bs= nb

    pivoting = allocated( piv )

    threads = OMP_GET_MAX_THREADS() / 2

    !$OMP PARALLEL PRIVATE(ik,fk) num_threads(threads)
    if( affi .and. threads.gt.1 ) call affinity( )
    do k=1, n     
      ik= (k-1)*bs+1
      fk= k*bs
      if (fk.gt.size(g1,1)) fk= size(g1,1)
      !$OMP SINGLE
      if( pivoting ) then
        CALL LDLt_Triang( g1(ik:fk), g2(ik:fk), lambda(ik:fk), diagC(ik:fk), L(k,k)%v, piv(ik:fk) )
      else
        CALL LDLt_Triang( g1(ik:fk), g2(ik:fk), lambda(ik:fk), diagC(ik:fk), L(k,k)%v )
      end if
      !$OMP END SINGLE
      !$OMP DO PRIVATE(ii,fi) schedule(runtime)
      do i = k+1, n
        ii= (i-1)*bs+1
        fi= i*bs
        if (fi.gt.size(g1,1)) fi= size(g1,1)
        CALL LDLt_Square( g1(ik:fk), g2(ik:fk), g1(ii:fi), g2(ii:fi), &
                          lambda(ik:fk), lambda(ii:fi), &
                          diagC(ik:fk), diagC(ii:fi), L(i,k)%v )      
      end do
      !$OMP END DO 
    end do
    !$OMP END PARALLEL 

  end subroutine LDLt_Blocks

  subroutine LDLt_Triang( g1, g2, lambda, diagC, L, piv ) 
    real(kind=r_kind), dimension(:,:), intent(inout) :: L
    real(kind=r_kind), dimension(:), intent(inout) :: diagC, g1, g2, lambda
    integer, dimension(:), intent(inout), optional :: piv

    real(kind=r_kind) :: Lik, d, lambdak_, g1k, g2k
    real(kind=r_kind) :: daux
    integer :: iaux
    integer :: n, i, j, k
    real(kind=r_kind), dimension(:), allocatable :: row
    integer, dimension(1) :: ind_piv
    n = size(L,1)

    if( present( piv ) ) then
      allocate( row( n ) )
      piv = (/ (i, i=1,n) /)
      do k = 1, n
         ind_piv = maxloc( abs( diagC(k:n) ) )
         j = ind_piv(1) + k - 1
         if( ind_piv(1).gt.1 ) then
           ! Swapping 
           iaux = piv(j)
           piv(j) = piv(k)
           piv(k) = iaux
           daux = g1(j)
           g1(j) = g1(k)
           g1(k) = daux
           daux = g2(j)
           g2(j) = g2(k)
           g2(k) = daux
           daux = lambda(j)
           lambda(j) = lambda(k)
           lambda(k) = daux
           daux = diagC(j)
           diagC(j) = diagC(k)
           diagC(k) = daux
           row(1:k-1) = L(j,1:k-1)
           L(j,1:k-1) = L(k,1:k-1)
           L(k,1:k-1) = row(1:k-1)
         end if
         d = diagC(k)
         lambdak_= lambda(k)
         g1k = g1(k)
         g2k = g2(k)
         L(k,k) = d
         L(k+1:n,k) = (-g2(k+1:n)*g1k + g1(k+1:n)*g2k) / ( (lambda(k+1:n) - lambdak_) * d )
         g1(k+1:n) = g1(k+1:n) - g1k * L(k+1:n,k)
         g2(k+1:n) = g2(k+1:n) - g2k * L(k+1:n,k)
         diagC(k+1:n) = diagC(k+1:n) - ( d * L(k+1:n,k)**2 )
      end do
    else
      do k=1, n
         d = diagC(k)
         lambdak_= lambda(k)
         g1k= g1(k)
         g2k= g2(k)
         L(k,k) = d
         do i=1+k, n
            Lik = ( (-g2(i)*g1k + g1(i)*g2k) / (lambda(i) - lambdak_) ) / d
            g1(i) = g1(i) - g1k * Lik
            g2(i) = g2(i) - g2k * Lik
            L(i,k) = Lik
            diagC(i) = diagC(i) - ( d * Lik**2 )
         end do 
      end do
    end if
    if( present( piv ) ) deallocate( row )
  end subroutine LDLt_Triang


  subroutine LDLt_Square( g1K, g2K, g1I, g2I, lambdaK, lambdaI, diagCK, diagCI, L ) 
    real(kind=r_kind), dimension(:), intent(in) :: lambdaK, lambdaI, g1K, g2K, diagCK
    real(kind=r_kind), dimension(:,:), intent(out) :: L
    real(kind=r_kind), dimension(:), intent(inout) :: diagCI, g1I, g2I
    real(kind=r_kind) :: Lik, d, lambdak_, g1k_, g2k_
    integer :: n1, n2, i, k
    n1 = size(L,1)
    n2 = size(L,2)

    do k=1, n2
       d = diagCK(k)
       lambdak_= lambdaK(k)
       g1k_= g1K(k)
       g2k_= g2K(k)
       do i=1, n1
          Lik = ( (-g2I(i)*g1k_ + g1I(i)*g2k_) / (lambdaI(i) - lambdak_) ) / d
          g1I(i) = g1I(i) - g1k_ * Lik
          g2I(i) = g2I(i) - g2k_ * Lik
          L(i,k) = Lik
          diagCI(i) = diagCI(i) - ( d * Lik**2 )
       end do 
    end do

  end subroutine LDLt_Square


  subroutine DTRSD_Blocks( LD, x, num_blocks, nb )
    type(block_type), dimension(:,:), intent(in) :: LD
    real(kind=r_kind), dimension(:), intent(inout):: x
    integer, intent(in) :: num_blocks, nb
    integer :: i, ii, idx, idx2, threads
    threads = OMP_GET_MAX_THREADS() / 2
    !$OMP PARALLEL PRIVATE(idx,ii,idx2) num_threads(threads)
    !$OMP DO 
    do i = 1, num_blocks
      idx = (i-1)*nb
      do ii = 1, LD(i,i)%x
        idx2 = ii+idx
        x(idx2) = x(idx2) / LD(i,i)%v(ii,ii)
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL 
  end subroutine DTRSD_Blocks


  subroutine allocate_LD(LD, num_blocks, m, nb)
    type(block_type), dimension(:,:), allocatable, intent(inout) :: LD
      integer, intent(in):: num_blocks, m, nb
    integer :: i, j
    allocate(LD(num_blocks,num_blocks))
    do j = 1, num_blocks
      do i = j, num_blocks
        LD(i,j)%x= nb
        LD(i,j)%y= nb
        if (mod(m,nb).ne.0) then
          if (i.eq.num_blocks) then
            LD(i,j)%x= mod(m,nb)
          end if
          if (j.eq.num_blocks) then
            LD(i,j)%y= mod(m,nb)
          end if
        end if
        allocate(LD(i,j)%v(LD(i,j)%x,LD(i,j)%y))
      end do
    end do
  end subroutine allocate_LD


  subroutine deallocate_LD(LD, num_blocks)
    type(block_type), dimension(:,:), allocatable, intent(inout) :: LD
    integer, intent(in):: num_blocks
    integer :: i, j
    do j = 1, num_blocks
        do i = j, num_blocks
        deallocate(LD(i,j)%v)
        end do
    end do
    deallocate(LD)
  end subroutine deallocate_LD  


  ! nrm1 operation
  !
  !> \brief Computes the 1-norm of a symmetric Toeplitz matrix T.
  !> \param[in]  T  First column (row) of matrix T.
  !> \return        1-norm of T.
  real(kind=r_kind) function sp_dts_nrm1(T)
    real(kind=r_kind), intent(in) :: T(:)

    integer :: n, i
    real(kind=r_kind), dimension(size(T)) :: a
    real(kind=r_kind) :: b, maximum = 0.0D+0

    n = size(T)
    a = T
    do i = 1, n/2+mod(n,2)
      b = sum(abs(a))
      if( b.gt.maximum ) maximum = b
      a = eoshift(a,-1,T(i+1))
    end do

    sp_dts_nrm1 = maximum
  end function sp_dts_nrm1


  ! SGEMV operation (not optimum)
  !
  !> \brief Performs the operation
  !>        b = b + alpha * T * x,
  !> where T is a symmetric Toeplitz matrix.
  !> \param[in]      alpha  Scalar alpha.
  !> \param[in]      T      First column (row) of matrix T.
  !> \param[in]      x      Vector x.
  !> \param[in,out]  b      Vector b. Changed on output.
  subroutine sp_dts_gemv(alpha, T, x, b)
    real(kind=r_kind), intent(in) :: alpha
    real(kind=r_kind), dimension(:), intent(in) :: T, x
    real(kind=r_kind), dimension(:), intent(inout) :: b

    integer :: n, i
    real(kind=r_kind), dimension(size(T)) :: a

    n = size(T)
    a = T
    do i = 1, n
      b = b + alpha * a * x(i)
      if( i.lt.n) a = eoshift(a,-1,T(i+1))
    end do

  end subroutine sp_dts_gemv


end module sp_dts
