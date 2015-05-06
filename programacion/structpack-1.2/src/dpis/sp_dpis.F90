!> \date April, 2011
!> \author StructPack Team
!> \brief Module for the solution of Tridiagonal Symmetric Toeplitz Linear Systems
module sp_dpis
  use iso_c_binding
  use sp_utils
  use sp_dst
  implicit none

  !public :: rojo_method, dst_method, ldlt_method, maxsvd, minsvd
  public :: sp_dpis_rojosv, sp_dpis_dstsv, &
    sp_dpis_ldltsv, sp_dpis_maxsvd, sp_dpis_minsvd, &
    sp_dpis_gemv

  private
  contains

  !> \brief Routine that calls the solver based on Rojo's method to return the solution of the linear system
  !>  T x = b,
  !> where T is tridiagonal symmetric Toeplitz.
  !> \param[in]         n       (integer) Size of matrix.
  !> \param[in]         t0      (real) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (real) Sub/Super--diagonal of the Toeplitz matrix.
  !> \param[in,out]     x       (real,dimension(:)) On input, the right hand side vector b. On output, the solution vector x.
  subroutine sp_dpis_rojosv(n, t0, t1, x)
    integer, intent(in) :: n
    real(kind=r_kind), intent(in) :: t0, t1
    real(kind=r_kind), intent(inout) :: x(n)
    !integer :: status

    integer :: ierr, i
    real(kind=r_kind) :: lambda, t, mmup, Cappa1, Cappa2, mu, cosmu, sinmu
    real(kind=r_kind), allocatable :: c(:), realpw(:), imagpw(:), realw(:), imagw(:), realp(:), imagp(:), realz(:), imagz(:)
    real(kind=r_kind), allocatable :: p(:), w(:), pw(:), z(:)

    lambda = -t0/t1
    allocate( c(n), stat=ierr )
    if( ierr.ne.0 ) stop "Error allocating memory"
    c      = -x/t1

    if( abs(lambda).lt.2.0_r_kind ) then
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !% Non Diagonal dominant case  %
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! print *,'Non diagonal dominant case'
        cosmu=lambda/2_r_kind 
        if( lambda.eq.0.0_r_kind ) then 
          sinmu=-1.0_r_kind
        else
          t=sign(sqrt( 4.0_r_kind/lambda**2 - 1.0_r_kind ),lambda) 
          sinmu=-cosmu*t 
        end if
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !% Solution of the 1-rank linear system
        !% 
        !% C z = c
        !%
        !% and of the linear system
        !% C p= e1
        !% 
        !% where e1=[1 0 0 ...0]
        !% where C = ( mu     -1 0 ... 0 )
        !%           ( -1 lambda -1 .. 0 )
        !%           (  0     -1 ........)
        !%
        !% mu is supposed to be |mu| = 1
        !% and mu + 1/mu = lambda.
        !% conj(mu)=cosmu+sinmu*j
        !% returns the real part of z.
        !% and the real and imaginary part of p 

        allocate( realp(n), imagp(n), realpw(n), imagpw(n), realz(n), imagz(n), realw(n), imagw(n), stat=ierr )
        if( ierr.ne.0 ) stop "Error allocating memory"

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%% Iterations %%%%%%%%%%%%%%%%%%%%%%%%%
        realpw(1) = 1.0_r_kind
        imagpw(1) = 0.0_r_kind 

        realw(1) = c(1) 
        imagw(1) = 0.0_r_kind 

        do i = 2,n
           realw(i) = c(i) + realw(i-1)*cosmu - imagw(i-1)*sinmu 
           imagw(i) =        realw(i-1)*sinmu + imagw(i-1)*cosmu 

           realpw(i) =        realpw(i-1)*cosmu - imagpw(i-1)*sinmu 
           imagpw(i) =        realpw(i-1)*sinmu + imagpw(i-1)*cosmu 

        end do

        realp(n) = realpw(n)*cosmu - imagpw(n)*sinmu 
        imagp(n) = imagpw(n)*cosmu + realpw(n)*sinmu 

        realz(n) = realw(n)*cosmu - imagw(n)*sinmu 
        imagz(n) = imagw(n)*cosmu + realw(n)*sinmu 

        do i=n-1,1,-1
           realp(i) = (realpw(i)+realp(i+1))*cosmu - (imagpw(i)+imagp(i+1))*sinmu 
           imagp(i) = (imagpw(i)+imagp(i+1))*cosmu + (realpw(i)+realp(i+1))*sinmu 

           realz(i) = (realw(i)+realz(i+1))*cosmu - (imagw(i)+imagz(i+1))*sinmu 
           imagz(i) = (imagw(i)+imagz(i+1))*cosmu + (realw(i)+realz(i+1))*sinmu 

        end do
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mmup = ((cosmu+realp(1))**2 + (-sinmu+imagp(1))**2) 

        Cappa1 = realz(1)*(cosmu+realp(1))+(imagz(1)*(-sinmu+imagp(1))) 
        Cappa1 = Cappa1/mmup 

        Cappa2 = -realz(1)*(-sinmu+imagp(1))+(imagz(1)*(cosmu+realp(1))) 
        Cappa2 = Cappa2/mmup 

        x = realz-Cappa1*realp+Cappa2*imagp 

        deallocate( realp, imagp, realpw, imagpw, realz, imagz, stat=ierr )
        if( ierr.ne.0 ) stop "Error deallocating memory"

    else 
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !% Diagonal dominant case  %
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !print *,'Non diagonal dominant case'

        mu = lambda/2_r_kind + sign(sqrt( lambda**2/4_r_kind - 1_r_kind ),lambda) 
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !% Solution of the 1-rank linear system
        !% 
        !% C z = c
        !%
        !% and of the linear system
        !% C p= e1
        !% 
        !% where e1=[1 0 0 ...0]
        !% where C = ( mu     -1 0 ... 0 )
        !%           ( -1 lambda -1 .. 0 )
        !%           (  0     -1 ........)
        !%

        allocate( w(n), p(n), pw(n), z(n), stat=ierr )
        if( ierr.ne.0 ) stop "Error allocating memory"

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !%%%%%%%%%%% Iterations %%%%%%%%%%%%%%%%%%%%%%%%%
        pw(1) = 1 
        w(1)  = c(1) 
        do i = 2,n
           w(i)  = c(i) + w(i-1)/mu 
           pw(i) = pw(i-1)/mu 
        end do
        z(n) = w(n)/mu 
        p(n) = pw(n)/mu 
        do i = n-1,1,-1
           z(i) = (w(i)+z(i+1))/mu 
           p(i) = (pw(i)+p(i+1))/mu  
        end do
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        mmup = z( 1 ) / ( mu + p( 1 ) ) 
        x = z - mmup * p 

        deallocate( w, p, pw, z, stat=ierr )
        if( ierr.ne.0 ) stop "Error deallocating memory"

    end if

  end subroutine sp_dpis_rojosv


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% DST method                                            %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !> \brief Routine that calls the solver based on the DST method to return the solution of the linear system
  !>  T x = b,
  !> where T is tridiagonal symmetric Toeplitz.
  !> \param[in]         n       (integer) Size of matrix.
  !> \param[in]         t0      (real) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (real) Sub/Super--diagonal of the Toeplitz matrix.
  !> \param[in,out]     x       (real,dimension(:)) On input, the right hand side vector b. On output, the solution vector x.
  subroutine sp_dpis_dstsv(n, t0, t1, x)
    integer, intent(in) :: n
    real(kind=r_kind), intent(in) :: t0, t1
    real(kind=r_kind), intent(inout) :: x(n)

    real(kind=r_kind), allocatable :: s(:), beta(:)
    real(kind=r_kind)              :: k
    integer :: m, o, ierr
    character(len=8) :: dst_type

    ! Choose dst_type: fftpack, mkl, auto
    dst_type = 'mkl'

    allocate( s(n), beta(n/2), stat=ierr )
    if( ierr.ne.0 ) stop "Error allocating memory"

    m = n/2
    o = n/2
    k = 1.0_r_kind / ( 2.0_r_kind * real(n+1) )

    s    = 0.0_r_kind
    s(1) = 1.0_r_kind
    call dst(s,dst_type)

    beta = s(2:n:2) / s(1:o)
    call dst(x,dst_type)
    x(1:m) = k * ( x(1:m) / ( t0 + t1 * beta ) )
    if( mod(n,2).eq.1 ) then
      x(m+1) = k * x(m+1) / t0
      m = m+1
    end if
    x(m+1:n) = k * ( x(m+1:n) / ( t0 - t1 * beta(o:1:-1) ) )
    call dst(x,dst_type)

    deallocate( s, beta, stat=ierr )
    if( ierr.ne.0 ) stop "Error deallocating memory"

  end subroutine sp_dpis_dstsv


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% LDLT method                                           %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !> \brief Routine that calls the solver based on the LDLt method to return the solution of the linear system
  !>  T x = b,
  !> where T is tridiagonal symmetric Toeplitz.
  !> \param[in]         n       (integer) Size of matrix.
  !> \param[in]         t0      (real) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (real) Sub/Super--diagonal of the Toeplitz matrix.
  !> \param[in,out]     x       (real,dimension(:)) On input, the right hand side vector b. On output, the solution vector x.
  !> \param[out]        status  (integer) Returns 0 on success.
  subroutine sp_dpis_ldltsv(n, t0, t1, x, status)
    real(kind=r_kind), intent(in) :: t0, t1
    integer, intent(in) :: n
    real(kind=r_kind), intent(inout) :: x(n)
    integer, intent(out) :: status

    integer :: i
    real(kind=r_kind), allocatable :: d(:), l(:)
    real(kind=r_kind) :: pivot

    allocate( d(n), l(n), stat=status )
    if( status.ne.0 ) stop "Error allocating memory"

    status = 0
    d( 1 ) = t0
    do i = 1,n-1
      pivot    = d( i )
      if( abs(pivot).lt.1.0e-15 ) then 
         status = -1
         !return 
      end if
      l( i   ) = t1 / pivot
      d( i+1 ) = t0 - t1 * l( i )
    end do
    do i = 2,n 
      x( i ) = x( i ) - L( i-1 ) * x( i-1 )
    end do
    x( n ) = x( n ) / D( n )
    do i = n-1,1,-1 
      x( i ) = x( i ) / D( i ) - L( i ) * x( i+1 )
    end do

    deallocate( d, l, stat=status )
    if( status.ne.0 ) stop "Error deallocating memory"

  end subroutine sp_dpis_ldltsv


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% SVD functions                                         %%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !> \brief Function that returns the maximum SVD of a symmetric tridiagonal Toeplitz matrix. 
  !> \param[in]         n       (integer) Size of matrix.
  !> \param[in]         t0      (real) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (real) Sub/Super--diagonal of the Toeplitz matrix.
  !> \return                    (real) Maximum SVD of the symmetric tridiagonal Toeplitz matrix.
  real(kind=r_kind) function sp_dpis_maxsvd(n, t0, t1)
    integer, intent(in) :: n
    real(kind=r_kind), intent(in) :: t0, t1

    real(kind=r_kind), parameter :: pi = 4.0_r_kind * atan( 1.0_r_kind )

    sp_dpis_maxsvd = t0 + t1 * 2.0_r_kind * cos( pi / (n+1) ) 

  end function sp_dpis_maxsvd


  !real(kind=r_kind) function minsvd( n, t0, t1 )
  !  integer, intent(in) :: n
  !  real(kind=r_kind), intent(in) :: t0, t1
  !
  !  real(kind=r_kind), allocatable, dimension(:) :: s
  !
  !  allocate( s(n) )
  !  call svd( n, t0, t1, s )
  !  minsvd = minabsvalue(s)
  !  deallocate( s )
  !
  !end function minsvd

  !real(kind=r_kind) function minsvd( n, t0, t1 )
  !  integer, intent(in)            :: n
  !  real(kind=r_kind), intent(in)  :: t0, t1
  !
  !  integer :: i
  !  real(kind=r_kind) :: pi = 4.0_r_kind * atan( 1.0_r_kind )
  !  real(kind=r_kind) :: tcheb0, tcheb1, tcheb, x, y
  !
  !  x = cos( pi / ( n+1 ) )
  !  tcheb0 = 1.0_r_kind
  !  tcheb1 = x
  !  minsvd = abs( t0 - 2.0_r_kind * tcheb1 * t1 )
  !  do i = 3,n/2+1
  !    tcheb = 2.0_r_kind * x * tcheb1 - tcheb0
  !    !y = abs( t0 - 2.0_r_kind * tcheb * t1 )
  !    y = abs( t0 - 2.0_r_kind * cos( (i-2)*pi/(n+1) ) * t1 )
  !    tcheb0 = tcheb1
  !    tcheb1 = tcheb
  !    if( y.lt.minsvd ) minsvd = y
  !  end do
  !
  !end function minsvd

  !real(kind=r_kind) function minsvd( n, t0, t1 )
  !  integer, intent(in)            :: n
  !  real(kind=r_kind), intent(in)  :: t0, t1
    !
  !  integer :: i
  !  real(kind=r_kind) :: pi = 4.0_r_kind * atan( 1.0_r_kind )
  !  real(kind=r_kind) :: c, x, alfa2, alfa3
  !
  !  c = t0/(2.0_r_kind*t1)
  !
  !  if (c < 1.0_r_kind) then
  !    x  = (n+1) * atan( sqrt((1.0_r_kind/c**2.0_r_kind)-1.0_r_kind) ) / pi
  !    alfa2  = abs(t0-2.0_r_kind*t1*cos(floor(x)*pi/(n+1)))
  !    alfa3  = abs(t0-2.0_r_kind*t1*cos(ceiling(x)*pi/(n+1)))
  !    minsvd = min(alfa2,alfa3)
  !  else 
  !    minsvd = abs(t0-2.0_r_kind*t1*cos(1.0_r_kind*pi/(n+1)))
  !  end if
  !
  !end function minsvd

  !> \brief Function that returns the minimum SVD of a symmetric tridiagonal Toeplitz matrix. 
  !> \param[in]         n       (integer) Size of matrix.
  !> \param[in]         t0      (real) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (real) Sub/Super--diagonal of the Toeplitz matrix.
  !> \return                    (real) Minimum SVD of the symmetric tridiagonal Toeplitz matrix.
  real(kind=r_kind) function sp_dpis_minsvd(n, t0, t1)
    integer, intent(in)            :: n
    real(kind=r_kind), intent(in)  :: t0, t1
  
    real(kind=r_kind), parameter :: pi = 4.0_r_kind * atan( 1.0_r_kind )
    real(kind=r_kind) :: c, x, alfa2, alfa3

    c = t0/(2.0_r_kind*t1)

    if (c < 1.0_r_kind) then
      x  = (n+1) * atan( sqrt((1.0_r_kind/c**2.0_r_kind)-1.0_r_kind) ) / pi
      alfa2  = abs(t0-2.0_r_kind*t1*cos(floor(x)*pi/(n+1)))
      alfa3  = abs(t0-2.0_r_kind*t1*cos(ceiling(x)*pi/(n+1)))
      sp_dpis_minsvd = min(alfa2,alfa3)
    else 
      sp_dpis_minsvd = abs(t0-2.0_r_kind*t1*cos(1.0_r_kind*pi/(n+1)))
    end if
  
  end function sp_dpis_minsvd


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !% SGEMV operation with tridiagonal Toeplitz (not optimum)%%
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !> \brief Routine that performs operation
  !>        b = b + alpha * T * x,
  !> where T is a symmetric tridiagonal Toeplitz matrix.
  !> \param[in]         n       (integer) Size of matrix.
  !> \param[in]         alpha   (real) Scalar alpha.
  !> \param[in]         t0      (real) Diagonal of the Toeplitz matrix.
  !> \param[in]         t1      (real) Sub/Super--diagonal of the Toeplitz matrix.
  !> \param[in]         x       (real,dimension(:)) Vector x.
  !> \param[in,out]     b       (real,dimension(:)) Vector b. Changed on output.
  subroutine sp_dpis_gemv(n, alpha, t0, t1, x, b)
    integer, intent(in) :: n
    real(kind=r_kind), intent(in) :: alpha, t0, t1
    real(kind=r_kind), intent(in) :: x(n)
    real(kind=r_kind), intent(inout) :: b(n)

    integer :: i

    b( 1 ) = b( 1 ) + alpha * ( t0 * x( 1 ) + t1 * x( 2 ) )
    do i = 2, n-1
      b( i ) = b( i ) + alpha * ( t1*x(i-1) + t0*x(i) + t1*x(i+1) )
    end do

    b( n ) = b( n ) + alpha * ( t1 * x( n-1 ) + t0 * x( n ) )

  end subroutine sp_dpis_gemv


  !Atencion: Mejorable
  subroutine svd( n, t0, t1, s )
    integer, intent(in)            :: n
    real(kind=r_kind), intent(in)  :: t0, t1
    real(kind=r_kind), intent(out) :: s(n)
  
    real(kind=r_kind), dimension(n/2) :: beta
    integer :: m, o
  
    m = n/2
    o = n/2
  
    s    = 0.0_r_kind
    s(1) = 1.0_r_kind
    call dst(s)
  
    beta = s(2:n:2) / s(1:o)
    s(1:m) =  t0 + t1 * beta 
    if( mod(n,2).eq.1 ) then
      s(m+1) = t0
      m = m+1
    end if
    s(m+1:n) = t0 - t1 * beta(o:1:-1)  
  
  end subroutine svd

  ! Auxiliar routines
  real(kind=r_kind) function minabsvalue( v )
    !integer, intent(in) :: n
    real(kind=r_kind), intent(in), dimension(:) :: v

    integer :: n, i

    n = size( v, 1 )
    minabsvalue = v( 1 )
    do i = 2,n
      if( abs(v(i)).lt.minabsvalue ) minabsvalue = v( i )
    end do

  end function minabsvalue


end module sp_dpis
