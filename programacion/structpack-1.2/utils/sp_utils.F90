module sp_utils
  implicit none

  ! precision
  !
  integer, parameter, public :: r_kind = kind(1d+0)
  real(kind=r_kind), parameter, public :: ZERO   =   0_r_kind, & ! 0.0d+0
                                          ONE    =   1_r_kind, & ! 1.0d+0
                                          TWO    =   2_r_kind, & ! 2.0d+0
                                          THREE  =   3_r_kind, & ! 3.0d+0
                                          FOUR   =   4_r_kind, & ! 4.0d+0
                                          PI     =   FOUR*atan(ONE) ! Number pi
 type,public:: block_type
      REAL(kind=r_kind), DIMENSION(:,:), ALLOCATABLE :: v
      INTEGER                                        :: x, y
 end type block_type

  ! cl_a#rgs
  !
  public :: cl_get , DTRSV_Blocks, DTRSM_Blocks
  private :: cl_get_int, cl_get_real, cl_get_char

  interface cl_get
    module procedure cl_get_int, cl_get_real, cl_get_char
  end interface

  ! threadsnumber
  !
  ! If Variable VarName is not defined the default value of num_threads will be used
  !
  character(len=100), parameter :: Num_threads_variable = "OMP_NUM_THREADS"
  integer, public :: num_threads = 1

contains

  ! cl_args
  !
  ! Function cl_get_int
  !
  integer function cl_get_int(tag, value)
    character(*), intent(in) :: tag
    integer, intent(out)     :: value

    integer         :: n_args, i, l_buffer
    character(1024) :: buffer

    n_args = command_argument_count()
    do i=1,n_args
      !call getarg(i, buffer, status)
      call get_command_argument(i, buffer)
      l_buffer = len_trim(buffer)
      if (l_buffer == 0) then
        cl_get_int = -1
        return
      end if

      if (buffer(:l_buffer) .eq. tag) then
        !call getarg(i+1, buffer, status)
        call get_command_argument(i+1, buffer)
        l_buffer = len_trim(buffer)
        if (l_buffer .lt. 0) then
          cl_get_int = -2
          return
        end if

        !read (buffer, '(i)', iostat=status) value
        !read (buffer, '(i)') value
        read (buffer, '(i9)') value
        !read (buffer, *) value
        !if (status .ne. 0) then
        !   cl_get_int = -3
        !   return
        !end if

        cl_get_int = 0
        return
      end if

    end do

    cl_get_int = -1
    return
  end function cl_get_int


  ! Function cl_get_real
  !
  integer function cl_get_real(tag, value)
    character(*), intent(in)       :: tag
    real(kind=r_kind), intent(out) :: value

    integer         :: n_args, i, l_buffer
    character(1024) :: buffer

    n_args = command_argument_count()
    do i=1,n_args
      !call getarg(i, buffer, status)
      call get_command_argument(i, buffer)
      l_buffer = len_trim(buffer)
      if (l_buffer == 0) then
        cl_get_real = -1
        return
      end if

      !if (buffer(:status) .eq. tag) then
      if (buffer(:l_buffer) .eq. tag) then
        !call getarg(i+1, buffer, status)
        call get_command_argument(i+1, buffer)
        l_buffer = len_trim(buffer)
        if (l_buffer .lt. 0) then
          cl_get_real = -2
          return
        end if

        !read (buffer, '(f20.0)', iostat=status) value
        read (buffer, *) value
        !if (status .ne. 0) then
        !   cl_get_real = -3
        !   return
        !end if

        cl_get_real = 0
        return
      end if

    end do

    cl_get_real = -1
    return
  end function cl_get_real


  ! Function cl_get_char
  !
  integer function cl_get_char(tag, value)
    character(*), intent(in)  :: tag
    character,    intent(out) :: value

    integer         :: n_args, i, l_buffer
    character(1024) :: buffer

    n_args = command_argument_count()
    do i=1,n_args
      !call getarg(i, buffer, status)
      call get_command_argument(i, buffer)
      l_buffer = len_trim(buffer)
      if (l_buffer .eq. 0) then
        cl_get_char = -1
        return
      end if

      if (buffer(:l_buffer) .eq. tag) then
        !call getarg(i+1, buffer, status)
        call get_command_argument(i+1, buffer)
        l_buffer = len_trim(buffer)
        if (l_buffer .lt. 0) then
          cl_get_char = -2
          return
        end if

        !read (buffer, '(a)', iostat=status) value
        read (buffer, *) value
        if (l_buffer .lt. 0) then
          cl_get_char = -3
          return
        end if

        cl_get_char = 0
        return
      end if

    end do

    cl_get_char = -1
    return
  end function cl_get_char


  ! ioutils
  !
  ! Function read_nat_env_var
  ! Reads a natural (>=0) environment variable; returns a negative value
  ! in case of error
  !
  integer function read_nat_env_var(VarName)
    character(len=100), intent(in) :: VarName

    character(len=100) :: VarValue
    integer            :: Stat

    ! Get VarValue from the environment variable VarName
    ! assuming that its value is not more than 100 characters.
    !
    !call get_environment_variable(name=VarName, value=VarValue, status=Stat), trim_name=.false. )
    call get_environment_variable(name=VarName, value=VarValue, status=Stat)
    if (Stat/=0) then
      !print *, "Environment variable ",VarName," not set!"
      read_nat_env_var = -1
    else
      !read ( VarValue, '(i)', iostat=Stat ) read_nat_env_var
      read ( VarValue, * ) read_nat_env_var
    end if
    return
  end function read_nat_env_var
  

    !> \brief Solves the UP/LO triangular linear system stored in blocks 
    !>        L * x = P * b,   (1)
    !> if L is lower triangular, or
    !>        L' * P * x = b,  (2)
    !> if L' is upper triangular.
    !> \param[in]      UPLO       Triangular factor is upper ('U') or lower ('L') stored. Warning: only checked for UPLO='L'.
    !> \param[in]      LD         Triangular factor stored in block format (block_type).
    !> \param[in,out]  x          On entry, the right hand side array b. On output, the solution array x.
    !> \param[in]      num_blocks Number of blocks of LD.
    !> \param[in]      bSize      Block size.            
    !> \param[in]      TRANS      If TRANS='N' solves (1), otherwise solves (2). 
    !> \param[in]      DIAG       If DIAG='U' diagonal entries are implicitly 1 and will not be referenced. Otherwise diagonal entries of LD will be referenced.
    !> \param[in]      threads    Number of active threads in the system. 
    !> \param[in]      piv        If piv exists and is allocated, it contains the permutations used for pivoting by routine LDLt_Blocks
    !>                            that computes the L factor so P*C*P' = LDL', being C a Cauchy-like matrix. Otherwise, P = I.
    subroutine DTRSV_Blocks( UPLO, LD, x, num_blocks, bSize, TRANS, DIAG, threads, piv)
        type(block_type),  dimension(:,:), intent( in  )  :: LD
        real(kind=r_kind), dimension(:),   intent(inout)  :: x
        integer,                           intent( in  )  :: num_blocks, bSize, threads
        character,                         intent( in  )  :: UPLO, TRANS, DIAG
        integer,           dimension(:),   intent(in), allocatable, optional :: piv
        
        LOGICAL LSAME
        EXTERNAL LSAME
        integer :: i, k, ii, ik, n, blocksize
        integer, dimension(bSize) :: pivtrans

        if( LSAME(UPLO,'U')) then
            if (LSAME(TRANS,'T')) then
                !$OMP PARALLEL num_threads(threads) private(ik,n)
                do k = 1, num_blocks     
                    ik = (k-1)*bSize+1
                    n  = ik + bSize - 1
                    !$OMP SINGLE
                    if( present(piv) .and. allocated( piv ) ) x(ik:n) = x( piv(ik:n)+ik-1 )
                    call DTRSV( 'U', TRANS, DIAG, bSize, LD(k,k)%v, bSize, x(ik), 1 )
                    !$OMP END SINGLE
                    !$OMP DO private(ii) schedule(runtime)
                    do i = k+1, num_blocks
                        ii = (i-1)*bSize+1
                        call DGEMV( TRANS, bSize,bSize, -1.0D+0, LD(k,i)%v, bSize, x(ik), 1, 1.0D+0, x(ii), 1 )
                    end do
                    !$OMP END DO
                end do
                !$OMP END PARALLEL 
            else
                !$OMP PARALLEL num_threads(threads) private(ik,n)
                do k = num_blocks, 1, -1
                    ik = (k-1)*bSize+1
                    n  = ik + bSize- 1
                    !$OMP SINGLE
                    call DTRSV( 'U', TRANS, DIAG, bSize, LD(k,k)%v, bSize, x(ik), 1 )
                    if( present(piv) .and. allocated( piv ) ) x(ik:n) = x( piv(ik:n)+ik-1 ) ! To revise
                    !$OMP END SINGLE
                    !$OMP DO private(ii) schedule(runtime)
                    do i = k-1, 1, -1
                        ii = (i-1)*bSize+1
                        call DGEMV( TRANS, bSize,bSize, -1.0D+0, LD(i,k)%v, bSize, x(ik), 1, 1.0D+0, x(ii), 1 )
                    end do
                    !$OMP END DO
                end do
                !$OMP END PARALLEL 
            end if
        else if( LSAME(UPLO, 'L')) then
             if( LSAME(TRANS,'N')) then
                !$OMP PARALLEL num_threads(threads) private(ik,n)
                do k = 1, num_blocks     
                    ik = (k-1)*bSize+1
                    n  = ik + LD(k,k)%x - 1
                    !$OMP SINGLE
                    if( present(piv) .and. allocated( piv ) ) x(ik:n) = x( piv(ik:n)+ik-1 )
                    call DTRSV( 'L', TRANS, DIAG, LD(k,k)%x, LD(k,k)%v, LD(k,k)%x, x(ik), 1 )
                    !$OMP END SINGLE
                    !$OMP DO private(ii) schedule(runtime)
                    do i = k+1, num_blocks
                        ii = (i-1)*bSize+1
                        call DGEMV( TRANS, LD(i,k)%x, LD(i,k)%y, -1.0D+0, LD(i,k)%v, LD(i,k)%x, x(ik), 1, 1.0D+0, x(ii), 1 )
                    end do
                    !$OMP END DO 
                end do
                !$OMP END PARALLEL 
            else
                !$OMP PARALLEL num_threads(threads) private(ik,n,blocksize,pivtrans)
                do k = num_blocks, 1, -1
                    blocksize = LD(k,k)%x
                    ik = (k-1)*bSize+1
                    n  = ik + blocksize - 1
                    !$OMP SINGLE
                    call DTRSV( 'L', TRANS, DIAG, blocksize, LD(k,k)%v, blocksize, x(ik), 1 )
                    if( present(piv) .and. allocated( piv ) ) then 
                      pivtrans(piv(ik:n)) = (/(i,i=1,blocksize)/)
                      x(ik:n) = x( pivtrans(1:blocksize)+ik-1 )
                    end if
                    !$OMP END SINGLE
                    !$OMP DO private(ii) schedule(runtime)
                    do i = k-1, 1, -1
                        ii = (i-1)*bSize+1
                        call DGEMV( TRANS, LD(k,i)%x, LD(k,i)%y, -1.0D+0, LD(k,i)%v, LD(k,i)%x, x(ik), 1, 1.0D+0, x(ii), 1 )
                    end do
                    !$OMP END DO
                end do
                !$OMP END PARALLEL 
            end if
        else 
            print *, "Invalid UPLO value"
        end if


    end subroutine DTRSV_Blocks
    

subroutine DTRSM_Blocks( LD, x, num_blocks, nb, TRANS, DIAG )
    type(block_type), dimension(:,:), intent(in) :: LD
    real(kind=r_kind), dimension(:,:), intent(inout):: x
    integer, intent(in) :: num_blocks, nb
    character, intent(in) :: TRANS, DIAG
    integer :: i, k, ii, ik, n
    LOGICAL LSAME
    EXTERNAL LSAME
    integer :: threads, OMP_GET_MAX_THREADS

    threads = OMP_GET_MAX_THREADS() 

    if (LSAME(TRANS,'T')) then
      do k = 1, num_blocks     
        ik = (k-1)*nb+1
        n  = ik + nb - 1
  
        CALL DTRSM( 'L', 'U', TRANS, DIAG, nb, size(x,2), 1.0d+0, LD(k,k)%v, nb, x(ik:n,:), nb )
        
        !$OMP PARALLEL DO num_threads(threads) private(i,ii)
        do i = k+1, num_blocks
          ii = (i-1)*nb+1
          CALL DGEMM(TRANS,'N', nb, size(x,2), nb, -1.0D+0, LD(k,i)%v, nb, x(ik:n,:), nb, 1.0d+0,x(ii:ii+nb-1,:),nb)
        end do
        !$OMP END PARALLEL DO 
        
      end do
    else
      do k = num_blocks, 1, -1
        ik = (k-1)*nb+1
        n  = ik + nb- 1
         CALL DTRSM( 'L', 'U', TRANS, DIAG, nb, size(x,2), 1.0d+0, LD(k,k)%v, nb, x(ik:n,:), nb )
        !$OMP PARALLEL DO num_threads(threads) private(i,ii)
        do i = k-1, 1, -1
          ii = (i-1)*nb+1
          CALL DGEMM(TRANS,'N', nb, size(x,2), nb, -1.0D+0, LD(i,k)%v, nb, x(ik:n,:), nb, 1.0d+0,x(ii:ii+nb-1,:),nb)

        end do
        !$OMP END PARALLEL DO
      end do
    end if
  end subroutine DTRSM_Blocks

end module sp_utils
