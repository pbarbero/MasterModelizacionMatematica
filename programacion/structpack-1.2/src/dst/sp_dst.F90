!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!  This module is an interface to the best dst routine    !
!  dsint (fftpack) or dst_mkl (mkl dst)                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module sp_dst
  use omp_lib
  use sp_utils
  use sp_ctimer
#ifdef USE_MKL
  use mkl_dfti  ! mkl dft
#endif
  implicit none

  public :: dst, dst_finalize, dst_time, dst_type

  real(kind=r_kind), save :: temps_dst=0.0_r_kind

  private

  ! List data types
  type ll_dst_node
     integer :: n
     logical :: fftpack = .true.
     real( kind=r_kind ), pointer, dimension(:) :: Wreal, v
#ifdef USE_MKL
     type(DFTI_DESCRIPTOR), pointer ::  Desc_Handle 
#endif
     type(ll_dst_node), pointer :: next => null()
  end type ll_dst_node

  type ll_dst_list
     type(ll_dst_node) :: head
  end type ll_dst_list

  type(ll_dst_list), save :: ll_init_sizes

  real( kind=r_kind ) :: tiempo_dst
  character(len=10) :: tipo_dst

contains
  
  subroutine dst( w, dst_type_ ) 
    
    !     dst performs a Discrete Sine Tranformation to a vector w. 
    !
    !     This routine tests in a first call whether it is better to use
    !     fftpack or mkl. 
    !
    !     The Transformation is NOT normalized. Two calls to this routine
    !     multiplies vector w by 2*(n+1).
    !
    !     The first call to this routine initializes some vectors and workspace
    !     used in the subsequent calls with the same size. 
    !     A call with a different size produces a new initialization.
    !
    !     This module must be compiled with openmp flag if will be use with
    !     concurrent threads.
    !
    !     dst_type_ is an optional argument that allows to force the method to use
    !     Values = { fftpack, mkl }
    !

    real(kind=r_kind), dimension(:), intent(inout), target :: w
    character(len=8), intent(in), optional :: dst_type_
    integer :: n, ierr

    tipo_dst = 'fftpack'
    if( present(dst_type_) ) then 
#ifdef USE_MKL
      if( dst_type_.eq.'mkl' .or. dst_type_.eq.'fftpack' ) then
        tipo_dst = dst_type_
      else
        tipo_dst = 'auto'
      end if
    else
      tipo_dst = 'auto'
#else
      tipo_dst = 'fftpack'
#endif
    end if

    n = size(w,1)
    
    !$OMP CRITICAL    
    if ( .not.ll_dst_exists_n(ll_init_sizes, n)) then
      ierr = init_dst(n)
      if (ierr .ne. 0) stop 'Error initializing dst()'
    end if

    call apply_dst(w)
    !$OMP END CRITICAL    

    return
  end subroutine dst
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  integer function init_dst(n) 
    integer, intent(in) :: n

    init_dst = ll_dst_add_n_dst(ll_init_sizes, n)

  end function init_dst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine apply_dst(w) 
    real(kind=r_kind), dimension(:), intent(inout) :: w

    integer :: n
#ifdef USE_MKL
    integer :: info
    real(kind=r_kind), allocatable :: v(:)
#endif

    type(ll_dst_node), pointer :: p_node
   
    ! Version no reentrante

    n = size(w,1)
    p_node => ll_dst_get_n_node(ll_init_sizes, n)
    if (.not. associated(p_node)) then
       print *, 'Error: trying to work with a non-initialized problem size'
       return
    end if

    if( p_node%fftpack ) then
      !print *,' '
      p_node%v(1:n) = w
      call dsint(n, p_node%v, p_node%Wreal)
      w = p_node%v(1:n)
#ifdef USE_MKL
    else
      ! No se sabe si el descriptor hay que protegerlo en una región crítica
      allocate(v(2*(n+1)), stat=info)
      if (info .ne. 0) then
        print *,'Error allocating memory in apply_dst'
      end if
      v(    1)       = 0.0_r_kind
      v(2:n+1)       = w
      v(  n+2)       = 0.0_r_kind
      v(n+3:2*(n+1)) = -w(n:1:-1)
      info = DftiComputeForward( p_node%Desc_Handle, v )
      call DftierrorManager( info )
      w = -v(3:2*(n+1):2)
      deallocate(v, stat=info)
      if (info .ne. 0) then
        print *,'Error allocating memory in apply_dst'
      end if
#endif
    end if

  end subroutine apply_dst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  integer function ll_dst_add_n_dst(ll, n)
    type(ll_dst_list), intent(inout) :: ll
    integer, intent(in) :: n
#ifdef USE_MKL
    real(kind=r_kind), allocatable :: v(:), w(:)
    real(kind=r_kind) :: t1, t2, tn, gtime, dsint_time, dst_mkl_time
#endif
    type(ll_dst_node), pointer :: node
    integer :: info
#ifdef USE_MKL
    integer :: i
#endif

    if (ll_dst_exists_n(ll, n)) then
      ll_dst_add_n_dst = 1
      return
    end if

    allocate(node, stat=info)
    if (info .ne. 0) then
      print *,'Error allocating memory in ll_dst_add_n_dst'
      ll_dst_add_n_dst = -1
    end if

    node%n = n

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method selected to fftpack                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( tipo_dst.eq.'fftpack' ) then
      allocate(node%v(n+1), node%Wreal(int(3*(n+1)+15)), stat=info)
      if (info .ne. 0) then
        print *,'Error allocating memory in ll_dst_add_n_dst'
        ll_dst_add_n_dst = -1
      end if
      node%fftpack = .true.
      call dsinti(n, node%Wreal)
    end if  ! tipo_dst = fftpack

#ifdef USE_MKL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method selected to mkl                     !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( tipo_dst.eq.'mkl' ) then
      info = DftiCreateDescriptor( node%Desc_Handle, DFTI_DOUBLE, DFTI_REAL, 1, 2*(node%n+1) )
      call DftierrorManager( info )
      info = DftiSetValue( node%Desc_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT )
      call DftierrorManager( info )
      info = DftiCommitDescriptor( node%Desc_Handle )
      call DftierrorManager( info )
      node%fftpack = .false.
    end if ! tipo_dst = mkl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method automatically chosen               !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( tipo_dst.eq.'auto' ) then
      allocate(node%v(n+1), node%Wreal(int(3*(n+1)+15)), stat=info)
      if (info .ne. 0) then
        print *,'Error allocating memory in ll_dst_add_n_dst'
        ll_dst_add_n_dst = -1
      end if

      allocate(v(n), w(2*(n+1)), stat=info)
      if (info .ne. 0) then
        print *,'Error allocating memory in ll_dst_add_n_dst'
        ll_dst_add_n_dst = -1
      end if

      v = [(i,i=1,n)] 

      call dsinti(n, node%Wreal)

      gtime = 0.0_r_kind
      call ctimer( t1, tn, tn, gtime )
      node%v(1:n) = v
      call dsint(n, node%v, node%Wreal)
      v = node%v(1:n) 
      call ctimer( t2, tn, tn, gtime )
      dsint_time = t2-t1

      info = DftiCreateDescriptor( node%Desc_Handle, DFTI_DOUBLE, DFTI_REAL, 1, 2*(node%n+1) )
      call DftierrorManager( info )
      info = DftiSetValue( node%Desc_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT )
      call DftierrorManager( info )
      info = DftiCommitDescriptor( node%Desc_Handle )
      call DftierrorManager( info )

      gtime = 0.0_r_kind
      call ctimer( t1, tn, tn, gtime )
      w(1) = 0.0_r_kind
      w(2:n+1) = v
      w(n+2) = 0.0_r_kind
      w(n+3:2*n+2) = -v(n:1:-1)
      info = DftiComputeForward( node%Desc_Handle, w )
      call DftierrorManager( info )
      v = -w(3:2*(n+1):2)
      call ctimer( t2, tn, tn, gtime )
      dst_mkl_time = t2-t1

      if( dsint_time.le.dst_mkl_time ) then
        info = DftiFreeDescriptor( node%Desc_Handle )
        call DftierrorManager( info )
        tiempo_dst = dsint_time
        tipo_dst = 'fftpack'
      else
        node%fftpack = .false.
        deallocate(node%v, node%Wreal, stat=info)
        if (info .ne. 0) then
          print *,'Error deallocating memory in ll_dst_add_n_dst'
          ll_dst_add_n_dst = -1
        end if
        tiempo_dst = dst_mkl_time
        tipo_dst = 'mkl'
      end if

      deallocate(v, w, stat=info)
      if (info .ne. 0) then
        print *,'Error deallocating memory in ll_dst_add_n_dst'
        ll_dst_add_n_dst = -1
      end if
    end if ! tipo_dst = auto
#endif

    node%next => ll%head%next
    ll%head%next => node

    ll_dst_add_n_dst = 0
  end function ll_dst_add_n_dst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function ll_dst_exists_n(ll, n)
    type(ll_dst_list), intent(in), target :: ll
    integer, intent(in) :: n

    type(ll_dst_node), pointer :: p_node
    
    ll_dst_exists_n = .false.

    p_node => ll%head
    do while (associated(p_node%next))
       p_node => p_node%next
       if (p_node%n .eq. n) then
          ll_dst_exists_n = .true.
          exit
       end if
    end do

  end function ll_dst_exists_n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function ll_dst_get_n_node(ll,n)
    type(ll_dst_node), pointer :: ll_dst_get_n_node
    type(ll_dst_list), intent(in), target :: ll
    integer, intent(in) :: n

    ll_dst_get_n_node => ll%head
    do while (associated(ll_dst_get_n_node%next))
       ll_dst_get_n_node => ll_dst_get_n_node%next
       if ( ll_dst_get_n_node%n .eq. n) then
          return
       end if
    end do

    nullify(ll_dst_get_n_node)

  end function ll_dst_get_n_node

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

#ifdef USE_MKL
  subroutine DftierrorManager( estado )
    integer, intent(in) :: estado

    if (estado .ne. 0) then 
      if( .not.DftiErrorClass( estado, DFTI_NO_ERROR ) ) then
        print *, 'Error: ',DftiErrorMessage( estado )
        stop
      end if
    end if
  end subroutine DftierrorManager
#endif 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  integer function dst_finalize()

    dst_finalize =  ll_dst_finalize(ll_init_sizes)
  end function dst_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer function ll_dst_finalize(ll)
    type(ll_dst_list), intent(inout), target :: ll
    integer :: i, n, ret
    
    n = ll_dst_length(ll)
    
    do i=1, n
       ret = ll_dst_remove_node(ll,1)
       if (ret .ne. 0) then
          ll_dst_finalize = -1
          return
       end if
    end do
    
    ll_dst_finalize = 0
  end function ll_dst_finalize
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ll_dst_length(ll)
    type(ll_dst_list), intent(in), target :: ll
    type(ll_dst_node), pointer :: p_node
    integer :: n

    p_node => ll%head
    n=0
    do while (associated(p_node%next))
       p_node => p_node%next
       n = n+1
    end do

    ll_dst_length = n
  end function ll_dst_length

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ll_dst_remove_node(ll, i)
    type(ll_dst_list), intent(inout), target :: ll
    integer, intent(in) :: i

    type(ll_dst_node), pointer :: p_node, p_ant
    integer :: j, info
    
    p_ant  => ll%head
    p_node => ll%head%next
    do j=1,(i-1)
       p_ant => p_node
       p_node => p_node%next
    end do

    p_ant%next => p_node%next

    deallocate(p_node, stat=info)
    if (info .ne. 0) then
       print *,'Error deallocating memory in ll_dst_remove_node'
       ll_dst_remove_node = -2
    end if

    ll_dst_remove_node = 0
  end function ll_dst_remove_node
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function dst_time ()
    real(kind=r_kind) :: dst_time

    dst_time=temps_dst
  end function dst_time
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function dst_type ()
    character(20) :: dst_type

    dst_type=tipo_dst
  end function dst_type

end module sp_dst
