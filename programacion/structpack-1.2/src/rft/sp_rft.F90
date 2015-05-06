!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
!  This module is an interface to the best real fft routine    !
!  dfftf (fftpack) or fft_mkl (mkl fft)                        !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module sp_rft
  use sp_utils   ! exports kind value r_kind
  use sp_ctimer   
#ifdef USE_MKL
  use mkl_dfti    ! mkl dft
#endif
  implicit none

  public :: rft, irft, rft_finalize

  private

  ! List data types
  type ll_rft_node
     integer :: n
     logical :: fftpack = .true.
     real( kind=r_kind ), pointer, dimension(:) :: Wreal
#ifdef USE_MKL
     type(DFTI_DESCRIPTOR), pointer ::  Desc_Handle 
#endif
     type(ll_rft_node), pointer :: next => null()
  end type ll_rft_node

  type ll_rft_list
     type(ll_rft_node) :: head
  end type ll_rft_list

  type(ll_rft_list), save :: ll_init_sizes


contains
  
  subroutine rft( w, rft_type_ ) 
    
    !     rft performs a Real Fast Fourier Tranformation to a vector w. 
    !
    !     This routine tests in a first call whether it is better to use
    !     fftpack or mkl. 
    !
    !     The Transformation is NOT normalized. A call to rft followed
    !     by a call to irft multiplies vector w by n
    !
    !     The first call to this routine initializes some vectors and workspace
    !     used in the subsequent calls with the same size. 
    !     A call with a different size produces a new initialization
    !
    !     This module must be compiled with openmp flag if will be use with
    !     conclurrent threads
    !
    !     rft_type_ is an optional argument that allows to force the method to use
    !     Values = { fftpack, mkl }
    !

    real(kind=r_kind), dimension(:), intent(inout), target :: w
    character(*), intent(in), optional :: rft_type_

    call rft_init( w, rft_type_ )
    call apply_rft(w,.true.)
    
  end subroutine rft
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine irft( w, rft_type_ ) 
       
    !     irft performs a Real Inverse Fast Fourier Tranformation to a vector w. 
    !   
    !     This routine tests in a first call whether it is better to use
    !     fftpack or mkl. 
    !   
    !     The Transformation is NOT normalized. A call to rft followed
    !     by a call to irft multiplies vector w by n
    !   
    !     The first call to this routine initializes some vectors and workspace
    !     used in the subsequent calls with the same size. 
    !     A call with a different size produces a new initialization
    !   
    !     This module must be compiled with openmp flag if will be use with
    !     conclurrent threads
    !   
    !     rft_type_ is an optional argument that allows to force the method to use
    !     Values = { fftpack, mkl }
    !   

    real(kind=r_kind), dimension(:), intent(inout), target :: w
    character(*), intent(in), optional :: rft_type_

    call rft_init( w, rft_type_ )
    call apply_rft(w,.false.)
            
  end subroutine irft 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine rft_init( w, rft_type_ )
    !
    !     rft_init performs the initialization, which is common for the two
    !     transf.
    !     Once a given size has been initialized with a method, it can't be
    !     reinitialized again with a different one.

    real(kind=r_kind), dimension(:), intent(inout), target :: w
    character(*), intent(in), optional :: rft_type_

    integer :: n, ierr
    character(len=8) :: tipo_rft
    
    n = size(w,1)
    
    !$OMP CRITICAL    
    if (.not. ll_rft_exists_n(ll_init_sizes, n)) then
      tipo_rft = 'fftpack'
      if( present(rft_type_) ) then
#ifdef USE_MKL
        if( rft_type_.eq.'mkl' .or. rft_type_.eq.'fftpack' ) then
          tipo_rft = rft_type_
        else
          tipo_rft = 'auto'
        end if
      else
        tipo_rft = 'auto'
#else
        tipo_rft = 'fftpack'
#endif
    end if
      ierr = init_rft(n,tipo_rft)
      if (ierr .ne. 0) stop 'Error initializing rft()'
    end if
    !$OMP END CRITICAL    
    
    return
    
  end subroutine rft_init
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  integer function init_rft(n,tipo_rft) 
    integer, intent(in) :: n
    character(len=8), intent(in) :: tipo_rft

    init_rft = ll_rft_add_n_rft(ll_init_sizes, n, tipo_rft )

  end function init_rft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine apply_rft(w,for_or_back) 
    real(kind=r_kind), dimension(:), intent(inout) :: w
    logical :: for_or_back

    integer :: n
#ifdef USE_MKL
    integer :: info
#endif

    type(ll_rft_node), pointer :: p_node
   
    n = size(w,1)
    !$OMP CRITICAL
    p_node => ll_rft_get_n_node(ll_init_sizes, n)
    !$OMP END CRITICAL
    if (.not. associated(p_node)) then
       print *, 'Error: trying to work with a non-initialized problem size'
       return
    end if

    if( p_node%fftpack ) then
      if( for_or_back ) then 
        call dfftf( n, w, p_node%Wreal )
      else 
        call dfftb( n, w, p_node%Wreal )
      end if
#ifdef USE_MKL
    else
      if( for_or_back ) then 
        info = DftiComputeForward( p_node%Desc_Handle, w )
      else 
        info = DftiComputeBackward( p_node%Desc_Handle, w )
      end if
      call DftierrorManager( info )
#endif
    end if

  end subroutine apply_rft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  integer function ll_rft_add_n_rft(ll, n, tipo_rft )
    type(ll_rft_list), intent(inout) :: ll
    integer, intent(in) :: n
    character(len=8), intent(in) :: tipo_rft

#ifdef USE_MKL
    integer :: i
    real(kind=r_kind), allocatable :: v(:)
    real(kind=r_kind) :: t1, t2, tn, gtime, dfftf_time, rft_mkl_time
#endif

    type(ll_rft_node), pointer :: node

    integer :: info

    if (ll_rft_exists_n(ll, n)) then
       ll_rft_add_n_rft = 1
       return
    end if

    allocate(node, stat=info)
    if (info .ne. 0) then
       print *,'Error allocating memory in ll_rft_add_n_rft'
       ll_rft_add_n_rft = -1
    end if

    node%n = n

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method selected to fftpack                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( tipo_rft.eq.'fftpack' ) then
      allocate(node%Wreal(int(2*n+15)), stat=info)
      if (info .ne. 0) then
         print *,'Error allocating memory in ll_rft_add_n_rft'
         ll_rft_add_n_rft = -1
      end if
      node%fftpack = .true.
      call dffti( n, node%Wreal )
    end if ! tipo_rft = fftpack

#ifdef USE_MKL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method selected to mkl                     !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(v(n), stat=info)
    if (info .ne. 0) then
       print *,'Error allocating memory in ll_rft_add_n_rft'
       ll_rft_add_n_rft = -1
    end if

    if( tipo_rft.eq.'mkl' ) then
      info = DftiCreateDescriptor( node%Desc_Handle, DFTI_DOUBLE, DFTI_REAL, 1, node%n )
      call DftierrorManager( info )
      info = DftiSetValue( node%Desc_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT )
      call DftierrorManager( info )
      info = DftiCommitDescriptor( node%Desc_Handle )
      call DftierrorManager( info )
      node%fftpack = .false.
    end if ! tipo_rft = mkl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method automatically chosen                !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( tipo_rft.eq.'auto' ) then
      allocate(node%Wreal(int(2*n+15)), stat=info)
      if (info .ne. 0) then
         print *,'Error allocating memory in ll_rft_add_n_rft'
         ll_rft_add_n_rft = -1
      end if

      call dffti(n, node%Wreal)

      ! Testing time for fftpack
      v = [(i,i=1,n)]
      gtime = 0.0_r_kind
      call ctimer( t1, tn, tn, gtime )
      call dfftf(n, v, node%Wreal)
      call ctimer( t2, tn, tn, gtime )
      dfftf_time = t2-t1

      ! Testing time for mkl
      info = DftiCreateDescriptor( node%Desc_Handle, DFTI_DOUBLE, DFTI_REAL, 1, node%n )
      call DftierrorManager( info )
      info = DftiSetValue( node%Desc_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT )
      call DftierrorManager( info )
      info = DftiCommitDescriptor( node%Desc_Handle )
      call DftierrorManager( info )

      v = [(i,i=1,n)]
      gtime = 0.0_r_kind
      call ctimer( t1, tn, tn, gtime )
      info   = DftiComputeForward( node%Desc_Handle, v )
      call ctimer( t2, tn, tn, gtime )
      rft_mkl_time = t2-t1

      if( dfftf_time.le.rft_mkl_time ) then
        info = DftiFreeDescriptor( node%Desc_Handle )
        call DftierrorManager( info )
      else
        node%fftpack = .false.
        deallocate(node%Wreal, stat=info)
        if (info .ne. 0) then
           print *,'Error deallocating memory in ll_rft_add_n_rft'
           ll_rft_add_n_rft = -1
        end if
      end if

    end if ! tipo_rft = auto
    deallocate(v, stat=info)
    if (info .ne. 0) then
       print *,'Error deallocating memory in ll_rft_add_n_rft'
       ll_rft_add_n_rft = -1
    end if
#endif

    node%next => ll%head%next
    ll%head%next => node

    ll_rft_add_n_rft = 0

  end function ll_rft_add_n_rft

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function ll_rft_exists_n(ll, n)
    type(ll_rft_list), intent(in), target :: ll
    integer, intent(in) :: n

    type(ll_rft_node), pointer :: p_node
    
    ll_rft_exists_n = .false.

    p_node => ll%head
    do while (associated(p_node%next))
       p_node => p_node%next
       if (p_node%n .eq. n) then
          ll_rft_exists_n = .true.
          exit
       end if
    end do

  end function ll_rft_exists_n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function ll_rft_get_n_node(ll,n)
    type(ll_rft_node), pointer :: ll_rft_get_n_node
    type(ll_rft_list), intent(in), target :: ll
    integer, intent(in) :: n

    ll_rft_get_n_node => ll%head
    do while (associated(ll_rft_get_n_node%next))
       ll_rft_get_n_node => ll_rft_get_n_node%next
       if ( ll_rft_get_n_node%n .eq. n) then
          return
       end if
    end do

    nullify(ll_rft_get_n_node)

  end function ll_rft_get_n_node

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

  integer function rft_finalize()

    rft_finalize =  ll_rft_finalize(ll_init_sizes)

  end function rft_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer function ll_rft_finalize(ll)
    type(ll_rft_list), intent(inout), target :: ll
    integer :: i, n, ret
    
    n = ll_rft_length(ll)
    
    do i=1, n
       ret = ll_rft_remove_node(ll,1)
       if (ret .ne. 0) then
          ll_rft_finalize = -1
          return
       end if
    end do
    
    ll_rft_finalize = 0
    
  end function ll_rft_finalize
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ll_rft_length(ll)
    type(ll_rft_list), intent(in), target :: ll
    type(ll_rft_node), pointer :: p_node
    integer :: n

    p_node => ll%head
    n=0
    do while (associated(p_node%next))
       p_node => p_node%next
       n = n+1
    end do

    ll_rft_length = n

  end function ll_rft_length

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ll_rft_remove_node(ll, i)
    type(ll_rft_list), intent(inout), target :: ll
    integer, intent(in) :: i

    type(ll_rft_node), pointer :: p_node, p_ant
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
       print *,'Error deallocating memory in ll_rft_remove_node'
       ll_rft_remove_node = -2
    end if

    ll_rft_remove_node = 0

  end function ll_rft_remove_node
  
end module sp_rft
