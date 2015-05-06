!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                           !
!  This module is an interface to the best dctii routines   !
!  dcosqf/dcosqb (fftpack) or a mkl-fft base implementation !
!                                                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module sp_dctii
  use sp_utils
  use sp_ctimer
#ifdef USE_MKL
  use mkl_dfti  ! mkl dft
#endif
  implicit none

  public :: dctii, idctii, dctii_finalize

  private
  ! List data types
  type ll_dctii_node
     integer :: n
     logical :: fftpack = .true.
     real( kind=r_kind ), pointer, dimension(:) :: w, Wreal, seno, coseno
#ifdef USE_MKL
     type(DFTI_DESCRIPTOR), pointer ::  Desc_Handle 
#endif
     type(ll_dctii_node), pointer :: next => null()
  end type ll_dctii_node

  type ll_dctii_list
     type(ll_dctii_node) :: head
  end type ll_dctii_list

  type(ll_dctii_list), save :: ll_init_sizes

  character(len=8) :: tipo_dctii

contains
  
  ! Direct Cosine II transformation
  !
  subroutine dctii( w, dctii_type_ ) 
    
    !     dctii performs a Discrete Cosine Tranformation II of a vector w. 
    !
    !     This routine tests in a first call whether it is better to use fftpack or mkl.
    !
    !     This transformation uses routine dcosqb of fftpack.
    !     It is NOT normalized. A call to this routine followed by a call to idctii multiplies vector w by 4*n.
    !
    !     The first call to this routine initializes some vectors and workspace
    !     used in the subsequent calls with the same size. 
    !     A call with a different size produces a new initialization.
    !
    !     This module must be compiled with the openmp flag to be used with concurrent threads.
    !
    !     dctii_type_ is an optional argument that allows to force the method to use
    !     Available values = { fftpack, mkl }
    !

    real(kind=r_kind), dimension(:), intent(inout), target :: w
    character(*), intent(in), optional :: dctii_type_

    call dctii_init( w, dctii_type_ ) 
    call apply_dctii(w)

  end subroutine dctii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  ! Inverse Cosine II transformation
  !
  subroutine idctii( w, dctii_type_ ) 
    
    !     idctii performs a Discrete Cosine Tranformation II of a vector w. 
    !
    !     This routine tests in a first call whether it is better to use fftpack or mkl. 
    !
    !     This transformation uses routine dcosqf of fftpack.
    !     It is NOT normalized. A call to this routine followed by a call to dctii multiplies vector w by 4*n.
    !
    !     The first call to this routine initializes some vectors and workspace
    !     used in the subsequent calls with the same size. 
    !     A call with a different size produces a new initialization.
    !
    !     This module must be compiled with the openmp flag to be used with concurrent threads.
    !
    !     dctii_type_ is an optional argument that allows to force the method to use
    !     Available values = { fftpack, mkl }
    !

    real(kind=r_kind), dimension(:), intent(inout), target :: w
    character(*), intent(in), optional :: dctii_type_

    call dctii_init( w, dctii_type_ ) 
    call apply_idctii(w)

  end subroutine idctii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine dctii_init( w, dctii_type_ ) 
    !
    !     dctii_init performs the initialization, which is common for the two transf.
    !     Once a given size has been initialized with a method, it can't be reinitialized again with a different one.
    !

    real(kind=r_kind), dimension(:), intent(inout), target :: w
    character(*), intent(in), optional :: dctii_type_

    integer :: n, ierr
    
    n = size(w,1)

    !$OMP CRITICAL    
    if ( .not.ll_dctii_exists_n(ll_init_sizes, n)) then
      tipo_dctii = 'fftpack'
      if( present(dctii_type_) ) then 
#ifdef USE_MKL
        if( dctii_type_.eq.'mkl' .or. dctii_type_.eq.'fftpack' ) then
          tipo_dctii = dctii_type_
        else
          tipo_dctii = 'auto'
        end if
      else
        tipo_dctii = 'auto'
#else
        tipo_dctii = 'fftpack'
#endif
      end if
    
      ierr = init_dctii(n)
      if (ierr .ne. 0) stop 'Error initializing dctii()'
    end if
    !$OMP END CRITICAL    

    return
    
  end subroutine dctii_init
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  integer function init_dctii(n) 
    integer, intent(in) :: n

    init_dctii = ll_dctii_add_n_dctii(ll_init_sizes, n)

  end function init_dctii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine apply_dctii(w) 
    real(kind=r_kind), dimension(:), intent(inout) :: w

    integer :: n
#ifdef USE_MKL
    integer :: i, j, info
#endif

    type(ll_dctii_node), pointer :: p_node
   
    !! Version no reentrante !!

    n = size(w,1)
    p_node => ll_dctii_get_n_node(ll_init_sizes, n)
    if (.not. associated(p_node)) then
       print *, 'Error: trying to work with a non-initialized problem size'
       return
    end if

    if( p_node%fftpack ) then
      call dcosqb(n, w, p_node%Wreal)
#ifdef USE_MKL
    else
      ! No se sabe si el descriptor hay que protegerlo en una región crítica
      p_node%w(1:n) = (/ w(1:n:2), w(n/2*2:2:-2) /)
      info   = DftiComputeForward( p_node%Desc_Handle, p_node%w )
      call DftierrorManager( info )
      w(1)   = TWO * p_node%w(1) * p_node%coseno(1)
      j = 2
      do i = 2, (n+1)/2
        w(i) = TWO * ( p_node%w(j) * p_node%coseno(i) - p_node%w(j+1) * p_node%seno(i) )
        j = j + 2
      end do
      if( mod(n,2).eq.0 ) then
        w(i) = TWO * ( p_node%w(n) * p_node%coseno(i) )
        i = i + 1
      end if
      j = j - 2
      do i = i, n
        w(i) = TWO * ( p_node%w(j) * p_node%coseno(i) + p_node%w(j+1) * p_node%seno(i) )
        j = j - 2
      end do
#endif
    end if

  end subroutine apply_dctii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  subroutine apply_idctii(w) 
    real(kind=r_kind), dimension(:), intent(inout) :: w

    integer :: n
#ifdef USE_MKL
    !integer :: i, info
    integer :: info
#endif

    type(ll_dctii_node), pointer :: p_node
   
    !! Version no reentrante !!

    n = size(w,1)
    p_node => ll_dctii_get_n_node(ll_init_sizes, n)
    if (.not. associated(p_node)) then
       print *, 'Error: trying to work with a non-initialized problem size'
       return
    end if

    if( p_node%fftpack ) then
      call dcosqf(n, w, p_node%Wreal)
#ifdef USE_MKL
    else
      ! No se sabe si el descriptor hay que protegerlo en una región crítica
      p_node%w( 1 )     = p_node%coseno( 1 ) * w( 1 ) * 0.5D+00
      p_node%w( 2:n:2 ) = ( p_node%coseno(2:n/2+1) * w( 2:n/2+1 ) - p_node%seno(2:n/2+1) * w( n:(n+3)/2:-1 ) ) * 0.5D+00
      p_node%w( 3:n:2 ) = ( - p_node%seno(2:(n+1)/2) * w( 2:(n+1)/2 ) - p_node%coseno(2:(n+1)/2) * w( n:(n+4)/2:-1 ) ) * 0.5D+00
      info = DftiComputeBackward( p_node%Desc_Handle, p_node%w )
      call DftierrorManager( info )
      w(1:n:2) = p_node%w(1:n/2) 
      w(2:n:2) = p_node%w(n:n/2+1:-1) 
#endif
    end if

  end subroutine apply_idctii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

  integer function ll_dctii_add_n_dctii(ll, n)
    type(ll_dctii_list), intent(inout) :: ll
    integer, intent(in) :: n

#ifdef USE_MKL
    real(kind=r_kind), allocatable :: v(:)
    real(kind=r_kind) :: t1, t2, tn, gtime, dcosqb_time, dctii_mkl_time
    !double complex :: twon
    complex(kind=r_kind) :: twon
    integer :: i, j
#endif

    type(ll_dctii_node), pointer :: node

    integer :: info

    if (ll_dctii_exists_n(ll, n)) then
       ll_dctii_add_n_dctii = 1
       return
    end if

    allocate(node, stat=info)
    if (info .ne. 0) then
       print *,'Error allocating memory in ll_dctii_add_n_dctii'
       ll_dctii_add_n_dctii = -1
    end if

    node%n = n

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method selected to fftpack                 !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( tipo_dctii.eq.'fftpack' ) then
      allocate(node%Wreal(int(3*n+15)), stat=info)
      if (info .ne. 0) then
         print *,'Error allocating memory in ll_dctii_add_n_dctii'
         ll_dctii_add_n_dctii = -1
      end if
      node%fftpack = .true.
      call dcosqi(n, node%Wreal)
    end if ! tipo_dctii = fftpack

#ifdef USE_MKL
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method selected to mkl                     !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(v(n), stat=info)
    if (info .ne. 0) then
       print *,'Error allocating memory in ll_dctii_add_n_dctii'
       ll_dctii_add_n_dctii = -1
    end if
    allocate(node%w(n), node%coseno(n), node%seno(n), stat=info)
    if (info .ne. 0) then
       print *,'Error allocating memory in ll_dctii_add_n_dctii'
       ll_dctii_add_n_dctii = -1
    end if

    if( tipo_dctii.eq.'mkl' ) then
      info = DftiCreateDescriptor( node%Desc_Handle, DFTI_DOUBLE, DFTI_REAL, 1, node%n )
      call DftierrorManager( info )
      info = DftiSetValue( node%Desc_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT )
      call DftierrorManager( info )
      info = DftiCommitDescriptor( node%Desc_Handle )
      call DftierrorManager( info )
      node%fftpack = .false.
      twon = TWO*n
      v = [(-i*PI/twon,i=1,n-1)] 
      node%coseno(1)   = TWO
      node%seno(1)     = ZERO
      node%coseno(2:n) = TWO*cos(v(1:n-1))
      node%seno(2:n)   = TWO*sin(v(1:n-1))
    end if ! tipo_dctii = mkl

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !! Method automatically chosen                !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if( tipo_dctii.eq.'auto' ) then
      allocate(node%Wreal(int(3*n+15)), stat=info)
      if (info .ne. 0) then
         print *,'Error allocating memory in ll_dctii_add_n_dctii'
         ll_dctii_add_n_dctii = -1
      end if

      call dcosqi(n, node%Wreal)

      ! Testing time for fftpack
      v = [(i,i=1,n)] 
      gtime = 0.0_r_kind
      call ctimer( t1, tn, tn, gtime )
      call dcosqb(n, v, node%Wreal)
      call ctimer( t2, tn, tn, gtime )
      dcosqb_time = t2-t1

      info = DftiCreateDescriptor( node%Desc_Handle, DFTI_DOUBLE, DFTI_REAL, 1, node%n )
      call DftierrorManager( info )
      info = DftiSetValue( node%Desc_Handle, DFTI_PACKED_FORMAT, DFTI_PACK_FORMAT )
      call DftierrorManager( info )
      info = DftiCommitDescriptor( node%Desc_Handle )
      call DftierrorManager( info )
      twon = TWO*n
      v = [(-i*PI/twon,i=1,n-1)] 
      node%coseno(1)   = TWO
      node%seno(1)     = ZERO
      node%coseno(2:n) = TWO*cos(v(1:n-1))
      node%seno(2:n)   = TWO*sin(v(1:n-1))

      ! Testing time for mkl
      v = [(i,i=1,n)] 
      gtime = 0.0_r_kind
      call ctimer( t1, tn, tn, gtime )
      node%w(1:n) = (/ v(1:n:2), v(n/2*2:2:-2) /)
      info   = DftiComputeForward( node%Desc_Handle, node%w )
      call DftierrorManager( info )
      v(1)   = TWO * node%w(1) * node%coseno(1)
      j = 2
      do i = 2, (n+1)/2
        v(i) = TWO * ( node%w(j) * node%coseno(i) - node%w(j+1) * node%seno(i) )
        j = j + 2
      end do
      if( mod(n,2).eq.0 ) then
        v(i) = TWO * ( node%w(n) * node%coseno(i) )
        i = i + 1
      end if
      j = j - 2
      do i = i, n
        v(i) = TWO * ( node%w(j) * node%coseno(i) + node%w(j+1) * node%seno(i) )
        j = j - 2
      end do
      call ctimer( t2, tn, tn, gtime )
      dctii_mkl_time = t2-t1

      if( dcosqb_time.le.dctii_mkl_time ) then
        info = DftiFreeDescriptor( node%Desc_Handle )
        call DftierrorManager( info )
        tipo_dctii = 'fftpack'
        deallocate(node%w, node%coseno, node%seno, stat=info)
        if (info .ne. 0) then
           print *,'Error deallocating memory in ll_dctii_add_n_dctii'
           ll_dctii_add_n_dctii = -1
        end if
      else
        node%fftpack = .false.
        deallocate(node%Wreal, stat=info)
        if (info .ne. 0) then
           print *,'Error deallocating memory in ll_dctii_add_n_dctii'
           ll_dctii_add_n_dctii = -1
        end if
        tipo_dctii = 'mkl'
      end if
      !print '(a,f10.5)','tiempo dcosqb_time =    ',dcosqb_time
      !print '(a,f10.5)','tiempo dctii_mkl_time = ',dctii_mkl_time
      !print *,'tipo dctii = ',tipo_dctii

    end if ! tipo_dctii = auto
    deallocate(v, stat=info)
    if (info .ne. 0) then
       print *,'Error deallocating memory in ll_dctii_add_n_dctii'
       ll_dctii_add_n_dctii = -1
    end if
#endif

    node%next => ll%head%next
    ll%head%next => node

    ll_dctii_add_n_dctii = 0

  end function ll_dctii_add_n_dctii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  logical function ll_dctii_exists_n(ll, n)
    type(ll_dctii_list), intent(in), target :: ll
    integer, intent(in) :: n

    type(ll_dctii_node), pointer :: p_node
    
    ll_dctii_exists_n = .false.

    p_node => ll%head
    do while (associated(p_node%next))
       p_node => p_node%next
       if (p_node%n .eq. n) then
          ll_dctii_exists_n = .true.
          exit
       end if
    end do

  end function ll_dctii_exists_n

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function ll_dctii_get_n_node(ll,n)
    type(ll_dctii_node), pointer :: ll_dctii_get_n_node
    type(ll_dctii_list), intent(in), target :: ll
    integer, intent(in) :: n

    ll_dctii_get_n_node => ll%head
    do while (associated(ll_dctii_get_n_node%next))
       ll_dctii_get_n_node => ll_dctii_get_n_node%next
       if ( ll_dctii_get_n_node%n .eq. n) then
          return
       end if
    end do

    nullify(ll_dctii_get_n_node)

  end function ll_dctii_get_n_node

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

  integer function dctii_finalize()

    dctii_finalize =  ll_dctii_finalize(ll_init_sizes)

  end function dctii_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  integer function ll_dctii_finalize(ll)
    type(ll_dctii_list), intent(inout), target :: ll
    integer :: i, n, ret
    
    n = ll_dctii_length(ll)
    
    do i=1, n
       ret = ll_dctii_remove_node(ll,1)
       if (ret .ne. 0) then
          ll_dctii_finalize = -1
          return
       end if
    end do
    
    ll_dctii_finalize = 0
    
  end function ll_dctii_finalize
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ll_dctii_length(ll)
    type(ll_dctii_list), intent(in), target :: ll
    type(ll_dctii_node), pointer :: p_node
    integer :: n

    p_node => ll%head
    n=0
    do while (associated(p_node%next))
       p_node => p_node%next
       n = n+1
    end do

    ll_dctii_length = n

  end function ll_dctii_length

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function ll_dctii_remove_node(ll, i)
    type(ll_dctii_list), intent(inout), target :: ll
    integer, intent(in) :: i

    type(ll_dctii_node), pointer :: p_node, p_ant
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
       print *,'Error deallocating memory in ll_dctii_remove_node'
       ll_dctii_remove_node = -2
    end if

    ll_dctii_remove_node = 0

  end function ll_dctii_remove_node
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module sp_dctii
