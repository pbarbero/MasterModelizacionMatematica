!> \brief Module for the computation of the product y = T*x with a convolution
!> \author StructPack team
!> \date April, 2012
module sp_conv
  use sp_utils   ! exports kind value r_kind
  use sp_rft
  
  implicit none
  
  public :: conv
  private

  contains
  
  !> \brief Performs the matrix-vector product
  !>              y = beta*y + alfa*T*x
  !> where T is a toeplitz matrix (currently should be square).
  !> \param[in]     alpha Scalar alpha. 
  !> \param[in]     t     First column of the Toeplitz matrix.
  !> \param[out]    x     Vector x.
  !> \param[in]     beta  Scalar beta. 
  !> \param[inout]  y     Vector y. 
  !> \param[in]     u_    (optional) If present, first row of the Toeplitz matrix.
  subroutine conv( alpha, t, x, beta, y, u_ ) 
    
    !     conv performs the product T*x by a convolution operation
    !
    real(kind=r_kind), dimension(:), intent(in)           :: t, x
    real(kind=r_kind), dimension(:), intent(in), optional :: u_
    real(kind=r_kind), dimension(:), intent(inout)        :: y
    real(kind=r_kind),               intent(in)           :: alpha, beta

    real(kind=r_kind), dimension(:), allocatable :: rft_t, rft_x

    integer :: m, n
    logical :: non_symmetric
    
    integer :: o, info, i
    real(kind=r_kind) :: r, s

    m = size(t,1)
    non_symmetric = present( u_ )
    
    if( non_symmetric ) then
      n = size(u_,1)
      if( size(x,1).ne.n .or. size(y,1).ne.m ) then
        print *,'conv: Dimensions error'
        return
      end if
    end if

    if( beta.ne.0.0_r_kind .and. beta.ne.1.0_r_kind ) then
      y = beta * y
    end if

    if( alpha.ne.0.0_r_kind ) then
      o = 2
        do while( o.lt.(2*m-1) ) 
        o = 2 * o
      end do
      allocate( rft_t(o), rft_x( o ), stat=info)
      if (info .ne. 0) then
         print *,'Error allocating memory in conv'
         return
      end if      
      if( non_symmetric ) then
        rft_t = (/t,(ZERO,i=0,o-2*m),u_(m:2:-1)/)
      else
        rft_t = (/t,(ZERO,i=0,o-2*m),t(m:2:-1)/)
      end if
      call rft( rft_t )
      rft_x = (/x,(ZERO,i=1,o-m)/)
      call rft( rft_x )
      !product of arrays
      if( non_symmetric ) then 
        rft_t(1) = rft_t(1) * rft_x(1)
        do i = 2, o-1, 2
          r = rft_t( i )
          s = rft_t( i+1 )
          rft_t( i )   = r * rft_x( i )   - s * rft_x( i+1 )
          rft_t( i+1 ) = r * rft_x( i+1 ) + s * rft_x( i )
        end do
        if( mod(o,2).eq.0 ) rft_t(o) = rft_t(o) * rft_x(o)
      else
        rft_t(3:o:2) = rft_t(2:o:2)
        rft_t = rft_t * rft_x
      end if
      call irft( rft_t )
      r = alpha / o
      if( beta.eq.0.0_r_kind ) then
        y = r * rft_t(1:m) 
      else
        y = y + r * rft_t(1:m) 
      end if
    end if

    deallocate( rft_t, rft_x, stat=info)

    return
    
  end subroutine conv
  
end module sp_conv
