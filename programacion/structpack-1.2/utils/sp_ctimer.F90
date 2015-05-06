module sp_ctimer
  use iso_c_binding, only: C_DOUBLE
  implicit none
  public :: ctimer
  interface
    subroutine ctimer(elapsed,ucpu,scpu,gtime) bind(C, name="ctimer_")
      use iso_c_binding, only: c_double
      real(kind=c_double) :: elapsed, ucpu, scpu, gtime
    end subroutine ctimer
  end interface
end module sp_ctimer
