module optical_properties_mod
  implicit none
  save
  !contains
  real*8, parameter :: wl_start=280.d0
  integer, parameter :: nwl=1!21 !121 uv wavelengths 
  integer, parameter :: nlayer=5

  real*8  :: l(nwl)
  integer :: n_pkt_wl(nwl)
  integer :: lcount(nlayer)

  real*8,dimension(nlayer,nwl):: u_s_all,u_a_all
  real*8,dimension(nwl):: g_all,n_all
end module optical_properties_mod
