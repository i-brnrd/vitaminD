module optical_properties_mod
  implicit none
  save
  !contains
  real*8, parameter :: wl_start=280.d0
  integer, parameter :: nwl=121
  integer, parameter :: nlayer=5

  real*8  :: l(nwl)
  integer :: n_phot_wl(nwl)
  integer :: lcount(nlayer)

  real*8,dimension(nlayer,nwl):: u_s,u_a
  real*8,dimension(nwl):: g_skin,n_skin
end module optical_properties_mod
