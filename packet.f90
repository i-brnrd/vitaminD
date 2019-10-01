module packet_mod
  implicit none
  save
  !contains
  integer, parameter :: nwl=121
  integer, parameter :: nlayer=5



  real*8  :: l(nwl)
  integer :: n_phot_wl(nwl)
  integer :: lcount(nlayer)


end module packet_mod
