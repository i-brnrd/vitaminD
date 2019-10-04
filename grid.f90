module grid_mod
  use optical_properties_mod, only:nwl,nlayer
  implicit none
  save
  !contains
  integer, parameter :: nxg=50,nyg=50,nzg=400
  !real*8, parameter ::

  real*8 :: xface(-3:nxg+4),yface(-3:nyg+4),zface(-3:nzg+4)
  real*8 :: MASK(nxg,nyg,nzg,nlayer)
  real*8 :: rhokap(nxg,nyg,nzg)

  REAL*8 :: PL_SUM(nwl,NXG,NYG,NZG)

  real*8 :: xmax,ymax,zmax
  real*8 :: delta

end module grid_mod
