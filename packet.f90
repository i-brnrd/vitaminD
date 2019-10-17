module packet_mod
use optical_properties_mod, only: nlayer 
  implicit none
  save

  !positional 
  real*8 xp,yp,zp !current position in grid
  real*8 nxp,nyp,nzp !current direction of packet 
  integer xcell,ycell,zcell

  !optical properties 
  real*8 ns,g,u_a(nlayer),u_s(nlayer)

  !HERE BE DRAGONS
  real*8 sint,cost,sinp,cosp,phi
  real*8 fi,fq,fu,fv


end module packet_mod
