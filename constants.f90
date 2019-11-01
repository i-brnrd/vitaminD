module constants_mod
  implicit none
  save

  real*8, parameter :: pi = 4.*atan(1.)
  real*8, parameter:: twopi=2.*pi
  real*8, parameter:: fourpi=4.*pi

  real*8, parameter :: h=6.626D-34 !Planck
  real*8, parameter :: c=3.D+8      !Speed of light
  real*8, parameter :: Na=6.022D+23  !Avogadro

end module constants_mod
