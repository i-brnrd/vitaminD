subroutine density(x,y,z,zmax,layer)

  implicit none

  real*8 x,y,z

  real*8 sine,pi
  real*8 a,w,b_epi,b_meli,b_basal,b_sc
  real*8 zmax
  real*8 um_to_cm
  integer layer

  !******CALCULATE DISTANCES TO SET UP DENSITY
  !all distances in CM
  !**'wavelength' for sinosoidal undulation 0.05mm
  a=0.003  !amplitude
  w=0.0015!

  um_to_cm=0.0001

  sine=a*(SIN(x/w)*COS(y/w))
  !bc: 0.002cm, bepi: 0.008cm, b_meli: 0.009cm, bbasal:0.01cm
  !b_sc=0.048
  !b_epi=sine + 0.042!base epidermis
  !b_meli=sine+0.041!base melanin
  !b_basal=sine+0.040  !base basal

  !b_sc=zmax-0.002
!  b_epi=sine+(zmax-0.008)
  !b_meli=sine+(zmax-0.009)
  !b_basal=sine+(zmax-0.01)

  ! Epiderman thinkness Jane Sandy Moller ! Sc = 14.8 cellular epidermis= 83.7



  b_sc=zmax-(0.00148)
  b_epi=sine+(zmax-0.00637)
  b_meli=sine+(zmax-0.00737)
  b_basal=sine+(zmax-0.00837)




  if (z.ge.b_sc)then
     layer=1
  elseif (z.ge.b_epi) then
     layer=2
  elseif (z.ge.b_meli) then
     layer=3
  elseif (z.ge.b_basal)then
     layer=4
  else
     layer=5
  endif


  return
end subroutine density
