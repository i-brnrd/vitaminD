!--------------------------------------------------------------------------------
  Program Estimators

  use grid_mod
  implicit none


  include 'photon.txt'
!  real*8,dimension(nxg,nyg,nxg)
   real*8 :: xmax,ymax,ZMAX
   real*8 fourpi
   real*8 lambdabins, n_spec, lum,
   real*8 :: nphotons, phcount







            read (10, file='plsum_test.dat', status='unknown', form='unformatted')
             write(10) PL_SUM
             close(10)
CALL PL_ESTIMATORS(ph_count,nphotons,fourpi,XMAX,YMAX,ZMAX,lambdabins,n_spec,lum,e_tot,e_uva,lumin,lambdabins)


   print*, 'hello', PL_SUM(1,100,100,100)
  end program Process_estimators
