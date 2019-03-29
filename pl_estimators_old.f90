module pl_estimators_mod
  implicit none
  save
contains
  !NONSENSE 
  
!--------------------------------------------------------------------------------
  SUBROUTINE PL_ESTIMATORS(ph_count,nphotons,fourpi,XMAX,wls,n_spec,L,e_tot,e_uva)

  use grid_mod
  implicit none
 
  integer, intent(in) :: ph_count,nphotons,n_spec
  real*8, intent(in) :: fourpi, xmax
  real*8, intent(in) :: wls(n_spec), L

  real*8 :: A,V

  real*8 :: hc 

  real*8, dimension(nxg,nyg,nzg) :: j_mean, fluence,energy
  real*8, dimension(nwl) :: e_basal,f_basal

  real*8 :: wl, sum, proportion
  real*8 :: l_norm,e_norm
  integer :: i,j,k,m
  integer :: n_basal

  real*8, intent(out) :: e_uva, e_tot 
  !Source Luminosity & Illumination Area 

  A=1.*(0.1**2) ! 1mm^2 =cm^2
  V = (2.*(XMAX)/NXG)**3
  
  print*,'Area',A
  print*,'Luminosity',L,'W/m^2'
  print*,'Volume of Voxel',v

  hc=(6.626*10.**(-34))*(3.*10.**8) !SI units

  ! initialise arrays & counters 
  DO I=1,NXG
     DO J=1,NYG
        DO K=1,NZG
           energy(i,j,k)=0.0
           j_mean(i,j,k)=0.0
        ENDDO
     ENDDO
  ENDDO
  do m=1,nwl
     e_basal(m)=0.0
     f_basal(m)=0.0
  enddo
  
  n_basal=0
  e_tot=0.
  e_uva=0.

  !write out sum over all wavelengths to array; and to basal layer only 
  DO I=1,NXG
     DO J=1,NYG
        DO K=1,NZG
           
           if (MASK(i,j,k,3)==1) then
              n_basal=n_basal+1
           endif
           
           do m=1,nwl
              j_mean(i,j,k)=j_mean(i,j,k)+PL_SUM(m,i,j,k)
              energy(i,j,k)=energy(i,j,k)+E_SUM(m,i,j,k)
              
              if (MASK(i,j,k,3)==1) then
                 e_basal(m)= e_basal(m) + (E_SUM(m,i,j,k))*0.0185
                 f_basal(m)= f_basal(m) + PL_SUM(m,i,j,k)
                 
              endif
           ENDDO


           
        ENDDO
     ENDDO
  ENDDO
  !Normalising

  l_norm=l*a/(nphotons*v)
  e_norm=l_norm*(10.**(-9)/hc)

  energy=e_norm*energy
  e_basal=e_norm*e_basal
  fluence=l_norm*j_mean
  f_basal=l_norm*f_basal


  do i=1,nwl
     l_norm=l*a/(n_phot_wl(i)*v)
     e_norm=l_norm*(10.**(-9)/hc)

     energy

     
  enddo


  
  
  !SLICE plots at mid-y axis
  open(10,file='./plots/flu3d.dat',status='replace')
  open(11,file='./plots/e3d.dat',status='replace')
  open(12,file='./plots/basal.dat',status='replace')

  j=50
  do i=1,nyg
     write(10,*)
     write(11,*)
     write(12,*)
     do k=1,nzg
        write(10,*) i,k,fluence(i,j,k)
        write(11,*) i,k,energy(i,j,k)
        if (MASK(i,j,k,3)==1) then
           write(12,*) i,k,energy(i,j,k)
        else
           write(12,*) i,k,0.0
        endif
     enddo
  enddo
  close(10)
  close(11)
  close(12)

  !AVERAGE FLUENCE at DEPTH
  OPEN(9,FILE="./plots/fl_depth.dat",status='replace')
  Do k=1,NZG
     sum = 0.
     DO I=1,NxG
        DO J=1,NYG
           sum= sum+fluence(i,j,k)
        enddo
     enddo
     sum = sum/((real(nyg))*(real(nxg))) !
     WRITE(9,*) k,sum
  ENDDO
  close (9)


  !average no UVA photons absorbed per  voxel...
  e_tot=0.
  e_uva=0.
  open(12,file='./plots/basal_e.dat',status='replace')
  open(13,file='./plots/basal_f.dat',status='replace')
  open(14,file='./plots/reproduced_spec.dat',status='unknown')
 

  do m=1,nwl
     wl = wls(1) + (real(m)*(wls(n_spec)-wls(1))/nwl)
     write(12,*) wl, e_basal(m)
     write(13,*) wl, f_basal(m)
     
     write(14,*) wl, n_phot_wl(m)
     
     e_tot=e_tot+e_basal(m)
     if (wl.gt.315) then
        e_uva=e_uva+e_basal(m)
     endif
  enddo

  close(12)
  close(13)
  close(14)
  
  print*,'number of basal cells', n_basal

  e_tot=(e_tot/real(n_basal))*(10.**(-12))
  !e_tot=e_tot*10.**(-15)
  !e_uva=e_uva*10.**(-15)
  e_uva=(e_uva/real(n_basal))*(10.**(-12))


  RETURN
END SUBROUTINE PL_ESTIMATORS

end module pl_estimators_mod




!!$  print*,'Total UVB Photons Per Basal cell',e_tot-e_uva
!!$  print*,'Total UVA Photons Per basal cell',e_uva
!!$  print*,'Total Photons Absorbed per cell',e_tot
!!$  
!!$  print*,'CPDs: UVB',(e_tot-e_uva)*0.05
!!$  print*,'CPDs: UVA',e_uva*0.0005
!!$  print*,'CPDs: total', e_uva*0.0005 + (e_tot-e_uva)*0.05 

!!$
!!$  OPEN(9,FILE="./plots/spec_depth.dat",status='replace')
!!$
!!$  do m=1,nwl
!!$     wl = lam(1) + (real(m)*(lam(n_spec)-lam(1))/nwl)
!!$     ii=0
!!$     DO K=101,77,-4
!!$        ii=ii+1
!!$        sum =0 
!!$        DO I=1,NXG
!!$           DO J=1,NYG
!!$              sum=sum+E_SUM(m,i,j,k)
!!$           ENDDO
!!$        ENDDO
!!$        sum = sum/((real(nyg))*(real(nxg)))
!!$        depth(ii)=sum
!!$     ENDDO
!!$     write(9,*) wl, depth(1), depth(2),depth(3),depth(4),depth(5),depth(6)
!!$  enddo
!!$
!!$  close (9)
