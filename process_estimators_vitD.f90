program process_estimators_vitD
use optical_properties_mod, only: nwl
use grid_mod, only: nxg,nyg,nzg, PL_SUM,zmax
use constants_mod
use vit_d_properties_mod
use vit_d_properties_init_mod

implicit none
real*8, dimension(nxg,nyg,nzg) :: j_mean, vitD
integer :: i,j,k,m
real*8 :: depth_val(nzg),rates
real*8 :: percentage_conc(n_species,nxg,nyg,nzg)

real*8 :: timestep

PL_SUM=0.d0
call vit_d_properties_init()
!read in the pathlengths file

open(10, file='pathlengths.dat',status='unknown',form='unformatted')
    read(10) PL_SUM
close(10)


 print*, c, h, Na

 percentage_conc = 0.d0 !molar cncentrations of each substance.....


 timestep = 1.d0

 !SET   of equations for the vitamin D formation stuff
 !spectral relative conversin rates???
 !substances are numbered as follows:
 !------------1: provitamin D, 2: previtamin d, 3: lumisterol,
 !------------4: tachysterol 5: vitamin D 6: any toxi/supra sterol or other decay product

 do m=1,nwl
     do i=1,nxg
      do j=1,nyg
        do k=1,nzg
          vitD(i,j,k)=vitD(i,j,k) + (PL_SUM(m,i,j,k))*rates(m)
         enddo
      enddo
     enddo
   enddo


  !   open(10,file='depths_from_postprocess.dat',status='replace')
!
     depth_val=0.0d0
    do k=1,nzg
       do j=1,nyg
         do i=1,nxg
           !print*, PL_SUM(1,i,j,k)
           depth_val(k)= depth_val(k) + vitD(i,j,k)
         enddo
       enddo
       depth_val(k)=depth_val(k)/(nxg*nyg)
       print*, depth_val(k)
!       write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)
      enddo




percentage_conc(1,:,:,:)=1.d0
  do m=1,nwl
    print*, quantum_yield(1,2,m), extinction(1,m)
    percentage_conc(2,:,:,:)=percentage_conc(2,:,:,:)+(extinction(1,m)*quantum_yield(1,2,m)*PL_SUM(m,:,:,:))
  enddo

  print*, sum(percentage_conc(2,:,:,:))


end program


function rates(wl)
  use vit_d_properties_mod

  !returns the distance along a given direction to a cell wall
  !given the positions of both cell walls
 implicit none
  real*8 :: rates
  integer,intent(in) :: wl

   real*8 :: AB,BC,BD
   AB=(extinction(2,wl)*quantum_yield(2,1,wl))/(extinction(1,wl)*quantum_yield(1,2,wl))
   BC=(extinction(2,wl)*quantum_yield(2,3,wl))/(extinction(3,wl)*quantum_yield(3,2,wl))
   BD=(extinction(2,wl)*quantum_yield(2,4,wl))/(extinction(4,wl)*quantum_yield(4,2,wl))

   rates=AB+BC+BD
 return
end function
!
! !verification thing shave been read in ok : learn how to do a loop in one line bitch
! j_mean=0.d0
!   do m=1,nwl
!     do i=1,nxg
!       do j=1,nyg
!         do k=1,nzg
!           j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))
!         enddo
!       enddo
!     enddo
!   enddo
!     open(10,file='depths_from_postprocess.dat',status='replace')
!
!     depth_val=0.0d0
!     do k=1,nzg
!       do j=1,nyg
!         do i=1,nxg
!           !print*, PL_SUM(1,i,j,k)
!           depth_val(k)= depth_val(k) + j_mean(i,j,k)
!         enddo
!       enddo
!       depth_val(k)=depth_val(k)/(nxg*nyg)
!       print*, depth_val(k)
!       write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)
!
!     enddo
