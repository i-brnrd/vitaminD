program process_estimators_vitD
use optical_properties_mod, only: nwl
use grid_mod, only: nxg,nyg,nzg, PL_SUM,zmax
use constants_mod
use vit_d_properties_mod
use vit_d_properties_init_mod

implicit none
real*8, dimension(nxg,nyg,nzg) :: j_mean, vitD
integer :: i,j,k,m
real*8 :: depth_val(nzg)
real*8 :: percentage_conc(n_species,nxg,nyg,nzg)

real*8 :: spec_depth(nwl)

real*8 :: pls,wl
real*8 :: timestep
real*8 ::AB,BA,BC,CB,BD,DB
real*8 :: r(6)

PL_SUM=0.d0
call vit_d_properties_init()
!read in the pathlengths file

open(10, file='pathlengths.dat',status='unknown',form='unformatted')
    read(10) PL_SUM
close(10)


print*, c, h, Na
print*, 1.d-9

 percentage_conc = 0.d0 !molar cncentrations of each substance.....


 timestep = 1.d0

 !SET   of equations for the vitamin D formation stuff
 !spectral relative conversin rates???
 !substances are numbered as follows:
 !------------1: provitamin D, 2: previtamin d, 3: lumisterol,
 !------------4: tachysterol 5: vitamin D 6: any toxi/supra sterol or other decay product

r=0.d0
!this is taking an horifficallu long time. why?

 do i=1,nxg

    do j=1,nyg

        do k=1,nzg
          r=0.d0
             do m=1,nwl
              wl = (real(m-1)+ 280.d0)*1.d-9
              pls=pl_sum(m,i,j,k)
              r(1)=r(1)+ (pls*wl/(h*c))* extinction(2,m)*quantum_yield(2,1,m)
              r(2)=r(2)+(pls*wl/(h*c)) * extinction(1,m)*quantum_yield(1,2,m)
              r(3)=r(3)+(pls*wl/(h*c)) * extinction(2,m)*quantum_yield(2,3,m)
              r(4)=r(4)+(pls*wl/(h*c)) * extinction(3,m)*quantum_yield(3,2,m)
              r(5)=r(5)+(pls*wl/(h*c)) * extinction(2,m)*quantum_yield(2,4,m)
              r(6)=r(6)+(pls*wl/(h*c)) * extinction(4,m)*quantum_yield(4,2,m)


       enddo
         vitD(i,j,k)=1.d0/(1.d0+(r(1)/r(2))+(r(3)/r(4))+(r(5)/r(6)))
      enddo
     enddo
   enddo
  !   open(10,file='depths_from_postprocess.dat',status='replace')
     depth_val=0.0d0
open(8, file='vitDdepth.txt',status='unknown')
    do k=1,nzg
       do j=1,nyg
         do i=1,nxg
           !print*, PL_SUM(1,i,j,k)
           depth_val(k)= depth_val(k) + vitD(i,j,k)
         enddo
       enddo

       depth_val(k)=depth_val(k)/(nxg*nyg)
       print*, depth_val(k)
        write(8,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)
        print*, k, (real(nzg-k)*(2.d0*zmax/real(nzg)))
!       write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)
enddo
close(8)
open(8, file='specdepth.txt',status='unknown')
spec_depth=0.d0
k=nxg-40
  do m=1,nwl
   do j=1,nyg
     do i=1,nxg
        spec_depth(m)= spec_depth(m) + PL_SUM(m,i,j,k)
      enddo
    enddo
          write(8,*) m+279, spec_depth(m)
  enddo


close(8)















end program

!
! function get_XYrate(m,pls,extinction,yield)
!   use constants_mod
!   !returns the thing done to the bloody thing
!  implicit none
!   real*8 :: get_XYrate
!   integer, intent(in):: m
!   real*8, intent(in)::pls(1,1,1,1), extinction, yield
!   real*8 :: wl
!
!   wl = (real(m-1)+ 280.d0)*1.d-9
!
!   get_XYrate=pls*wl/(h*c))*extinction*yield
!
!  return
! end function
