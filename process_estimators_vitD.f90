program process_estimators_vitD
use optical_properties_mod, only: nwl
use grid_mod, only: nxg,nyg,nzg, PL_SUM,zmax
use constants_mod
use vit_d_properties_mod
use vit_d_properties_init_mod

implicit none
real*8, dimension(nxg,nyg,nzg) :: j_mean, vitD
integer :: i,j,k,m
real*8 :: depth_val(nzg),rates, get_XYrate
real*8 :: percentage_conc(n_species,nxg,nyg,nzg)

real*8 :: timestep
real*8 ::AB,BA,BC,CB,BD,DB

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

  Ab=0
!this is taking an horifficallu long time. why?

       do i=1,nxg
       print*, 'on cell x....',i
      do j=1,nyg
           print*, 'on cell y....',j

        do k=1,nzg
             print*, 'on cell k....',k
          do m=1,nwl
          BA=BA+get_XYrate(m,extinction(2,m),quantum_yield(2,1,m))

          AB=AB+get_XYrate(m,extinction(1,m),quantum_yield(1,2,m))
          !BC=BC+get_XYrate(m,extinction(2,m),quantum_yield(2,3,m))
        !  CB=CB+get_XYrate(m,extinction(3,m),quantum_yield(3,2,m))
          !BD=BD+get_XYrate(m,extinction(2,m),quantum_yield(2,4,m))
        !  DB=DB+get_XYrate(m,extinction(4,m),quantum_yield(4,2,m))
         enddo
         vitD(i,j,k)=1.d0+(BA/AB)
      enddo
     enddo
   enddo
  !   open(10,file='depths_from_postprocess.dat',status='replace')
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


end program




function get_XYrate(m,extinction,yield)
  use constants_mod
  use grid_mod, only: pl_sum
  !returns the thing done to the bloody thing
 implicit none
  real*8 :: get_XYrate
  integer, intent(in):: m
  real*8, intent(in)::extinction, yield
  real*8 :: wl
  integer:: i,j,k



  wl = (real(m-1)+ 280.d0)*1.d-9


  do i=1,NXG
    do j=1,NYG
      do k=1,NZG

        get_XYrate=pl_sum(m,i,j,k)!*wl/(h*c))*extinction*yield

      enddo
    enddo
  enddo

 return
end function
