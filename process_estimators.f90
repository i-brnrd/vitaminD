program process_estimators
use optical_properties_mod, only: nwl
use grid_mod, only: nxg,nyg,nzg, PL_SUM,zmax
use constants_mod
use vit_d_properties_mod
use vit_d_properties_init_mod

implicit none
real*8, dimension(nxg,nyg,nzg) :: j_mean
integer :: i,j,k,m
real*8 :: depth_val(nzg)



PL_SUM=0.d0



call vit_d_properties_init()


stop

!read in the pathlengths file

    open(10, file='pathlengths.dat',status='unknown',form='unformatted')
     read(10) PL_SUM
    close(10)


    do i=1,nxg
      do j=1,nyg
        do k=1,nzg
          j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))
        enddo
      enddo
    enddo




    open(10,file='depths_from_postprocess.dat',status='replace')

    depth_val=0.0d0
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !print*, PL_SUM(1,i,j,k)
          depth_val(k)= depth_val(k) + j_mean(i,j,k)
        enddo
      enddo
      depth_val(k)=depth_val(k)/(nxg*nyg)
      write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)

    enddo

    close(10)
end program
