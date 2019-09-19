!--------------------------------------------------------------------------------
  Program Process_estimators

  use grid_mod
  use pl_estimators_mod
  implicit none


  ! Idea here is to post process from J mean....

  include 'photon.txt'
!  real*8,dimension(nxg,nyg,nxg)
   real*8 :: xmax,ymax,ZMAX
   real*8 lambdabins, n_spec, l, lumin(NWL),lum
   real*8 :: nphotons
   character*30 fname,fname2
   integer :: ph_count

!Import jmean.....

    print*, 'hello'

    open(10, file='nb_uvb/j_mean.dat', status='unknown', form='unformatted')
    read(10) PL_SUM

        open(11,file='j_mean.dat',status='replace')
        open(12,file='energy.dat',status='replace')

        read(11,*) j_mean
        read(12,*) energy

        close(11)
        close(12)
        print*,l_sum

! already calulated the jmean path length estimator.
! Check hhave i only importred jmean?


! print out values at each layer
open(10,file='./plots/depths.dat',status='replace')
depth_val=0.0d0
do k=1,nzg
   do j=1,nyg
      do i=1,nxg
          !print*, PL_SUM(1,i,j,k)
          depth_val(k)= depth_val(k) + j_mean(i,j,k)
      enddo
   enddo
   depth_val(k)=depth_val(k)/(nxg*nyg)
   write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg)))/(1.E-4), depth_val(k)/L_sum
   !print*, (real(nzg-k)*(2.d0*zmax/real(nzg)))/(1.E-4), depth_val(k)/L_sum
enddo

close(10)



    close(10)








  end program Process_estimators







    ! print out values at each layer
   open(10,file='./plots/depths.dat',status='replace')
    depth_val=0.0d0
    do k=1,nzg
       do j=1,nyg
          do i=1,nxg
              !print*, PL_SUM(1,i,j,k)
              depth_val(k)= depth_val(k) + j_mean(i,j,k)
          enddo
       enddo
       depth_val(k)=depth_val(k)/(nxg*nyg)
       write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg)))/(1.E-4), depth_val(k)/L_sum
       !print*, (real(nzg-k)*(2.d0*zmax/real(nzg)))/(1.E-4), depth_val(k)/L_sum
    enddo

   close(10)



    stop
   !open(11,file='./plots/borehole.txt',status='replace')
  !write(11,*) 0.00, 1.0

    do i=1,nzg
      print*, zface(i), depth_val(i)
      !write(11,*) zface(i),depth_val(i)/l
    enddo
    stop
    close(11)

    !SLICE plots at mid-y axisnal (at the very least for some sound advice).
    open(10,file='./plots/flu3d.dat',status='replace')
    open(11,file='./plots/e3d.dat',status='replace')
    open(12,file='./plots/basal.dat',status='replace')

    j=50
    do i=1,nxg
       write(10,*)
       write(11,*)
       write(12,*)
       do k=1,nzg
          write(10,*) i,k,fluence(i,j,k)
          write(11,*) i,k,energy(i,j,k)
          if (MASK(i,j,k,basal)==1) then
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




  close (9)
    RETURN
  END SUBROUTINE PL_ESTIMATORS

  end module pl_estimators_mod
