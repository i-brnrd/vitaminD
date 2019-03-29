module pl_estimators_mod
  implicit none
  save
contains

!--------------------------------------------------------------------------------
  SUBROUTINE PL_ESTIMATORS(ph_count,nphotons,XMAX,YMAX,ZMAX,wls,L,lumin)

    use grid_mod
    implicit none

    integer, intent(in) :: ph_count,nphotons

    real*8, intent(in) :: xmax,ymax,zmax
    real*8, intent(in) :: wls(nwl)
    real*8 :: lumin(nwl)
    real*8 :: A,V
    real*8 :: hc, L
    real*8, dimension(nxg,nyg,nzg) :: j_mean, fluence,energy
    real*8, dimension(nwl) :: e_basal,f_basal,toplayer
    real*8 :: wl, sum, proportion
    real*8 :: l_norm,e_norm
    real*8 :: l_sum
    integer :: i,j,k,m
    integer :: n_basal
    real*8 :: depth(6)
    integer:: kk
    real*8 :: fluence_depth(nwl)
    real *8 :: flu_per(2)
    real*8 :: depth_val(nzg)
    integer:: basal
    real*8:: e_tot, e_uva
    !Source Luminosity & Illumination Area


    print*, nxg,nyg,nzg


    print*, l


    print*, '****'
    !print*, lumin
    A=(2.d0*(xmax))*(2.d0*(ymax))

  !  V = (2.*(XMAX)/NXG)**3  !Each cell has the same volume
    V= (2.d0*(XMAX/NXG))*(2.d0*(YMAX/NYG))*(2.d0*(ZMAX/NZG))
    print*,'Area',A
    !L=L/(100.d0*100.d0) !W/cm
    print*,'Luminosity',L,'W/m^2'
    print*,'Volume of Voxel',v
    print*, 'cells', nxg,nyg,nzg

    !print*, 'wavelength dependant luminosity', lumin

    !lumin=lumin/(100.d0*100.d0)

    hc=(6.626d0*10.d0**(-34))*(3.d0*10.d0**8) !SI units

    basal=4

    ! initialise arrays & counters
    n_basal=0
    e_tot=0.
    e_uva=0.
    flu_per=0.
    l_sum=0.
    DO I=1,NXG
       DO J=1,NYG
          DO K=1,NZG
             energy(i,j,k)=0.0
             j_mean(i,j,k)=0.0
             if (MASK(i,j,k,basal)==1) then
                n_basal=n_basal+1
             endif
          ENDDO
       ENDDO
    ENDDO
    do m=1,nwl
       e_basal(m)=0.0
       f_basal(m)=0.0
       toplayer(m)=0.0
    enddo

    !start path length counters
    print*, nwl

    do m=1,nwl
      l_sum=l_sum+lumin(m)
       if (n_phot_wl(m).eq.0) cycle
       l_norm=lumin(m)*a/(real(n_phot_wl(m))*v)

       !print*, l_norm ,lumin(m),a,n_phot_wl(m),v
       e_norm=l_norm*(10.**(-9)/hc)
       do i=1,nxg
          do j=1,nyg
             do k=1,nzg
                j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))*l_norm
                energy(i,j,k)=energy(i,j,k) + (E_SUM(m,i,j,k))*e_norm
                if (k==96) then
                  fluence_depth(m)=fluence_depth(m) + (PL_SUM(m,i,j,k))*l_norm
                end if
                if(k==101) THEN
                  toplayer(m)=toplayer(m)+(PL_SUM(m,i,j,k))*l_norm
                endif
                if (MASK(i,j,k,basal)==1) then
                   e_basal(m)= e_basal(m) + (E_SUM(m,i,j,k)*e_norm)
                 f_basal(m)= f_basal(m) + PL_SUM(m,i,j,k)*l_norm
              endif
             enddo
          enddo
       enddo

    enddo

    print*,l_sum


    ! print out values at each layer
   open(10,file='./plots/depths.dat')
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
       print*, (real(nzg-k)*(2.d0*zmax/real(nzg)))/(1.E-4), depth_val(k)/L_sum
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
