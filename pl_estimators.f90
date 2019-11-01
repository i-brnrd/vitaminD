module pl_estimators_mod
  implicit none
  save
contains


!what I want from this module
!obviously a hugmungoes tidy
! I want OPTION to write out full pl_ grid of wls.

!ando ok - idea here; yes can write out the integrated depth plot too just as an itional checker
!--------------------------------------------------------------------------------
  SUBROUTINE PL_ESTIMATORS(pkt_count,incident_spec_irr)
    use optical_properties_mod, only: n_pkt_wl, l, nwl
    use grid_mod, only: nxg,nyg,nzg,xmax,ymax,zmax, PL_SUM
    use constants_mod
    implicit none

    integer, intent(in) :: pkt_count
    real*8, intent(in) :: incident_spec_irr(nwl)

    real*8 :: Area

    !unnecessary
    real*8 :: wls(nwl)
    real*8 :: lumin(nwl)

    real*8 :: V
    real*8, dimension(nxg,nyg,nzg) :: j_mean
    real*8 :: l_norm, l_sum
    integer :: i,j,k,m
    real*8 :: depth_val(nzg)



    real*8 :: answer



    print*, nxg,nyg,nzg


    wls=l
    lumin=incident_spec_irr

    l_sum=sum(incident_spec_irr)

    print*, '****'
    !print*, lumin
    Area=(2.d0*(xmax))*(2.d0*(ymax))
    V= (2.d0*(XMAX/NXG))*(2.d0*(YMAX/NYG))*(2.d0*(ZMAX/NZG))

    print*,'Area',Area, ' in cm**2'
    !L=L/(100.d0*100.d0) !W/cm
    !gprint*,'Luminosity',L,'W/m^2'
    print*,'Volume of Voxel',v
    print*, 'cells', nxg,nyg,nzg

    j_mean=0.d0
    do m=1,nwl
      if (n_pkt_wl(m).eq.0) cycle
      l_norm=lumin(m)*area/(real(n_pkt_wl(m))*v)
      do i=1,nxg
        do j=1,nyg
          do k=1,nzg
            j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))*l_norm
          enddo
        enddo
      enddo
    enddo

!I would like to write out a dataset that is pl_sum* l_norm for each wavelength i think; lets try it
!And then it would be tis dataset I interrogate and load up and do stuff with
!two options: one is write aout a humongously long file.
!another: qrite out the unformattted array BUT each pl_sum entry has already been normed
!lets try that one

!   do i=1,nxg
!     do j=1,nyg
!       do k=1,nzg
!         do m=1,nwl
!           l_norm=lumin(m)*area/(real(n_pkt_wl(m))*v)
!           if (n_pkt_wl(m).eq.0) then
!             answer=0.d0
!           else
!             answer=(PL_SUM(m,i,j,k))*l_norm
!           write()
!
!         j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))*l_norm
!       enddo
!     enddo
!   enddo
! enddo

    j_mean=0.d0
    do m=1,nwl
      if (n_pkt_wl(m).eq.0) cycle
      l_norm=lumin(m)*area/(real(n_pkt_wl(m))*v)
      do i=1,nxg
        do j=1,nyg
          do k=1,nzg
            j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))*l_norm
          enddo
        enddo
      enddo
    enddo

    open(10,file='./plots/depths1.dat',status='replace')

    depth_val=0.0d0
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !print*, PL_SUM(1,i,j,k)
          depth_val(k)= depth_val(k) + j_mean(i,j,k)
        enddo
      enddo
      depth_val(k)=depth_val(k)/(nxg*nyg)
      write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)/L_sum
    !print*, (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)/L_sum
    enddo

    close(10)



      do m=1,nwl
        if (n_pkt_wl(m).eq.0) then
          PL_SUM(m,:,:,:)=0.d0
        else
          l_norm=lumin(m)*area/(real(n_pkt_wl(m))*v)
          PL_SUM(m,:,:,:)=l_norm*PL_SUM(m,:,:,:)
        endif
      enddo



          j_mean=0.d0
          do m=1,nwl
            if (n_pkt_wl(m).eq.0) cycle
            do i=1,nxg
              do j=1,nyg
                do k=1,nzg
                  j_mean(i,j,k)=j_mean(i,j,k) + (PL_SUM(m,i,j,k))
                enddo
              enddo
            enddo
          enddo





      print*, '...ABOUT TO WRITE TO THE UFORMATTED ARRAY......'
      open(10,file='./pathlengths.dat',status='replace',form='unformatted')
        write(10) PL_SUM
      close(10)
      print*, 'written it out!!!!'


    !I would still like to do this.....
    ! print out values at each depth
    open(10,file='./plots/depths2.dat',status='replace')

    depth_val=0.0d0
    do k=1,nzg
      do j=1,nyg
        do i=1,nxg
          !print*, PL_SUM(1,i,j,k)
          depth_val(k)= depth_val(k) + j_mean(i,j,k)
        enddo
      enddo
      depth_val(k)=depth_val(k)/(nxg*nyg)
      write(10,*) (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)/L_sum
    !print*, (real(nzg-k)*(2.d0*zmax/real(nzg))), depth_val(k)/L_sum
    enddo

    close(10)

    RETURN
  END SUBROUTINE PL_ESTIMATORS

end module pl_estimators_mod
