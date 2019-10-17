module gridset_mod
  implicit none
  save
contains


  subroutine gridset(kappa)!(xface,yface,zface,MASK,xmax,ymax,zmax,kappa)
    use optical_properties_mod
    use grid_mod, ONLY: nxg,nyg,nzg,xface,yface,zface,MASK,xmax,ymax,zmax,delta,grid_max

    implicit none

    real*8 kappa
    real*8 x,y,z,rho,taueq,taupole

    integer i,j,k,p,layer


    delta=1.e-6*(2.*zmax/nzg) !small roundoff error number for when close to cell walls

    print *, 'Setting up density grid....'

    !**********  Linear Cartesian grid. Set up grid faces.. ****************
    !do i=1,nxg+1
      do i=-3,nxg+4
       xface(i)=(i-1)*2.d0*xmax/nxg
    end do
    !do i=1,nyg+1
      do i=-3,nyg+4
       yface(i)=(i-1)*2.d0*ymax/nyg
    end do
    !do i=1,nzg+1
      do i=-3,nzg+4
       zface(i)=(i-1)*2.d0*zmax/nzg
    end do
    !FIX xmax, ymax,ymas to be exactlt the same as relevaent face... 
    grid_max(1)=xface(nxg+1)
    grid_max(2)=yface(nyg+1)
    grid_max(3)=zface(nzg+1)

    !**********  Linear Cartesian grid. Set up grid faces ****************

    !**************  Loop through x, y, and z to set up grid density.  ****
    do i=1,nxg
       do j=1,nyg
          do k=1,nzg
             x=xface(i)-xmax+xmax/nxg
             y=yface(j)-ymax+ymax/nyg
             z=zface(k)-zmax+zmax/nzg


             !**********************Call density setup subroutine
             call density(x,y,z,zmax,layer)

             do p=1,nlayer
                if (layer.eq.p) then
                   MASK(i,j,k,p)=1 !mask
                   lcount(p)=lcount(p)+1 !counter fornumber of cells in a layer
                else
                   MASK(i,j,k,p)=0
                endif
             enddo
          end do
       end do
    end do


    return
  end subroutine gridset
end module gridset_mod
