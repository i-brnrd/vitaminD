module gridset_mod
  implicit none
  save
contains


  subroutine gridset(kappa)!(xface,yface,zface,MASK,xmax,ymax,zmax,kappa)
    use optical_properties_mod
    use grid_mod, ONLY: nxg,nyg,nzg,xface,yface,zface,MASK,xmax,ymax,zmax,delta

    implicit none

    real*8 kappa
    real*8 x,y,z,rho,taueq,taupole

    integer i,j,k,p,layer


    delta=1.e-6*(2.*zmax/nzg) !small roundoff error number for when close to cell walls

    print *, 'Setting up density grid....'

    !**********  Linear Cartesian grid. Set up grid faces ****************
    !do i=1,nxg+1
      do i=-3,nxg+4
       xface(i)=(i-1)*2.*xmax/nxg
    end do
    !do i=1,nyg+1
      do i=-3,nyg+4
       yface(i)=(i-1)*2.*ymax/nyg
    end do
    !do i=1,nzg+1
      do i=-3,nzg+4
       zface(i)=(i-1)*2.*zmax/nzg
    end do
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



!!$
!!$
!!$****************** Calculate equatorial and polar optical depths ****
!!$      taueq=0.
!!$      taupole=0.
!!$     do i=1,nxg
!!$         taueq=taueq+rhokap(i,nyg/2,nzg/2)
!!$      enddo
!!$      do i=1,nzg
!!$        taupole=taupole+rhokap(nxg/2,nyg/2,i)
!!$      enddo
!!$      taueq=taueq*2.*xmax/nxg
!!$      taupole=taupole*2.*zmax/nzg
!!$      print *,'taueq = ',taueq,'  taupole = ',taupole
!!$
!!$      ************** Write out density grid as unformatted array
!!$      open(10,file='density.dat',status='unknown',form='unformatted')
!!$           write(10) rhokap
!!$      close(10)
!!$
!!$     open(10,file='density_slice.dat',status='unknown')
!!$      do k=1,4
!!$     do i=1,nxg
!!$       write(10,*) (rhokap(i,j,101,k), j=1,nyg)
!!$     end do
!!$     enddo
!!$
!!$      do i=1,nxg
!!$       do j=1,nyg
!!$           x=xface(i)-xmax+xmax/nxg
!!$           y=yface(j)-ymax+ymax/nyg
!!$         write(10,*) x,y,rhokap(i,j,101)
!!$
!!$       end do
!!$      end do

    return
  end subroutine gridset
end module gridset_mod
