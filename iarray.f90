module iarray_mod
implicit none
  save
contains

subroutine iarray!(xface,yface,zface,rhokap,MASK,PL_SUM)
  use packet_mod
  use grid_mod
  implicit none


  integer i,j,k,p,m

  !**** Initialize array values to be zero

  do i=1,nxg+1
     xface(i)=0.
  end do
  do i=1,nyg+1
     yface(i)=0.
  end do
  do i=1,nzg+1
     zface(i)=0.
  end do

  do i=1,nxg
     do j=1,nyg
        do k=1,nzg
           do p=1,4
              MASK(i,j,k,p)=0.
           enddo
           rhokap(i,j,k)=0.
           do m=1,nwl
              PL_SUM(m,I,J,K)=0.
              e_sum(m,i,j,k)=0.
              !e_dna_sum(m,i,j,k)=0.
           enddo

        end do
     end do
  end do

  do m=1,nwl
     n_phot_wl(m)=0
     l(m)=0.
  enddo
  do i=1,nlayer
     lcount(i)=0
  enddo


  return
end subroutine iarray
end module iarray_mod
