! module tauint4_mod
!   implicit none
!   save
! contains

program tauint4!(j,xp,yp,zp,nxp,nyp,nzp,&
  !   xcell,ycell,zcell,tflag,u_s,u_a,b_wl,seg_flag)

!use grid_mod, only: xmax,ymax,zmax,nxg,nyg,nzg 
  implicit none


!INTENTS 
real*8 :: xp,yp,zp,nxp,nyp,nzp
integer:: xcell,ycell,zcell 
real*8:: u_s, u_a

real*8 :: ran 
real*8:: tau,taurun
real*8:: d, dcell

real*8 :: pos(3),dir(3)
integer :: cell(3),celli,cellj,cellk

logical :: face_flag(3)
!INPUTS FROM PROGRAM/ MODULES
xp=1.
yp=1.
zp=1.
nxp=1.
nyp=0.
nzp=0.
xcell=1.
ycell=1.
zcell=1. 

!Switch coordinate systems 
pos(1)=xp
pos(2)=yp
pos(3)=zp

cell(1)=xcell
cell(2)=ycell
cell(3)=zcell 

dir(1)=nxp
dir(2)=nyp
dir(3)=nzp


!Set cumumlative ditances to 0 
taurun=0.
d=0. 

!get random optical depth 
call random_number(ran)
tau=-log(ran)



do
	face_flag=(/.false.,.false.,.false./) !LEwis has this 
	print*, face_flag
	!find ditance to each wall : try it with LEwis's first then try with yours 


	!take shortest once as d_cell 


exit

enddo
!Switch coordinate systems 
xp=pos(1)
yp=pos(2)
zp=pos(3)

xcell=cell(1)
ycell=cell(2)
zcell=cell(3)


end program
!end module

!want to be able to pass ALL 3 into this as a nice vecotr 
subroutine dir_test(direction,current_pos,face,grid_max,cell,dist)
  !dx is the distace long the firection of movement to the next xface
  !UNLESS dx < delta. if it is we just move it TO the next face.
  !like it's already so close to the next face just move it to the next face?
  !THIS ALSO CHANGES ci though
  !can a FUNCTION change a modularised variable???? 
  !use grid_mod, ONLY: delta
 implicit none
 	real*8 delta
   real*8 direction, current_pos, face(3),grid_max,dist
   integer cell
   !real*8 di_test



   if(direction.gt.0.) then
    print*, 'gt 0'
      dist=(face(3)-current_pos)/direction
      if(dist.lt.delta) then
         current_pos=face(3)
         cell=cell+1
        dist=(face(3)-current_pos)/direction
      endif
   elseif(direction.lt.0.) then
    print*, 'lt 0'
      dist=(face(2)-current_pos)/direction
      if(dist.lt.delta) then
         current_pos=face(2)
         dist=(face(1)-current_pos)/direction
        cell=cell-1
      endif
  elseif(direction.eq.0.) then
    print*, 'eq 0'
      dist=1.e2*grid_max
   endif


 return
end subroutine dir_test

! subroutine dist_to_cell_wall(celli, cellj, cellk, xcur, ycur, zcur, next_wall,wall_dist)
!    !funtion that returns distant to nearest wall and which wall that is (x,y or z)
!    !
!    !
!       use iarray,      only : xface, yface, zface
!       use photon_vars, only : nxp, nyp, nzp

!       implicit none

!       real,    intent(INOUT) :: xcur, ycur, zcur
!       logical, intent(INOUT) :: dir(:)
!       integer, intent(INOUT) :: celli, cellj, cellk
!       real                   :: dx, dy, dz


!       if(nxp > 0.)then
!          dx = (xface(celli+1) - xcur)/nxp
!       elseif(nxp < 0.)then
!          dx = (xface(celli) - xcur)/nxp
!       elseif(nxp == 0.)then
!          dx = 100000.
!       end if

!       if(nyp > 0.)then
!          dy = (yface(cellj+1) - ycur)/nyp
!       elseif(nyp < 0.)then
!          dy = (yface(cellj) - ycur)/nyp
!       elseif(nyp == 0.)then
!          dy = 100000.
!       end if

!       if(nzp > 0.)then
!          dz = (zface(cellk+1) - zcur)/nzp
!       elseif(nzp < 0.)then
!          dz = (zface(cellk) - zcur)/nzp
!       elseif(nzp == 0.)then
!          dz = 100000.
!       end if

!       wall_dist = min(dx, dy, dz)
!       if(wall_dist < 0.)print'(A,7F9.5)','dcell < 0.0 warning! ',wall_dist,dx,dy,dz,nxp,nyp,nzp
!       if(wall_dist == dx)next_wall=(/.TRUE., .FALSE., .FALSE./)
!       if(wall_dist == dy)next_wall=(/.FALSE., .TRUE., .FALSE./)
!       if(wall_dist == dz)next_wall=(/.FALSE., .FALSE., .TRUE./)
!       if(.not.dir(1) .and. .not.dir(2) .and. .not.dir(3))print*,'Error in dir flag'
!       return

!    end subroutine dist_to_cell_wall 

