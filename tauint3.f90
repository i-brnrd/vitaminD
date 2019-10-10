 module tauint3_mod
   implicit none
   save
 contains
!PUT DELTAS INT
subroutine tauint3(j,xp,yp,zp,nxp,nyp,nzp,&
     xcell,ycell,zcell,tflag,u_s,u_a,b_wl,seg_flag)
use optical_properties_mod, only: nlayer 
use grid_mod, only: xmax,ymax,zmax,xface,yface,zface, MASK

implicit none

!INTENTS 
integer, intent(in) :: j, b_wl
real*8, intent(inout) :: xp,yp,zp,nxp,nyp,nzp,u_s(nlayer), u_a(nlayer)
integer, intent(inout):: xcell,ycell,zcell 

!mcrt, path lengths & optical 
real*8:: ran !random random_number
real*8:: tau,taurun, taucell !optical depths 
real*8:: d, d_cell  !actual distance travelled 
real*8 :: rk_tot,u_a_tot !total optical properties over all species 


real*8 :: pos(3),dir(3),grid_max(3),faces(3,3)
integer :: cell(3)


integer :: seg_flag, tflag
logical :: face_flag(3),edge_flag(3),complete
integer :: i 

!SWITCH COORDINATE SYSTEMS: Sandboxes tauint2 to an extent
!PUT DELTAS IN THERE
pos(1)=xp + xmax
pos(2)=yp + ymax
pos(3)=zp + zmax

cell(1)=xcell
cell(2)=ycell
cell(3)=zcell 

dir(1)=nxp
dir(2)=nyp
dir(3)=nzp

grid_max(1)=xmax
grid_max(2)=ymax
grid_max(3)=zmax

!Intitalise edge flags and cumulative distances 
edge_flag=(/.false.,.false.,.false./)
face_flag=(/.false.,.false.,.false./)
taurun=0.
d=0. 
!get random optical depth 
call random_number(ran)
tau=-log(ran)


do
	!find the next cell wall and distance to it 
	call find_next_cell(cell,pos,dir,grid_max,d_cell,face_flag,faces)

	!set optical properties for current cell 
  rk_tot=0. 
  u_a_tot=0.
  do i=1,nlayer
    rk_tot=rk_tot+ (MASK(cell(1),cell(2),cell(3),i)*(u_s(i) + u_a(i)))
    u_a_tot=u_a_tot+ (MASK(cell(1),cell(2),cell(3),i)*u_a(i))  
  enddo
  taucell=d_cell*rk_tot

  if (taurun + taucell .lt.tau) then !packet passes through current voxel 
    !Update path length counters, and position, and return to main program 
    taurun=taurun + taucell
    d=d+d_cell
    complete=.false.
    !PL_SUM(b_wl,CI,CJ,CK)=PL_SUM(b_wl,CI,CJ,CK)+d_cell 

    call update_position(complete,d_cell,pos,dir,cell,faces,face_flag,edge_flag)
    print*, 'EDGE FLAG', edge_flag,'FACE FLAG',face_flag,'CELL', cell 
    if (any(edge_flag)) then
      print*, 'MProg out of bounds WHOOOO', cell
      call boundary_behaviour(edge_flag,complete,cell,grid_max)
      stop
    endif

  else !within the voxel, we reach the interaction location
    !Update path length counters, and position, and return to main program 
    d_cell=((tau-taurun)/rk_tot)
    d=d+d_cell
    complete=.true. 
    print*, 'LEAVING',j
    !	PL_SUM(b_wl,CI,CJ,CK)=PL_SUM(b_wl,CI,CJ,CK)+d_cell 
    call update_position(complete,d_cell,pos,dir,cell,faces,face_flag) !with an exit flag 
    EXIT
  endif 
enddo

!complete  is not the SAME as t_flag
!complete is tau_reached 
!Switch coordinate systems back to those used in main modules 
xp=pos(1) -xmax
yp=pos(2) -ymax
zp=pos(3) -zmax

xcell=cell(1)
ycell=cell(2)
zcell=cell(3)

end subroutine
end module


subroutine boundary_behaviour(edge_flag,complete,cell,grid_max,pos)
use grid_mod, only: nxg,nyg,nzg,delta
implicit none 
real*8, intent(inout) :: pos(3)
real*8,intent(in) :: grid_max(3)
logical, intent(inout) :: edge_flag(3), complete
real*8 ::deltas(3)
integer :: cell(3),ng(3)
integer :: i, n_grid(3)
print*, 'IN BOUnDARY CONDItionS', edge_flag

deltas=delta
ng=(/nxg,nyg,nzg/)

do i=1,2
if (edge_flag(i)) then
  print*, 'x-y repeats nedeed',edge_flag
    if (cell(i).lt.1) then
      print*,'move to max', i, pos(i) 
      pos(i)=(2.d0*grid_max(i))-deltas(i)
      print*,'moved to max', i, pos(i) 
    elseif (cell(i).gt.ng(i)) then
     print*, 'move to min',i,pos(i)
     print*, deltas(i)
     pos(i)=deltas(i)
     print*,'moved to min', i, pos(i) 
    else 
     print*,'Error in boundary_behaviour',edge_flag,cell
    endif 
 endif   
enddo

if(edge_flag(3)) then
  print*, 'z behaviour needed!!'
  if (cell(3).lt.1) then
    print*,'bottom of grid...check we exit'
  elseif(cell(3).gt.ng(3))then
    print*,'need to do some fresnel!! ' 
  else 
    print*,'Error in boundary_behaviour',edge_flag,cell
  endif 
endif 

edge_flag=(/.false.,.false.,.false./)

end subroutine


subroutine update_position(complete,d_cell,current_position,dir,cell,faces,face_flag,edge_flag)!(cell,pos,dir,grid_max,dist_min,face_flag))
use grid_mod, only: delta,xface,yface,zface 
implicit none

logical :: face_flag(3),edge_flag(3)
real*8,intent(in) :: faces(3,3)
integer, intent(inout):: cell(3)
real*8, intent(inout) :: current_position(3)
real*8, intent(in) :: dir(3),d_cell
real*8 :: deltas(3)
logical, intent(in) :: complete
integer::i

integer:: low, up 

deltas=delta 

call get_cells(current_position,cell,edge_flag)

if(complete)then
          do i=1,3
          current_position(i)=current_position(i)+dir(i)*d_cell 
        enddo
else
  do i=1,3
   if(face_flag(i))then!use faces
            if(dir(i).gt.0.) then
              current_position(i)= faces(3,i) + deltas(i)
            elseif(dir(i).lt.0.) then
               current_position(i) = faces(2,i) - deltas(i)
            else
               print*,'Error in x dir in update_pos', dir,face_flag
            end if
     else 
     current_position(i) = current_position(i)+ dir(i)*d_cell   
   endif 
  enddo
  call get_cells(current_position,cell,edge_flag)
endif
end subroutine update_position

subroutine get_cells(current_position, cell,edge_flag)
!updates the current voxel based upon position
       use grid_mod, only : xface, yface, zface, nxg,nyg,nzg,xmax,ymax,zmax
       use search_bisec_mod
       implicit none

       real*8,    intent(IN):: current_position(3)
       integer, intent(OUT) :: cell(3)
       logical, intent(OUT) :: edge_flag(3)
       integer :: low,up
  
        call search_bisec(current_position(1),xface,low,up)
        cell(1)=low + lbound(xface,1) -1 


        call search_bisec(current_position(2),xface,low,up)
        cell(2)=low + lbound(yface,1) -1 

        call search_bisec(current_position(3),zface,low,up)
        cell(3)=low + lbound(zface,1) -1


        if ((cell(1).lt.1).or.(cell(1).gt.nxg))edge_flag(1)=.true.
        if ((cell(2).lt.1).or.(cell(2).gt.nyg))edge_flag(2)=.true.
        if ((cell(3).lt.1).or.(cell(3).gt.nzg))edge_flag(3)=.true.


        if ((cell(1).lt.1).or.(cell(1).gt.nxg)) print*, 'XXXXXXXXXXXXXX'
        if ((cell(2).lt.1).or.(cell(2).gt.nyg)) print*, 'YYYYYYYYYYYYY'
        if ((cell(3).lt.1).or.(cell(3).gt.nzg)) print*, 'ZZZZZZZZZZZZZZ'

   end subroutine get_cells


!!!!!!!--------------------------
subroutine find_next_cell(cell,pos,dir,grid_max,dist_min,face_flag,faces)
	!difference between lewis & Kenny: we DO NOT update ci,cj,ck in here. 
use grid_mod, only: xface,yface,zface
implicit none
real*8, intent(in) ::pos(3),dir(3),grid_max(3) !MODULE
integer, intent(inout) :: cell(3)
logical, intent(out) :: face_flag(3)
real*8 :: faces(3,3), distance(3), dist_min
integer :: i

logical:: edge_flag(3)

!function
real*8 :: dist_to_cell_wall
	!initita;ise face)flag and distance 
	face_flag=(/.false.,.false.,.false./) !LEwis has this 
	!if (.not.any(face_flag))  print*,'all are false good'
	distance=0.d0


	faces(:,1)=(/xface(cell(1)-1),xface(cell(1)),xface(cell(1)+1)/)
	faces(:,2)=(/yface(cell(2)-1),yface(cell(2)), yface(cell(2)+1)/)
	faces(:,3)=(/zface(cell(3)-1),zface(cell(3)),zface(cell(3)+1)/)

	do i=1,3
		distance(i)=dist_to_cell_wall(dir(i),pos(i),faces(:,i),grid_max(i))
     if (distance(i).lt.0.) then
      print*,'xory',i,'cell', cell(i),'direct',dir(i),'positio',pos(i),'3 faces',faces(:,i),grid_max(i)
     call get_cells(pos,cell,edge_flag)
      print*, cell
     endif 
	enddo


	dist_min=minval(distance(:))
	if(dist_min < 0.)print*,'dcell < 0.0 warning!! ',dist_min,distance,dir
  if(dist_min < 0.)stop
	face_flag(minloc(distance(:)))=.true.
	print*, dist_min,face_flag
	if (.not.any(face_flag))  print*,'Next Face Hit Flag ERROR'
end subroutine find_next_cell

function dist_to_cell_wall(direction,current_position,face_position,grid_max)
 implicit none
  real*8 :: dist_to_cell_wall
   real*8,intent(in) :: direction, current_position, face_position(3),grid_max

   if(direction.gt.0.) then
      dist_to_cell_wall=(face_position(3)-current_position)/direction
   elseif(direction.lt.0.) then
      dist_to_cell_wall=(face_position(2)-current_position)/direction
  elseif(direction.eq.0.) then
      dist_to_cell_wall=1.e2*grid_max
   endif
 return
end function







