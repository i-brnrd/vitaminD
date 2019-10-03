module tauint2_mod
  implicit none
  save
contains

subroutine tauint2(j,xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
     xcell,ycell,zcell,tflag,iseed,delta,u_s,u_a,b_wl,seg_flag)

  use grid_mod
  implicit none

  integer, intent(in)::j
  integer tflag,iseed,xcell,ycell,zcell,count,seg_flag
  real*8 xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax
  real*8 ran
  real*8 gsmax


  integer ci,cj,ck
  real*8 tau,taurun,taucell,d,d1,dcell,xcur,ycur,zcur,dsx,dsy,dsz
  real*8 dx,dy,dz,smax,delta

  real*8 rk, uat !rhokappa,energy_absorbed
  real*8 u_s(nlayer),u_a(nlayer)
  real*8 wl
  integer  i, b_wl

  real*8 ns,ca,ct1,ct2,st1,st2,rf !snell & fresnel components (ct1:costheta1 etc.)


  !*****tflag=0 means photon is in envelope
  tflag=0
  !**** generate random optical depth tau
  call random_number(ran)
  tau=-log(ran)

  !**********!SETUP REFRACTIVE INDEX PARAMETERS
  NS=1.38                   !REFRACTIVE INDEX OF SKIN
  CA=ASIN(1./NS)            !CRITICAL ANGLE

  !*****set the cumulative distance and optical depth (d and taurun)
  !*****along the photon path to zero.  set the current photon coordinates.
  !*****note that the origin of the (xcur,ycur,zcur) system is at the
  !*****bottom corner of the grid.
  taurun=0.
  d=0
  xcur=xp+xmax
  ycur=yp+ymax
  zcur=zp+zmax

  ci=xcell
  cj=ycell
  ck=zcell
  smax=gsmax(nxp,nyp,nzp,xmax,ymax,zmax,xcur,ycur,zcur)

  !*****integrate through grid

  if(smax.lt.delta) then
     tflag=1
     return
  endif


  do

     !*****find distance to next x, y, and z cell walls.
     !*****note that dx is not the x-distance, but the actual distance along
     !*****the direction of travel to the next x-face, and likewise for dy and dz.

     if(nxp.gt.0.) then
        dx=(xface(ci+1)-xcur)/nxp
        if(dx.lt.delta) then
           xcur=xface(ci+1)
           ci=ci+1
           dx=(xface(ci+1)-xcur)/nxp
        endif
     elseif(nxp.lt.0.) then
        dx=(xface(ci)-xcur)/nxp
        if(dx.lt.delta) then
           xcur=xface(ci)
           dx=(xface(ci-1)-xcur)/nxp
           ci=ci-1
        endif
     elseif(nxp.eq.0.) then
        dx=1.e2*xmax
     endif

     if(nyp.gt.0.) then
        dy=(yface(cj+1)-ycur)/nyp
        if(dy.lt.delta) then
           ycur=yface(cj+1)
           cj=cj+1
           dy=(yface(cj+1)-ycur)/nyp
        endif
     elseif(nyp.lt.0.) then
        dy=(yface(cj)-ycur)/nyp
        if(dy.lt.delta) then
           ycur=yface(cj)
           dy=(yface(cj-1)-ycur)/nyp
           cj=cj-1
        endif
     elseif(nyp.eq.0.) then
        dy=1.e2*ymax
     endif

     if(nzp.gt.0.) then
        dz=(zface(ck+1)-zcur)/nzp
        if(dz.lt.delta) then
           zcur=zface(ck+1)
           ck=ck+1
           dz=(zface(ck+1)-zcur)/nzp
        endif
     elseif(nzp.lt.0.) then
        dz=(zface(ck)-zcur)/nzp
        if(dz.lt.delta) then
           zcur=zface(ck)
           dz=(zface(ck-1)-zcur)/nzp
           ck=ck-1
        endif
     elseif(nzp.eq.0.) then
        dz=1.e2*zmax
     endif

     !this is a bit of a cheat: see screenshots in ./debug
     !*******THIS USED TO BE (dx.eq.0.)

     if( (dx.eq.0.) .or. ((abs(dx)).lt.(delta)) ) dx=1.e2*xmax
     if( (dy.eq.0.) .or. ((abs(dy)).lt.(delta)) ) dy=1.e2*ymax
     if( (dz.eq.0.) .or. ((abs(dz)).lt.(delta)) ) dz=1.e2*zmax

     !*****find distance to next cell wall -- minimum of dx, dy, and dz

     dcell=amin1(dx,dy,dz)

     if(dcell.le.0.) then
!!$        print *,'tauint2: dcell < 0',dx,dy,dz,nxp,nyp,nzp
!!$        print *,xcur,ycur,zcur,ci,cj,ck
        print*,'dcell < 0 fault @ at pkt',j,dcell,dx,dy,dz
        seg_flag=1
        return
     endif
     if(dx.lt.0.) dcell=amin1(dy,dz)
     if(dy.lt.0.) dcell=amin1(dx,dz)
     if(dz.lt.0.) dcell=amin1(dx,dy)

!!$
!!$     if((ci.gt.nxg).or.(ci.lt.0))  print*,j,'celli',ci,xcur
!!$     if((cj.gt.nyg).or.(cj.lt.0)) print*,j,'cellj',cj,ycur
!!$     if((ck.gt.nyg).or.(ck.lt.0)) print*,j,'cellk',ck,xcur
      !print*, ci,cj,ck
     if((ci.gt.nxg).or.(ci.lt.1).or.(cj.gt.nyg).or.(cj.lt.1).or.(ck.gt.nzg).or.(ck.lt.1)) then
        print*,'cell out of bounds @ pkt',j, ci,cj,ck
        seg_flag=1
        return
     endif




     !*****distances are only zero if photon is on cell wall.  if it is
     !*****on cell wall then set to arbitrary large distance, since we will
     !*****in fact hit another wall

     !*****optical depth to next cell wall is
     !*****taucell= (distance to cell)*(opacity of current cell)

     rk=0.
     uat=0.
     do i=1,nlayer
        rk=(MASK(ci,cj,ck,i)*(u_s(i) + u_a(i)))+rk
        uat= (MASK(ci,cj,ck,i)*u_a(i)) + uat
     enddo
     !rhokap(ci,cj,ck)=rk
     !print*,'RK',rhokap(ci,cj,ck)

     !taucell=dcell*rhokap(ci,cj,ck)
     taucell=dcell*rk


     !*****if taurun+taucell>tau then scatter at distance d+d1.
     !*****update photon position and cell.
     !*****if taurun+taucell<tau then photon moves distance dcell
     !*****(i.e. ends up on next cell wall) and update photon position
     !*****and cell.

     !*****this is PL Counterbgf

     if((taurun+taucell).ge.tau) then

        d1=(tau-taurun)/rk
        d=d+d1
        taurun=taurun+taucell
        xcur=xcur+d1*nxp
        ycur=ycur+d1*nyp
        zcur=zcur+d1*nzp

        PL_SUM(b_wl,CI,CJ,CK)=PL_SUM(b_wl,CI,CJ,CK)+D1



        !***************Linear Grid ************************
        ci=int(nxg*xcur/(2.*xmax))+1
        cj=int(nyg*ycur/(2.*ymax))+1
        ck=int(nzg*zcur/(2.*zmax))+1
        !****************************************************
        EXIT                !Exit Do loop

     ELSE

        d=d+dcell
        taurun=taurun+taucell

        if (d.gt.(smax*0.999)) then
           !     WHICH EDGE WILL IT HIT?
           if ((dx.lt.dy).AND.(dx.lt.dz)) then

              PL_SUM(b_wl,CI,CJ,CK)=PL_SUM(b_wl,CI,CJ,CK)+DCELL
            !E_SUM(b_wl,ci,cj,ck)=E_SUM(b_wl,ci,cj,ck)+(dcell*ea)

              xcur=xcur+dcell*nxp
              ycur=ycur+dcell*nyp
              zcur=zcur+dcell*nzp

              xcur=(2.*xmax)-xcur !flip

              !***************Linear Grid ************************
              ci=int(nxg*xcur/(2.*xmax))+1
              cj=int(nyg*ycur/(2.*ymax))+1
              ck=int(nzg*zcur/(2.*zmax))+1
              !****************************************************

              smax=gsmax(nxp,nyp,nzp,xmax,ymax,zmax,xcur,ycur,zcur)
              d=0.
              CYCLE

           elseif ((dy.lt.dx).and.(dy.lt.dz)) then

              !     IF GOING TO HIT YEDGE, LOOP

              PL_SUM(b_wl,CI,CJ,CK)=PL_SUM(b_wl,CI,CJ,CK)+DCELL
            !  E_SUM(b_wl,ci,cj,ck)=E_SUM(b_wl,ci,cj,ck)+(dcell*ea)

              xcur=xcur+dcell*nxp
              ycur=ycur+dcell*nyp
              zcur=zcur+dcell*nzp

              ycur=(2.*ymax)-ycur
              !***************Linear Grid ************************
              ci=int(nxg*xcur/(2.*xmax))+1
              cj=int(nyg*ycur/(2.*ymax))+1
              ck=int(nzg*zcur/(2.*zmax))+1
              !****************************************************

              smax=gsmax(nxp,nyp,nzp,xmax,ymax,zmax,xcur,ycur,zcur)
              d=0.
              CYCLE

           ELSE          ! HIT Z
              IF(nzp.le.0) then ! hitting bottom edge; let it go
                 tflag =1
                 EXIT
              ELSE       !hitting skin surface... reflect or transmit...

                 CT1=NZP

                 IF (ACOS(CT1).GT.CA) THEN !TIR

                    PL_SUM(b_wl,CI,CJ,CK)=PL_SUM(b_wl,CI,CJ,CK)+DCELL
                    !E_SUM(b_wl,ci,cj,ck)=E_SUM(b_wl,ci,cj,ck)+(dcell*ea)

                    xcur=xcur+dcell*nxp
                    ycur=ycur+dcell*nyp
                    zcur=zcur+dcell*nzp
                    !***************Linear Grid ************************
                    ci=int(nxg*xcur/(2.*xmax))+1
                    cj=int(nyg*ycur/(2.*ymax))+1
                    ck=int(nzg*zcur/(2.*zmax))+1
                    !****************************************************

                    NZP = -NZP !CONTINUE WITH INTEGRATION
                    CT1=NZP

                    smax=gsmax(nxp,nyp,nzp,xmax,ymax,zmax,xcur,ycur,zcur)
                    d=0.
                    CYCLE

                 ELSE    !CALCULATE SNELL& FRESNEL EQU. COMPONENTS

                    ST1=SIN(ACOS(CT1))
                    ST2=ST1*NS
                    CT2=COS(ASIN(ST2))

                    RF = ((NS*CT1-CT2)/(NS*CT1 + CT2))**2
                    RF = RF + ((NS*CT2-CT1)/(NS*CT2 + CT1))**2
                    RF = RF/2.

                    call random_number(ran)
                    IF (RAN.LT.RF) THEN !REFLECT PHOTON

                       PL_SUM(b_wl,CI,CJ,CK)=PL_SUM(b_wl,CI,CJ,CK)+DCELL
                       !E_SUM(b_wl,ci,cj,ck)=E_SUM(b_wl,ci,cj,ck)+dcell*ea
                       xcur=xcur+dcell*nxp
                       ycur=ycur+dcell*nyp
                       zcur=zcur+dcell*nzp
                       !***************Linear Grid ************************
                       ci=int(nxg*xcur/(2.*xmax))+1
                       cj=int(nyg*ycur/(2.*ymax))+1
                       ck=int(nzg*zcur/(2.*zmax))+1
                       !***************************************************

                       nzp = -nzp
                       ct1=nzp

                       smax=gsmax(nxp,nyp,nzp,xmax,ymax,zmax,xcur,ycur,zcur)
                       d=0.
                       CYCLE

                    else ! I don't know
                       tflag=1
                       EXIT
                    ENDIF
                 endif
              endif
           endif

        else

           PL_SUM(b_wl,CI,CJ,CK)=PL_SUM(b_wl,CI,CJ,CK)+DCELL
           !E_SUM(b_wl,ci,cj,ck)=E_SUM(b_wl,ci,cj,ck)+dcell*ea

           xcur=xcur+dcell*nxp
           ycur=ycur+dcell*nyp
           zcur=zcur+dcell*nzp
           !***************Linear Grid ************************
           ci=int(nxg*xcur/(2.*xmax))+1
           cj=int(nyg*ycur/(2.*ymax))+1
           ck=int(nzg*zcur/(2.*zmax))+1
           !****************************************************
        endif


     endif
  end do

  !*****calculate photon final position.  if it escapes envelope then
  !*****set tflag=1.  if photon doesn't escape leave tflag=0 and update
  !*****photon position.
  if((d.ge.(smax*0.999))) then
     tflag=1
  else

     xp=xcur-xmax
     yp=ycur-ymax
     zp=zcur-zmax
     xcell=ci
     ycell=cj
     zcell=ck

  endif

  return
end subroutine

end module tauint2_mod

FUNCTION  gsmax(nxp,nyp,nzp,xmax,ymax,zmax,xcur,ycur,zcur)
  implicit none
  real*8 gsmax
  real*8 nxp,nyp,nzp,xmax,ymax,zmax
  real*8 xcur,ycur,zcur,dsx,dsy,dsz

  if(nxp.gt.0.) then
     dsx=(2.*xmax-xcur)/nxp
  elseif(nxp.lt.0.) then
     dsx=-xcur/nxp
  elseif(nxp.eq.0.) then
     dsx=1.e2*xmax
  endif

  if(nyp.gt.0.) then
     dsy=(2.*ymax-ycur)/nyp
  elseif(nyp.lt.0.) then
     dsy=-ycur/nyp
  elseif(nyp.eq.0.) then
     dsy=1.e2*ymax
  endif

  if(nzp.gt.0.) then
     dsz=(2.*zmax-zcur)/nzp
  elseif(nzp.lt.0.) then
     dsz=-zcur/nzp
  elseif(nzp.eq.0.) then
     dsz=1.e2*zmax
  endif

  gsmax=amin1(dsx,dsy,dsz)
  return

end function gsmax
