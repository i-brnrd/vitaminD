!************************************************************
program mcpolar
  include 'modules_full_list.txt'
  implicit none


  !*****Parameter declarations ****************************************
  !input params file will include IN SOME WAY
  integer nphotons   !number of photon packets sent into simulation
  !real*8 xmax,ymax,zmax !this is a grid parameter

  real*8 diff !diffuse fractions of incident light
  character*30 :: fname_incident_irradiation
  real*8, parameter :: pi = 4.*atan(1.)
  real*8, parameter:: twopi=2.*pi
  real*8, parameter:: fourpi=4.*pi

  !input params but explicitley filled up by spectra (so like filled up variables)
  !isla that is a variable you plum
  real*8 incident_spec_irr(nwl)
  real*8 tot_irr
  real*8 :: cdf(nwl)

  integer:: pkt_count,scatter_count   !Counts packets that actually travel through medium (not reflected)
  real*8 ran
  real*8 albedo
!properties of a PACKET/ al the packets....
  !include 'photon.txt'  !why !why !why
  real*8 hgg,g2  !Henyey Greenstein phase function variables
  real*8 u_s_old(nlayer),u_a_old(nlayer),g
  real*8 ns ! refractive index of skin
  integer b_wl
  integer xcell,ycell,zcell

 ! (loop) indices, counters and seeds &
 integer j  !
 INTEGER I !GET RID OF THIS it needs to be placed into a sub
 !*****! ummmm what now
 real*8 kappa,pl,pc,sc
 integer tflag, seg_flag, r_flag

 logical:: diffuse_flag


 print*,'**********************'
 print*,'Simulation of MC-UVRT through upper layers of skin'

 include 'parameters.txt' !MAKE ANOTE OF WHAT THIS INCLUDES !JAQUES VERIFICATION ON OFF TRUE FALSE MAYBE
 diff=1.d0
 fname_incident_irradiation='./spectra/bb_uva.csv'

 print*, 'model size(cm)', xmax,ymax,zmax
 print*, 'sim in',nlayer, ' layers'
 print*, 'using photons', nphotons

call iarray !Initialise the arrays shared by the modules to 0

! Initialise counter :
pkt_count=0
!Skin parameters/ Optical Properties parameters/ Tissue parameters
ns=1.38   !CHECK THIS i think its in tauint2: it is, what a mess. Thing to do: hunt it down
!WHAT ARE THESE FFS
kappa=1.  !GET RID OF AS I DON@T USE IT
pl=0.
pc=0.
sc=1.

!INITIALISE GRID
call gridset(kappa)

!INITIALISE PACKETS OPTICAL PROPERTIES
if (nwl.gt.1) then
  call random_seed()
  !**Initialise Spectra**
  print*, wl_start
  do i=1,nwl
    l(i)= wl_start + real(i) -1.
  enddo
  call optical_properties_set(u_a,u_s,g_skin)
  call load_spec2(fname_incident_irradiation,incident_spec_irr,1.d0)
  call get_cdf(cdf,l,incident_spec_irr,nwl)
  tot_irr=sum(incident_spec_irr)
  !___________________
  print*,'Incident spectral irradiance loaded from ', fname_incident_irradiation
  print*,'Total irradiance:',tot_irr
  print*,'Diffuse Fraction',diff

  !JAQUES VERIFICATION: DEFO into SUBROUTINE BUT WHEN modules are wrangled.
else

  print*, 'Jaques verification mode on nwls:', nwl
  Diff = 0.d0
  tot_irr=1.d0
  l(1)= 630.d0
  incident_spec_irr(1)=1.d0
  print *, tot_irr, l
  print*, 'diffuse fraction', diff
  !print*,'Total Luminolssity:',Lum
  call random_seed()
  g=0.9
  do i=1,nlayer
    u_a_old(i)=0.23
    u_s_old(i)=(21.)/(1.-g)
  enddo

endif

!---------------------------------------------------------------------------------------------------------------
!MCRT
!-------------------------------------------------------------------------------------------------------------------------

     !**** Loop over nph photons from each source *************************
scatter_count=0

do j=1,20
    if(mod(j,100000).eq.0)then
      print *, j,' scattered photons completed'
    end if

        !*****Release photon from point source *******************************
    call sourceph(twopi,xcell,ycell,zcell,cdf,diff,b_wl,diffuse_flag)
        if (nwl.eq.1) b_wl=1

        if (diffuse_flag) then
           call reflect(1.d0,ns,cost,pi,r_flag) !fresnel reflection
           if (r_flag.eq.1) then
              CYCLE
           endif
        endif

        pkt_count=pkt_count+1 !counts actual packets reaching the medium
        seg_flag=0
        !add one to wavelength bin
        n_phot_wl(b_wl)= n_phot_wl(b_wl) + 1
! obtain all the optical properties for the photon
        u_a_old=u_a(:,b_wl)
        u_s_old=u_s(:,b_wl)
        g=g_skin(b_wl)
        hgg=g
        g2=hgg**2      ! Henyey-Greenstein parameter, hgg^2

        !******Find FIRST scattering location
       ! if((xcell.lt.1).or.(ycell.lt.1).or.(zcell.lt.1)) EXIT
       !IMPORTANT NOTEy its arguably faster for a single opprop set to be given tp tauint2. i don't know.
        call tauint3(j,xp,yp,zp,nxp,nyp,nzp,&
             xcell,ycell,zcell,tflag,u_s_old,u_a_old,b_wl,seg_flag)
       
        if (seg_flag.eq.1) CYCLE  !why
        !********Photon scatters in grid until it exits (tflag=1)
        tflag=0
        seg_flag=0
        do while(tflag.eq.0)
          albedo=0.d0
           do i =1,nlayer
              albedo=albedo+u_s_old(i)/(u_s_old(i)+u_a_old(i))*MASK(xcell,ycell,zcell,i)
           enddo
           call random_number(ran)
           if(ran.lt.albedo) then
              cost=nzp
        
              !************Scatter photon into new direction and update Stokes parameters
              call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,&
                   fi,fq,fu,fv,pl,pc,sc,hgg,g2,pi,twopi)
              scatter_count=scatter_count+1
           else
  
              EXIT
           endif

         
           !************Find next scattering location
            !if((xcell.lt.1).or.(ycell.lt.1).or.(zcell.lt.1)) EXIT
           call tauint3(j,xp,yp,zp,nxp,nyp,nzp,&
                xcell,ycell,zcell,tflag,u_s_old,u_a_old,b_wl,seg_flag)
           !print*,xcell,ycell,zcell
           if (seg_flag.eq.1) EXIT
        end do
       
     end do ! end loop over nph photons
     stop
     
     print*, 'total number scatterings',scatter_count, 'per packet', real(scatter_count)/real(nphotons)
     !--------------------------------------------------------------------------------
     !     CONVERT PATH LENGTH COUNTER(S) INTO PATH LENGTH ESIMATORS
     print*,'photon count', pkt_count, nphotons, pkt_count/nphotons
     CALL PL_ESTIMATORS(pkt_count,nphotons,l,tot_irr,incident_spec_irr)
     print*,'Average number of scatterings = ',(scatter_count/nphotons)

endprogram mcpolar
