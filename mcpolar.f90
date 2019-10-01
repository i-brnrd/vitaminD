!************************************************************
program mcpolar
  include 'modules_full_list.txt'
  implicit none
  include 'photon.txt'  !why !why !why

  !*****Parameter declarations ****************************************
  !****! known and understood
  !input params file will include IN SOME WAY
  integer nphotons   !number of photon packets sent into simulation
  real*8 xmax,ymax,zmax
  real*8 diff !diffuse fractions of incident light
  character*30 :: fname_incident_irradiation
  real*8, parameter :: pi = 4.*atan(1.)
  real*8, parameter:: twopi=2.*pi
  real*8, parameter:: fourpi=4.*pi

  !input params but explicitley filled up by spectra (so like filled up variables)
  !isla that is a variable you plum
  !change fname, fname2....argh....
  real*8 incident_spec_irr(nwl)
  real*8 tot_irr
  real*8 :: c(nwl)
  real*8, dimension(:,:), allocatable :: epi_s,sc_s,stratc,epi,eumel,phmel,dna,ohb,dhb !Action Spectra

  character*30 fname,fname2, filename

  integer pkt_count,scatter_count   !Counts packets that actually travel through medium (not reflected)

  !verbose countin integers
  !prefer i_wl for the waveengths
  integer i_wl

  real*8,dimension(nlayer,nwl):: u_s_test,u_a_test
  real*8:: g_skin_test(nwl)
  real*8 :: spec(nwl)
! not cecked not sorted.
  real*8 hgg,g2  !Henyey Greenstein phase function variables
  real*8 delta
  real*8 U_S(nlayer),U_A(nlayer),g
  real*8 ns ! refractive index of skin

 ! (loop) indices, counters and seeds & FLAGS: SET FLAGS TO 0???
  integer iseed
  integer j
  integer xcell,ycell,zcell
  integer tflag, seg_flag, r_flag
  INTEGER I,k,p
  integer d_flag
  integer sz
  integer b_wl
  !*****! ummmm what now
  real*8 kappa,albedo,pl,pc,sc
  real*8 e

  real*8 WL
  real*8 dnaval
  real*8 ran
  integer  A,B
  integer CASE
  real*8 :: tot(2),uva(2)

  print*,'**********************'
  print*,'Simulation of MC-UVRT through upper layers of skin'

  include 'parameters.txt' !MAKE ANOTE OF WHAT THIS INCLUDES !JAQUES VERIFICATION ON OFF TRUE FALSE MAYBE
  diff=1.d0
  fname_incident_irradiation='./spectra/bb_uva.csv'

  print*, 'model size(cm)', xmax,ymax,zmax
  print*, 'sim in',nlayer, ' layers'
  print*, 'using photons', nphotons

    call iarray !Initialise the arrays shared by the modules to 0
  ! Initialise counter :why? what is this and where is it used? withn MCRT, and within PL counters
    pkt_count=0
  !Skin parameters/ Optical Properties parameters/ Tissue parameters
    ns=1.38   !CHECK THIS i think its in tauint2: it is, what a mess. Thing to do: hunt it down
    !WHAT ARE THESE FFS
    kappa=1.  !GET RID OF AS I DON@T USE IT
    pl=0.
    pc=0.
    sc=1.

!INITIALISE GRID
    call gridset(xmax,ymax,zmax,kappa)
!INITIALISE PACKETS OPTICAL PROPERTIES :
!everything relies on the wavelengths now being these, eeek.
print*, wl_start
do i=1,nwl
  l(i)= wl_start + real(i) -1.
  print*, l(i)
enddo

 if (nwl.gt.1) then
   include 'include_old_op_props.txt'
  !**Initialise Spectra**
    call random_seed()

!Check sum, luminosity, and CDF all work ok. Lumin, lum, stupid names. integrated_lum useful? Idk
!only actually USE 'lum' as in total luminosity much much later (pl estimators) AND i include
!that once w ehave got the CDF I o't need lumin
! certainly don't need so many repeats fo the same dataset.
!what wsa it- incident_irradiance? lamp? source_irradiance? source_spec_irr?
    call load_spec2(fname_incident_irradiation,incident_spec_irr,1.d0)
    call get_cdf(c,l,incident_spec_irr,nwl)
    tot_irr=sum(incident_spec_irr)
    print*,'Incident spectral irradiance loaded from ', fname_incident_irradiation
    print*,'Total irradiance:',tot_irr
    print*,'Diffuse Fraction',diff

    !    now can do



!leave this point when sum, lum, diff are loaded in ok and results the SAME USING THE NEW SPECTRUM

  !JAQUES VERIFICATION: DEFO into SUBROUTINE
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
     wl=630.d0
   do i=1,nlayer
     u_a(i)=0.23
     u_s(i)=(21.)/(1.-g)
  enddo

 endif


! I ARRAY Should set all arrays in all modules to zero maybe; anddo it a bit more elegantly??
!

    do i=1,nwl
         l(i)=280.+real(i) -1.
    enddo

        !ONCE THIS IS FIXED: shove the other stuff into a text file and leave as 'include . txt'; the other optical properties loading up; so that the optical prperties lookup still works ok
        !then we need to make the BIG change to the code regarding intro: nex committ will be pre adoption of hte packet mod/ object/
    call op_prop_set(u_a_test,u_s_test,g_skin_test)

    print*,'layer 2 cell count',lcount(2),'layer 3 cell count',lcount(3)
!---------------------------------------------------------------------------------------------------------------


!-------------------------------------------------------------------------------------------------------------------------

     !*****Set small distance for use in optical depth integration routines
     !*****for roundoff effects when crossing cell walls
     delta=1.e-6*(2.*zmax/nzg)
!!$
     !**** Loop over nph photons from each source *************************
     scatter_count=0

     do j=1,nphotons
        if(mod(j,100000).eq.0)then
           print *, j,' scattered photons completed'
        end if

        !*****Release photon from point source *******************************
        call sourceph(xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi,&
             fi,fq,fu,fv,xmax,ymax,zmax,twopi,xcell,ycell,zcell,iseed,L,&
             C,DELTA,WL,diff,d_flag)

             !ns=1.39-((wl-279)*0.0001)

        if (d_flag.eq.1) then  !if diffuse photon, check fresnel
           !direct photons ~2% misdiagnosed as diffuse ...
           call reflect(1.d0,ns,cost,pi,r_flag) !fresnel reflection
           if (r_flag.eq.1) then
              CYCLE
           endif
        endif

        pkt_count=pkt_count+1

        seg_flag=0
        e=1./wl

        !bin wavelength

        if (nwl.gt.1) then
        b_wl = ceiling( ( wl - l(1) ) * nwl / ( l(nwl) - l(1) ) )
      else
          b_wl = 1
        endif
        !add one to wavelength bin
        n_phot_wl(b_wl)= n_phot_wl(b_wl) + 1

! obtain all the optical properties for the photon
        call op_props(wl,u_a,u_s,g,epi,eumel,phmel,dna,ohb,dhb,stratc,epi_s,sc_s)
        !do i=1,nlayer
        !  u_a(i)=0.23
        !  u_s(i)=(21.)/(1.-g)
      ! enddo

        dnaval=u_a(4)/0.0185       ! dnaval : what is this?
!!$     !Jaques Verification :PUT THIS IN A SUBROUTINE

        hgg=g
        g2=hgg**2      ! Henyey-Greenstein parameter, hgg^2

        !******Find FIRST scattering location
       ! if((xcell.lt.1).or.(ycell.lt.1).or.(zcell.lt.1)) EXIT
        call tauint2(j,xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
             xcell,ycell,zcell,tflag,iseed,delta,u_s,u_a,b_wl,e,seg_flag,dnaval)



          ! print*,'after tauint2',xp,yp,zp, xcell,ycell,zcell,seg_flag

        if (seg_flag.eq.1) EXIT
        !********Photon scatters in grid until it exits (tflag=1)
        tflag=0
        seg_flag=0
        do while(tflag.eq.0)
          albedo=0.d0
           do i =1,nlayer
              albedo=albedo+u_s(i)/(u_s(i)+u_a(i))*MASK(xcell,ycell,zcell,i)
           enddo
           call random_number(ran)
           if(ran.lt.albedo) then
              cost=nzp
              !************Scatter photon into new direction and update Stokes parameters

              call stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,&
                   fi,fq,fu,fv,pl,pc,sc,hgg,g2,pi,twopi,iseed)
              scatter_count=scatter_count+1

           else

              EXIT
           endif

           !************Find next scattering location
            !if((xcell.lt.1).or.(ycell.lt.1).or.(zcell.lt.1)) EXIT
           call tauint2(j,xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
                xcell,ycell,zcell,tflag,iseed,delta,u_s,u_a,b_wl,e,seg_flag,dnaval)
           !print*,xcell,ycell,zcell
           if (seg_flag.eq.1) EXIT
        end do



     end do ! end loop over nph photons

     print*, 'total number scatterings',scatter_count, 'per packet', real(scatter_count)/real(nphotons)

     !--------------------------------------------------------------------------------
     !     CONVERT PATH LENGTH COUNTER(S) INTO PATH LENGTH ESIMATORS


     print*,'photon count', pkt_count, nphotons, pkt_count/nphotons


     CALL PL_ESTIMATORS(pkt_count,nphotons,XMAX,YMAX,ZMAX,l,tot_irr,incident_spec_irr,fname)


     print*,'Average number of scatterings = ',(scatter_count/nphotons)





endprogram mcpolar
!-------------------------------------------------------------------------------
