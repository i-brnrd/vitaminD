
!************************************************************
program mcpolar
  include 'modules_full_list.txt'
  implicit none

  include 'photon.txt'

  !*****Parameter declarations ****************************************
  ! known and understood
  integer nphotons   !number of photon packets sent into simulation
  real*8 hgg,g2  !Henyey Greenstein phase function variables
  real*8 delta
  real*8 xmax,ymax,zmax
  real*8,parameter :: pi = 4.*atan(1.)
  real*8, parameter:: twopi=2.*pi
  real*8, parameter:: fourpi=4.*pi
  real*8 U_S(nlayer),U_A(nlayer),g
  real*8 ns ! refractive index of skin
  character*30 fname,fname2


 ! (loop) indices, counters and seeds & FLAGS: SET FLAGS TO 0???
  real*8 nscatt
  integer iseed
  integer j
  integer xcell,ycell,zcell
  integer tflag, seg_flag, r_flag
  INTEGER I,k,p


  ! ummmm what now
  real*8 kappa,albedo,pl,pc,sc
  real*8 e
  integer sz
  integer ph_count, d_flag
  real*8 diff !diffuse fractions of incident light
  integer b_wl
  real*8 WL
  real*8 dnaval
  real*8 lumin(nwl), lum
  real*8 :: l(nwl), c(nwl)
  real*8 ran
  integer  A,B
  integer CASE
  real*8 :: tot(2),uva(2)
  real*8, allocatable :: input_spec(:,:) !input spectrum
  real*8, dimension(:,:), allocatable :: epi_s,sc_s,stratc,epi,eumel,phmel,dna,ohb,dhb !Action Spectra
  real*8, allocatable :: dummy(:,:) !dummy allocatable spectrum

  real*8,allocatable:: ref_spec(:,:)
  real*8, allocatable :: trans_spec(:,:)

  real*8 :: sum

  real*8 :: spectrum(2,nwl)

  real*8 :: interpol

  ! Initialise Variables
  ph_count=0

  !include 'parameters.txt' !MAKE ANOTE OF WHAT THIS INCLUDES
    nphotons = 100000000
    xmax=0.025
    ymax=0.025
    zmax=0.02
     print*, xmax,ymax,zmax
  !Skin parameters/ Optical Properties parameters/ Tissue parameters
    ns=1.38
    !WHAT ARE THESE FFS
    kappa=1.
    pl=0.
    pc=0.
    sc=1.

   print*, nphotons


! have a params file for optical properties with the number of layers

  print*,'Simulation of MC-UVRT through upper layers of skin stratified into', nlayer, ' layer/s'

!!! so each layer has a filename and a set of properties
  !****Initialise action spectra for the optical properties


 if (nwl.gt.1) then
     !Stratum Corneum
     fname='./spectra/sc_abs_2.txt' !epi_as.txt'
     call load_spec(fname,sz,stratc)
     fname='./spectra/sc_scatt_2.txt' !epi_as.txt'fname
     call load_spec(fname,sz,sc_s)
     stratc(2,:)=stratc(2,:)*10.
     sc_s(2,:)=sc_s(2,:)*10.
     !EPIDERMIS
     fname='./spectra/epi_abs_2.txt'
     call load_spec(fname,sz,epi)
     fname='./spectra/epi_scatt_2.txt'
     call  load_spec(fname,sz,epi_s)
     epi(2,:)=epi(2,:)*10.  !datasets are per mm, need per cm
     epi_s(2,:)=epi_s(2,:)*10. !datasets are per mm, need per cm
     !Melanin
     fname='./spectra/eumel_as.txt'
     call load_spec(fname,sz,eumel)
     fname='./spectra/phmel_as.txt'
     call load_spec(fname,sz,phmel)
     !DNA
     !fname='./spectra/dna2.txt' !out by 1000 (stupid matlab)
     fname='./spectra/b924712b-f1.txt' !out by 1000 (stupid matlab)
      call load_spec(fname,sz,dna)
      !dna(2,:)=dna(2,:)
     !Heamoblobin
     fname='./spectra/ohb_as.txt'
     call load_spec(fname,sz,ohb)
     fname='./spectra/dhb_as.txt'
     call load_spec(fname,sz,dhb)

  !**Initialise Spectra**

     do i=1,nwl
       l(i)=280.+real(i) -1.
     enddo

        call random_seed()
        diff=1.d0
        !lum=50.1/(100.**2)
        !fname='./spectra/solar_spec.txt'
        fname='./spectra/bb_uva.csv'
        !fname='./spectra/mh_uva1.csv'
        !fname='./spectra/nb_uvb.csv'

        !fname='./spectra/f_uva1.csv'

     call load_spec(fname,sz,input_spec)


!interpolating input spectrum via what method?

     sum =0.0
      do i=1,nwl
       call interpolate(l(i),input_spec,interpol)
        if (interpol.lt.0.00) then
          interpol=0.00
        end if
        spectrum(1,i)=l(i)
        spectrum(2,i)=interpol
        sum = sum + interpol
        enddo
      call cdf(c,spectrum,nwl)

      lum=sum
      lumin=spectrum(2,:)


          print*,'spectra loaded'
          print*,'input spectrum:', fname
          print*,'Total Luminosity:',Lum
          print*,'Diffuse Fraction',diff

          !**VERIFICATION
          !Verify Optical Properties

        !  do i=1,nwl
          !  print*,wl
          !   wl=l(i)
          !   call op_props(wl,u_a,u_s,g,epi,eumel,phmel,dna,ohb,dhb,stratc,epi_s,sc_s)
          !!   write(12,*) wl,u_a(1),u_a(2),u_a(3),u_a(4),u_a(5)
          !   write(13,*) wl,u_s(1),u_s(2),u_s(5), g

        !  enddo
        !  close(12)
        !  close(13)


        !  enddo
        !  close(12)
        !  close(13)





   else

     print*, 'Jaques verification mode on nwls:', nwl
     Diff = 0.d0
     lum=1.d0
     l(1)= 630.d0
     lumin(1)=1.d0
     print *, lum, l
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



    call iarray
    call gridset(xmax,ymax,zmax,kappa)


    ! try to writeout gridmasks?
      ! just call density









    print*,'layer 2 cell count',lcount(2),'layer 3 cell count',lcount(3)
    !print*, lcount


                      ! wl= 315.d0

              !  call op_props(wl,u_a,u_s,g,epi,eumel,phmel,dna,ohb,dhb,stratc,epi_s,sc_s)
                !      print*, wl,u_a(1),u_a(2)
                !      print*, wl,u_s(1),u_s(2)
                !      wl=365.d0
              !call op_props(wl,u_a,u_s,g,epi,eumel,phmel,dna,ohb,dhb,stratc,epi_s,sc_s)
                !      print*, wl,u_a(1),u_a(2)
                !      print*, wl,u_s(1),u_s(2)


     print*, 'here'

     !*****Set small distance for use in optical depth integration routines
     !*****for roundoff effects when crossing cell walls
     delta=1.e-6*(2.*zmax/nzg)
!!$
     !**** Loop over nph photons from each source *************************
     nscatt=0

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

        ph_count=ph_count+1

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
              nscatt=nscatt+1

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

     print*, 'total number scatterings',NSCATT, 'per packet', real(NSCATT)/real(nphotons)

     !--------------------------------------------------------------------------------
     !     CONVERT PATH LENGTH COUNTER(S) INTO PATH LENGTH ESIMATORS


     print*,'photon count', ph_count, nphotons, ph_count/nphotons


     CALL PL_ESTIMATORS(ph_count,nphotons,XMAX,YMAX,ZMAX,l,lum,lumin,fname)


     print*,'Average number of scatterings = ',(nscatt/nphotons)





endprogram mcpolar
!-------------------------------------------------------------------------------
