
!************************************************************
program mcpolar

  use grid_mod
  use iarray_mod
  use gridset_mod
  use get_dim_mod
  use load_mod
  use load_spec_mod
  use search_bisec_mod
  use interpolate_mod
  use op_props_mod
  use source_ph_mod
  use tauint2_mod
  use reflect_mod
  use pl_estimators_mod

  implicit none 

  include 'photon.txt'

  !*****Parameter declarations ****************************************
  integer nphotons,iseed,j,xcell,ycell,zcell,tflag, seg_flag, r_flag
  real*8 nscatt
  real*8 kappa,albedo,hgg,pl,pc,sc,xmax,ymax,zmax

  !real*8 sz(3) !(x,y,z)=(1,2,3)

  real*8 e
  real*8 pi,twopi,fourpi,g2,delta

  integer n_spec, sz

  integer ph_count, d_flag

  real*8 ns ! refractive index of skin

  character*30 fname
  real*8 WL
  real*8 U_S(4),U_A(4),g

  real*8 diff,lum !diffuse fractions of incident light

  integer b_wl
  INTEGER I,k
  integer ecount,mcount

  real*8 ran
  integer  A,B

  !comparison help
  integer CASE
  real*8 :: e_uva,e_tot
  real*8 :: tot(2),uva(2)


  real*8, allocatable :: input_spec(:,:) !input spectrum
  real*8, dimension(:,:), allocatable :: epi,eumel,phmel,dna,ohb,dhb !Action Spectra
  real*8, dimension(:),allocatable :: L, C !wavelengths & CDF

  real*8, allocatable :: dummy(:,:) !dummy allocatable spectrum 

  !*****Set up constants, pi and 2*pi  ********************************
  pi=4.*atan(1.)
  twopi=2.*pi
  fourpi=4.*pi
  !Input Parameters
  kappa=1.
  pl=0.
  pc=0.
  sc=1.
  xmax=0.05
  ymax=0.05
  zmax=0.05

  ph_count=0

  ns=1.38
  !**Initialise Spectra**
  !INITIALISE THE SOLAR/SUNBED SPECTRUM CDF

 
        call random_seed()
        Diff=0.13
        lum=50.1/(100.**2)
        fname='./spectra/solar_spec.txt'
!!$
!!$     diff=0.13
!!$     Lum=49.11/(100.**2)
!!$     fname='./spectra/ss_uva.txt'
   
!!$
        
        diff=1. !0.!1.
        lum=283.65/(100.**2)
        fname='./spectra/sunbed_spec.txt'

!!$        Lum=277.01/(100.**2)! w/m^2->w/cm^2
!!$        fname='./spectra/sb_uva.txt'

     print*,'input spectrum:', fname
     print*,'Luminosity:',Lum
     print*,'Diffuse Fraction',diff

     call load_spec(fname,n_spec,input_spec)
     allocate(l(n_spec),c(n_spec))
     l(:)=input_spec(1,:)
     call cdf(c,input_spec,n_spec)

     !****Initialise action spectra for the optical properties
     !Epidermis
     fname='./spectra/epi_as.txt' !_as_the.txt'
     call load_spec(fname,sz,epi)

     !Melanin 
     fname='./spectra/eumel_as.txt'
     call load_spec(fname,sz,eumel)
     fname='./spectra/phmel_as.txt'
     call load_spec(fname,sz,phmel)

     !DNA
     fname='./spectra/dna_as.txt' !out by 1000 (stupid matlab)
     call load_spec(fname,sz,dna)
     dna(1,:)=dna(1,:)*1000.
     dna(2,:)=dna(2,:)*1000.

     !Heamoblobin 
     fname='./spectra/ohb_as.txt'
     call load_spec(fname,sz,ohb)
     fname='./spectra/dhb_as.txt'
     call load_spec(fname,sz,dhb)

     print*,'spectra loaded'


     call iarray
     call gridset(xmax,ymax,zmax,kappa)

     !**finding out total number of cells in EPIDERMIS
     ecount=0
     mcount=0
     do i=1,nxg

        do  j=1,nyg
           do k=1,nzg
              if (MASK(i,j,k,1)==1) then
                 ecount=ecount+1
              else if (Mask(i,j,k,4)==1) then
                 mcount=mcount+1
              endif
           enddo
        enddo
     enddo
     print*,ecount,mcount



     !**VERIFICATION
     !Verify Optical Properties
     open(12,file='./plots/abs.dat',status='replace')
     open(13,file='./plots/ops.dat',status='replace')

     do i=1,n_spec

        wl=l(i)
        call op_props(wl,u_a,u_s,g,epi,eumel,phmel,dna,ohb,dhb,ecount,mcount)
        write(12,*) wl,u_a(1),u_a(2),u_a(3),u_a(4)
        write(13,*) wl, u_s(1), g

     enddo
     close(12)
     close(13)

     !*****Set small distance for use in optical depth integration routines 
     !*****for roundoff effects when crossing cell walls
     delta=1.e-6*(2.*xmax/nxg)
!!$
     !**** Loop over nph photons from each source *************************
!!$  ! nscatt=0
     nphotons=1000000
     do j=1,nphotons
       
        if(mod(j,1000000).eq.0)then
           print *, j,' scattered photons completed'
        end if

        !*****Release photon from point source *******************************
        call sourceph(xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi,&
             fi,fq,fu,fv,xmax,ymax,zmax,twopi,xcell,ycell,zcell,nxg,nyg,nzg,iseed,L,&
             C,DELTA,WL,n_spec,diff,d_flag)

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
        b_wl = ceiling( ( wl - l(1) ) * nwl / ( l(n_spec) - l(1) ) )
        n_phot_wl(b_wl)=n_phot_wl(b_wl) +1

        !* *****COMPUTE HG PARAMETER FROM WL:WHY NOT LAYER DEP

        call op_props(wl,u_a,u_s,g,epi,eumel,phmel,dna,ohb,dhb,ecount,mcount)
!!$     !Jaques Verification 
!!$     g=0.9
!!$     do i=1,4        
!!$        u_a(i)=0.23
!!$        u_s(i)=(21.)/(1.-g)
!!$     enddo

        hgg=g       
        g2=hgg**2      ! Henyey-Greenstein parameter, hgg^2

        !******Find FIRST scattering location
       ! if((xcell.lt.1).or.(ycell.lt.1).or.(zcell.lt.1)) EXIT
        call tauint2(j,xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,&
             xcell,ycell,zcell,tflag,iseed,delta,u_s,u_a,b_wl,e,seg_flag)
        if (seg_flag.eq.1) EXIT
        !********Photon scatters in grid until it exits (tflag=1) 
        tflag=0
        seg_flag=0
        do while(tflag.eq.0)

           albedo=0.0
           do i =1,4 
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
                xcell,ycell,zcell,tflag,iseed,delta,u_s,u_a,b_wl,e,seg_flag)
           !print*,xcell,ycell,zcell
           if (seg_flag.eq.1) EXIT
        end do
     end do ! end loop over nph photons

     !--------------------------------------------------------------------------------
     !     CONVERT PATH LENGTH COUNTER(S) INTO PATH LENGTH ESIMATORS

     print*,'photon count', ph_count
     CALL PL_ESTIMATORS(ph_count,nphotons,fourpi,XMAX,l,n_spec,lum,e_tot,e_uva)
     print*,'Average number of scatterings = ',(nscatt/nphotons)
!!$
!!$
!!$  print*,e_tot, e_uva
!!$
!!$  print*,'Total UVB Photons Per Basal cell',e_tot-e_uva
!!$  print*,'Total UVA Photons Per basal cell',e_uva
!!$  print*,'Total Photons Absorbed per cell',e_tot
!!$
!!$  print*,'CPDs: UVB',(e_tot-e_uva)*0.05
!!$  print*,'CPDs: UVA',e_uva*0.0005
!!$  print*,'CPDs: total', e_uva*0.0005 + (e_tot-e_uva)*0.05

     tot(1)=e_tot
     uva(1)=e_uva

     print*,'******'
     print*,'******'
     deallocate(l,c)


endprogram mcpolar
!-------------------------------------------------------------------------------
