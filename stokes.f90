      subroutine stokes(nxp,nyp,nzp,sint,cost,sinp,cosp,phi,&
                   fi,fq,fu,fv,pl,pc,sc,hgg,g2,pi,twopi,iseed)

      implicit none

      integer iseed
      real*8 nxp,nyp,nzp,sint,cost,sinp,cosp,phi,fi,fq,fu,fv
      real*8 pl,pc,sc,hgg,g2,pi,twopi,wght


      real*8 p1,p2,p3,p4,a11,a12,a13,a21,a22,a23,a24,a31,a32,a33,a34
      real*8 a42,a43,a44,a,rprob,si,sq,su,sv
      real*8 costp,sintp,phip
      real*8 bmu,b,ri1,ri3,cosi3,sini3,cosb2,sinbt,sini2,bott,cosdph
      real*8 cosi2,sin2i3,sin2i2,cos2i3,cos2i2,sin2,cos2,sin2cos1
      real*8 cos2sin1,cosi1,sini1,sin2i1,cos2i1

      real*8 ran

      wght=fi
      fi=fi/wght
      fq=fq/wght
      fu=fu/wght
      fv=fv/wght


!***** isotropic scattering ******************************
if(hgg.eq.0.) then
   call random_number(ran)
        cost=2.*ran-1.
        sint=(1.-cost*cost)
        if(sint.le.0.)then
          sint=0.
        else
          sint=sqrt(sint)
        endif
   call random_number(ran)
        phi=twopi*ran
        sinp=sin(phi)
        cosp=cos(phi)

        nxp=sint*cosp
        nyp=sint*sinp
        nzp=cost

        goto 100
      endif

      costp=cost
      sintp=sint
      phip=phi

!***** electron scattering ********************************
10 continue

!!$c      bmu=1.-2.*ran2(iseed)
!!$c      cosb2=bmu*bmu
!!$c      b=cosb2-1.
!!$c      p1=1.+cosb2
!!$c      p2=b
!!$c      p3=2.*bmu
!!$c      p4=0.
!!$c***********************************************************

      !***** dust scattering ********************************
      call random_number(ran)
      bmu=((1.+g2)-((1.-g2)/(1.-hgg+2.*hgg*ran))**2)/(2.*hgg)
      cosb2=bmu**2
      b=cosb2-1.
      call dustmat(p1,p2,p3,p4,bmu,cosb2,pl,pc,sc,hgg,g2,pi)
      a=p1
!***********************************************************

      if(abs(bmu).gt.1.) then
         if(bmu.gt.1.) then
            bmu=1.
            cosb2=1.
            b=0.
         else
            bmu=-1.
            cosb2=1.
            b=0.
         end if
      end if
      sinbt=sqrt(1.-cosb2)
      call random_number(ran)
      ri1=twopi*ran

      if(ri1.gt.pi) then
         ri3=twopi-ri1
         cosi3=cos(ri3)
         sini3=sin(ri3)
         sin2i3=2.*sini3*cosi3
         cos2i3=2.*cosi3*cosi3-1.
         a11=p1
         a12=p2*cos2i3
         a13=p2*sin2i3

!!$c****** for electron scattering **********
!!$c         rprob=a11*fi+a12*fq+a13*fu
!!$c         if((2.*ran2(iseed)).gt.rprob) goto 10
!!$c         a=rprob
!!$c******************************************

         if(bmu.eq.1.) then
            goto 100
         else
            if(bmu.eq.-1.) then
               fu=-fu
               goto 100
            end if
         end if

         cost=costp*bmu+sintp*sinbt*cosi3
         if(abs(cost).lt.1.) then
            sint=abs(sqrt(1.-cost*cost))
            sini2=sini3*sintp/sint
            bott=sint*sinbt
            cosi2=costp/bott-cost*bmu/bott
         else
            sint=0.
            sini2=0.
            if(cost.ge.1.)  cosi2=-1.
            if(cost.le.-1.) cosi2=1.
         end if

         cosdph=-cosi2*cosi3+sini2*sini3*bmu
         if(abs(cosdph).gt.1.) then
            if(cosdph.gt.1.) then
               cosdph=1.
            else
               cosdph=-1.
            end if
         end if
         phi=phip+acos(cosdph)
         if(phi.gt.twopi) phi=phi-twopi
         if(phi.lt.0.)    phi=phi+twopi

         sin2i2=2.*sini2*cosi2
         cos2i2=2.*cosi2*cosi2-1.
         sin2=sin2i2*sin2i3
         cos2=cos2i2*cos2i3
         sin2cos1=sin2i2*cos2i3
         cos2sin1=cos2i2*sin2i3

         a21=p2*cos2i2
         a22=p1*cos2-p3*sin2
         a23=p1*cos2sin1+p3*sin2cos1
         a24=-p4*sin2i2
         a31=-p2*sin2i2
         a32=-p1*sin2cos1-p3*cos2sin1
         a33=-p1*sin2+p3*cos2
         a34=-p4*cos2i2
         a42=-p4*sin2i3
         a43=p4*cos2i3
         a44=p3

!      elseif(ri1.le.pi) then
      else  

         cosi1=cos(ri1)
         sini1=sin(ri1)
         sin2i1=2.*sini1*cosi1
         cos2i1=2.*cosi1*cosi1-1.
         a11=p1
         a12=p2*cos2i1
         a13=-p2*sin2i1
!!$
!!$c************* for electron scattering ****************
!!$c         rprob=a11*fi+a12*fq+a13*fu
!!$c         if((2.*ran2(iseed)).gt.rprob) goto 10
!!$c         a=rprob
!!$c*******************************************************

         if(bmu.eq.1.) then
            goto 100
         else
            if(bmu.eq.-1.) then
               fu=-fu
               goto 100
            end if
         end if

         cost=costp*bmu+sintp*sinbt*cosi1
         if(abs(cost).lt.1.) then
            sint=abs(sqrt(1.-cost*cost))
            sini2=sini1*sintp/sint
            bott=sint*sinbt
            cosi2=costp/bott-cost*bmu/bott
         else
            sint=0.
            sini2=0.
            if(cost.ge.1.)  cosi2=-1.
            if(cost.le.-1.) cosi2=1.
         end if

         cosdph=-cosi1*cosi2+sini1*sini2*bmu
         if(abs(cosdph).gt.1.) then
            if(cosdph.gt.1.) then
               cosdph=1.
            else
               cosdph=-1.
            end if
         end if
         phi=phip-acos(cosdph)
         if(phi.gt.twopi) phi=phi-twopi
         if(phi.lt.0.)    phi=phi+twopi

         sin2i2=2.*sini2*cosi2
         cos2i2=2.*cosi2*cosi2-1.
         sin2=sin2i2*sin2i1
         cos2=cos2i2*cos2i1
         sin2cos1=sin2i2*cos2i1
         cos2sin1=cos2i2*sin2i1

         a21=p2*cos2i2
         a22=p1*cos2-p3*sin2
         a23=-p1*cos2sin1-p3*sin2cos1
         a24=p4*sin2i2
         a31=p2*sin2i2
         a32=p1*sin2cos1+p3*cos2sin1
         a33=-p1*sin2+p3*cos2
         a34=-p4*cos2i2
         a42=p4*sin2i1
         a43=p4*cos2i1
         a44=p3

      end if

      si=(a11*fi+a12*fq+a13*fu)/a
      sq=(a21*fi+a22*fq+a23*fu+a24*fv)/a
      su=(a31*fi+a32*fq+a33*fu+a34*fv)/a
      sv=(a42*fq+a43*fu+a44*fv)/a

      fi=si*wght
      fq=sq*wght
      fu=su*wght
      fv=sv*wght

      cosp=cos(phi)
      sinp=sin(phi)

      nxp=sint*cosp
      nyp=sint*sinp
      nzp=cost

100   continue

      return
      end subroutine


	subroutine dustmat(p1,p2,p3,p4,cost,cost2,pl,pc,sc,hgg,g2,pi)
!!$
!!$c				a.d. code october 25, 1989
!!$c				revised baw apr 10, 1990
!!$c	**********************************************************************
!!$c
!!$c	this program calculates the elements of the phase matrix for a
!!$c	simple representation of the mrn dust mixture using the algorithms
!!$c	for the ultraviolet region due to richard l. white ap.j. 229, 954,
!!$c	1979.
!!$c
!!$c	***********************************************************************
!!$c
!!$c	    cost = cos(angle) of scattering (i.e. angle between incident 
!!$c	        photon and scattered photon)
!!$c	    g = mean value of cosine of scattering angle (henyey-greenstein)
!!$c	    pl = peak linear polarization
!!$c	    pc = peak value of linear to circular conversion
!!$c	    sc = asymmetry of the circular polarization.
!!$c	    p1 = intensity phase function
!!$c	    p2 = polarization function
!!$c	    p3 = skew polarization
!!$c	    p4 = circular polarization
!!$c
!!$c	    the scattering matrix for (i,q,u,v) is of the form
!!$c
!!$c		 p1    p2    0	  0
!!$c		 p2    p1    0	  0
!!$c		 0     0     p3  -p4
!!$c		 0     0     p4   p3
!!$c
!!$c	**********************************************************************

      implicit none

      real*8 p1,p2,p3,p4,cost,cost2,pl,pc,sc,hgg,g2,pi
      real*8 phi,f,f2,c

      p1 = (1 - g2)/(1+g2-2*hgg*cost)**1.5
      p2 = -pl*p1*(1-cost2)/(1+cost2)
      p3 = p1*2*cost/(1+cost2)
      if(abs(cost).gt.1.) then 
        print *, 'in dustmat, cost.gt.1',cost
        if(cost.gt.1) then     
           cost=1.
        else    
           cost=-1.     
        end if 
      end if
!     angle in degrees!
      phi=acos(cost)*180./pi
      f=3.13*phi*exp(-7.0*phi/180.)
!     now convert to radians
      f2=(phi+sc*f)*pi/180.
!     fphi= (1+3.13*sc*exp(-7.0*phi/pi))*phi
      c=(cos(f2))**2
      p4 = -pc*p1*(1-c)/(1+c)
!!$
!!$c******* Isotropic scattering
!!$c      p1 = 1.0
!!$c      p3 = 1.0
!!$c      p2 = 0.0
!!$c      p4 = 0.0  

	return
	end subroutine


