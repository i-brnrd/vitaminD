module fresnel_mod
  implicit none
  save
contains
  subroutine fresnel(n1,n2,ct1,reflect_flag)
    use constants_mod, only: pi
    implicit none
    real*8,intent(in) :: n1,n2!refractive indices 
    real*8,intent(in) :: ct1 !cos theta 1
    real*8 :: ca, ti, tt, si, ci, st, ct
    real*8 :: rs, rp, R
    real*8 :: ran
    logical :: reflect_flag

    real*8 :: st1,st2,ct2,rf 

    reflect_flag=.false.
    CA=ASIN(1./N1)

    !check for TIR 
    if (n1.gt.n2) then
       ca=asin(1./n1)
       if (acos(ct1).gt.ca) then
          print*,'reflect:TIR: doesnt do anything yet!no else or if or anything but'
          stop
       endif
    endif
    
    ! !Fresnel Reflection 
    ! ST1=SIN(ACOS(CT1))
    ! ST2=ST1*N1
    ! CT2=COS(ASIN(ST2))

    !  RF = ((N1*CT1-CT2)/(N1*CT1 + CT2))**2
    !  RF = RF + ((N1*CT2-CT1)/(N1*CT2 + CT1))**2
    !  RF = RF/2.

    ti=acos(ct1)
    ti=pi-ti
    
    si=sin(ti)
    ci=cos(ti)
    
    st= (n1/n2)*si
    tt=asin(st)
    ct=cos(tt)

    rs= ((n1*ci - n2*ct)/(n1*ci + n2*ct))**2
    rp= ((n1*ct - n2*ci)/(n1*ct + n2*ci))**2

   R=0.5*(rs+rp)

   print*, 'reflect fresnetl prob', R, RF
       
       call random_number(ran)
       if (ran.lt.R) then
          reflect_flag=.true.
       endif

       return
          
  end subroutine fresnel

  subroutine reflect()
  use packet_mod, only: cost, nzp 
  implicit none 
  nzp=-nzp
  cost=nzp

  print*, 'REFLECT IS FINE.... '
  end subroutine

end module fresnel_mod
