module reflect_mod
  implicit none
  save
contains
  subroutine reflect(n1,n2,ct1,pi,r_flag)
    implicit none
    real*8,intent(in) :: n1,n2,pi !refractive index of  skin
    real*8,intent(in) :: ct1 !cos theta 1
    real*8 :: ca, ti, tt, si, ci, st, ct
    real*8 :: rs, rp, R
    real*8 :: ran
    integer r_flag

    r_flag=0
    CA=ASIN(1./N1)

    if (n1.gt.n2) then
       ca=asin(1./n1)
       if (acos(ct1).gt.ca) then
          print*,'reflect:TIR: doesnt do anything yet!no else or if or anything'
       endif
    endif
    
    !Fresnel Reflection 

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
       
       call random_number(ran)
       if (ran.lt.R) then
          r_flag=1
       endif

       return
          
  end subroutine reflect
   end module reflect_mod
