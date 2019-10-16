module n_interface_mod 
  implicit none
  save
contains
!this module is caluclating Frensel reflectio probabilities
!ONLY for a flat surface of the skin.....
!this is incredibly un abstract 
  subroutine n_interface(n1,n2,ct1,diffuse_flag,reflect_flag)
    use constants_mod, only: pi
    implicit none
  real*8,intent(in) :: n1,n2 !refractive indices 
  real*8,intent(in) :: ct1 !cos theta 1
  logical, intent(in) :: diffuse_flag
  logical, intent(out):: reflect_flag
  real*8 :: Ca
  real*8 :: ti, tt, si, ci, st, ct
  real*8 :: rs, rp, R
  real*8 :: ran
  real*8 :: st1,st2,ct2,rf 

  real*8 :: theta1, theta2, rrs, rrp, rrr

  logical :: TIR 
    reflect_flag=.false.

    TIR=.false.


  if (ct1.lt.0.d0) then !it means pkt is entering simulation 
    if (.not.diffuse_flag) then 
      print*, 'not a diffuse packet, fresnel reducing to...'
      R=((n1-n2)/(n1+n2))**2
      print*, R, 'do i need to do anything to reflect flag?? '

    else !diffuse flag is true; 

      print*, 'is a diffuse packet, coming from source need to do a frensel reflection'
      theta1=dacos(ct1)
      theta1=pi-theta1 
      !do a fresnel with the alterest theta 
      call fresnel(n1,n2,theta1,r) 
      print*,r 
    endif 

  else 
      print*,'internal reflaction is happening '
      !@ct1>0; bouncing back UP, 
     ! do a fresnel  or a TI
     if (n2.lt.n1) then
      print*, 'TIR n1:',n1,'n2',n2,'just check right way round. '
      ca=dasin(n2/n1)
      print*, 'critical angle', ca, 'costheta',dacos(ct1)
      if (acos(ct1).gt.ca) then
        print*,'reflect:TIR!!!!what do I need to do to the flag??'
        TIR=.true.
      else 
        theta1=dacos(ct1)
        print*, 'nternal, not TIR,do a fresnel'
        call fresnel(n1,n2,theta1,r)
        print*,r
      endif 
  endif
    
 endif



  print*, 'reflect fresnetl prob', R, 'TIR FLAG', TIR 
  if (TIR) then
    reflect_flag=.true.
  else
    call random_number(ran)
    if (ran.lt.R) then
      reflect_flag=.true.
    endif
  endif 

    return
          
  end subroutine n_interface

  subroutine fresnel(n1,n2,theta1,r)
  implicit none 
  real*8, intent(in) :: n1,n2,theta1

  real*8 :: rrs,rrp,rrr,r, theta2

    theta2=(n1/n2)*SIN(theta1)

    rrs = (n1*cos(theta1) - n2*cos(theta2))/(n1*cos(theta1)+ n2*cos(theta2))
    rrs=rrs**2 

    rrp=(n1*cos(theta2) - n2*cos(theta1))/(n1*cos(theta2) + n2*cos(theta1))
    rrp=rrp**2 
    print*,'rrs,rrp',  rrs,rrp
    rrr=0.5*(rrs+rrp)

    print*, rrr

    r=rrr

    if (r.gt.1)print*, 'R greaater than 1 in fresnel huge error argh'

  end subroutine fresnel



  subroutine reflect(z_dir)
  use packet_mod, only: cost
  implicit none 
  real*8, intent(inout) :: z_dir
  z_dir= -z_dir
  cost=z_dir
  end subroutine

end module n_interface_mod
