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
  real*8,intent(in) :: ct1 !cos theta 1, or nzp, same thing
  logical, intent(in) :: diffuse_flag
  logical, intent(out):: reflect_flag

  real*8 :: ca, theta1
  real*8 :: R
  real*8 :: ran
  logical :: TIR

  !initialise decision vars to 0
  reflect_flag=.false.
  TIR=.false.
  R=0.00

  if (ct1.lt.0.d0) then !it means pkt is entering medium
    if (.not.diffuse_flag) then !if packet is not diffuse, used reduced fresnel eqn.
      R=((n1-n2)/(n1+n2))**2
    else !pkt is difffuse, use full fresnel. swap quadrant of theta to fix direction frame of packet
      theta1=dacos(ct1)
      theta1=pi-theta1
      call fresnel(n1,n2,theta1,r)
    endif
  else !packet is within medium- need to check if it leaves or reflects off skin surface
     if (n2.lt.n1) then
      ca=dasin(n2/n1)
      if (acos(ct1).gt.ca) then !if angle of incidence > crit angle, Total Internal Reflection
        TIR=.true.
      else  !else call fresnel
        theta1=dacos(ct1)
        call fresnel(n1,n2,theta1,r)
      endif
  endif

 endif
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
  !check probability of reflection for an incident packet
  implicit none
  real*8, intent(in) :: n1,n2,theta1
  real*8 :: theta2
  real*8 :: rs,rp,r

    theta2=(n1/n2)*SIN(theta1)

    rs = (n1*cos(theta1) - n2*cos(theta2))/(n1*cos(theta1)+ n2*cos(theta2))
    rs=rs**2

    rp=(n1*cos(theta2) - n2*cos(theta1))/(n1*cos(theta2) + n2*cos(theta1))
    rp=rp**2

    r=0.5*(rs+rp)


    if (r.gt.1) print*, 'R greater than 1 in fresnel: please investigate broken physics R:',r,'theta1:',theta1,'theta2',theta2

  end subroutine fresnel



  subroutine reflect(z_dir)
    !WARNING: takes in single direction(z) and uses to reset cost
  use packet_mod, only: cost
  implicit none
  real*8, intent(inout) :: z_dir
  z_dir= -z_dir
  cost=z_dir
  end subroutine

end module n_interface_mod
