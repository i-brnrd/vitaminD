module interpolate_mod
  implicit none
  save
contains

  subroutine interpolate(wl,array,val)
    use search_bisec_mod
    implicit none

    real*8, intent(in) :: wl
    real*8, intent(in) :: array(:,:)

    real*8, intent(out) :: val
    integer :: l,u 


    ! search by bisection- get x bounds
    call search_bisec(wl,array(1,:),l,u)

    
    !use interpolation to find corresponding y
    val=(wl-array(1,l))/(array(1,u)-array(1,l))
    val=(array(2,u)-array(2,l))*val
    val=array(2,l) + val


    return
  end subroutine interpolate
end module interpolate_mod
