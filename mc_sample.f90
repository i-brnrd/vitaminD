module mc_sample_mod
  implicit none
  save
contains
 

  subroutine mc_sample(C,L,WL)
    use search_bisec_mod
  !****---------------------------
  ! ****Returns sampled WL in nm
  !***----------------------------
  implicit none

  real*8, intent(in) :: l(:),c(:)
  real*8, intent(out) :: wl
  real*8 :: ran
  integer:: low, up


  call random_number(ran)
  call search_bisec(ran,c,low,up)

  !use interpolation to find corresponding wavelength
  wl = (l(up) - l(low)) * (ran-c(low)) / (c(up)-c(low))
  wl = l(low) + wl

  return 
end subroutine mc_sample

end module 
