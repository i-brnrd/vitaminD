module iarray_mod
implicit none
  save
contains

subroutine iarray
  use optical_properties_mod
  use grid_mod
  implicit none

  !**** Initialize array values to be zero
  xface =0.
  yface=0.
  zface=0.

  mask=0.
  rhokap=0.
  pl_sum=0.

  n_pkt_wl=0
  l=0.
  lcount=0

  u_s_all=0.
  u_a_all=0.
  g_all=0.
  n_all=0.



  return
end subroutine iarray
end module iarray_mod
