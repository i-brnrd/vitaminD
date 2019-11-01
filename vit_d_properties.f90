module vit_d_properties_mod
  use optical_properties_mod, only: nwl
  implicit none
  save
  !contains
  !following van dijk et al


  integer, parameter :: n_species = 6 ! substances are numbered as follows:
  !------------1: provitamin D, 2: previtamin d, 3: lumisterol,
  !------------4: tachysterol 5: vitamin D 6: any toxi/supra sterol or other decay product
  real*8  :: extinction(n_species,nwl)

  real*8 :: quantum_yield(n_species,n_species,nwl)

end module vit_d_properties_mod
