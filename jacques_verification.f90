module jacques_verification_mod
  implicit none
  save
contains

subroutine jacques_verification(diffuse_fraction, tot_irr,l,incident_spec_irr)
use optical_properties_mod, only: nwl,nlayer,u_s_all,u_a_all,g_all,n_all 
implicit none
real*8, intent(out) :: diffuse_fraction, tot_irr
real*8, dimension(nwl) :: l,incident_spec_irr
integer:: i 

print*, 'Jaques verification mode ON: nwls:', nwl
print*, 'Setting optical properties....'
  
diffuse_fraction = 0.d0

l(1)= 630.d0
incident_spec_irr(1)=1.d0
tot_irr=1.d0
  


g_all(1)=0.9
n_all(1)=1.38

u_a_all(:,1)=0.23
u_s_all(:,1)=(21.)/(1.-g_all(1))


print*, 'diffuse fraction', diffuse_fraction
print*,'total irradiance',tot_irr
print*,'spectral irradiance', incident_spec_irr


end subroutine
end module