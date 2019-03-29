module op_props_mod
  implicit none
  save
contains

  subroutine op_props(wl,u_a,u_s,g,epi,eumel,phmel,dna,ohb,dhb,stratc,epi_s,sc_s)
    use grid_mod
    use search_bisec_mod
    use interpolate_mod

    implicit none

    real*8, intent(in):: wl
    real*8, intent(out):: u_a(:),u_s(:),g

    real*8, dimension(:,:), intent(in):: epi,eumel,phmel,dna,ohb,dhb,stratc,epi_s,sc_s


    real*8 u_a_base, u_s_base !base values from omlc

    integer sz

    real*8 e_sc
    real*8 e_epi !extinction coefficients
    real*8 e_eu, e_ph, e_mel, V_m, v_m_epi
    real*8 e_dna
    real*8 e_dhb, e_ohb, e_hgb
    real*8 v_b,s02

    real*8 sc_scatt
    real*8 epi_scatt

    real*8 n_mel !'how many' cells

    real*8 lg10



    !HG FUCNTION VALUE
    g=0.62+0.29*WL*(10.**(-3))
    !'Base' Values
    u_a_base= 7.84*(10.**8)*(wl**(-3.255))
    u_s_base=1.752*(10.**8)*(wl**(-2.33)) + 134.67*(wl**(-0.494))
    !u_s_base=(u_s_base)/(1.-g)

    lg10 = Log(10.)


    !SRATUM CORNEUM
    call interpolate(wl,stratc,e_sc)
    call interpolate(wl,sc_s,sc_scatt)
    !EPIDERMIS
    call interpolate(wl, epi, e_epi)
    call interpolate(wl,epi_s,epi_scatt)



    !MELANINLAYER
    V_m=0.02 !p. 135 patricks thesis; skin type 1 only basal? very low ..
    v_m_epi=0.04
!!$    call interpolate(wl,eumel,e_eu)
!!$    e_eu=e_eu * 80. !80 g/L is concentration
!!$    call interpolate(wl,phmel,e_ph)
!!$    e_ph=e_ph * 12. !12 g/L is concentration
!!$    e_mel=e_eu + e_ph
!!$    e_mel=e_mel*lg10


    e_mel=6.6*(10.**11)*(wl**(-3.33))

    !DNA
    call interpolate(wl,dna,e_dna)

    !DERMIS
    v_b=0.02   !blood volume only 2%?!?really?
    S02=0.75
    call interpolate(wl,ohb,e_ohb)
    call interpolate(wl,dhb,e_dhb)
    e_ohb=(e_ohb/66500.)*150.
    e_dhb=(e_dhb/66500.)*150.
    e_hgb=e_ohb*S02+e_dhb*(1.-s02)
    e_hgb=e_hgb*lg10

    !ASSIGN PROPERTIES

    n_mel=(lcount(2)*v_m) ! number of 'melanin voxels' in epidermis


    !STRATUM CORNEUM
    u_a(1)=e_sc
    u_s(1)=sc_scatt
    !u_s(1)=u_s_base

    !EPIDERMIS
    !u_a(2)=e_sc
    u_a(2)=e_epi - v_m_epi*e_mel
    u_s(2)=epi_scatt
    !MELANIN LAYER
    u_a(3)=e_epi + (n_mel/lcount(3))*e_mel

    !u_a(3)=e_mel*v_m + u_a_base*(1.-V_m)!0.28 + u_a(1)*(1.-0.28) !u_a_base*(1.-V_m)
    !u_a(2)=(n_mel/mcount)*e_mel  + (1-(n_mel/mcount))*u_a(1)
    u_s(3)=epi_scatt

    !BASAL LAYER
    u_a(4)=e_dna*lg10*0.0185
    u_s(4)=u_s_base!0.0!u_s_base

    !DERMIS
    u_a(5)=e_hgb*V_b + u_a_base*(1.-v_b)
    u_s(5)=u_s_base




    return
  END SUBROUTINE op_props
end module op_props_mod
