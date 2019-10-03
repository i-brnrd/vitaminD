module op_prop_set_mod
  implicit none
  save
contains
  !!!!!!N
   !!!!ns=1.39-((wl-279)*0.0001)
  subroutine op_prop_set(u_a,u_s,g_skin)
    use packet_mod, ONLY: nwl, nlayer,lcount,l
    use load_spec2_mod
    use search_bisec_mod
    use interpolate_mod

    implicit none

    real*8,dimension(nlayer,nwl):: u_s,u_a,n_skin
    real*8 :: g_skin(nwl)
    real*8,dimension (nwl):: u_a_base, u_s_base,u_a_mel, u_a_epidermis,oxy_hgb,deoxy_hgb !base values from omlc "CHECK THIS CHECK THIS urgent****"
    real*8 spec(nwl),wl
    integer :: i_wl, layer
    real*8 :: jaques_basal

    character*30 :: filename
    real*8 lg10
    real*8 :: V_m, v_m_epi   ! boake
    real*8 n_mel !'how many' cells
    real*8 v_b,sO2, hgb

    print*, 'URGENT: op props need source/cited'
    !CAN IARRAY JUST SET TO ZERO WITHOUT THE COUNTING NONSENSE: yes probably
        u_s=0.
        u_a=0.
        n_skin=0.
        g_skin=0.
          lg10 = Log(10.)
          V_m=0.02 !p. 135 patricks thesis; skin type 1 only basal? very low ..
          v_m_epi=0.04 !eeeerrrrrgg
          n_mel=(lcount(2)*v_m) ! number of 'melanin voxels' in epidermis
          v_b=0.02   !blood volume only 2%?!?really?
          SO2=0.75
          hgb=150./66500.
          !INITIALISE THIS ELSEWHERE IT IS USED EVERYWHERE

 !THIS NEEDS CHECKED BALANCED AND TIDIED
  !isla everything here needs, seriously, to be checked and cited some of this might be wrong

      do i_wl=1,nwl
        wl=l(i_wl)
        g_skin(i_wl)=0.62+0.29*wl*(10.**(-3))!????????????????
        u_a_base(i_wl)=7.84*(10.**8)*(wl**(-3.255))  !??'Base' Values; Jaques et al 1998 skin optics OMLC News
        u_s_base(i_wl)=1.752*(10.**8)*(wl**(-2.33)) + 134.67*(wl**(-0.494))  !??'Base' Values; Jaques et al 1998 skin optics OMLC News
        u_a_mel(i_wl)=6.6*(10.**11)*(wl**(-3.33))  !what is this its omlc bu twhat
      enddo

    layer=1 !STRATUM CORNEUM
    !Absorption
    filename='./spectra/sc_abs_2.txt'
    call load_spec2(filename,spectrum=spec,unit_to_cm=10.d0)
    u_a(layer,:)=spec
    !scattering
    filename='./spectra/sc_scatt_2.txt'  !isla this looks horrific, really steppy and weird
    call load_spec2(filename,spec,10.d0)
    u_s(layer,:)=spec
    !--------------------------------------------------
    layer =2 !EPIDERMIS
    !Absorption
    filename='./spectra/epi_abs_2.txt'
    call load_spec2(filename,u_a_epidermis,10.d0)
    u_a(layer,:)=u_a_epidermis-(v_m_epi*u_a_mel)
    !Scattering
    filename='./spectra/epi_scatt_2.txt'
    call load_spec2(filename,spec,10.d0)
    u_s(layer,:)=spec
    !-------------------------------------------------------
    layer=3 !MELANIN LAYER
    !Absorption
    u_a(layer,:)= u_a_epidermis + (n_mel/real(lcount(3))) * u_a_mel
    !Scattering
    u_s(layer,:)=u_s(2,:)
    !-------------------------------------------------------
    layer=4 !DNA LAYER DNA: units?? is it out by 1000???
    !Absorption
    filename='./spectra/b924712b-f1.txt'
    call load_spec2(filename,spec,1.d0)
    u_a(layer,:)=lg10*0.0185*spec
    !Scattering
    u_s(layer,:)=u_s(2,:)
    !-------------------------------------------------------
    layer = 5 !DERMIS
    !Absorption
    filename='./spectra/ohb_as.txt'
    call load_spec2(filename,spec,1.d0)
    oxy_hgb=spec

    filename='./spectra/dhb_as.txt'
    call load_spec2(filename,spec,1.d0)
    deoxy_hgb=spec
    u_a(layer,:)=(v_b*hgb*(SO2*oxy_hgb + (1.d0-sO2)*deoxy_hgb)*lg10)+u_a_base*(1.d0-v_b)

    !scattering
    u_s(layer,:)=u_s_base

    return
  END SUBROUTINE op_prop_set

end module op_prop_set_mod
