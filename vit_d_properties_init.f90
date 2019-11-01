module vit_d_properties_init_mod
  implicit none
  save
contains
  !!!!!!N
   !!!!ns=1.39-((wl-279)*0.0001)
  subroutine vit_d_properties_init()
    use vit_d_properties_mod
    use optical_properties_mod, only : nwl
    use load_spec2_mod
    use search_bisec_mod
    use interpolate_mod
    implicit none

    integer :: species
    character*70 :: filename
    real*8 :: spec(nwl)


    extinction=0.d0
    quantum_yield=0.d0
    !substances are numbered as follows:
    !------------1: provitamin D, 2: previtamin d, 3: lumisterol,
    !------------4: tachysterol 5: vitamin D 6: any toxi/supra sterol or other decay product
    species=1 !provitamin D3
      !Extinction
      filename='./vitamin_D_spectra/proD3_molar_absorptivity.txt'
      call load_spec2(filename,spectrum=spec,unit_to_cm=1.d0)
      extinction(species,:)=spec
      !quantum_yield FROM 1 to 2
      filename='./vitamin_D_spectra/pro-pre_quantum_yield.txt'
      call load_spec2(filename,spectrum=spec,unit_to_cm=1.d0)
      quantum_yield(1,2,:)=spec

    species=2!Previtamin D
      !Extincton
      filename='./vitamin_D_spectra/preD3_molar_absorptivity.txt'
      call load_spec2(filename,spec,1.d0)
      extinction(species,:)=spec
      !quantum yield FROM 2 to 1
      filename='./vitamin_D_spectra/pre-pro_quantum_yield.txt'
      call load_spec2(filename,spec,1.d0)
      quantum_yield(2,1,:)=spec
      !quantum yield FROM 2 to 3
      filename='./vitamin_D_spectra/pre-lum_quantum_yield.txt'
      call load_spec2(filename,spec,1.d0)
      quantum_yield(2,3,:)=spec
      !quantum yield FROM 2 to 4
      filename='./vitamin_D_spectra/pre-tachy_quantum_yield.txt'
      call load_spec2(filename,spec,1.d0)
      quantum_yield(2,4,:)=spec
      !quantum yield FROM 2 to 6
      filename='./vitamin_D_spectra/pre-tox_quantum_yield.txt'
      call load_spec2(filename,spec,1.d0)
      quantum_yield(2,6,:)=spec

    species=3 !lumisterol
      !Extinction
      filename='./vitamin_D_spectra/lumisterol_molar_absorptivity.txt'
      call load_spec2(filename,spec,1.d0)
      extinction(species,:)=spec
      !quantum yield FROM 3 to 2
      filename='./vitamin_D_spectra/lum-pre_quantum_yield.txt'
      call load_spec2(filename,spec,1.d0)
      quantum_yield(3,2,:)=spec
      !quantum yield FROM 3 to 6 : In the equation but not listed
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    species=4 !tachysterol
      !Extinction
      filename='./vitamin_D_spectra/tachysterol_molar_absorbivity.txt'
      call load_spec2(filename,spec,1.d0)
      extinction(species,:)=spec
      !quantum_yield from 4 to 2
      filename='./vitamin_D_spectra/tachy-pre_quantum_yield.txt'
      call load_spec2(filename,spec,1.d0)
      quantum_yield(4,2,:)=spec
      !quantum yield FROM 4 to 6 !In the equation biut not listed
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !-------------------------------------------------------
    species=5 !vitamin D
    !Extinction
    filename='./vitamin_D_spectra/vitaminD3_molar_absorptivity.txt'
    call load_spec2(filename,spec,1.d0)
    extinction(species,:)=spec
    !quantum_yield from 5 to 6
    filename='./vitamin_D_spectra/d3-tox_quantum_yield.txt'
    call load_spec2(filename,spec,1.d0)
    quantum_yield(5,6,:)=spec
    !-------------------------------------------------------


    return
  END SUBROUTINE vit_d_properties_init

end module vit_d_properties_init_mod
