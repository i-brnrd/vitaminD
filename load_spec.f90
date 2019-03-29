module load_spec_mod
  implicit none
  save
  contains

    subroutine load_spec(fname,sz,spectrum)
      use get_dim_mod
      use load_mod
      !--------------------------------------------------------------------
      ! loads a spectrum 
      !-------------------------------------------------------------------
      implicit none
      character*25, intent(in):: fname
      integer, intent(out)::sz
      real*8, allocatable, intent(out)::spectrum(:,:)

      
      !deallocate(spectrum)
      call get_dim(fname,sz)
      allocate(spectrum(2,sz))
      call load(fname,sz,spectrum)

      return 
    end subroutine load_spec
end module load_spec_mod 

      
