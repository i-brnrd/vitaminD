module load_spec2_mod
  implicit none
  save
  contains

    subroutine load_spec2(filename,spectrum,unit_to_cm)
      use get_dim_mod
      use load_mod
      use interpolate_mod
      use grid_mod, ONLY: nwl
      !--------------------------------------------------------------------
      ! loads a spectrum and interpolates BETWEEN 280-nwl in bins of 1nm
      ! if reading from spec is less than zero we return 0
      ! OPTIONAL ARGUMENT convert_to_cm.....
      !nwl, l(i) needs to be moduled
      !-------------------------------------------------------------------
      implicit none

      character*70, intent(in):: filename
      real*8, intent(in):: unit_to_cm
      integer :: length
      real*8, allocatable ::raw_data(:,:)
      real*8, intent(out) :: spectrum(nwl)
      real*8 :: l(nwl)
      real*8 :: interpol
      integer :: i


      !load in the raw data from the file
      call get_dim(filename,length)
      allocate(raw_data(2,length))
      call load(filename,length,raw_data)

      do i=1,nwl
           l(i)=280.+real(i) -1.
      enddo

      do i=1,nwl
         call interpolate(l(i),raw_data,interpol)
            if (interpol.lt.0.00) then
              interpol=0.00
              end if
              spectrum(i)=interpol*unit_to_cm

              enddo

      !      !deallocate(spectrum) ??
      return
    end subroutine load_spec2
end module load_spec2_mod
