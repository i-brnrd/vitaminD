module load_mod
  implicit none
  save
contains

  subroutine load(fname,sz,spectrum)
    !--------------------------------------------------------------------
    ! loads a spectrum
    !-------------------------------------------------------------------

    implicit none
    character*70, intent(in):: fname
    integer, intent(in)::sz
    real*8, intent(inout)::spectrum(2,sz)

    integer i

    !initialise allocated array
    do i=1,sz
       spectrum(1,i)=0.
       spectrum(2,i)=0.
    enddo

    open(8, file=fname,status="old")

    read(8,*)                 !line 1 is header
    do i=1,sz
       read(8,*) spectrum(1,i), spectrum(2,i)
    enddo
    close(8)

    return
  end subroutine load


end module load_mod
