module get_dim_mod
  implicit none
  save
contains


  subroutine get_dim(name,sz)

    !*****************************
    !get dim of input file
    !ALWAYS assuming 1 header line
    !*****************************

    implicit none 

    integer, intent(out):: sz
    character*25:: name
    integer :: i,io
    !*****************************

    sz=0
    open(8, file=name,status="old")
    read(8,*,IOSTAT=io)  !line 1 is header
    do 
       read(8,*,IOSTAT=io)
       if (io/=0) exit
       sz=sz+1
    enddo
    close(8)

  end subroutine get_dim


end module get_dim_mod
  
