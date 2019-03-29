module search_bisec_mod
  implicit none
  save

 contains
      subroutine search_bisec(v,x,lower,upper)
!     -------------------------------------------------------------
!     given value v; search an (**ascending order**)  array x(n) 
!     by bisection, returning the indices that bracket value v
!     -------------------------------------------------------------
      implicit none
      real*8, intent(in):: v, x(:)
      integer :: n,mid
      integer, intent(out)::lower,upper

      n=size(x)
      upper=n
      lower=1
      
      do while ((upper - lower) .gt. 1)
         mid= (upper + lower)/2
         if ( v .gt. x(mid)) then
            upper=upper
            lower=mid
         else
            upper=mid
            lower=lower
         endif
      enddo

      
      return
      end subroutine 
!***********************************************************************
end module search_bisec_mod
