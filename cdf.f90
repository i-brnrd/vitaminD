      
subroutine cdf(c,spec,ndim)
  !-------------------------------------------------
  !given a spectrum, returns (CDF)
  !uses trapezoidal integration 
  !-------------------------------------------------
  implicit none
  
  integer, intent(in):: ndim
  real*8, intent(in):: spec(2,ndim)
  real*8, intent(inout):: c(ndim)
  real*8 :: f(ndim),l(ndim)

  integer :: i

  !initialise each array
  do i=1,ndim
     l(i)=0.
     f(i) = 0.
     c(i) = 0.
  enddo

  l(:)=spec(1,:)
  f(:)=spec(2,:)
 
  !numerical integration to produce cumulative distribution
  do i=1,ndim-1
     c(i+1)= c(i) + 0.5*( f(i+1) + f(i) ) * ( l(i+1) - l(i) )
  enddo
 
  !normalise cumulative distribution to get cdf
  c=c/c(ndim)
 
  return 
end subroutine cdf
