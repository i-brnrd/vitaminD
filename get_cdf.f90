
subroutine get_cdf(cdf,x_vals,y_vals,ndim)
  !-------------------------------------------------
  !given x&y values returns (CDF)
  !uses trapezoidal integration
  !-------------------------------------------------
  implicit none

  integer, intent(in):: ndim
  real*8, intent(in):: x_vals(ndim),y_vals(ndim)
  real*8, intent(inout):: cdf(ndim)
  integer :: i

  !initialise each array
  cdf=0.d0

  !numerical integration to produce cumulative distribution
  do i=1,ndim-1
     cdf(i+1)= cdf(i) + 0.5*( y_vals(i+1) + y_vals(i) ) * ( x_vals(i+1) - x_vals(i) )
  enddo

  !normalise cumulative distribution to get cdf
  cdf=cdf/cdf(ndim)

  return
end subroutine get_cdf
