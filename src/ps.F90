#include "config_F90.h"

module ps

#if defined(ONE_D)
#  include "ps1D.F90"
#elif defined(THREE_D)
#  include "ps3D.F90"
#endif

subroutine derivate_in_log_grid(a, b, nrval, f, dfdr)
  real(r8), intent(in) :: a, b
  integer,  intent(in) :: nrval
  real(r8), intent(IN) :: f(nrval)
  real(r8), intent(out) :: dfdr(nrval)

  real(r8) :: x,y
  integer :: i

  x = 1.0_r8 - exp(-2*a)
  y = 1.0_r8 - exp(-a)

  dfdr(1) = (1/(y*b))*exp(-a)*(f(2)-f(1))  
  do i = 2, nrval-1
    dfdr(i) = (1/(x*b))*exp(-i*a)*(f(i+1)-f(i-1))
  enddo
  dfdr(nrval) = (1/(y*b))*exp(-(nrval-1)*a)*(f(nrval)-f(nrval-1))

  return
end subroutine derivate_in_log_grid

end module ps
