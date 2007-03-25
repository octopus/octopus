#include "global.h"


subroutine doperate(np, nn, w, opi, fi, fo)
  implicit none

  integer, intent(in) :: np
  integer, intent(in) :: nn
  FLOAT,   intent(in) :: w(1:nn)
  integer, intent(in) :: opi(1:nn, 1:np)
  FLOAT,   intent(in) :: fi(1:np)
  FLOAT,   intent(out):: fo(1:np) 

  integer :: ii

  do ii = 1, np
    fo(ii) = sum(w(1:nn)  * fi(opi(1:nn, ii)))
  end do

end subroutine doperate


subroutine zoperate(np, nn, w, opi, fi, fo)
  implicit none

  integer, intent(in) :: np
  integer, intent(in) :: nn
  FLOAT,   intent(in) :: w(1:nn)
  integer, intent(in) :: opi(1:nn, 1:np)
  CMPLX,   intent(in) :: fi(1:np)
  CMPLX,   intent(out):: fo(1:np) 

  integer :: ii

  do ii = 1, np
    fo(ii) = sum(w(1:nn)  * fi(opi(1:nn, ii)))
  end do
  
end subroutine zoperate

