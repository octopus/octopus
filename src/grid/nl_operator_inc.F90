!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

! The next two routines are wrappers around the operator sum.
! The reason for those wrappers is that, by help of the parameter
! list, compilers can generate more efficient code because of they then
! know the arrays not to be shaped.
!
! ---------------------------------------------------------
subroutine X(operate)(np, np_part, nn, w, opi, fi, fo)
  integer, intent(in) :: np
  integer, intent(in) :: np_part
  integer, intent(in) :: nn
  FLOAT,   intent(in) :: w(1:nn)
  integer, intent(in) :: opi(1:nn, 1:np)
  R_TYPE,   intent(in) :: fi(1:np_part)
  R_TYPE,   intent(out):: fo(1:np) 
  
  integer :: ii
  !$omp parallel do
  do ii = 1, np
    fo(ii) = sum(w(1:nn)  * fi(opi(1:nn, ii)))
  end do
  !$omp end parallel do
  
end subroutine X(operate)

! ---------------------------------------------------------
subroutine X(operate_olu)(np, np_part, nn, w, opi, fi, fo)
  integer, intent(in) :: np
  integer, intent(in) :: np_part
  integer, intent(in) :: nn
  FLOAT,   intent(in) :: w(1:nn)
  integer, intent(in) :: opi(1:nn, 1:np)
  R_TYPE,  intent(in) :: fi(1:np_part)
  R_TYPE,  intent(out):: fo(1:np) 
  
  integer :: ii, jj, kk
  R_TYPE :: a0, a1, a2, a3, a4, a5, a6, a7, a8, a9, aa, ab

#ifdef R_TCOMPLEX
#define DEPTH 6
#else
#define DEPTH 12
#endif

  !$omp parallel do private(kk, a0, a1, a2, a3, a4, a5 a6, a7, a8, a9, aa, ab)
  do ii = 1, (np - DEPTH + 1), DEPTH
    
    a0 = w(1) * fi(opi(1, ii  ))
    a1 = w(1) * fi(opi(1, ii+1))
    a2 = w(1) * fi(opi(1, ii+2))
    a3 = w(1) * fi(opi(1, ii+3))
    a4 = w(1) * fi(opi(1, ii+4))
    a5 = w(1) * fi(opi(1, ii+5))
#if DEPTH > 6
    a6 = w(1) * fi(opi(1, ii+6))
    a7 = w(1) * fi(opi(1, ii+7))
    a8 = w(1) * fi(opi(1, ii+8))
    a9 = w(1) * fi(opi(1, ii+9))
    aa = w(1) * fi(opi(1, ii+10))
    ab = w(1) * fi(opi(1, ii+11))
#endif
    
    do kk= 2, nn
      a0 = a0 + w(kk) * fi(opi(kk, ii  ))
      a1 = a1 + w(kk) * fi(opi(kk, ii+1))
      a2 = a2 + w(kk) * fi(opi(kk, ii+2))
      a3 = a3 + w(kk) * fi(opi(kk, ii+3))
      a4 = a4 + w(kk) * fi(opi(kk, ii+4))
      a5 = a5 + w(kk) * fi(opi(kk, ii+5))
#if DEPTH > 6
      a6 = a6 + w(kk) * fi(opi(kk, ii+6))
      a7 = a7 + w(kk) * fi(opi(kk, ii+7))
      a8 = a8 + w(kk) * fi(opi(kk, ii+8))
      a9 = a9 + w(kk) * fi(opi(kk, ii+9))
      aa = aa + w(kk) * fi(opi(kk, ii+10))
      ab = ab + w(kk) * fi(opi(kk, ii+11))
#endif
    end do

      fo(ii)   = a0
      fo(ii+1) = a1
      fo(ii+2) = a2
      fo(ii+3) = a3
      fo(ii+4) = a4
      fo(ii+5) = a5
#if DEPTH > 6
      fo(ii+6) = a6
      fo(ii+7) = a7
      fo(ii+8) = a8
      fo(ii+9) = a9
      fo(ii+10) = aa
      fo(ii+11) = ab
#endif

  end do
  !$omp end parallel do
  
  do jj = ii, np
    fo(jj) = sum(w(1:nn)  * fi(opi(1:nn, jj)))
  end do

#undef DEPTH
end subroutine X(operate_olu)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
