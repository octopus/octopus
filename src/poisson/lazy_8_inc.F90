!! Copyright (C) 2011 X. Andrade
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

!!****h* BigDFT/lazy_8
!! NAME
!!   lazy_8
!!
!! FUNCTION
!!   Filters for interpolating scaling functions (order 8)
!!
!! SOURCE
!!
integer, parameter :: m=10
real(kind=8), dimension(-m:m) ::  ch,cg,cht,cgt

!******** coefficients for wavelet transform *********************
do i=-m,m
  ch(i)=0.d0
  cht(i)=0.d0
  cg(i)=0.d0
  cgt(i)=0.d0
end do

! The normalization is chosen such that a constant function remains the same constant
! on each level of the transform

ch(-7)=-5.d0/2048.d0
ch(-6)=0.d0
ch(-5)=49.d0/2048.d0
ch(-4)=0.d0
ch(-3)=-245.d0/2048.d0
ch(-2)=0.d0
ch(-1)=1225.d0/2048.d0
ch( 0)=1.d0
ch( 1)=1225.d0/2048.d0
ch( 2)=0.d0
ch( 3)=-245.d0/2048.d0
ch( 4)=0.d0
ch( 5)=49.d0/2048.d0
ch( 6)=0.d0
ch( 7)=-5.d0/2048.d0
!
cht( 0)=1.d0

! g coefficients from h coefficients
do i=-m,m-1
  cg(i+1)=cht(-i)*(-1)**(i+1)
  cgt(i+1)=ch(-i)*(-1)**(i+1)
end do
!!***

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
