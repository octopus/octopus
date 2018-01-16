!> @file
!!    Modulefile for handling basic definitions for the PSolver
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module PSbase
  use f_precisions
  implicit none
  private
  ! General precision, density and the potential types, to be moved in a low-levle module
  integer, parameter, public :: gp=f_double  !< general-type precision
  integer, parameter, public :: dp=f_double  !< density-type precision
end module PSbase
