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
!! $Id: opt_control.F90 3099 2007-07-23 14:21:35Z lorenzen $

#include "global.h"

module opt_control_constants_m

  implicit none

  integer, parameter ::  &
    oct_is_groundstate   = 1,      &
    oct_is_excited       = 2,      &
    oct_is_superposition = 3,      &
    oct_is_userdefined   = 4         
    
  integer, parameter ::  &
    oct_tg_groundstate   = 1,      &
    oct_tg_excited       = 2,      &
    oct_tg_superposition = 3,      &
    oct_tg_userdefined   = 4,      &
    oct_tg_local         = 5,      &        
    oct_tg_td_local      = 6   

  integer, parameter ::  &
    oct_algorithm_zbr98 = 1,       &
    oct_algorithm_zr98  = 2,       &
    oct_algorithm_wg05  = 3       
 
  integer, parameter ::  &
    oct_targetmode_static = 0,     &
    oct_targetmode_td     = 1
  
  integer, parameter ::  &
    oct_tgtype_state = 1,          &
    oct_tgtype_local = 2
  
  integer, parameter ::  &
    oct_ftype_gauss = 1,           &
    oct_ftype_step  = 2

end module opt_control_constants_m
