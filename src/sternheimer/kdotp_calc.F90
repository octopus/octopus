!!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: kdotp_calc.F90 2548 2006-11-06 21:42:27Z xavier $

#include "global.h"

module kdotp_calc_m
  use datasets_m
  use elf_m
  use functions_m
  use grid_m
  use global_m
  use hamiltonian_m
  use loct_parser_m
  use linear_response_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use pert_m
  use profiling_m
  use states_m
  use sternheimer_m
  use system_m
  use xc_m

  implicit none

  private
  public ::                        &
     kdotp_wfs_tag,                &
     kdotp_rho_tag

contains

  character(len=100) function kdotp_rho_tag(dir) result(str)
    integer, intent(in) :: dir

    !this function has to be consistent with oct_search_file_lr in liboct/oct_f.c

    call push_sub('kdotp_calc.kdotp_rho_tag')

    write(str, '(a,i1)') 'rho_', dir

    call pop_sub()

  end function kdotp_rho_tag
  
  character(len=100) function kdotp_wfs_tag(dir) result(str)
    integer, intent(in) :: dir 

    call push_sub('kdotp_calc.kdotp_wfs_tag')

    write(str, '(a,i1)') "wfs_", dir

    call pop_sub()

  end function kdotp_wfs_tag
  
end module kdotp_calc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
