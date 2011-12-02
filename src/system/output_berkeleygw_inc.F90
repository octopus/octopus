!! Copyright (C) 2011 D. Strubbe
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
!! $Id: output_etsf_inc.F90 5880 2009-09-03 23:44:44Z dstrubbe $

! ---------------------------------------------------------

subroutine output_berkeleygw_init(nst, bgw)
  integer, intent(in) :: nst
  type(output_bgw_t), intent(out) :: bgw

  PUSH_SUB(output_berkeleygw_init)

  !%Variable BerkeleyGW_NumberBands
  !%Type integer
  !%Default all states
  !%Section Output::BerkeleyGW
  !%Description
  !% Wavefunctions for bands up to this number will be output.
  !%End
  call parse_integer(datasets_check('BerkeleyGW_NumberBands'), nst, bgw%nbands)

  !%Variable BerkeleyGW_Vxc_diag_nmin
  !%Type integer
  !%Default 1
  !%Section Output::BerkeleyGW
  !%Description
  !% Lowest band for which to write diagonal exchange-correlation matrix elements.
  !%End
  call parse_integer(datasets_check('BerkeleyGW_Vxc_diag_nmin'), 1, bgw%vxc_diag_nmin)
  
  !%Variable BerkeleyGW_Vxc_diag_nmax
  !%Type integer
  !%Default nst
  !%Section Output::BerkeleyGW
  !%Description
  !% Highest band for which to write diagonal exchange-correlation matrix elements.
  !%End
  call parse_integer(datasets_check('BerkeleyGW_Vxc_diag_nmax'), nst, bgw%vxc_diag_nmax)
  
  !%Variable BerkeleyGW_Vxc_offdiag_nmin
  !%Type integer
  !%Default 1
  !%Section Output::BerkeleyGW
  !%Description
  !% Lowest band for which to write off-diagonal exchange-correlation matrix elements.
  !%End
  call parse_integer(datasets_check('BerkeleyGW_Vxc_offdiag_nmin'), 1, bgw%vxc_offdiag_nmin)
  
  !%Variable BerkeleyGW_Vxc_offdiag_nmax
  !%Type integer
  !%Default nst
  !%Section Output::BerkeleyGW
  !%Description
  !% Highest band for which to write off-diagonal exchange-correlation matrix elements.
  !%End
  call parse_integer(datasets_check('BerkeleyGW_Vxc_offdiag_nmax'), nst, bgw%vxc_offdiag_nmax)
  
  !%Variable BerkeleyGW_Complex
  !%Type logical
  !%Default false
  !%Section Output::BerkeleyGW
  !%Description
  !% Even when wavefunctions, density, and XC potential could be real in reciprocal space,
  !% they will be output as complex.
  !%End
  call parse_logical(datasets_check('BerkeleyGW_Complex'), .false., bgw%complex)
  
  !%Variable BerkeleyGW_WFN_filename
  !%Type string
  !%Default WFN
  !%Section Output::BerkeleyGW
  !%Description
  !% Filename for the wavefunctions.
  !%End
  call parse_string(datasets_check('BerkeleyGW_WFN_filename'), 'WFN', bgw%wfn_filename)

  POP_SUB(output_berkeleygw_init)
end subroutine output_berkeleygw_init

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
