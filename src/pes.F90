!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module PES_m
#if !defined(DISABLE_PES) && defined(HAVE_FFT)
  use global_m
  use lib_oct_parser_m
  use io_m
  use units_m
  use fft_m
  use mesh_m
  use states_m
  use simul_box_m

  implicit none

  type PES_rc_t
    integer          :: npoints   ! how many points we store the wf
    integer, pointer :: points(:) ! which points to use
    character(len=30), pointer :: filenames(:) ! filenames
    complex(r4), pointer :: wf(:,:,:,:,:)
  end type PES_rc_t

  type PES_mask_t
    CMPLX, pointer :: k(:,:,:,:,:,:) ! masked wf in momentum space
    FLOAT, pointer :: r(:,:,:,:,:)   ! summed masked density in real space

    type(fft_t) :: fft
  end type PES_mask_t

  type PES_t
    logical :: calc_rc
    type(PES_rc_t)   :: rc

    logical :: calc_mask
    type(PES_mask_t) :: mask
  end type PES_t

contains

  ! ---------------------------------------------------------
  subroutine PES_init(p, m, sb, st, ab, save_iter)
    type(pes_t),    intent(out)   :: p
    type(mesh_t),   intent(inout) :: m
    type(simul_box_t), intent(in) :: sb
    type(states_t),    intent(in) :: st
    integer,           intent(in) :: ab, save_iter

    !%Variable CalcPES_rc
    !%Type logical
    !%Default no
    !%Section Time Dependent::PES
    !%Description
    !%  If <tt>true</tt>, store the wave functions in specific points in order to 
    !% calculate the photo-electron spectrum in a point far in the box as proposed in 
    !% A. Pohl, P.-G. Reinhard, and E. Suraud Phys. Rev. Lett. <b>84</b>, 5090 (2000).
    !%End
    call loct_parse_logical(check_inp('CalcPES_rc'), .false., p%calc_rc)
    if(p%calc_rc) then
      p%calc_rc = .true.
      call PES_rc_init(p%rc, m, st, save_iter)
    end if

    p%calc_mask = .false.
    ! have the mask, and we are working in the velocity gauge
    if(ab == 2) then
      !%Variable CalcPES_mask
      !%Type logical
      !%Default no
      !%Section Time Dependent::PES
      !%Description
      !% If <tt>true</tt>, calculate the photo-electron spectrum using the mask method
      !% (M. Marques, D. Varsano, H. Appel, E.K.U. Gross and A. Rubio to be submitted). 
      !% In order for this to work, masking boundaries are necessary 
      !% (<tt>AbsorbingBoundaries == 2</tt>).
      !%End
      call loct_parse_logical(check_inp('CalcPES_Mask'), .false., p%calc_mask)
      if(p%calc_mask) then
        call PES_mask_init(p%mask, m, sb, st)
      end if
    end if

  end subroutine PES_init


  ! ---------------------------------------------------------
  subroutine PES_end(p)
    type(PES_t), intent(inout) :: p

    if(p%calc_rc)   call PES_rc_end  (p%rc)
    if(p%calc_mask) call PES_mask_end(p%mask)

  end subroutine PES_end


  ! ---------------------------------------------------------
  subroutine PES_doit(p, m, st, ii, dt, mask)
    type(PES_t),    intent(inout) :: p
    type(mesh_t),   intent(in)    :: m
    type(states_t), intent(in)    :: st
    FLOAT,          intent(in)    :: dt
    FLOAT,          pointer       :: mask(:)
    integer,        intent(in)    :: ii

    if(p%calc_rc)   call PES_rc_doit  (p%rc, st, ii)
    if(p%calc_mask) call PES_mask_doit(p%mask, m, st, dt, mask)

  end subroutine PES_doit


  ! ---------------------------------------------------------
  subroutine PES_output(p, m, st, iter, save_iter, dt)
    type(PES_t),    intent(in) :: p
    type(mesh_t),   intent(in) :: m
    type(states_t), intent(in) :: st
    integer,        intent(in) :: iter, save_iter
    FLOAT,          intent(in) :: dt

    if(p%calc_rc)   call PES_rc_output   (p%rc, st, iter, save_iter, dt)
    if(p%calc_mask) call PES_mask_output (p%mask, m, st, "PES")

  end subroutine PES_output

#include "pes_rc.F90"
#include "pes_mask.F90"

#endif
end module PES_m

