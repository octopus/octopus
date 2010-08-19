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

#include "global.h"

module PES_m
  use datasets_m
  use fft_m
  use global_m
  use io_m
  use mesh_m
  use messages_m
  use parser_m
  use profiling_m
  use simul_box_m
  use states_m
  use unit_m
  use unit_system_m

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
  subroutine PES_init(pes, mesh, sb, st, ab, save_iter)
    type(pes_t),    intent(out)   :: pes
    type(mesh_t),   intent(inout) :: mesh
    type(simul_box_t), intent(in) :: sb
    type(states_t),    intent(in) :: st
    integer,           intent(in) :: ab, save_iter

    PUSH_SUB(PES_init)

    !%Variable CalcPES_rc
    !%Type logical
    !%Default no
    !%Section Time-Dependent::PES
    !%Description
    !%  If <tt>true</tt>, store the wavefunctions at specific points in order to 
    !% calculate the photoelectron spectrum at a point far in the box as proposed in 
    !% A. Pohl, P.-G. Reinhard, and E. Suraud, <i>Phys. Rev. Lett.</i> <b>84</b>, 5090 (2000).
    !%End
    call parse_logical(datasets_check('CalcPES_rc'), .false., pes%calc_rc)
    if(pes%calc_rc) then
      pes%calc_rc = .true.
      call PES_rc_init(pes%rc, mesh, st, save_iter)
    end if

    pes%calc_mask = .false.
    ! have the mask, and we are working in the velocity gauge
    if(ab == 2) then
      !%Variable CalcPES_mask
      !%Type logical
      !%Default no
      !%Section Time-Dependent::PES
      !%Description
      !% If <tt>true</tt>, calculate the photo-electron spectrum using the mask method
      !% (M. Marques, D. Varsano, H. Appel, E.K.U. Gross and A. Rubio, to be submitted). 
      !% For this to work, masking boundaries are necessary (<tt>AbsorbingBoundaries == 2</tt>).
      !%End
      call parse_logical(datasets_check('CalcPES_Mask'), .false., pes%calc_mask)
      if(pes%calc_mask) then
        call PES_mask_init(pes%mask, mesh, sb, st)
      end if
    end if

    POP_SUB(PES_init)
  end subroutine PES_init


  ! ---------------------------------------------------------
  subroutine PES_end(pes)
    type(PES_t), intent(inout) :: pes

    PUSH_SUB(PES_end)

    if(pes%calc_rc)   call PES_rc_end  (pes%rc)
    if(pes%calc_mask) call PES_mask_end(pes%mask)

    POP_SUB(PES_end)
  end subroutine PES_end


  ! ---------------------------------------------------------
  subroutine PES_calc(pes, mesh, st, ii, dt, mask)
    type(PES_t),    intent(inout) :: pes
    type(mesh_t),   intent(in)    :: mesh
    type(states_t), intent(in)    :: st
    FLOAT,          intent(in)    :: dt
    FLOAT,          pointer       :: mask(:)
    integer,        intent(in)    :: ii

    PUSH_SUB(PES_calc)

    if(pes%calc_rc)   call PES_rc_calc  (pes%rc, st, ii)
    if(pes%calc_mask) call PES_mask_calc(pes%mask, mesh, st, dt, mask)

    POP_SUB(PES_calc)
  end subroutine PES_calc


  ! ---------------------------------------------------------
  subroutine PES_output(pes, mesh, st, iter, save_iter, dt)
    type(PES_t),    intent(in) :: pes
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(in) :: st
    integer,        intent(in) :: iter, save_iter
    FLOAT,          intent(in) :: dt

    PUSH_SUB(PES_output)

    if(pes%calc_rc)   call PES_rc_output   (pes%rc, st, iter, save_iter, dt)
    if(pes%calc_mask) call PES_mask_output (pes%mask, mesh, st, "PES")

    POP_SUB(PES_output)
  end subroutine PES_output

#include "pes_rc_inc.F90"
#include "pes_mask_inc.F90"

end module PES_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
