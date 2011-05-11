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
  use io_binary_m
  use mesh_m
  use index_m
  use messages_m
  use parser_m
  use profiling_m
  use simul_box_m
  use states_m
  use unit_m
  use unit_system_m
  use mpi_m
  use hamiltonian_m
  use geometry_m
  use lasers_m
  use varinfo_m

  implicit none

  type PES_rc_t
    integer          :: npoints   ! how many points we store the wf
    integer, pointer :: points(:) ! which points to use
    character(len=30), pointer :: filenames(:) ! filenames
    complex(r8), pointer :: wf(:,:,:,:,:)
    integer, pointer ::rankmin(:)  !partion of the mesh containing the points
  end type PES_rc_t

  type PES_mask_t
    CMPLX, pointer :: k(:,:,:,:,:,:) ! masked wf in momentum space
    FLOAT, pointer :: r(:,:,:,:,:)   ! summed masked density in real space
    
    FLOAT, pointer :: vec_pot(:,:)   ! vector potential from the laser
    type(fft_t) :: fft
    
    FLOAT :: energyMax 
    FLOAT :: energyStep 

  end type PES_mask_t

  type PES_t
    logical :: calc_rc
    type(PES_rc_t)   :: rc

    logical :: calc_mask
    type(PES_mask_t) :: mask
    
  end type PES_t

  integer, parameter ::     &
    PHOTOELECTRON_NONE = 0, &
    PHOTOELECTRON_RC   = 2, &
    PHOTOELECTRON_MASK = 4

contains

  ! ---------------------------------------------------------
  subroutine PES_init(pes, mesh, sb, st, ab, save_iter,hm, max_iter,dt)
    type(pes_t),         intent(out)     :: pes
    type(mesh_t),        intent(inout)   :: mesh
    type(simul_box_t),   intent(in)      :: sb
    type(states_t),      intent(in)      :: st
    integer,             intent(in)      :: ab, save_iter
    type(hamiltonian_t), intent(in)      :: hm
    integer,             intent(in)      :: max_iter
    FLOAT,               intent(in)      :: dt

    integer :: photoelectron_flags

    PUSH_SUB(PES_init)

    call messages_obsolete_variable('CalcPES_rc', 'PhotoElectronSpectrum')
    call messages_obsolete_variable('CalcPES_mask', 'PhotoElectronSpectrum')

    !%Variable PhotoElectronSpectrum
    !%Type flag
    !%Default no
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%This variable controls the method used for the calculation of
    !%the photoelectron spectrum. You can specify more than one value
    !%by giving them as a sum, for example:
    !% <tt>PhotoElectronSpectrum = pes_rc + pes_mask</tt>
    !%Option none 0
    !% The photoelectron spectrum is not calculated. This is the default.
    !%Option pes_rc 2
    !% Store the wavefunctions at specific points in order to 
    !% calculate the photoelectron spectrum at a point far in the box as proposed in 
    !% A. Pohl, P.-G. Reinhard, and E. Suraud, <i>Phys. Rev. Lett.</i> <b>84</b>, 5090 (2000).
    !%Option pes_mask 4
    !% Calculate the photo-electron spectrum using the mask method
    !% (M. Marques, D. Varsano, H. Appel, E.K.U. Gross and A. Rubio, to be submitted). 
    !% For this to work, masking boundaries are necessary (<tt>AbsorbingBoundaries == 2</tt>).
    !%End
    call parse_integer(datasets_check('PhotoElectronSpectrum'), PHOTOELECTRON_NONE, photoelectron_flags)
    if(.not.varinfo_valid_option('PhotoElectronSpectrum', photoelectron_flags, is_flag=.true.)) then
      call input_error('PhotoElectronSpectrum')
    end if
    
    pes%calc_rc = iand(photoelectron_flags, PHOTOELECTRON_RC) /= 0
    pes%calc_mask = iand(photoelectron_flags, PHOTOELECTRON_MASK) /= 0

    if(pes%calc_mask .and. ab /= MASK_ABSORBING) then
      message(1) = 'PhotoElectronSpectrum = pes_mask requires AbsorbingBoundaries = mask'
      call messages_fatal(1)
    end if
    
    if(pes%calc_rc) call PES_rc_init(pes%rc, mesh, st, save_iter)
    if(pes%calc_mask) call PES_mask_init(pes%mask, mesh, sb, st,hm,max_iter,dt)

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
  subroutine PES_calc(pes, mesh, st, ii, dt, mask,hm,geo,iter)
    type(PES_t),    intent(inout) :: pes
    type(mesh_t),   intent(in)    :: mesh
    type(states_t), intent(in)    :: st
    FLOAT,          intent(in)    :: dt
    FLOAT,          pointer       :: mask(:)
    integer,        intent(in)    :: ii
    integer,        intent(in)    :: iter
    type(hamiltonian_t), intent(in)    :: hm
    type(geometry_t), intent(in)    :: geo

    PUSH_SUB(PES_calc)

    if(pes%calc_rc)   call PES_rc_calc  (pes%rc, st,mesh, ii)
    if(pes%calc_mask) call PES_mask_calc(pes%mask, mesh, st, dt, mask&
         ,hm,geo,iter)

    POP_SUB(PES_calc)
  end subroutine PES_calc


  ! ---------------------------------------------------------
  subroutine PES_output(pes, mesh, st, iter, save_iter, dt)
    type(PES_t),    intent(inout) :: pes
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(in) :: st
    integer,        intent(in) :: iter, save_iter
    FLOAT,          intent(in) :: dt

    PUSH_SUB(PES_output)
    
    if(pes%calc_mask .AND. st%parallel_in_states) call PES_mask_collect(pes%mask, st,mesh)

    if(mpi_grp_is_root(mpi_world)) then

      if(pes%calc_rc)   call PES_rc_output   (pes%rc, st, iter, save_iter, dt)
      if(pes%calc_mask) call PES_mask_output (pes%mask, mesh, st, "td.general/PES")

    endif

    POP_SUB(PES_output)
  end subroutine PES_output

  ! ---------------------------------------------------------
  subroutine PES_restart_write(pes, mesh, st)
    type(PES_t),    intent(in) :: pes
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(in) :: st

    PUSH_SUB(PES_restart_write)

    if(mpi_grp_is_root(mpi_world)) then

      if(pes%calc_mask) call PES_mask_restart_write (pes%mask, mesh, st)

    endif

    POP_SUB(PES_restart_write)
  end subroutine PES_restart_write

  ! ---------------------------------------------------------
  subroutine PES_restart_read(pes, mesh, st)
    type(PES_t),    intent(inout) :: pes
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(in) :: st

    PUSH_SUB(PES_restart_read)

      if(pes%calc_mask) call PES_mask_restart_read (pes%mask, mesh, st)


    POP_SUB(PES_restart_read)
  end subroutine PES_restart_read


  ! ---------------------------------------------------------
  subroutine PES_init_write(pes, mesh, st)
    type(PES_t),    intent(in)  :: pes
    type(mesh_t),   intent(in)  :: mesh
    type(states_t), intent(in)  :: st


    PUSH_SUB(PES_init_write)

    if(mpi_grp_is_root(mpi_world)) then

      if(pes%calc_rc)   call PES_rc_init_write (pes%rc, mesh, st)

    endif

    POP_SUB(PES_init_write)
  end subroutine PES_init_write


#include "pes_rc_inc.F90"
#include "pes_mask_inc.F90"

end module PES_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
