!! Copyright (C) 2006-2011 M. Marques, U. De Giovannini
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
!! $Id$

#include "global.h"

module pes_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use mesh_m
  use messages_m
  use mpi_m
  use output_m
  use parser_m
  use pes_mask_m
  use pes_rc_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use varinfo_m
    
  implicit none

  private

  public ::                             &
    pes_t,                              &
    pes_end,                            &
    pes_init,                           &
    pes_init_write,                     &
    pes_calc,                           &
    pes_output,                         &
    pes_load,                           &
    pes_dump

  type pes_t
    logical :: calc_rc
    type(pes_rc_t) :: rc

    logical :: calc_mask
    type(pes_mask_t) :: mask
    
  end type pes_t

  integer, parameter ::     &
    PHOTOELECTRON_NONE = 0, &
    PHOTOELECTRON_RC   = 2, &
    PHOTOELECTRON_MASK = 4

contains

  ! ---------------------------------------------------------
  !elemental (PUSH/POP)_SUB are not PURE.
  subroutine pes_rc_nullify(this)
    type(pes_rc_t), intent(out) :: this
    !
    PUSH_SUB(pes_rc_nullify)
    !this%npoints  = 0
    this%points    =>null()
    this%filenames =>null()
    this%wf        =>null()
    this%rankmin   =>null()
    POP_SUB(pes_rc_nullify)
    return
  end subroutine pes_rc_nullify

  ! ---------------------------------------------------------
  !elemental (PUSH/POP)_SUB are not PURE.
  subroutine pes_mask_nullify(this)
    type(pes_mask_t), intent(out) :: this
    !
    PUSH_SUB(pes_mask_nullify)
    this%k=>null()
    !this%ll=0
    !this%np=0
    !this%spacing=M_ZERO
    this%mesh=>null()
    !call cube_nullify(this%cube)
    this%ext_pot=>null()
    this%Mk=>null()
    !call cube_function_nullify(this%cM)
    this%mask_R=>null()
    !this%shape=0
    this%Lk=>null()
    !this%enlarge
    !thisenlarge_nfft
    !this%llr
    !this%energyMax 
    !this%energyStep 
    !this%sw_evolve
    !this%back_action
    !this%add_psia
    !this%interpolate_out
    !this%filter_k
    !this%mode
    !this%pw_map_how
    !call fft_nullify(this%fft)
#if defined(HAVE_NFFT) 
    !call nfft_nullify(this%nfft)
#endif
    !call tdpsf_nullify(this%psf)
    POP_SUB(pes_mask_nullify)
    return
  end subroutine pes_mask_nullify

  ! ---------------------------------------------------------
  !elemental (PUSH/POP)_SUB are not PURE.
  subroutine pes_nullify(this)
    type(pes_t), intent(out) :: this
    !
    PUSH_SUB(pes_nullify)
    !this%calc_rc=.false.
    call pes_rc_nullify(this%rc)
    !this%calc_mask=.false.
    call pes_mask_nullify(this%mask)
    POP_SUB(pes_nullify)
    return
  end subroutine pes_nullify

  ! ---------------------------------------------------------
  subroutine pes_init(pes, mesh, sb, st, save_iter, hm, max_iter, dt)
    type(pes_t),         intent(out)   :: pes
    type(mesh_t),        intent(inout) :: mesh
    type(simul_box_t),   intent(in)    :: sb
    type(states_t),      intent(in)    :: st
    integer,             intent(in)    :: save_iter
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: max_iter
    FLOAT,               intent(in)    :: dt

    character(len=50)    :: str
    integer :: photoelectron_flags

    PUSH_SUB(pes_init)
    
    call pes_nullify(pes)

    call messages_obsolete_variable('CalcPES_rc', 'PhotoElectronSpectrum')
    call messages_obsolete_variable('CalcPES_mask', 'PhotoElectronSpectrum')

    !%Variable PhotoElectronSpectrum
    !%Type flag
    !%Default none
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
    !% Calculate the photo-electron spectrum using the mask method.
    !% U. De Giovannini, D. Varsano, M. A. L. Marques, H. Appel, E. K. U. Gross, and A. Rubio,
    !% <i>Phys. Rev. A</i> <b>85</b>, 062515 (2012).
    !%End

    call parse_integer(datasets_check('PhotoElectronSpectrum'), PHOTOELECTRON_NONE, photoelectron_flags)
    if(.not.varinfo_valid_option('PhotoElectronSpectrum', photoelectron_flags, is_flag = .true.)) then
      call input_error('PhotoElectronSpectrum')
    end if
    
    pes%calc_rc = iand(photoelectron_flags, PHOTOELECTRON_RC) /= 0
    pes%calc_mask = iand(photoelectron_flags, PHOTOELECTRON_MASK) /= 0

    !Header Photoelectron info
    if(pes%calc_rc .or. pes%calc_mask) then 
      write(str, '(a,i5)') 'Photoelectron'
      call messages_print_stress(stdout, trim(str))
    end if 

    
    if(pes%calc_rc) call pes_rc_init(pes%rc, mesh, st, save_iter)
    if(pes%calc_mask) call pes_mask_init(pes%mask, mesh, sb, st,hm,max_iter,dt)


    !Footer Photoelectron info
    if(pes%calc_rc .or. pes%calc_mask) then 
      call messages_print_stress(stdout)
    end if 

    POP_SUB(pes_init)
  end subroutine pes_init


  ! ---------------------------------------------------------
  subroutine pes_end(pes)
    type(pes_t), intent(inout) :: pes

    PUSH_SUB(pes_end)

    if(pes%calc_rc)   call pes_rc_end  (pes%rc)
    if(pes%calc_mask) call pes_mask_end(pes%mask)

    POP_SUB(pes_end)
  end subroutine pes_end


  ! ---------------------------------------------------------
  subroutine pes_calc(pes, mesh, st, ii, dt, iter)
    type(pes_t),         intent(inout) :: pes
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st
    integer,             intent(in)    :: ii
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter

    PUSH_SUB(pes_calc)

    if(pes%calc_rc)   call pes_rc_calc  (pes%rc, st, mesh, ii)
    if(pes%calc_mask) call pes_mask_calc(pes%mask, mesh, st, dt, iter)

    POP_SUB(pes_calc)
  end subroutine pes_calc


  ! ---------------------------------------------------------
  subroutine pes_output(pes, mesh, st, iter, outp, dt, gr, geo)
    type(pes_t),      intent(inout) :: pes
    type(mesh_t),     intent(in)    :: mesh
    type(states_t),   intent(in)    :: st
    integer,          intent(in)    :: iter
    type(output_t),   intent(in)    :: outp
    FLOAT,            intent(in)    :: dt
    type(grid_t),     intent(inout) :: gr
    type(geometry_t), intent(in)    :: geo

    PUSH_SUB(pes_output)
    
    if(mpi_grp_is_root(mpi_world)) then
      if(pes%calc_rc) call pes_rc_output(pes%rc, st, iter, outp%output_interval, dt)
    endif

    if(pes%calc_mask) call pes_mask_output (pes%mask, mesh, st,outp, "td.general/PESM", gr, geo,iter)

    POP_SUB(pes_output)
  end subroutine pes_output


  ! ---------------------------------------------------------
  subroutine pes_dump(restart, pes, st)
    type(restart_t), intent(in) :: restart
    type(pes_t),     intent(in) :: pes
    type(states_t),  intent(in) :: st

    PUSH_SUB(pes_dump)

    if (.not. restart_skip(restart) .and. pes%calc_mask) then
      message(1) = "Info: Writing PES mask."
      call messages_info(1)

      call pes_mask_dump(restart, pes%mask, st)

      message(1) = "Info: Finished writing PES mask."
      call messages_info(1)
    end if

    POP_SUB(pes_dump)
  end subroutine pes_dump


  ! ---------------------------------------------------------
  subroutine pes_load(restart, pes, st)
    type(restart_t), intent(in)    :: restart
    type(pes_t),     intent(inout) :: pes
    type(states_t),  intent(inout) :: st

    PUSH_SUB(pes_load)

    if (.not. restart_skip(restart) .and. pes%calc_mask) then
      message(1) = "Info: Reading PES mask."
      call messages_info(1)

      call pes_mask_load(restart, pes%mask, st)

      message(1) = "Info: Finished reading PES mask."
      call messages_info(1)
    end if

    POP_SUB(pes_load)
  end subroutine pes_load


  ! ---------------------------------------------------------
  subroutine pes_init_write(pes, mesh, st)
    type(pes_t),    intent(in)  :: pes
    type(mesh_t),   intent(in)  :: mesh
    type(states_t), intent(in)  :: st


    PUSH_SUB(pes_init_write)

    if(mpi_grp_is_root(mpi_world)) then

      if(pes%calc_rc)   call pes_rc_init_write (pes%rc, mesh, st)

    endif

    POP_SUB(pes_init_write)
  end subroutine pes_init_write



end module pes_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
