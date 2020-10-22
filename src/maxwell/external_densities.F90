!! Copyright (C) 2020 N. Tancogne-Dejean
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

#include "global.h"

module external_densities_oct_m
  use global_oct_m
  use lalg_basic_oct_m
  use math_oct_m
  use states_mxll_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use string_oct_m
  use tdfunction_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m  
  use varinfo_oct_m
  use wfs_elec_oct_m
  use xc_oct_m
  
  implicit none

  private

  public ::                               &
    get_rs_density_ext,                   &
    external_current_init,                &
    external_current_calculation

  integer, parameter, public ::           &
    EXTERNAL_CURRENT_PARSER      = 0,     &
    EXTERNAL_CURRENT_TD_FUNCTION = 1

contains

  !----------------------------------------------------------
  subroutine get_rs_density_ext(st, mesh, time, rs_current_density_ext)
    type(states_mxll_t), intent(inout) :: st
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: time
    CMPLX,     optional, intent(inout) :: rs_current_density_ext(:,:)

    FLOAT, allocatable :: current(:,:,:)
    integer :: idim

    PUSH_SUB(get_rs_density_ext)

    SAFE_ALLOCATE(current(1:mesh%np, 1:mesh%sb%dim, 1))  !< The 1 in the last column is a dummy to use batch routines

    call external_current_calculation(st, mesh, time, current(:, :, 1))
    call build_rs_current_state(current(:, :, 1), mesh, rs_current_density_ext(:, :), st%ep(:), mesh%np)
    do idim = 1, mesh%sb%dim
      call lalg_scal(mesh%np, -M_ONE, rs_current_density_ext(:, idim))
    end do

    SAFE_DEALLOCATE_A(current)

    POP_SUB(get_rs_density_ext)
  end subroutine get_rs_density_ext


  !----------------------------------------------------------
  subroutine external_current_init(st, namespace, mesh)
    type(states_mxll_t), intent(inout) :: st
    type(mesh_t),        intent(in)    :: mesh
    type(namespace_t),   intent(in)    :: namespace

    type(block_t)        :: blk
    integer              :: ip, il, nlines, ncols, idir, ierr
    FLOAT                :: j_vector(MAX_DIM), dummy(MAX_DIM), xx(MAX_DIM), rr, omega
    character(len=1024)  :: tdf_expression, phase_expression

    type(profile_t), save :: prof

    PUSH_SUB(external_current_init)

    call profiling_in(prof, 'EXTERNAL_CURRENT_INIT')

    !%Variable UserDefinedMaxwellExternalCurrent
    !%Type block
    !%Section MaxwellStates
    !%Description
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedMaxwellExternalCurrent
    !% <br>&nbsp;&nbsp; current_parser      | "expression_x_dir1" | "expression_y_dir1" | "expression_z_dir1"
    !% <br>&nbsp;&nbsp; current_parser      | "expression_x_dir2" | "expression_y_dir2" | "expression_z_dir2"
    !% <br>&nbsp;&nbsp; current_td_function | "amplitude_j0_x"    | "amplitude_j0_y"    | "amplitude_j0_z"    | omega   | envelope_td_function_name | phase
    !% <br>%</tt>
    !%
    !% Description about UserDefinedMaxwellExternalCurrent follows
    !%
    !%Option current_parser 0
    !% description follows
    !%Option current_td_function 1
    !% description follows
    !%End

    if(parse_block(namespace, 'UserDefinedMaxwellExternalCurrent', blk) == 0) then

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)

      st%external_current_number = nlines
      SAFE_ALLOCATE(st%external_current_modus(nlines))
      SAFE_ALLOCATE(st%external_current_string(MAX_DIM, nlines))
      SAFE_ALLOCATE(st%external_current_amplitude(1:mesh%np, MAX_DIM, nlines))
      SAFE_ALLOCATE(st%external_current_td_function(nlines))
      SAFE_ALLOCATE(st%external_current_omega(nlines))
      SAFE_ALLOCATE(st%external_current_td_phase(nlines))

      ! read all lines
      do il = 1, nlines
        ! Check that number of columns is four, five, six or seven.
        ncols = parse_block_cols(blk, il - 1)
        if((ncols  /=  4) .and. (ncols /= 5) .and. (ncols /= 6) .and. (ncols /= 7)) then
          message(1) = 'Each line in the MaxwellExternalCurrent block must have'
          message(2) = 'four, five, six or or seven columns.'
          call messages_fatal(2, namespace=namespace)
        end if

        call parse_block_integer(blk, il - 1, 0, st%external_current_modus(il))

        if (st%external_current_modus(il) == EXTERNAL_CURRENT_PARSER) then
          ! parse formula string
          do idir = 1, st%dim
            call parse_block_string(blk, il - 1, idir, st%external_current_string(idir, il))
            call conv_to_C_string(st%external_current_string(idir, il))
          end do
        else if (st%external_current_modus(il) == EXTERNAL_CURRENT_TD_FUNCTION) then
          do ip = 1, mesh%np
            call mesh_r(mesh, ip, rr, coords = xx)
            do idir = 1, st%dim
              call parse_block_string(blk, il - 1, idir, st%external_current_string(idir, il))
              call conv_to_C_string(st%external_current_string(idir, il))
              call parse_expression(j_vector(idir), dummy(idir), st%dim, xx, rr, M_ZERO, &
                st%external_current_string(idir, il))
              j_vector(idir) = units_to_atomic(units_inp%energy/(units_inp%length**2), j_vector(idir))
            end do
            st%external_current_amplitude(ip, 1:st%dim, il) = j_vector(1:st%dim)
          end do
          call parse_block_float(blk, il-1, 4, omega, unit_one/units_inp%time)
          st%external_current_omega(il) = omega
          call parse_block_string(blk, il-1, 5, tdf_expression)
          call tdf_read(st%external_current_td_function(il), namespace, trim(tdf_expression), ierr)
          if(parse_block_cols(blk, il-1) > 6) then
            call parse_block_string(blk, il-1, 6, phase_expression)
            call tdf_read(st%external_current_td_phase(il), namespace, trim(phase_expression), ierr)
            if (ierr /= 0) then
              write(message(1),'(3A)') 'Error in the "', trim(tdf_expression), '" field defined in the TDExternalFields block:'
              write(message(2),'(3A)') 'Time-dependent phase function "', trim(phase_expression), '" not found.'
              call messages_warning(2, namespace=namespace)
            end if
          else
            call tdf_init(st%external_current_td_phase(il))
          end if
        end if
      end do
      call parse_block_end(blk)
    end if

    call profiling_out(prof)

    POP_SUB(external_current_init)
  end subroutine external_current_init

  !----------------------------------------------------------
  subroutine external_current_calculation(st, mesh, time, current)
    type(states_mxll_t), intent(inout) :: st
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: time
    FLOAT,               intent(inout) :: current(:,:)

    integer :: ip, jn, idir
    FLOAT   :: xx(MAX_DIM), rr, tt, j_vector(MAX_DIM), dummy(MAX_DIM), amp(MAX_DIM)
    CMPLX   :: exp_arg
    type(profile_t), save :: prof
    FLOAT   :: tmp_amp, phase

    PUSH_SUB(external_current_calculation)

    call profiling_in(prof, "EXTERNAL_CURRENT_CALC")

    current(:,:) = M_ZERO
    do jn = 1, st%external_current_number
      if (st%external_current_modus(jn) == EXTERNAL_CURRENT_PARSER) then
        do ip = 1, mesh%np
          call mesh_r(mesh, ip, rr, coords = xx)
          do idir = 1, st%dim
            tt = time
            call parse_expression(j_vector(idir), dummy(idir), st%dim, xx, rr, tt, &
              & trim(st%external_current_string(idir,jn)))
            j_vector(idir) = units_to_atomic(units_inp%energy/(units_inp%length**2), j_vector(idir))
          end do
          current(ip, 1:st%dim) = current(ip, 1:st%dim) + j_vector(1:st%dim)
        end do

      else if(st%external_current_modus(jn) == EXTERNAL_CURRENT_TD_FUNCTION) then
        exp_arg = st%external_current_omega(jn) * time + tdf(st%external_current_td_phase(jn),time)
        phase = TOFLOAT(exp(-M_zI*exp_arg))
        tmp_amp = tdf(st%external_current_td_function(jn), time)
        do ip = 1, mesh%np
          amp(1:st%dim) = st%external_current_amplitude(ip, 1:st%dim, jn) * tmp_amp
          current(ip, 1:st%dim) = current(ip, 1:st%dim) + amp(1:st%dim) * phase
        end do
      end if
    end do

    call profiling_out(prof)

    POP_SUB(external_current_calculation)
  end subroutine external_current_calculation
  
end module external_densities_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
