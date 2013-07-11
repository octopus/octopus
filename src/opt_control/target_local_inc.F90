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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id: target_local_inc.F90 $


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_local(gr, tg)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg

    integer             :: ip
    FLOAT               :: xx(MAX_DIM), rr, psi_re, psi_im
    character(len=1024) :: expression
    PUSH_SUB(target_init_local)

    !%Variable OCTLocalTarget
    !%Type string
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If <tt>OCTTargetOperator = oct_tg_local</tt>, then one must supply a function
    !% that defines the target. This should be done by defining it through a string, using 
    !% the variable <tt>OCTLocalTarget</tt>.
    !%End
    if(parse_isdef('OCTLocalTarget') /= 0) then
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))
      tg%rho = M_ZERO
      call parse_string(datasets_check('OCTLocalTarget'), "0", expression)
      call conv_to_C_string(expression)
      do ip = 1, gr%mesh%np
        call mesh_r(gr%mesh, ip, rr, coords = xx)
        ! parse user-defined expression
        call parse_expression(psi_re, psi_im, gr%sb%dim, xx, rr, M_ZERO, expression)
        tg%rho(ip) = psi_re
      end do
    else
      message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
      message(2) = 'of this target in variable "OCTLocalTarget".'
      call messages_fatal(2)
    end if

    POP_SUB(target_init_local)
  end subroutine target_init_local


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_local(tg)
    type(target_t),   intent(inout) :: tg
    PUSH_SUB(target_end_local)
    SAFE_DEALLOCATE_P(tg%rho)
    POP_SUB(target_end_local)
  end subroutine target_end_local


  ! ----------------------------------------------------------------------
  subroutine target_output_local(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    integer :: ierr
    PUSH_SUB(target_output_local)
    
    call loct_mkdir(trim(dir))
    if(outp%how /= 0) then
      call dio_function_output(outp%how, trim(dir), 'local_target', gr%mesh, &
        tg%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
    end if

    POP_SUB(target_output_local)
  end subroutine target_output_local
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_local(gr, tg, psi) result(j1)
    type(grid_t),     intent(inout) :: gr
    type(target_t),   intent(inout) :: tg
    type(states_t),   intent(inout) :: psi

    integer :: is
    PUSH_SUB(target_j1_local)

    j1 = M_ZERO
    do is = 1, psi%d%spin_channels
      j1 = j1 + dmf_dotp(gr%mesh, tg%rho, psi%rho(:, is))
    end do

    POP_SUB(target_j1_local)
  end function target_j1_local


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_local(tg, gr, psi_in, chi_out)
    type(target_t),    intent(inout) :: tg
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out

    integer :: ik, idim, ist, ip
    PUSH_SUB(target_chi_local)

    do ik = 1, psi_in%d%nik
      do idim = 1, psi_in%d%dim
        do ist = psi_in%st_start, psi_in%st_end
          do ip = 1, gr%mesh%np
            chi_out%zpsi(ip, idim, ist, ik) = psi_in%occ(ist, ik) * tg%rho(ip) * psi_in%zpsi(ip, idim, ist, ik)
          end do
        end do
      end do
    end do

    POP_SUB(target_chi_local)
  end subroutine target_chi_local



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
