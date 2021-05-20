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


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_spin(tg, namespace)
    type(target_t),    intent(inout) :: tg
    type(namespace_t), intent(in)    :: namespace
        

    type(block_t)       :: blk
    integer :: jst
    CMPLX :: alpha(3)

    PUSH_SUB(target_init_spin)


    message(1) =  'Info: Using a spin target'
    call messages_info(1)

    !%Variable OCTTargetSpin
    !%Type block
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% (Experimental) Specify the targeted spin as a 3-component vector. It will be normalized.
    !%End
    if(parse_is_defined(namespace, 'OCTTargetSpin')) then
      if(parse_block(namespace, 'OCTTargetSpin', blk) == 0) then
        alpha = M_z0
        do jst = 1, parse_block_cols(blk, 0)
          call parse_block_cmplx(blk, 0, jst - 1, alpha(jst))
        end do
        call parse_block_end(blk)
        
        alpha = alpha/sqrt(dot_product(alpha, alpha))

        tg%spin_matrix(1, 1) = alpha(3)
        tg%spin_matrix(2, 2) = -alpha(3)
        tg%spin_matrix(1, 2) = alpha(1) - M_zI * alpha(2)
        tg%spin_matrix(2, 1) = alpha(1) + M_zI * alpha(2)
        
      else
        call messages_variable_is_block(namespace, 'OCTTargetSpin')
      end if

    else
      message(1) = 'If "OCTTargetOperator = oct_tg_spin", then you must'
      message(2) = 'supply a "OCTTargetSpin" block.'
      call messages_fatal(2)
    end if


    POP_SUB(target_init_spin)
  end subroutine target_init_spin



  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_spin(tg, gr, psi) result(j1)
    type(target_t),      intent(in) :: tg
    type(grid_t),        intent(in) :: gr
    type(states_elec_t), intent(in) :: psi
    
    integer :: i, j
    CMPLX, allocatable :: zpsi(:, :)
    
    PUSH_SUB(target_j1_spin)

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:tg%st%d%dim))
    
    call states_elec_get_state(psi, gr%mesh, 1, 1, zpsi)
    
    j1 = M_ZERO
     do i = 1, 2
       do j = 1, 2
         j1 = j1 + real(tg%spin_matrix(i,j)*zmf_dotp(gr%mesh, zpsi(:, i), zpsi(:, j)))
      end do
    end do

    SAFE_DEALLOCATE_A(zpsi)
    
    POP_SUB(target_j1_spin)
  end function target_j1_spin



  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_spin(tg, gr, psi_in, chi_out)
    type(target_t),       intent(in)    :: tg
    type(grid_t),         intent(in)    :: gr
    type(states_elec_t),  intent(in)    :: psi_in
    type(states_elec_t),  intent(inout) :: chi_out
    
    integer :: i, j
    CMPLX, allocatable :: zpsi(:, :), zchi(:, :)
    
    PUSH_SUB(target_chi_spin)

    SAFE_ALLOCATE(zpsi(1:gr%mesh%np, 1:tg%st%d%dim))
    SAFE_ALLOCATE(zchi(1:gr%mesh%np, 1:tg%st%d%dim))

    call states_elec_get_state(psi_in, gr%mesh, 1, 1, zpsi)
    
    zchi(1:gr%mesh%np, 1:tg%st%d%dim) = CNST(0.0)

    do i = 1, 2
      do j = 1, 2
        zchi(1:gr%mesh%np, i) = zchi(1:gr%mesh%np, i) + tg%spin_matrix(i, j)*zpsi(1:gr%mesh%np, j)
      end do
    end do

    call states_elec_set_state(chi_out, gr%mesh, 1, 1, zchi)
    
    SAFE_DEALLOCATE_A(zpsi)
    SAFE_DEALLOCATE_A(zchi)
    
    POP_SUB(target_chi_spin)
  end subroutine target_chi_spin


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
