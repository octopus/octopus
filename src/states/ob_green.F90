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
!! $Id: h.F90 4037 2008-04-03 13:30:00Z xavier $

! Lead-related routines for open boundaries.

#include "global.h"

module ob_green_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use math_m
  use messages_m
  use ob_interface_m
  use profiling_m
  use simul_box_m

  implicit none
  private

  public ::         &
    lead_self_energy

contains

  ! calculate the lead self-energy
  ! prefer the fast method if possible, but if not converged check also other method
  subroutine lead_self_energy(energy, diag, offdiag, intf, self_energy, h_is_real)
    FLOAT,             intent(in)  :: energy        ! Energy to calculate Green`s function for.
    CMPLX,             intent(in)  :: diag(:, :)    ! Diagonal block of lead Hamiltonian.
    CMPLX,             intent(in)  :: offdiag(:, :) ! Off-diagonal block of lead Hamiltonian.
    type(interface_t), intent(in)  :: intf          ! => gr%intf(il)
    CMPLX,             intent(inout) :: self_energy(:, :) ! The calculated self_energy.
    logical,           intent(in)  :: h_is_real     ! Is the Hamiltonian real? (no vector potential)

    integer            :: np, np_uc
    FLOAT              :: threshold, eta, residual, residual2, eps
    CMPLX              :: c_energy
    CMPLX, allocatable :: green2(:, :)

    PUSH_SUB(lead_self_energy)

    np = intf%np_intf
    np_uc = intf%np_uc
    
    eta = CNST(1e-7) ! FIXME: read from input. relative eta
    threshold = CNST(1e-14) ! FIXME: read from input.
    eps = CNST(1e-5) ! FIXME: read from input.
    c_energy = energy + eta*M_zI

    ! if we cannot reduce the unit cell
    ! we have only one way of calculating the Green`s function: with the Sancho method
    if(.not.intf%reducible) then ! use large matrices (rank = np_uc)
      call lead_green_sancho(c_energy, diag, offdiag, np_uc, self_energy, threshold, h_is_real)
      residual = calc_residual_green(energy, self_energy, diag, offdiag, intf)
      if(in_debug_mode) then ! write info
        write(message(1), '(a)') 'Non-reducible unit cell'
        write(message(2), '(a,e10.3)') 'Sancho-Residual = ', residual
        call messages_info(2)
      end if
    else
      ! 1. calculate with the fastest method
      call lead_green_umerski(c_energy, diag, offdiag, intf, self_energy)
      ! 2. check if converged
      residual = calc_residual_green(energy, self_energy, diag, offdiag, intf)
      if(in_debug_mode) then ! write info
        write(message(1), '(a,e10.3)') 'Umerski-Residual = ', residual
        call messages_info(1)
      end if
      if(residual.gt.eps) then
        ! 3. if not converged try the other method while saving the old
        SAFE_ALLOCATE(green2(1:np, 1:np))
        call lead_green_sancho(c_energy, diag, offdiag, np, green2, threshold, h_is_real)
        residual2 = calc_residual_green(energy, green2, diag, offdiag, intf)
        if(in_debug_mode) then ! write info
          write(message(1), '(a,e10.3)') 'Sancho-Residual = ', residual2
          call messages_info(1)
        end if
        ! 4. if also not converged choose the best one and give warning
        if(residual2.lt.eps) then ! finally converged
          self_energy = green2
        else ! write warning
          if(residual2.lt.residual) then
            self_energy = green2
          end if
          message(1) = "The surface Green's function did not converge properly"
          message(2) = 'with either the decimation technique or the closed form!'
          write(message(3), '(a,e10.3)') 'Umerski-Residual = ', residual
          write(message(4), '(a,e10.3)') 'Sancho-Residual  = ', residual2
          message(5) = 'The better converged version is taken.'
          call messages_warning(5)
        end if
        ! cleanup
        SAFE_DEALLOCATE_A(green2)
      end if
    end if
    ! now calculate the self-energy
    if(intf%il.eq.LEFT) then
      call lalg_trmm(np, np, 'U', 'N', 'L', M_z1, offdiag, self_energy)
      call lalg_trmm(np, np, 'U', 'T', 'R', M_z1, offdiag, self_energy)
    else
      call lalg_trmm(np, np, 'L', 'N', 'L', M_z1, offdiag, self_energy)
      call lalg_trmm(np, np, 'L', 'T', 'R', M_z1, offdiag, self_energy)
    end if

    POP_SUB(lead_self_energy)
  end subroutine lead_self_energy


  ! check if the Green`s function gives the correct density of states
  ! if not compute its Hermitian conjugate
  subroutine fix_green(np, green)
    integer, intent(in)    :: np
    CMPLX,   intent(inout) :: green(:, :)

    integer  :: j
    FLOAT    :: dos

    PUSH_SUB(fix_green)

    ! calculate DOS=-Tr(imag(green)) and check if negative
    dos = -aimag(green(1, 1))
    do j = 2, np
      dos = dos - aimag(green(j, j))
    end do

    if(in_debug_mode) then ! write info
      write(message(1), '(a,e10.3)') 'density of states = ', abs(dos)
      call messages_info(1)
    end if

    if(dos.lt.M_ZERO) then
      green(:, :) = transpose(conjg(green(:, :)))
      if(in_debug_mode) then ! write info
        message(1) = "surface Green's function changed to its hermitian conjugate"
        call messages_info(1)
      end if
    end if

    POP_SUB(fix_green)
  end subroutine fix_green


  ! calculate the residual of the Green`s function
  FLOAT function calc_residual_green(energy, green, diag, offdiag, intf) result(residual)
    FLOAT,   intent(in) :: energy
    CMPLX,   intent(in) :: green(:, :)
    CMPLX,   intent(in) :: diag(:, :)
    CMPLX,   intent(in) :: offdiag(:, :)
    type(interface_t), intent(in)  :: intf ! => gr%intf(il)

    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :)
    FLOAT              :: det
    integer            :: j, np, np_uc

    PUSH_SUB(calc_residual_green)

    np = intf%np_intf
    np_uc = intf%np_uc
    
    SAFE_ALLOCATE(tmp1(1:np, 1:np))
    SAFE_ALLOCATE(tmp2(1:np, 1:np))

    ! 1. calculate the residual InfNorm(Inverse(energy-h-offdiag*green*offdiag^T)-green)
    tmp1(1:np, 1:np) = -diag(1:np, 1:np)
    forall (j = 1:np) tmp1(j, j) = tmp1(j, j) + energy

    tmp2(1:np, 1:np) = green(1:np, 1:np)
    if(mod(intf%il+1,2)+1.eq.1) then
      call lalg_trmm(np, np, 'U', 'N', 'L', M_z1, offdiag, tmp2)
      call lalg_trmm(np, np, 'U', 'T', 'R', M_z1, offdiag, tmp2)
    else
      call lalg_trmm(np, np, 'L', 'N', 'L', M_z1, offdiag, tmp2)
      call lalg_trmm(np, np, 'L', 'T', 'R', M_z1, offdiag, tmp2)
    end if

    tmp1(1:np, 1:np) = tmp1(1:np, 1:np) - tmp2(1:np, 1:np)
    det = lalg_inverter(np, tmp1, invert = .true.)
!    tmp1(1:np, 1:np) = tmp1(1:np, 1:np) - green(1:np, 1:np)
    residual = abs(M_ONE - infinity_norm(tmp1)/infinity_norm(green))

    SAFE_DEALLOCATE_A(tmp1)
    SAFE_DEALLOCATE_A(tmp2)

    POP_SUB(calc_residual_green)
  end function calc_residual_green


  ! ---------------------------------------------------------
  ! Calculate head of the semi-infinite surface Green`s function with the
  ! algorithm from the paper
  ! Highly convergent schemes for the calculation of bulk and surface Green`s functions
  ! M. P. Lopez Sancho, J. M. Lopez Sancho, and J. Rubio (1984)
  ! J. Phys. F: Met. Phys. 15 (1985) 851-858
  subroutine lead_green_sancho(energy, diag, offdiag, np, green, threshold, h_is_real)
    CMPLX,   intent(in)  :: energy        ! Energy to calculate Green`s function for (already shifted).
    CMPLX,   intent(in)  :: diag(:, :)    ! Diagonal block of lead Hamiltonian.
    CMPLX,   intent(in)  :: offdiag(:, :) ! Off-diagonal block of lead Hamiltonian.
    integer, intent(in)  :: np            ! Number of interface points.
    CMPLX,   intent(out) :: green(:, :)   ! The calculated Green`s function.
    FLOAT,   intent(in)  :: threshold     ! this defines convergence
    logical, intent(in)  :: h_is_real     ! Is the Hamiltonian real? (no vector potential)

    CMPLX, allocatable :: e(:, :), es(:, :), a(:, :), b(:, :), inv(:, :)
    CMPLX, allocatable :: tmp1(:, :), tmp2(:, :), tmp3(:, :)
    integer            :: i, j
    FLOAT              :: det, old_norm, norm, res

    PUSH_SUB(lead_green_sancho)

    SAFE_ALLOCATE(   e(1:np, 1:np))
    SAFE_ALLOCATE(  es(1:np, 1:np))
    SAFE_ALLOCATE(   a(1:np, 1:np))
    SAFE_ALLOCATE(   b(1:np, 1:np))
    SAFE_ALLOCATE( inv(1:np, 1:np))
    SAFE_ALLOCATE(tmp1(1:np, 1:np))
    SAFE_ALLOCATE(tmp2(1:np, 1:np))
    SAFE_ALLOCATE(tmp3(1:np, 1:np))

    tmp1 = M_ZERO
    tmp2 = M_ZERO
    tmp3 = M_ZERO

    ! Fill with start values.
    e(1:np, 1:np) = diag(1:np, 1:np)
    a(1:np, 1:np) = offdiag(1:np, 1:np)
    b(1:np, 1:np) = transpose(conjg(offdiag(1:np, 1:np)))
    
    es(1:np, 1:np) = diag(1:np, 1:np)
    old_norm = M_ZERO

    do i = 1, 1000 ! FIXME: read from input. 2^1000 efective layers
      inv(1:np, 1:np) = - e(1:np, 1:np)

      forall (j = 1:np) inv(j, j) = inv(j, j) + energy
      det = lalg_inverter(np, inv, invert=.true.)

      call lalg_gemm(np, np, np, M_z1, a, inv, M_z0, tmp2)
      call lalg_gemm(np, np, np, M_z1, tmp2, b, M_z0, tmp1)

      es(1:np, 1:np) = es(1:np, 1:np) + tmp1(1:np, 1:np)

      norm = infinity_norm(es)
      res  = abs(old_norm/norm-M_ONE)

      if(res.lt.threshold) exit

      old_norm = norm

      e(1:np, 1:np) = e(1:np, 1:np) + tmp1(1:np, 1:np)
      call lalg_gemm(np, np, np, M_z1, b, inv, M_Z0, tmp1)
      call lalg_gemm(np, np, np, M_z1, tmp1, a, M_z1, e)
      tmp3(1:np, 1:np) = a(1:np, 1:np)
      call lalg_gemm(np, np, np, M_z1, tmp2, tmp3, M_z0, a)
      tmp3(1:np, 1:np) = b(1:np, 1:np)
      call lalg_gemm(np, np, np, M_z1, tmp1, tmp3, m_z0, b)
    end do

    green(1:np, 1:np) = - es(1:np, 1:np)

    forall (j = 1:np) green(j, j) = green(j, j) + energy
    det = lalg_inverter(np, green, invert = .true.)
    if (h_is_real) then ! the Green`s function is complex symmetric
      call matrix_symmetric_average(green, np)
    end if ! otherwise it is general complex

    if(in_debug_mode) then ! write some info
      message(1) = "surface Green's function: iterative algorithm"
      call messages_info(1)
    end if

    call fix_green(np, green)

    SAFE_DEALLOCATE_A(e)
    SAFE_DEALLOCATE_A(es)
    SAFE_DEALLOCATE_A(a)
    SAFE_DEALLOCATE_A(b)
    SAFE_DEALLOCATE_A(inv)
    SAFE_DEALLOCATE_A(tmp1)
    SAFE_DEALLOCATE_A(tmp2)
    SAFE_DEALLOCATE_A(tmp3)

    POP_SUB(lead_green_sancho)
  end subroutine lead_green_sancho


  ! create the Moebius transformation matrix which is to be diagonalized
  ! (or used to create the matrices for the final transformation matrix)
  subroutine create_moeb_trans_matrix(intf, energy, diag, offdiag, matrix)
    type(interface_t), intent(in)  :: intf          ! => gr%intf(il)
    CMPLX,             intent(in)  :: energy
    CMPLX,             intent(in)  :: diag(1:intf%np_intf, 1:intf%np_intf)
    CMPLX,             intent(in)  :: offdiag(1:intf%np_intf, 1:intf%np_intf)
    CMPLX,             intent(out) :: matrix(1:2*intf%np_intf, 1:2*intf%np_intf)

    integer  :: ip, np, np2
    CMPLX, allocatable   :: x1(:,:), x2(:,:)

    PUSH_SUB(create_moeb_trans_matrix)

    np  = intf%np_intf
    np2 = 2*np
    SAFE_ALLOCATE( x1(1:np, 1:np) )
    SAFE_ALLOCATE( x2(1:np, 1:np) )

    ! 1. create matrix x ( x = {{0,offdiag^(-1)},{-offdiag^H,(energy-diag)*offdiag^(-1)}} )
    forall(ip = 1:np) matrix(1:np, ip) = M_z0

    ! use o2 as tmp variable
    forall(ip = 1:np) x1(1:np, ip) = offdiag(1:np, ip)

    if (mod(intf%il+1,2)+1.eq.1) then
      call lalg_invert_upper_triangular(np, x1)
    else
      call lalg_invert_lower_triangular(np, x1)
    end if

    forall(ip = 1:np) matrix(ip, np+1:np2) = x1(ip, 1:np)
    forall(ip = 1:np) matrix(np+1:np2, ip) = -conjg(offdiag(ip, 1:np))

    ! use o4 as tmp variable
    forall(ip = 1:np) x2(1:np, ip) = -diag(1:np, ip)
    forall (ip = 1:np) x2(ip, ip) = x2(ip, ip) + energy

    if (mod(intf%il+1,2)+1.eq.1) then
      call lalg_trmm(np, np, 'U', 'N', 'R', M_z1, x1, x2)
    else
      call lalg_trmm(np, np, 'L', 'N', 'R', M_z1, x1, x2)
    end if
    
    forall (ip = 1:np) matrix(np+1:np2, np+ip) = x2(1:np, ip)

    SAFE_DEALLOCATE_A(x1)
    SAFE_DEALLOCATE_A(x2)

    POP_SUB(create_moeb_trans_matrix)
  end subroutine create_moeb_trans_matrix

  
  ! extract from the periodic Hamiltonian the submatrices of the size
  ! of the interface
  subroutine extract_sub_matrices(intf, diag, offdiag, h, v)
    type(interface_t), intent(in)  :: intf          ! => gr%intf(il)
    CMPLX,   intent(in)    :: diag(1:intf%np_uc, 1:intf%np_uc)   ! the lead diagonal (of the unit cell)
    CMPLX,   intent(in)    :: offdiag(1:intf%np_uc, 1:intf%np_uc)! the lead off-diagonal (of the unit cell)
    CMPLX,   intent(out)   :: h(:, :, :) ! the diagonal sub blocks
    CMPLX,   intent(out)   :: v(:, :, :) ! the off-diagonal sub blocks

    integer :: ni, iblock, np

    PUSH_SUB(extract_sub_matrices)

    np = intf%np_intf
    ! extract sub diagonals
    do iblock = 1, intf%nblocks-1
      ni = (iblock-1)*np
      h(1:np, 1:np, iblock) = diag(ni+1:ni+np, ni+1:ni+np)
      if(mod(intf%il+1,2)+1.eq.1) then
        ! lower blocks
        v(1:np, 1:np, iblock) = diag(ni+1:ni+np, ni+np+1:ni+2*np)
      else
        ! upper blocks
        v(1:np, 1:np, iblock) = diag(ni+np+1:ni+2*np, ni+1:ni+np)
      end if
    end do
    ! now the last block
    ni = (intf%nblocks-1)*np
    h(1:np, 1:np, intf%nblocks) = diag(ni+1:ni+np, ni+1:ni+np)
    if(mod(intf%il+1,2)+1.eq.1) then
      ! upper block of offdiag
      v(1:np, 1:np, intf%nblocks) = offdiag(ni+1:ni+np, 1:np)
    else
      ! lower block of offdiag
      v(1:np, 1:np, intf%nblocks) = offdiag(1:np, ni+1:ni+np)
    end if

    POP_SUB(extract_sub_matrices)
  end subroutine extract_sub_matrices


  ! compute the semi-infinite surface Green`s function
  ! Algorithm taken from A. Umerski, Closed-form solutions to surface Green`s functions
  ! http://www.city.ac.uk/sems/dps/mathematics/research/nanostructures/prb55_5266.pdf
  subroutine lead_green_umerski(energy, diag, offdiag, intf, green)
    CMPLX,     intent(in)  :: energy
    CMPLX,     intent(in)  :: diag(:, :)
    CMPLX,     intent(in)  :: offdiag(:, :)
    type(interface_t), intent(in)  :: intf ! => gr%intf(il)
    CMPLX,     intent(out) :: green(:, :)

    integer              :: ib, np, np2
    CMPLX, allocatable   :: x(:,:), x2(:,:), o2(:,:), o4(:,:), d(:), h(:, :, :), v(:, :, :)
    FLOAT                :: det

    PUSH_SUB(lead_green_umerski)

    ! check if intf%np_uc is a integer multiple of the intf%np_intf
    ASSERT(mod(intf%np_uc, intf%np_intf).eq.0)

    np = intf%np_intf
    np2 = 2*np
    SAFE_ALLOCATE( x(1:np2, 1:np2) )
    SAFE_ALLOCATE( x2(1:np2, 1:np2) )
    SAFE_ALLOCATE( o2(1:np, 1:np) )
    SAFE_ALLOCATE( o4(1:np, 1:np) )
    SAFE_ALLOCATE( d(1:np2) )
    SAFE_ALLOCATE( h(1:np, 1:np, 1:intf%nblocks) )
    SAFE_ALLOCATE( v(1:np, 1:np, 1:intf%nblocks) )

    ! 1. Reduce any super-structure to the size of the normal surface Green`s function.
    ! 1.(a) Extract all block-diagonal and hopping matrices with the size of the
    ! interface matrix [np x np], that makes intf%nblocks of these matrices for both h and v
    call extract_sub_matrices(intf, diag, offdiag, h, v)
    ! 1.(b) Create the Moebius transformation matrices for each sub-system and multiply all
    ! x = x(m)*...*x(2)*x(1)
    do ib=1, intf%nblocks
      call create_moeb_trans_matrix(intf, energy, h(:, :, ib), v(:, :, ib), x2)
      if (ib.eq.1) then
        x = x2
      else
        x = matmul(x2, x)
      endif
    end do
    
    ! 2. compute diagonalization matrix s, s^(-1)*x*s = d
    call lalg_eigensolve_nonh(np2, x, d)
    ! the eigenvalues and the corresponding eigenvectors have to be sorted in the following way:
    ! let d = diag(d(1),d(2),d(3),...,d(2n))
    ! with |d(1)| < |d(2)| < ... < |d(2n)|
    call matrix_sort(np2, x, d)
    ! 3. extract submatrices ( o2 and o4 of S = {{o1,o2},{o3,o4}};)
    o2(1:np, 1:np) = x(1:np, np+1:np2)
    o4(1:np, 1:np) = x(np+1:np2, np+1:np2)
    ! 4. calculate green = o2*o4^(-1)
    det = lalg_inverter(np, o4, invert = .true.)
    call lalg_gemm(np, np, np, M_z1, o2, o4, M_z0, green)
    
    if(in_debug_mode) then ! write some info
      message(1) = "surface Green's function: direct algorithm (Moebius transformation)"
      call messages_info(1)
    end if

    call fix_green(np, green)

    SAFE_DEALLOCATE_A(x)
    SAFE_DEALLOCATE_A(x2)
    SAFE_DEALLOCATE_A(o2)
    SAFE_DEALLOCATE_A(o4)
    SAFE_DEALLOCATE_A(d)
    SAFE_DEALLOCATE_A(h)
    SAFE_DEALLOCATE_A(v)

    POP_SUB(lead_green_umerski)
  end subroutine lead_green_umerski

end module ob_green_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
