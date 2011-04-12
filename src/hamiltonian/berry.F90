!! Copyright (C) 2010 D. Strubbe
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

module berry_m
  use global_m
  use grid_m
  use lalg_adv_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use profiling_m
  use simul_box_m
  use smear_m
  use states_m
  use states_dim_m

  implicit none

  private
  public ::                 &
    berry_dipole,           &
    berry_phase_det,        &
    berry_potential,        &
    berry_energy_correction

contains

  ! ---------------------------------------------------------
  ! Uses the single-point Berry`s phase method to calculate dipole moment in a periodic system
  ! This is only accurate in the limit of a large supercell.
  ! It is implemented only for an orthogonal unit cell.
  ! mu = - eL/2*pi Im ln <Psi|exp(-i(2*pi/L)x)|Psi>
  ! E Yaschenko, L Fu, L Resca, R Resta, Phys. Rev. B 58, 1222-1229 (1998)
  ! Single-point Berry`s phase method for dipole should not be used when there is more than one k-point.
  ! in this case, finite differences should be used to construct derivatives with respect to k
  FLOAT function berry_dipole(st, mesh, dir) result(dipole)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: mesh
    integer,        intent(in) :: dir

    integer :: ik
    CMPLX :: det

    PUSH_SUB(berry_dipole)

    if(.not. smear_is_semiconducting(st%smear)) then
      message(1) = "Warning: single-point Berry's phase dipole calculation not correct without integer occupations."
      call messages_warning(1)
    endif

    dipole = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end ! determinants for different spins and k-points multiply since matrix is block-diagonal
      det  = berry_phase_det(st, mesh, dir, ik)
      dipole = dipole + aimag(log(det))
    enddo

    dipole = -(mesh%sb%lsize(dir) / M_PI) * dipole

    POP_SUB(berry_dipole)
  end function berry_dipole


  ! ---------------------------------------------------------
  ! Uses the single-point Berry`s phase method to calculate dipole moment in a periodic system
  ! This is only accurate in the limit of a large supercell.
  ! It is implemented only for an orthogonal unit cell.
  ! mu = - eL/2*pi Im ln <Psi|exp(-i(2*pi/L)x)|Psi>
  ! E Yaschenko, L Fu, L Resca, R Resta, Phys. Rev. B 58, 1222-1229 (1998)
  ! Single-point Berry`s phase method for dipole should not be used when there is more than one k-point.
  ! in this case, finite differences should be used to construct derivatives with respect to k
  CMPLX function berry_phase_det(st, mesh, dir, ik) result(det)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: mesh
    integer,        intent(in) :: dir
    integer,        intent(in) :: ik

    integer ist, ist2, idim, ip
    CMPLX, allocatable :: matrix(:, :), tmp(:)

    PUSH_SUB(berry_phase_det)

    SAFE_ALLOCATE(matrix(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(tmp(1:mesh%np))

    do ist = 1, st%nst
      do ist2 = 1, st%nst
        matrix(ist, ist2) = M_Z0
        do idim = 1, st%d%dim ! spinor components
            
          if(states_are_complex(st)) then
            forall(ip = 1:mesh%np)
              tmp(ip) = conjg(st%zpsi(ip, idim, ist, ik)) * &
                exp(-M_zI * (M_PI / mesh%sb%lsize(dir)) * mesh%x(ip, dir)) * st%zpsi(ip, idim, ist2, ik)
              ! factor of two removed from exp since actual lattice vector is 2 * lsize
            end forall
          else
            forall(ip = 1:mesh%np)
              tmp(ip) = st%dpsi(ip, idim, ist, ik) * &
                exp(-M_zI * (M_PI / mesh%sb%lsize(dir)) * mesh%x(ip, dir)) * st%dpsi(ip, idim, ist2, ik)
            end forall
          end if
          
          matrix(ist, ist2) = matrix(ist, ist2) + zmf_integrate(mesh, tmp)
        end do
      enddo !ist2
    enddo !ist
      
    det = lalg_determinant(st%nst, matrix(1:st%nst, 1:st%nst), invert = .false.) ** st%smear%el_per_state

    SAFE_DEALLOCATE_A(matrix)
    SAFE_DEALLOCATE_A(tmp)

    POP_SUB(berry_phase_det)
  end function berry_phase_det

 
  ! ---------------------------------------------------------
  ! local potential for electric enthalpy of uniform field in single-point Berry phase
  ! P Umari et al., Phys Rev Lett 95, 207602 (2005) eqs (3), (7)
  ! E * (e L / 2 pi) Im e^(i 2 pi r / L) / z  
  subroutine berry_potential(st, mesh, E_field, pot)
    type(states_t), intent(in)  :: st
    type(mesh_t),   intent(in)  :: mesh
    FLOAT,          intent(in)  :: E_field(:) ! mesh%sb%dim
    FLOAT,          intent(out) :: pot(:,:)   ! mesh%np, st%d%nspin

    integer :: ispin, idir
    CMPLX :: factor, det

    PUSH_SUB(berry_potential)
    
    pot(1:mesh%np, 1:st%d%nspin) = M_ZERO

    do ispin = 1, st%d%nspin
      do idir = 1, mesh%sb%periodic_dim
        if(abs(E_field(idir)) > M_EPSILON) then
          ! calculate the ip-independent part first
          det = berry_phase_det(st, mesh, idir, ispin)
          if(abs(det) .gt. M_EPSILON) then
            factor = E_field(idir) * (mesh%sb%lsize(idir) / M_PI) / det
          else
            ! If det = 0, mu = -infinity, so this condition should never be reached
            ! if things are working properly.
            write(message(1),*) "Divide by zero: dir = ", idir, " Berry-phase determinant = ", det
            call messages_fatal(1)
          endif
          pot(1:mesh%np, ispin) = pot(1:mesh%np, ispin) + &
            aimag(factor * exp(M_PI * M_zI * mesh%x(1:mesh%np, idir) / mesh%sb%lsize(idir)))
        endif
      enddo
    enddo

    POP_SUB(berry_potential)
  end subroutine berry_potential


  ! ---------------------------------------------------------
  FLOAT function berry_energy_correction(st, mesh, E_field, vberry) result(delta)
    type(states_t), intent(in) :: st
    type(mesh_t),   intent(in) :: mesh
    FLOAT,          intent(in) :: E_field(:)  ! mesh%sb%periodic_dim
    FLOAT,          intent(in) :: vberry(:,:) ! mesh%np, st%d%nspin

    integer :: ispin, ist, ip, idir, idim
    FLOAT, allocatable :: vpsi(:,:)

    PUSH_SUB(berry_energy_correction)

    SAFE_ALLOCATE(vpsi(1:mesh%np, 1:st%d%dim))

    ! first we calculate expectation value of Berry potential, to subtract off
    delta = M_ZERO
    do ispin = 1, st%d%nspin
      do ist = 1, st%nst
        forall(ip = 1:mesh%np, idim = 1:st%d%dim) vpsi(ip, idim) = st%dpsi(ip, idim, ist, ispin) * vberry(ip, ispin)
        delta = delta + dmf_dotp(mesh, st%d%dim, st%dpsi(1:mesh%np, 1:st%d%dim, ist, ispin), vpsi(1:mesh%np, 1:st%d%dim))
      enddo
    enddo
    delta = -delta * st%smear%el_per_state

    ! the real energy contribution is -mu.E
    do idir = 1, mesh%sb%periodic_dim
      delta = delta - berry_dipole(st, mesh, idir) * E_field(idir)
    enddo

    SAFE_DEALLOCATE_A(vpsi)
    POP_SUB(berry_energy_correction)
  end function berry_energy_correction

end module berry_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
