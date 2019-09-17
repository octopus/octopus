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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module berry_oct_m
  use global_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m

  implicit none

  private
  public ::                 &
    berry_dipole,           &
    berry_phase_det,        &
    berry_phase_matrix,     &
    berry_potential,        &
    berry_energy_correction

contains

  ! ---------------------------------------------------------
  !> Uses the single-point Berry`s phase method to calculate dipole moment in a periodic system
  !!
  !! This is only accurate in the limit of a large supercell.
  !! It is implemented only for an orthogonal unit cell.
  !! \f[
  !! \mu = - eL/2*\pi Im ln <\Psi|exp(-i(2*\pi/L)x)|\Psi>
  !! \f]
  !! E Yaschenko, L Fu, L Resca, R Resta, Phys. Rev. B 58, 1222-1229 (1998)
  !! Single-point Berry`s phase method for dipole should not be used when there is more than one k-point.
  !! in this case, finite differences should be used to construct derivatives with respect to k
  FLOAT function berry_dipole(st, mesh, dir) result(dipole)
    type(states_elec_t), intent(in) :: st
    type(mesh_t),        intent(in) :: mesh
    integer,             intent(in) :: dir

    integer :: ik
    CMPLX :: det

    PUSH_SUB(berry_dipole)

    dipole = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end ! determinants for different spins and k-points multiply since matrix is block-diagonal
      det  = berry_phase_det(st, mesh, dir, ik)
      dipole = dipole + aimag(log(det))
    end do

    dipole = -(mesh%sb%lsize(dir) / M_PI) * dipole

    POP_SUB(berry_dipole)
  end function berry_dipole


  ! ---------------------------------------------------------
  !! \f[
  !! det <\psi_i|exp(-i(2*\pi/L)x)|\psi_j>
  !! \f]
  !! E Yaschenko, L Fu, L Resca, R Resta, Phys. Rev. B 58, 1222-1229 (1998)
  CMPLX function berry_phase_det(st, mesh, dir, ik) result(det)
    type(states_elec_t), intent(in) :: st
    type(mesh_t),        intent(in) :: mesh
    integer,             intent(in) :: dir
    integer,             intent(in) :: ik

    integer :: ist, noccst, gvector(3)
    CMPLX, allocatable :: matrix(:, :), tmp(:), phase(:)
    CMPLX, allocatable :: psi(:, :), psi2(:, :)

    PUSH_SUB(berry_phase_det)

    ! find how many states are occupied on this k-point. Formalism only works if semiconducting anyway.
    noccst = 0
    do ist = 1, st%nst
      if(st%occ(ist, ik) > M_EPSILON) noccst = ist
    end do

    SAFE_ALLOCATE(matrix(1:noccst, 1:noccst))

    gvector(:) = 0
    gvector(dir) = 1
    call berry_phase_matrix(st, mesh, noccst, ik, ik, gvector, matrix)
      
    if(noccst > 0) then
      det = lalg_determinant(noccst, matrix(1:noccst, 1:noccst), invert = .false.) ** st%smear%el_per_state
    else
      det = M_ONE
    end if

    SAFE_DEALLOCATE_A(matrix)
    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(phase)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(psi2)

    POP_SUB(berry_phase_det)
  end function berry_phase_det


  ! ---------------------------------------------------------
  subroutine berry_phase_matrix(st, mesh, nst, ik, ik2, gvector, matrix)
    type(states_elec_t), intent(in)  :: st
    type(mesh_t),        intent(in)  :: mesh
    integer,             intent(in)  :: nst
    integer,             intent(in)  :: ik
    integer,             intent(in)  :: ik2
    integer,             intent(in)  :: gvector(:) !< (3)
    CMPLX,               intent(out) :: matrix(:,:) !< (nst, nst)

    integer :: ist, ist2, idim, ip, idir
    CMPLX, allocatable :: tmp(:), phase(:)
    CMPLX, allocatable :: psi(:, :), psi2(:, :)
    ! FIXME: real/cplx versions

    PUSH_SUB(berry_phase_matrix)

    if(st%parallel_in_states) then
      call messages_not_implemented("Berry phase parallel in states")
    end if

    SAFE_ALLOCATE(tmp(1:mesh%np))
    SAFE_ALLOCATE(phase(1:mesh%np))
    SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(psi2(1:mesh%np, 1:st%d%dim))

    phase = M_ZERO
    do idir = 1, mesh%sb%dim
      if(gvector(idir) == 0) cycle
      forall(ip = 1:mesh%np)
        phase(ip) = phase(ip) + exp(-M_zI*gvector(idir)*(M_PI/mesh%sb%lsize(idir))*mesh%x(ip, idir))
        ! factor of two removed from exp since actual lattice vector is 2*lsize
      end forall
    end do

    do ist = 1, nst
      call states_elec_get_state(st, mesh, ist, ik, psi)
      do ist2 = 1, nst
        call states_elec_get_state(st, mesh, ist2, ik2, psi2)
        matrix(ist, ist2) = M_Z0
        do idim = 1, st%d%dim ! spinor components
            
          forall(ip = 1:mesh%np)
            tmp(ip) = conjg(psi(ip, idim))*phase(ip)*psi2(ip, idim)
          end forall
          
          matrix(ist, ist2) = matrix(ist, ist2) + zmf_integrate(mesh, tmp)
        end do
      end do !ist2
    end do !ist
      
    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(phase)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(psi2)

    POP_SUB(berry_phase_matrix)
  end subroutine berry_phase_matrix

 
  ! ---------------------------------------------------------
  !> local potential for electric enthalpy of uniform field in single-point Berry phase
  !!
  !! P Umari et al., Phys Rev Lett 95, 207602 (2005) eqs (3), (7)
  !! \f[
  !! E * (e L / 2 \pi) Im e^{i 2 \pi r / L} / z  
  !! \f]
  subroutine berry_potential(st, mesh, E_field, pot)
    type(states_elec_t), intent(in)  :: st
    type(mesh_t),        intent(in)  :: mesh
    FLOAT,               intent(in)  :: E_field(:) !< (mesh%sb%dim)
    FLOAT,               intent(out) :: pot(:,:)   !< (mesh%np, st%d%nspin)

    integer :: ispin, idir
    CMPLX :: factor, det

    PUSH_SUB(berry_potential)
    
    pot(1:mesh%np, 1:st%d%nspin) = M_ZERO

    do ispin = 1, st%d%nspin
      do idir = 1, mesh%sb%periodic_dim
        if(abs(E_field(idir)) > M_EPSILON) then
          ! calculate the ip-independent part first
          det = berry_phase_det(st, mesh, idir, ispin)
          if(abs(det) > M_EPSILON) then
            factor = E_field(idir) * (mesh%sb%lsize(idir) / M_PI) / det
          else
            ! If det = 0, mu = -infinity, so this condition should never be reached
            ! if things are working properly.
            write(message(1),*) "Divide by zero: dir = ", idir, " Berry-phase determinant = ", det
            call messages_fatal(1, namespace=st%namespace)
          end if
          pot(1:mesh%np, ispin) = pot(1:mesh%np, ispin) + &
            aimag(factor * exp(M_PI * M_zI * mesh%x(1:mesh%np, idir) / mesh%sb%lsize(idir)))
        end if
      end do
    end do

    POP_SUB(berry_potential)
  end subroutine berry_potential


  ! ---------------------------------------------------------
  FLOAT function berry_energy_correction(st, mesh, E_field, vberry) result(delta)
    type(states_elec_t), intent(in) :: st
    type(mesh_t),        intent(in) :: mesh
    FLOAT,               intent(in) :: E_field(:)  !< (mesh%sb%periodic_dim)
    FLOAT,               intent(in) :: vberry(:,:) !< (mesh%np, st%d%nspin)

    integer :: ispin, idir

    PUSH_SUB(berry_energy_correction)

    ! first we calculate expectation value of Berry potential, to subtract off
    delta = M_ZERO
    do ispin = 1, st%d%nspin
      delta = delta - dmf_dotp(mesh, st%rho(1:mesh%np, ispin), vberry(1:mesh%np, ispin))
    end do

    ! the real energy contribution is -mu.E
    do idir = 1, mesh%sb%periodic_dim
      delta = delta - berry_dipole(st, mesh, idir) * E_field(idir)
    end do

    POP_SUB(berry_energy_correction)
  end function berry_energy_correction

end module berry_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
