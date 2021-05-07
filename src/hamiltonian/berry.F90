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
  use eigensolver_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use parser_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use space_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use v_ks_oct_m

  implicit none

  private
  public ::                 &
    berry_t,                &
    berry_init,             &
    berry_perform_internal_scf, &
    calc_dipole,            &
    berry_dipole,           &
    berry_potential,        &
    berry_energy_correction

  type berry_t
    private

    integer :: max_iter_berry  !< max number of electronic iterations before updating density, for Berr  y potential
  end type berry_t
contains

  ! ---------------------------------------------------------
  subroutine berry_init(this, namespace)
    type(berry_t),     intent(inout) :: this
    type(namespace_t), intent(in)    :: namespace

    PUSH_SUB(berry_init)

    !%Variable MaximumIterBerry
    !%Type integer
    !%Default 10
    !%Section SCF::Convergence
    !%Description
    !% Maximum number of iterations for the Berry potential, within each SCF iteration.
    !% Only applies if a <tt>StaticElectricField</tt> is applied in a periodic direction.
    !% The code will move on to the next SCF iteration even if convergence
    !% has not been achieved. -1 means unlimited.
    !%End
    call parse_variable(namespace, 'MaximumIterBerry', 10, this%max_iter_berry)
    if(this%max_iter_berry < 0) this%max_iter_berry = huge(this%max_iter_berry)


    POP_SUB(berry_init)
  end subroutine berry_init

  ! ---------------------------------------------------------
  subroutine berry_perform_internal_scf(this, namespace, space, eigens, gr, st, hm, iter, ks, ions)
    type(berry_t),            intent(in)    :: this
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(eigensolver_t),      intent(inout) :: eigens
    type(grid_t),             intent(in)    :: gr
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    integer,                  intent(in)    :: iter
    type(v_ks_t),             intent(inout) :: ks
    type(ions_t),             intent(in)    :: ions

    integer :: iberry, idir
    logical :: berry_conv
    FLOAT :: dipole_prev(1:MAX_DIM), dipole(1:MAX_DIM)
    FLOAT, parameter :: tol = CNST(1e-5)

    PUSH_SUB(berry_perform_internal_scf)

    ASSERT(allocated(hm%vberry))

    if(st%parallel_in_states) then
      call messages_not_implemented("Berry phase parallel in states", namespace=namespace)
    end if

    call calc_dipole(dipole, space, gr%mesh, st, ions)

    do iberry = 1, this%max_iter_berry
      eigens%converged = 0
      call eigensolver_run(eigens, namespace, gr, st, hm, iter)

      !Calculation of the Berry potential
      call berry_potential(st, namespace, space, gr%mesh, hm%ep%E_field, hm%vberry)

      !Calculation of the corresponding energy 
      hm%energy%berry = berry_energy_correction(st, space, gr%mesh, &
         hm%ep%E_field(1:space%periodic_dim), hm%vberry(1:gr%mesh%np, 1:hm%d%nspin))
  
      !We recompute the KS potential
      call v_ks_calc(ks, namespace, space, hm, st, ions, calc_current=.false.)

      dipole_prev = dipole
      call calc_dipole(dipole, space, gr%mesh, st, ions)
      write(message(1),'(a,9f12.6)') 'Dipole = ', dipole(1:space%dim)
      call messages_info(1)
  
      berry_conv = .true.
      do idir = 1, space%periodic_dim
        if(abs(dipole_prev(idir)) > CNST(1e-10)) then
          berry_conv = berry_conv .and. (abs(dipole(idir) - dipole_prev(idir)) < tol &
            .or.(abs((dipole(idir) - dipole_prev(idir)) / dipole_prev(idir)) < tol))
        else
          berry_conv = berry_conv .and. (abs(dipole(idir) - dipole_prev(idir)) < tol)
        end if
      end do

      if(berry_conv) exit
    end do

    POP_SUB(berry_perform_internal_scf)
  end subroutine berry_perform_internal_scf

  ! ---------------------------------------------------------
  !TODO: This should be a method of the electronic system directly
  ! that can be exposed
  subroutine calc_dipole(dipole, space, mesh, st, ions)
    FLOAT,                 intent(out)   :: dipole(:)
    type(space_t),         intent(in)    :: space
    type(mesh_t),          intent(in)    :: mesh
    type(states_elec_t),   intent(in)    :: st
    type(ions_t),          intent(in)    :: ions

    integer :: ispin, idir
    FLOAT :: e_dip(space%dim + 1, st%d%nspin), n_dip(space%dim), nquantumpol

    PUSH_SUB(calc_dipole)

    ASSERT(.not. ions%latt%nonorthogonal)

    dipole(1:space%dim) = M_ZERO

    do ispin = 1, st%d%nspin
      call dmf_multipoles(mesh, st%rho(:, ispin), 1, e_dip(:, ispin))
    end do

    n_dip = ions%dipole()

    do idir = 1, space%dim
      ! in periodic directions use single-point Berry`s phase calculation
      if(idir  <=  space%periodic_dim) then
        dipole(idir) = -n_dip(idir) - berry_dipole(st, mesh, idir)

        ! use quantum of polarization to reduce to smallest possible magnitude
        nquantumpol = nint(dipole(idir)/norm2(ions%latt%rlattice(:, idir)))
        dipole(idir) = dipole(idir) - nquantumpol * norm2(ions%latt%rlattice(:, idir))
       ! in aperiodic directions use normal dipole formula

      else
        e_dip(idir + 1, 1) = sum(e_dip(idir + 1, :))
        dipole(idir) = -n_dip(idir) - e_dip(idir + 1, 1)
      end if
    end do

    POP_SUB(calc_dipole)
  end subroutine calc_dipole


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

    !There is a missing reduction here. 
    !However, as we are limited to a single k-point at the moment, this is not a problem.
    ASSERT(st%d%nik == 1)

    dipole = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end ! determinants for different spins and k-points multiply since matrix is block-diagonal
      det  = berry_phase_det(st, mesh, dir, ik)
      dipole = dipole + aimag(log(det))
    end do

    dipole = -(norm2(mesh%sb%latt%rlattice(:,dir)) / (M_TWO*M_PI)) * dipole

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
      det = lalg_determinant(noccst, matrix(1:noccst, 1:noccst), preserve_mat=.false.) ** st%smear%el_per_state
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

    ASSERT(.not. st%parallel_in_states)
    ASSERT(.not. mesh%sb%latt%nonorthogonal)

    SAFE_ALLOCATE(tmp(1:mesh%np))
    SAFE_ALLOCATE(phase(1:mesh%np))
    SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(psi2(1:mesh%np, 1:st%d%dim))

    phase = M_ZERO
    do idir = 1, mesh%sb%dim
      if(gvector(idir) == 0) cycle
      do ip = 1, mesh%np
        phase(ip) = phase(ip) + exp(-M_zI*gvector(idir)*(M_TWO*M_PI/norm2(mesh%sb%latt%rlattice(:,idir)))*mesh%x(ip, idir))
      end do
    end do

    do ist = 1, nst
      call states_elec_get_state(st, mesh, ist, ik, psi)
      do ist2 = 1, nst
        call states_elec_get_state(st, mesh, ist2, ik2, psi2)
        matrix(ist, ist2) = M_Z0
        do idim = 1, st%d%dim ! spinor components

          do ip = 1, mesh%np
            tmp(ip) = conjg(psi(ip, idim))*phase(ip)*psi2(ip, idim)
          end do

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
  subroutine berry_potential(st, namespace, space, mesh, E_field, pot)
    type(states_elec_t), intent(in)  :: st
    type(namespace_t),   intent(in)  :: namespace
    type(space_t),       intent(in)  :: space
    type(mesh_t),        intent(in)  :: mesh
    FLOAT,               intent(in)  :: E_field(:) !< (space%dim)
    FLOAT,               intent(out) :: pot(:,:)   !< (mesh%np, st%d%nspin)

    integer :: ispin, idir
    CMPLX :: factor, det

    PUSH_SUB(berry_potential)

    if(mesh%sb%latt%nonorthogonal) then
      call messages_not_implemented("Berry phase for non-orthogonal cells.")
    end if
  
    if(st%d%nik > 1) then
      call messages_not_implemented("Berry phase with k-points.")
    end if

    pot(1:mesh%np, 1:st%d%nspin) = M_ZERO

    do ispin = 1, st%d%nspin
      do idir = 1, space%periodic_dim
        if(abs(E_field(idir)) > M_EPSILON) then
          ! calculate the ip-independent part first
          det = berry_phase_det(st, mesh, idir, ispin)
          if(abs(det) > M_EPSILON) then
            factor = E_field(idir) * (norm2(mesh%sb%latt%rlattice(:,idir)) / (M_TWO*M_PI)) / det
          else
            ! If det = 0, mu = -infinity, so this condition should never be reached
            ! if things are working properly.
            write(message(1),*) "Divide by zero: dir = ", idir, " Berry-phase determinant = ", det
            call messages_fatal(1, namespace=namespace)
          end if
          pot(1:mesh%np, ispin) = pot(1:mesh%np, ispin) + &
            aimag(factor * exp(M_TWO * M_PI * M_zI * mesh%x(1:mesh%np, idir) / norm2(mesh%sb%latt%rlattice(:,idir))))
        end if
      end do
    end do

    POP_SUB(berry_potential)
  end subroutine berry_potential


  ! ---------------------------------------------------------
  FLOAT function berry_energy_correction(st, space, mesh, E_field, vberry) result(delta)
    type(states_elec_t), intent(in) :: st
    type(space_t),       intent(in) :: space
    type(mesh_t),        intent(in) :: mesh
    FLOAT,               intent(in) :: E_field(:)  !< (space%periodic_dim)
    FLOAT,               intent(in) :: vberry(:,:) !< (mesh%np, st%d%nspin)

    integer :: ispin, idir

    PUSH_SUB(berry_energy_correction)

    ! first we calculate expectation value of Berry potential, to subtract off
    delta = M_ZERO
    do ispin = 1, st%d%nspin
      delta = delta - dmf_dotp(mesh, st%rho(1:mesh%np, ispin), vberry(1:mesh%np, ispin))
    end do

    ! the real energy contribution is -mu.E
    do idir = 1, space%periodic_dim
      delta = delta - berry_dipole(st, mesh, idir) * E_field(idir)
    end do

    POP_SUB(berry_energy_correction)
  end function berry_energy_correction

end module berry_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
