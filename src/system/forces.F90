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
!! $Id$

#include "global.h"

module forces_m
  use batch_m
  use batch_ops_m
  use born_charges_m
#ifdef HAVE_OPENCL
  use cl
#endif
  use comm_m
  use derivatives_m
  use epot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use hamiltonian_base_m
  use index_m
  use io_m
  use kpoints_m
  use lalg_basic_m
  use lasers_m
  use linear_response_m
  use loct_math_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use profiling_m
  use projector_m
  use octcl_kernel_m
  use opencl_m
  use simul_box_m
  use species_m
  use species_pot_m
  use states_m
  use states_dim_m
  use symm_op_m
  use symmetrizer_m
  use types_m

  implicit none

  private
  public ::                    &
    forces_calculate,          &
    dforces_from_potential,    &
    zforces_from_potential,    &
    dforces_derivative,        &
    zforces_derivative,        &
    dforces_born_charges,      &
    zforces_born_charges,      &
    total_force_calculate,     &
    forces_costate_calculate

  type(profile_t), save :: prof_comm


  type(geometry_t), pointer :: geo_
  type(grid_t), pointer :: gr_
  type(hamiltonian_t), pointer :: hm_
  type(states_t), pointer :: psi_
  type(states_t), pointer :: chi_
  integer, pointer :: i_, j_, ist_, ik_, iatom_
  CMPLX, allocatable :: derpsi_(:, :, :)

contains

  ! ---------------------------------------------------------
  !> This computes the total forces on the ions created by the electrons
  !! (it excludes the force due to possible time-dependent external fields).
  subroutine total_force_calculate(gr, geo, ep, st, x)
    type(grid_t),     intent(inout) :: gr
    type(geometry_t), intent(in)    :: geo
    type(epot_t),     intent(inout) :: ep
    type(states_t),   intent(inout) :: st
    FLOAT, intent(inout)            :: x(MAX_DIM)

    type(profile_t), save :: forces_prof

    call profiling_in(forces_prof, "FORCES")
    PUSH_SUB(total_force_calculate)

    x = M_ZERO
    if (states_are_real(st) ) then 
      call dtotal_force_from_potential(gr, geo, ep, st, x)
    else
      call ztotal_force_from_potential(gr, geo, ep, st, x)
    end if

    POP_SUB(total_force_calculate)
    call profiling_out(forces_prof)
  end subroutine total_force_calculate


  subroutine forces_costate_calculate(gr, geo, hm, psi, chi, f, q)
    type(grid_t), target, intent(inout) :: gr
    type(geometry_t), target, intent(inout) :: geo
    type(hamiltonian_t), target, intent(inout) :: hm
    type(states_t), target, intent(inout) :: psi
    type(states_t), target, intent(inout) :: chi
    FLOAT,            intent(inout) :: f(:, :)
    FLOAT,            intent(in)    :: q(:, :)

    integer :: jatom, idim, jdim
    integer, target :: i, j, ist, ik, iatom
    FLOAT :: r, w2r_, w1r_, xx(MAX_DIM), gq, fq, dq, abserr
    type(profile_t), save :: forces_prof

    call profiling_in(forces_prof, "FORCES")
    PUSH_SUB(forces_costate_calculate)

    f = M_ZERO
    do iatom = 1, geo%natoms
      do jatom = 1, geo%natoms
        if(jatom == iatom) cycle
        xx(1:gr%sb%dim) = geo%atom(jatom)%x(1:gr%sb%dim) - geo%atom(iatom)%x(1:gr%sb%dim)
        r = sqrt( sum( xx(1:gr%sb%dim)**2 ) )
        w2r_ = w2r(geo%atom(iatom)%spec, geo%atom(jatom)%spec, r)
        w1r_ = w1r(geo%atom(iatom)%spec, geo%atom(jatom)%spec, r)
        do idim = 1, gr%sb%dim
          do jdim = 1, gr%sb%dim
            f(iatom, idim) = f(iatom, idim) + (q(jatom, jdim) - q(iatom, jdim)) * w2r_ * (M_ONE/r**2) * xx(idim) * xx(jdim)
            f(iatom, idim) = f(iatom, idim) - (q(jatom, jdim) - q(iatom, jdim)) * w1r_ * (M_ONE/r**3) * xx(idim) * xx(jdim)
            if(jdim == idim) then
              f(iatom, idim) = f(iatom, idim) + (q(jatom, jdim) - q(iatom, jdim)) * w1r_ * (M_ONE/r)
            end if
          end do
        end do
      end do
    end do

    SAFE_ALLOCATE(derpsi_(1:gr%mesh%np_part, 1:gr%sb%dim, 1:psi%d%dim))

    dq = CNST(0.005)

    geo_ => geo
    gr_ => gr
    hm_ => hm
    i_ => i
    j_ => j
    ist_ => ist
    ik_ => ik
    iatom_ => iatom
    psi_ => psi
    chi_ => chi
    
    do ist = 1, psi%nst
      do ik = 1, psi%d%nik
        derpsi_ = M_z0
        call zderivatives_grad(gr%der, psi%zpsi(:, 1, ist, ik), derpsi_(:, :, 1))
        do iatom = 1, geo%natoms
          do i = 1, gr%sb%dim
            do j = 1, gr%sb%dim
              call loct_numerical_derivative(geo%atom(iatom)%x(j), dq, gq, abserr, gofq)
              f(iatom, i) = f(iatom, i) - M_TWO * psi%occ(ist, ik) * q(iatom, j) * gq
            end do
            call loct_numerical_derivative(geo%atom(iatom)%x(i), dq, fq, abserr, fofq)
            f(iatom, i) = f(iatom, i) + M_TWO * fq
          end do
        end do
      end do
    end do
    SAFE_DEALLOCATE_A(derpsi_)

    POP_SUB(forces_costate_calculate)
    call profiling_out(forces_prof)

  contains

    FLOAT function wr(speca, specb, r)
      type(species_t), intent(in) :: speca, specb
      FLOAT, intent(in) :: r
      wr = species_zval(speca) * species_zval(specb) / r
    end function wr
    
    FLOAT function w1r(speca, specb, r)
      type(species_t), intent(in) :: speca, specb
      FLOAT, intent(in) :: r
      w1r = - species_zval(speca) * species_zval(specb) / r**2
    end function w1r
    
    FLOAT function w2r(speca, specb, r)
      type(species_t), intent(in) :: speca, specb
      FLOAT, intent(in) :: r
      w2r = M_TWO * species_zval(speca) * species_zval(specb) / r**3
    end function w2r
  end subroutine forces_costate_calculate
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine gofq(q, res)
    FLOAT, intent(in) :: q
    FLOAT, intent(inout) :: res

    FLOAT :: qold
    CMPLX, allocatable :: viapsi(:, :)
    qold = geo_%atom(iatom_)%x(j_)
    geo_%atom(iatom_)%x(j_) = q
    SAFE_ALLOCATE(viapsi(1:gr_%mesh%np_part, 1:psi_%d%dim))
    viapsi = M_z0
    call zhamiltonian_apply_atom (hm_, geo_, gr_, iatom_, psi_%zpsi(:, :, ist_, ik_), viapsi)
    res = real( zmf_dotp(gr_%mesh, viapsi(:, 1), derpsi_(:, i_, 1)) , REAL_PRECISION)
    geo_%atom(iatom_)%x(j_) = qold
    SAFE_DEALLOCATE_A(viapsi)
  end subroutine gofq
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine fofq(q, res)
    FLOAT, intent(in) :: q
    FLOAT, intent(inout) :: res

    FLOAT :: qold
    CMPLX, allocatable :: viapsi(:, :)
    qold = geo_%atom(iatom_)%x(i_)
    geo_%atom(iatom_)%x(i_) = q
    SAFE_ALLOCATE(viapsi(1:gr_%mesh%np_part, 1:psi_%d%dim))
    viapsi = M_z0
    call zhamiltonian_apply_atom (hm_, geo_, gr_, iatom_, psi_%zpsi(:, :, ist_, ik_), viapsi)
    res = real(M_zI * zmf_dotp(gr_%mesh, chi_%zpsi(:, 1, ist_, ik_), viapsi(:, 1)), REAL_PRECISION)
    geo_%atom(iatom_)%x(i_) = qold
    SAFE_DEALLOCATE_A(viapsi)
  end subroutine fofq
  ! ---------------------------------------------------------


  ! ---------------------------------------------------------
  subroutine forces_calculate(gr, geo, hm, st, t, dt)
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(hamiltonian_t), intent(inout) :: hm
    type(states_t),      intent(inout) :: st
    FLOAT,     optional, intent(in)    :: t
    FLOAT,     optional, intent(in)    :: dt

    integer :: i, j, iatom, idir
    FLOAT :: x(MAX_DIM), time
    FLOAT, allocatable :: force(:, :), total_force(:)
    type(profile_t), save :: forces_prof

    call profiling_in(forces_prof, "FORCES")
    PUSH_SUB(forces_calculate)

    x(:) = M_ZERO
    time = M_ZERO
    if(present(t)) time = t

    ! the ion-ion term is already calculated
    do i = 1, geo%natoms
      geo%atom(i)%f(1:gr%sb%dim) = hm%ep%fii(1:gr%sb%dim, i)
    end do

    SAFE_ALLOCATE(force(1:gr%mesh%sb%dim, 1:geo%natoms))
    
    if (states_are_real(st) ) then 
      call dforces_from_potential(gr, geo, hm, st, force)
    else
      call zforces_from_potential(gr, geo, hm, st, force)
    end if

    if(hm%ep%force_total_enforce) call forces_set_total_to_zero(geo, force)

    do iatom = 1, geo%natoms
      do idir = 1, gr%mesh%sb%dim
        geo%atom(iatom)%f(idir) = geo%atom(iatom)%f(idir) + force(idir, iatom)
      end do
    end do

    SAFE_DEALLOCATE_A(force)
    
    !\todo forces due to the magnetic fields (static and time-dependent)
    if(present(t)) then
      do j = 1, hm%ep%no_lasers
        select case(laser_kind(hm%ep%lasers(j)))
        case(E_FIELD_ELECTRIC)
          x(1:gr%sb%dim) = M_ZERO
          call laser_field(hm%ep%lasers(j), x(1:gr%sb%dim), t)
          do i = 1, geo%natoms
            ! Here the proton charge is +1, since the electric field has the usual sign.
            geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) &
             + species_zval(geo%atom(i)%spec)*x(1:gr%mesh%sb%dim)
          end do
    
        case(E_FIELD_VECTOR_POTENTIAL)
          ! Forces are correctly calculated only if the time-dependent
          ! vector potential has no spatial dependence.
          ! The full force taking account of the spatial dependence of A should be:
          ! F = q [- dA/dt + v x \nabla x A]

          x(1:gr%sb%dim) = M_ZERO
          call laser_electric_field(hm%ep%lasers(j), x(1:gr%sb%dim), t, dt) !convert in E field (E = -dA/ c dt)
          do i = 1, geo%natoms
            ! Also here the proton charge is +1
            geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) &
             + species_zval(geo%atom(i)%spec)*x(1:gr%mesh%sb%dim)
          end do

        case(E_FIELD_MAGNETIC, E_FIELD_SCALAR_POTENTIAL)
          write(message(1),'(a)') 'The forces are currently not properly calculated if time-dependent'
          write(message(2),'(a)') 'magnetic fields are present.'
          call messages_fatal(2)
        end select
      end do
    end if

    if(associated(hm%ep%E_field)) then
      do i = 1, geo%natoms
        ! Here the proton charge is +1, since the electric field has the usual sign.
        geo%atom(i)%f(1:gr%mesh%sb%dim) = geo%atom(i)%f(1:gr%mesh%sb%dim) &
          + species_zval(geo%atom(i)%spec)*hm%ep%E_field(1:gr%mesh%sb%dim)
      end do
    end if
    
    POP_SUB(forces_calculate)
    call profiling_out(forces_prof)

  end subroutine forces_calculate

  ! ----------------------------------------------------------------------

  subroutine forces_set_total_to_zero(geo, force)
    type(geometry_t),    intent(in)    :: geo
    FLOAT,               intent(inout) :: force(:, :)

    FLOAT, allocatable :: total_force(:)
    integer :: iatom

    SAFE_ALLOCATE(total_force(1:geo%space%dim))

    total_force(1:geo%space%dim) = CNST(0.0)
    do iatom = 1, geo%natoms
      total_force(1:geo%space%dim) = total_force(1:geo%space%dim) + force(1:geo%space%dim, iatom)/geo%natoms
    end do

    do iatom = 1, geo%natoms
      force(1:geo%space%dim, iatom) = force(1:geo%space%dim, iatom) - total_force(1:geo%space%dim)
    end do

    SAFE_DEALLOCATE_A(total_force)

  end subroutine forces_set_total_to_zero

#include "undef.F90"
#include "real.F90"
#include "forces_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "forces_inc.F90"

end module forces_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
