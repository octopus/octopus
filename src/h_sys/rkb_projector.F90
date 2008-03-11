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
!! $Id: epot.F90 2648 2007-01-09 19:08:10Z lorenzen $

#include "global.h"

module rkb_projector_m
  use global_m
  use grid_m
  use lalg_basic_m
  use mesh_m
  use messages_m
  use simul_box_m
  use submesh_m
  use ps_m
  use specie_m
  use specie_pot_m
  use geometry_m
  use mpi_m
  use mpi_debug_m

  implicit none

  private
  public :: &
       rkb_projector_t,    &
       rkb_projector_null, &
       rkb_projector_init, &
       rkb_project,        &
       rkb_project_bra,    &
       rkb_project_ket,    &
       rkb_projector_end

  ! The rkb_projector data type holds the KB projectors build with total angular
  ! momentum eigenfunctions.
  ! This way the spin-orbit coupling in straighforwardly included.
  type rkb_projector_t
    private
    integer          :: n_s ! number of points inside the sphere
    CMPLX,   pointer :: bra(:, :)
    CMPLX,   pointer :: ket(:, :, :, :)
    FLOAT            :: f(2, 2, 2)

    ! The following variables are only used to compute the forces,
    ! because in that case we will ignore the spin-orbit coupling.
    FLOAT,   pointer :: p(:, :)     ! projectors
    FLOAT            :: e(2)        ! KB energies
  end type rkb_projector_t


contains

  ! ---------------------------------------------------------
  subroutine rkb_projector_null(rkb_p)
    type(rkb_projector_t), intent(out) :: rkb_p

    nullify(rkb_p%bra)
    nullify(rkb_p%ket)
    nullify(rkb_p%p)

  end subroutine rkb_projector_null

  ! ---------------------------------------------------------
  subroutine rkb_projector_init(rkb_p, sm, gr, a, l, lm)
    type(rkb_projector_t), intent(inout) :: rkb_p
    type(submesh_t),       intent(in)    :: sm
    type(grid_t),          intent(in)    :: gr
    type(atom_t),          intent(in)    :: a
    integer,               intent(in)    :: l, lm
 
    integer :: is, i
    FLOAT :: v, dv(3), x(MAX_DIM)
    CMPLX :: zv



    rkb_p%n_s = sm%ns

    !Allocate memory
    ALLOCATE(rkb_p%bra(rkb_p%n_s, 2),        rkb_p%n_s*2)
    ALLOCATE(rkb_p%ket(rkb_p%n_s, 2, 2, 2),  rkb_p%n_s*2*2*2)
    ALLOCATE(rkb_p%p(rkb_p%n_s, 2),          rkb_p%n_s*2)

    !Build projectors
    do is = 1, rkb_p%n_s
      x(1:MAX_DIM) = sm%x(is, 1:MAX_DIM)

      ! i runs over j=l+1/2 and j=l-1/2
      do i = 1, 2
        call specie_nl_projector(a%spec, x, l, lm, i, zv)
        rkb_p%bra(is, i) = conjg(zv)
        
        rkb_p%ket(is, i, 1, 1) = zv
        rkb_p%ket(is, i, 2, 2) = zv
        if (lm /= l) then
          call specie_nl_projector(a%spec, x, l, lm+1, i, zv)
          rkb_p%ket(is, i, 2, 1) = zv
        else
          rkb_p%ket(is, i, 2, 1) = M_z0
        end if
        if (lm /= -l) then
          call specie_nl_projector(a%spec, x, l, lm-1, i, zv)
          rkb_p%ket(is, i, 1, 2) = zv
        else
          rkb_p%ket(is, i, 1, 2) = M_z0
        end if
        
        call specie_real_nl_projector(a%spec, x, l, lm, i, v, dv)
        rkb_p%p(is, i) = v
      end do
    end do
    
    ! The l and m dependent prefactors are included in the KB energies
    rkb_p%f(1, 1, 1) = real(l + lm + 1, REAL_PRECISION)
    rkb_p%f(1, 2, 1) = sqrt(real((l + lm + 1)*(l - lm), REAL_PRECISION))
    rkb_p%f(1, 1, 2) = sqrt(real((l - lm + 1)*(l + lm), REAL_PRECISION))
    rkb_p%f(1, 2, 2) = real(l - lm + 1, REAL_PRECISION)
    rkb_p%f(2, 1, 1) = real(l - lm, REAL_PRECISION)
    rkb_p%f(2, 2, 1) = -sqrt(real((l + lm + 1)*(l - lm), REAL_PRECISION))
    rkb_p%f(2, 1, 2) = -sqrt(real((l - lm + 1)*(l + lm), REAL_PRECISION))
    rkb_p%f(2, 2, 2) = real(l + lm, REAL_PRECISION)
    rkb_p%f = rkb_p%f/real(2*l + 1, REAL_PRECISION)

    rkb_p%f(1, :, :) = rkb_p%f(1, :, :)*a%spec%ps%h(l, 1, 1)
    rkb_p%f(2, :, :) = rkb_p%f(2, :, :)*a%spec%ps%h(l, 2, 2)

    ! The projectors used to compute the forces should be averaged. 
    ! The weights are included in the KB energies
    rkb_p%e(1) = a%spec%ps%h(l, 1, 1)*real(l+1, REAL_PRECISION)/real(2*l+1, REAL_PRECISION)
    rkb_p%e(2) = a%spec%ps%h(l, 2, 2)*real(l,   REAL_PRECISION)/real(2*l+1, REAL_PRECISION)

  end subroutine rkb_projector_init

  ! ---------------------------------------------------------
  subroutine rkb_projector_end(rkb_p)
    type(rkb_projector_t), intent(inout) :: rkb_p

    if (associated(rkb_p%bra)) deallocate(rkb_p%bra)
    if (associated(rkb_p%ket)) deallocate(rkb_p%ket)
    if (associated(rkb_p%p))   deallocate(rkb_p%p)

  end subroutine rkb_projector_end

  ! ---------------------------------------------------------
  subroutine rkb_project(mesh, sm, rkb_p, psi, ppsi)
    type(mesh_t),          intent(in)  :: mesh
    type(submesh_t),       intent(in)  :: sm
    type(rkb_projector_t), intent(in)  :: rkb_p
    CMPLX,                 intent(in)  :: psi(:, :)  ! psi(kb%n_s, 2)
    CMPLX,                 intent(out) :: ppsi(:, :) ! ppsi(kb%n_s, 2)

    integer :: j, idim, jdim, n_s
    CMPLX :: uvpsi

    call push_sub('rkb_projector.rkb_project')

    n_s = rkb_p%n_s
    ppsi = M_z0

    do idim = 1, 2   
      do j = 1, 2

        if (all(rkb_p%f(j, :, idim) == M_ZERO)) cycle

        uvpsi = zzsm_integrate_prod(mesh, sm, psi(1:n_s, idim), rkb_p%bra(1:n_s, j))

        do jdim = 1, 2
          if (rkb_p%f(j, jdim, idim) == M_ZERO) cycle
          call lalg_axpy(n_s, rkb_p%f(j, jdim, idim)*uvpsi, rkb_p%ket(:, j, jdim, idim), ppsi(:, jdim))
        end do

      end do
    end do
    
    call pop_sub()
  end subroutine rkb_project


  ! ---------------------------------------------------------
  subroutine rkb_project_bra(mesh, sm, rkb_p, psi, uvpsi, phase)
    type(mesh_t),          intent(in)  :: mesh
    type(submesh_t),       intent(in)  :: sm
    type(rkb_projector_t), intent(in)  :: rkb_p
    CMPLX,                 intent(in)  :: psi(:, :)
    CMPLX,                 intent(out) :: uvpsi(1:2, 1:2)    
    CMPLX, optional,       intent(in)  :: phase(:)

    integer :: idim, n_s, ip, is

    CMPLX, allocatable :: bra(:, :)

    call push_sub('rkb_projector.rkb_project_bra')

    uvpsi = M_ZERO

    n_s = rkb_p%n_s

    ALLOCATE(bra(1:n_s, 2), n_s*2)

    if(mesh%use_curvlinear) then
      bra(1:n_s, 1) = rkb_p%bra(1:n_s, 1)*mesh%vol_pp(sm%jxyz(1:n_s))
      bra(1:n_s, 2) = rkb_p%bra(1:n_s, 2)*mesh%vol_pp(sm%jxyz(1:n_s))
    else
      bra(1:n_s, 1:2) = rkb_p%bra(1:n_s, 1:2)*mesh%vol_pp(1)
    end if

    if(present(phase)) then
      bra(1:n_s, 1) = bra(1:n_s, 1)*phase(1:n_s)
      bra(1:n_s, 2) = bra(1:n_s, 2)*phase(1:n_s)
    end if

    do idim = 1, 2
      do is = 1, n_s
        ip = sm%jxyz(is)
        uvpsi(idim, 1) = uvpsi(idim, 1) + psi(ip, idim)*bra(is, 1)
        uvpsi(idim, 2) = uvpsi(idim, 2) + psi(ip, idim)*bra(is, 2)
      end do
    end do

    deallocate(bra)
  end subroutine rkb_project_bra

  ! ---------------------------------------------------------
  subroutine rkb_project_ket(mesh, sm, rkb_p, uvpsi, psi, phase)
    type(mesh_t),          intent(in)    :: mesh
    type(submesh_t),       intent(in)    :: sm
    type(rkb_projector_t), intent(in)    :: rkb_p
    CMPLX,                 intent(in)    :: uvpsi(1:2, 1:2)
    CMPLX,                 intent(inout) :: psi(:, :)
    CMPLX, optional,       intent(in)    :: phase(:)

    integer :: idim, jdim, n_s, ip, is
    CMPLX :: aa

    call push_sub('rkb_projector.rkb_project_bra')

    n_s = rkb_p%n_s

    if(.not. present(phase)) then
      
      do jdim = 1, 2
        do is = 1, n_s
          ip = sm%jxyz(is)
          aa = M_z0
          do idim = 1, 2
            aa = aa + rkb_p%f(1, jdim, idim)*uvpsi(idim, 1)*rkb_p%ket(is, 1, jdim, idim)
            aa = aa + rkb_p%f(2, jdim, idim)*uvpsi(idim, 2)*rkb_p%ket(is, 2, jdim, idim)
          end do
          psi(ip, jdim) = psi(ip, jdim) + aa
        end do
      end do
      
    else
      
      do jdim = 1, 2
        do is = 1, n_s
          ip = sm%jxyz(is)
          aa = M_z0
          do idim = 1, 2
            aa = aa + rkb_p%f(1, jdim, idim)*uvpsi(idim, 1)*rkb_p%ket(is, 1, jdim, idim)
            aa = aa + rkb_p%f(2, jdim, idim)*uvpsi(idim, 2)*rkb_p%ket(is, 2, jdim, idim)
          end do
          psi(ip, jdim) = psi(ip, jdim) + aa*conjg(phase(is))
        end do
      end do

    end if
    
    call pop_sub()
  end subroutine rkb_project_ket
  
end module rkb_projector_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
