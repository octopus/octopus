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
  use profiling_m
  use ps_m
  use species_m
  use geometry_m
  use mpi_m

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

  !> The rkb_projector data type holds the KB projectors build with total angular
  !! momentum eigenfunctions.
  !! This way the spin-orbit coupling in straighforwardly included.
  type rkb_projector_t
    private
    integer          :: n_s !< number of points inside the sphere
    CMPLX,   pointer :: bra(:, :)
    CMPLX,   pointer :: ket(:, :, :, :)
    FLOAT            :: f(2, 2, 2)
  end type rkb_projector_t


contains

  ! ---------------------------------------------------------
  subroutine rkb_projector_null(rkb_p)
    type(rkb_projector_t), intent(out) :: rkb_p

    PUSH_SUB(rkb_projector_null)

    nullify(rkb_p%bra)
    nullify(rkb_p%ket)

    POP_SUB(rkb_projector_null)
  end subroutine rkb_projector_null

  ! ---------------------------------------------------------
  subroutine rkb_projector_init(rkb_p, sm, a, l, lm, so_strength)
    type(rkb_projector_t), intent(inout) :: rkb_p
    type(submesh_t),       intent(in)    :: sm
    type(atom_t),          intent(in)    :: a
    integer,               intent(in)    :: l, lm
    FLOAT,                 intent(in)    :: so_strength
 
    integer :: is, i
    FLOAT :: x(MAX_DIM)
    CMPLX :: zv
    type(ps_t), pointer :: ps

    PUSH_SUB(rkb_projector_init)

    rkb_p%n_s = sm%np

    !Allocate memory
    SAFE_ALLOCATE(rkb_p%bra(1:rkb_p%n_s, 1:2))
    SAFE_ALLOCATE(rkb_p%ket(1:rkb_p%n_s, 1:2, 1:2, 1:2))

    !Build projectors
    do is = 1, rkb_p%n_s
      x(1:ubound(sm%x, 2)) = sm%x(is, 1:ubound(sm%x, 2))

      ! i runs over j=l+1/2 and j=l-1/2
      do i = 1, 2
        call species_nl_projector(a%spec, x, l, lm, i, zv)
        rkb_p%bra(is, i) = conjg(zv)
        
        rkb_p%ket(is, i, 1, 1) = zv
        rkb_p%ket(is, i, 2, 2) = zv
        if (lm /= l) then
          call species_nl_projector(a%spec, x, l, lm+1, i, zv)
          rkb_p%ket(is, i, 2, 1) = zv
        else
          rkb_p%ket(is, i, 2, 1) = M_z0
        end if
        if (lm /= -l) then
          call species_nl_projector(a%spec, x, l, lm-1, i, zv)
          rkb_p%ket(is, i, 1, 2) = zv
        else
          rkb_p%ket(is, i, 1, 2) = M_z0
        end if
      end do
    end do
    
    ! The l and m dependent prefactors are included in the KB energies
    rkb_p%f(1, 1, 1) = real(l + so_strength*lm + 1, REAL_PRECISION)
    rkb_p%f(1, 2, 1) = so_strength*sqrt(real((l + lm + 1)*(l - lm), REAL_PRECISION))
    rkb_p%f(1, 1, 2) = so_strength*sqrt(real((l - lm + 1)*(l + lm), REAL_PRECISION))
    rkb_p%f(1, 2, 2) = real(l - so_strength*lm + 1, REAL_PRECISION)
    rkb_p%f(2, 1, 1) = real(l - so_strength*lm, REAL_PRECISION)
    rkb_p%f(2, 2, 1) = -so_strength*sqrt(real((l + lm + 1)*(l - lm), REAL_PRECISION))
    rkb_p%f(2, 1, 2) = -so_strength*sqrt(real((l - lm + 1)*(l + lm), REAL_PRECISION))
    rkb_p%f(2, 2, 2) = real(l + so_strength*lm, REAL_PRECISION)
    rkb_p%f = rkb_p%f/real(2*l + 1, REAL_PRECISION)

    ps => species_ps(a%spec)
    rkb_p%f(1, :, :) = rkb_p%f(1, :, :) * ps%h(l, 1, 1)
    rkb_p%f(2, :, :) = rkb_p%f(2, :, :) * ps%h(l, 2, 2)
    nullify(ps)

    POP_SUB(rkb_projector_init)
  end subroutine rkb_projector_init

  ! ---------------------------------------------------------
  subroutine rkb_projector_end(rkb_p)
    type(rkb_projector_t), intent(inout) :: rkb_p

    PUSH_SUB(rkb_projector_end)

    SAFE_DEALLOCATE_P(rkb_p%bra)
    SAFE_DEALLOCATE_P(rkb_p%ket)

    POP_SUB(rkb_projector_end)
  end subroutine rkb_projector_end

  ! ---------------------------------------------------------
  subroutine rkb_project(mesh, sm, rkb_p, psi, ppsi)
    type(mesh_t),          intent(in)    :: mesh
    type(submesh_t),       intent(in)    :: sm
    type(rkb_projector_t), intent(in)    :: rkb_p
    CMPLX,                 intent(in)    :: psi(:, :)  !< psi(kb%n_s, 2)
    CMPLX,                 intent(inout) :: ppsi(:, :) !< ppsi(kb%n_s, 2)

    CMPLX :: uvpsi(1:2, 1:2)
#ifdef HAVE_MPI
    CMPLX :: uvpsi_tmp(1:2, 1:2)
#endif

    PUSH_SUB(rkb_project)

    call rkb_project_bra(mesh, sm, rkb_p, psi, uvpsi)

#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) then
      call MPI_Allreduce(uvpsi, uvpsi_tmp, 4, MPI_CMPLX, MPI_SUM, mesh%vp%comm, mpi_err)
      uvpsi = uvpsi_tmp
    end if
#endif

    call rkb_project_ket(rkb_p, uvpsi, ppsi)

    POP_SUB(rkb_project)
  end subroutine rkb_project

  ! ---------------------------------------------------------
  !> THREADSAFE
  subroutine rkb_project_bra(mesh, sm, rkb_p, psi, uvpsi)
    type(mesh_t),          intent(in)  :: mesh
    type(submesh_t),       intent(in)  :: sm
    type(rkb_projector_t), intent(in)  :: rkb_p
    CMPLX,                 intent(in)  :: psi(:, :)
    CMPLX,                 intent(out) :: uvpsi(1:2, 1:2)    

    integer :: idim, n_s, is

    CMPLX, allocatable :: bra(:, :)
#ifndef HAVE_OPENMP
    PUSH_SUB(rkb_project_bra)
#endif
    uvpsi = M_ZERO

    n_s = rkb_p%n_s

    SAFE_ALLOCATE(bra(1:n_s, 1:2))

    if(mesh%use_curvilinear) then
      bra(1:n_s, 1) = rkb_p%bra(1:n_s, 1)*mesh%vol_pp(sm%map(1:n_s))
      bra(1:n_s, 2) = rkb_p%bra(1:n_s, 2)*mesh%vol_pp(sm%map(1:n_s))
    else
      bra(1:n_s, 1:2) = rkb_p%bra(1:n_s, 1:2)*mesh%vol_pp(1)
    end if

    do idim = 1, 2
      do is = 1, n_s
        uvpsi(idim, 1) = uvpsi(idim, 1) + psi(is, idim)*bra(is, 1)
        uvpsi(idim, 2) = uvpsi(idim, 2) + psi(is, idim)*bra(is, 2)
      end do
    end do

    SAFE_DEALLOCATE_A(bra)
#ifndef HAVE_OPENMP
    POP_SUB(rkb_project_bra)
#endif
  end subroutine rkb_project_bra

  ! ---------------------------------------------------------
  !> THREADSAFE
  subroutine rkb_project_ket(rkb_p, uvpsi, psi)
    type(rkb_projector_t), intent(in)    :: rkb_p
    CMPLX,                 intent(in)    :: uvpsi(1:2, 1:2)
    CMPLX,                 intent(inout) :: psi(:, :)

    integer :: idim, jdim, n_s, is
    CMPLX :: aa
#ifndef HAVE_OPENMP
    PUSH_SUB(rkb_project_ket)
#endif
    n_s = rkb_p%n_s

    do jdim = 1, 2
      do is = 1, n_s
        aa = M_z0
        do idim = 1, 2
          aa = aa + rkb_p%f(1, jdim, idim)*uvpsi(idim, 1)*rkb_p%ket(is, 1, jdim, idim)
          aa = aa + rkb_p%f(2, jdim, idim)*uvpsi(idim, 2)*rkb_p%ket(is, 2, jdim, idim)
        end do
        psi(is, jdim) = psi(is, jdim) + aa
      end do
    end do
#ifndef HAVE_OPENMP
    POP_SUB(rkb_project_ket)
#endif
  end subroutine rkb_project_ket
  
end module rkb_projector_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
