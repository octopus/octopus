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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: epot.F90 2648 2007-01-09 19:08:10Z lorenzen $

#include "global.h"

module rkb_projector_m
  use global_m
  use grid_m
  use mesh_m
  use messages_m
  use simul_box_m
  use ps_m
  use specie_m
  use specie_pot_m
  use geometry_m
  use mpi_m
  use mpi_debug_m

  private
  public :: &
       rkb_projector_t,    &
       rkb_projector_null, &
       rkb_projector_init, &
       rkb_project,        &
       rkb_dproject,       &
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
    FLOAT,   pointer :: dp(:, :, :) ! projectors derivatives
    FLOAT            :: e(2)        ! KB energies
  end type rkb_projector_t


contains

  ! ---------------------------------------------------------
  subroutine rkb_projector_null(rkb_p)
    type(rkb_projector_t), intent(out) :: rkb_p

    nullify(rkb_p%bra)
    nullify(rkb_p%ket)
    nullify(rkb_p%p)
    nullify(rkb_p%dp)

  end subroutine rkb_projector_null

  ! ---------------------------------------------------------
  subroutine rkb_projector_init(rkb_p, n_s, jxyz, gr, a, l, lm)
    type(rkb_projector_t), intent(inout) :: rkb_p
    integer,               intent(in)    :: n_s
    integer,               intent(in)    :: jxyz(:)
    type(grid_t),          intent(in)    :: gr
    type(atom_t),          intent(in)    :: a
    integer,               intent(in)    :: l, lm
 
    integer :: j, k, i
    FLOAT :: v, dv(3), r, x(3), x_in(3)
    CMPLX :: zv

    rkb_p%n_s = n_s

    !Allocate memory
    ALLOCATE(rkb_p%bra(n_s, 2),        n_s*2)
    ALLOCATE(rkb_p%ket(n_s, 2, 2, 2),  n_s*2*2*2)
    ALLOCATE(rkb_p%p(n_s, 2),          n_s*2)
    ALLOCATE(rkb_p%dp(n_s, 3, 2),      n_s*3*2)

    !Build projectors
    do j = 1, rkb_p%n_s

      do k = 1, 3**gr%sb%periodic_dim
        x_in(:) = gr%m%x(jxyz(j), :) - gr%sb%shift(k,:)
        x(:) = x_in(:) - a%x
        r = sqrt(sum(x*x))
        if (r > a%spec%ps%rc_max + gr%m%h(1)) cycle

        ! i runs over j=l+1/2 and j=l-1/2
        do i = 1, 2
          call specie_nl_projector(a%spec, a%x, x_in, l, lm, i, zv)
          rkb_p%bra(j, i) = conjg(zv)

          rkb_p%ket(j, i, 1, 1) = zv
          rkb_p%ket(j, i, 2, 2) = zv
          if (lm /= l) then
            call specie_nl_projector(a%spec, a%x, x_in, l, lm+1, i, zv)
            rkb_p%ket(j, i, 2, 1) = zv
          else
            rkb_p%ket(j, i, 2, 1) = M_z0
          end if
          if (lm /= -l) then
            call specie_nl_projector(a%spec, a%x, x_in, l, lm-1, i, zv)
            rkb_p%ket(j, i, 1, 2) = zv
          else
            rkb_p%ket(j, i, 1, 2) = M_z0
          end if

          call specie_real_nl_projector(a%spec, a%x, x_in, l, lm, i, v, dv)
          rkb_p%p(j, i) = v
          rkb_p%dp(j, :, i) = dv
        end do
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
    if (associated(rkb_p%dp))  deallocate(rkb_p%dp)

  end subroutine rkb_projector_end

  ! ---------------------------------------------------------
  subroutine rkb_project(mesh, rkb_p, psi, ppsi, phases)
    type(mesh_t),          intent(in)  :: mesh
    type(rkb_projector_t), intent(in)  :: rkb_p
    CMPLX,                 intent(in)  :: psi(:, :)  ! psi(kb%n_s, 2)
    CMPLX,                 intent(out) :: ppsi(:, :) ! ppsi(kb%n_s, 2)
    CMPLX, optional,       intent(in)  :: phases(:)

    integer :: j, idim, jdim, n_s
    CMPLX :: uvpsi
#if defined(HAVE_MPI)
    CMPLX :: tmp
#endif

    n_s = rkb_p%n_s
    ppsi = M_z0

    do idim = 1, 2   
      do j = 1, 2
        if (all(rkb_p%f(j, :, idim) == M_ZERO)) cycle
        uvpsi = sum(psi(1:n_s, idim)*rkb_p%bra(1:n_s, j))

#if defined(HAVE_MPI)
        if(mesh%parallel_in_domains) then
          call MPI_Allreduce(uvpsi, tmp, 1, MPI_CMPLX, MPI_SUM, mesh%vp%comm, mpi_err)
          uvpsi = tmp
        end if
#endif

        do jdim = 1, 2
          if (rkb_p%f(j, jdim, idim) == M_ZERO) cycle
          if (present(phases)) then
            ppsi(1:n_s, jdim) = ppsi(1:n_s, jdim) + &
                 rkb_p%f(j, jdim, idim) * uvpsi * rkb_p%ket(1:n_s, j, jdim, idim) * conjg(phases(1:n_s))
          else
            ppsi(1:n_s, jdim) = ppsi(1:n_s, jdim) + &
                 rkb_p%f(j, jdim, idim) * uvpsi * rkb_p%ket(1:n_s, j, jdim, idim)
          end if
        end do
      end do
    end do
    
  end subroutine rkb_project

  !------------------------------------------------------------------------------
  ! X(kb_dproject) calculates:
  !  \sum_{i}^3\sum{k}^3 p%e(i) <psi|rkb_p%dp(:, k, i)><rkb_p%p(:, i)|psi>
  !------------------------------------------------------------------------------
  function rkb_dproject(mesh, rkb_p, psi, phases) result(res)
    type(mesh_t),          intent(in) :: mesh
    type(rkb_projector_t), intent(in) :: rkb_p
    CMPLX,                 intent(in) :: psi(:, :) ! psi(rkb%n_s, 2)
    CMPLX, optional,       intent(in) :: phases(:)
    CMPLX :: res(3)

    integer :: n_s, i, k, idim
    CMPLX :: uvpsi
    CMPLX, allocatable :: ppsi(:)
#if defined(HAVE_MPI)
    CMPLX :: tmp
    CMPLX :: tmp2(3)
#endif

    res = M_z0
    n_s = rkb_p%n_s
    ALLOCATE(ppsi(n_s), n_s)

    do idim = 1, 2
      do i = 1, 2
        if (rkb_p%e(i) == M_ZERO) cycle
        uvpsi = sum(psi(1:n_s, idim)*rkb_p%p(1:n_s, i))
#if defined(HAVE_MPI)
        if(mesh%parallel_in_domains) then
          call MPI_Allreduce(uvpsi, tmp, 1, MPI_CMPLX, MPI_SUM, mesh%vp%comm, mpi_err)
          uvpsi = tmp
        end if
#endif
        
        do k = 1, 3
          if (present(phases)) then
            ppsi(1:n_s) = rkb_p%e(i) * uvpsi * rkb_p%dp(1:n_s, k, i) * conjg(phases(1:n_s))
          else
            ppsi(1:n_s) = rkb_p%e(i) * uvpsi * rkb_p%dp(1:n_s, k, i)
          end if
          res(k) = res(k) + sum(conjg(psi(1:n_s, idim)) * ppsi(1:n_s))
        end do
      end do
    end do

    deallocate(ppsi)

#if defined(HAVE_MPI)
    if(mesh%parallel_in_domains) then
      call MPI_Allreduce(res(1), tmp2(1), 3, MPI_CMPLX, MPI_SUM, mesh%vp%comm, mpi_err)
      res = tmp2
    end if
#endif

  end function rkb_dproject

end module rkb_projector_m
