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
!! $Id$


! ---------------------------------------------------------
! Electronic acceleration (to calculate harmonic spectrum...)
! It is calculated as:
!
! d2<x>/dt2 = d<p>/dt + i<[H,[V_nl,x]]> =
!           = i<[V_l,p]> + i<[V_nl,p]> - E(t)N + i<[H,[V_nl,x]]>
!
! WARNING: This subroutine only works if ions are not
!          allowed to move
! ---------------------------------------------------------
subroutine td_calc_tacc(gr, geo, st, h, acc, t)
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(inout) :: geo
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: h
  FLOAT,               intent(in)    :: t
  FLOAT,               intent(out)   :: acc(MAX_DIM)

  FLOAT :: field(MAX_DIM), x(MAX_DIM)
  CMPLX, allocatable :: hzpsi(:,:), hhzpsi(:,:), xzpsi(:,:,:), vnl_xzpsi(:,:), conj(:)
  integer  :: j, k, i, ik, ist, idim

#if defined(HAVE_MPI)
  FLOAT   :: y(MAX_DIM)
#endif

  call push_sub('td_calc.td_calc_tacc')

  ! The term i<[V_l,p]> + i<[V_nl,p]> may be considered as equal but opposite to the
  ! force exerted by the electrons on the ions. COMMENT: This has to be thought about.
  ! Maybe we are forgetting something....
  x = M_ZERO
  call epot_forces(gr, geo, h%ep, st)
  do i = 1, geo%natoms
    x = x - geo%atom(i)%f
  end do
  acc = x

  ! Adds the laser contribution : i<[V_laser, p]>
  if(h%ep%no_lasers > 0) then
    call epot_laser_field(gr%sb, h%ep, t, field)
    acc(1:3) = acc(1:3) - st%qtot*field(1:3)
  end if

  if(h%ep%nvnl <= 0) then
    call pop_sub()
    return
  end if

  ! And now, i<[H,[V_nl,x]]>
  x = M_ZERO
  ALLOCATE(hzpsi (NP, st%d%dim), NP*st%d%dim)
  ALLOCATE(hhzpsi(3, NP),        3*NP)
  ALLOCATE(conj  (NP),           NP)

  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end

      call zhpsi(h, gr, st%zpsi(:, :, ist, ik), hzpsi(:,:), ik)
      do k = 1, st%d%dim
        hzpsi(:, k) = hzpsi(:, k) + epot_laser_scalar_pot(gr%m%np, gr, h%ep, t)*st%zpsi(:, k, ist, ik)
      end do

      ALLOCATE(xzpsi    (NP, st%d%dim, 3), NP*st%d%dim*3)
      ALLOCATE(vnl_xzpsi(NP, st%d%dim),    NP*st%d%dim)
      xzpsi = M_z0
      do k = 1, NP
        do j = 1, NDIM
          xzpsi(k, 1:st%d%dim, j) = gr%m%x(j, k) * st%zpsi(k, 1:st%d%dim, ist, ik)
        end do
      end do

      do j = 1, NDIM
        vnl_xzpsi = M_z0
        call zvnlpsi(h, gr%m, xzpsi(:,:, j), vnl_xzpsi(:,:), ik)

        do idim = 1, st%d%dim
          conj = conjg(hzpsi(:, idim))
          x(j) = x(j) - 2*st%occ(ist, ik)*zmf_dotp(gr%m, conj, vnl_xzpsi(:, idim) )
        end do
      end do

      xzpsi = M_z0
      do k = 1, NP
        do j = 1, NDIM
          xzpsi(k, 1:st%d%dim, j) = gr%m%x(j, k) * hzpsi(k, 1:st%d%dim)
        end do
      end do

      do j = 1, NDIM
        vnl_xzpsi = M_z0
        call zvnlpsi(h, gr%m, xzpsi(:,:, j), vnl_xzpsi(:,:), ik)
        do idim = 1, st%d%dim
          conj = conjg(st%zpsi(:, idim, ist, ik))
          x(j) = x(j) + 2*st%occ(ist, ik)* &
            zmf_dotp(gr%m, conj, vnl_xzpsi(:, idim) )
        end do
      end do
      deallocate(xzpsi, vnl_xzpsi)

    end do
  end do
  deallocate(hzpsi, hhzpsi, conj)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(x(1), y(1), NDIM, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    x = y
  end if
#endif
  acc = acc + x

  call pop_sub()
end subroutine td_calc_tacc
