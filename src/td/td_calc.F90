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

module td_calc_m
  use c_pointer_m 
  use forces_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_base_m
  use hamiltonian_m
  use lasers_m
  use loct_math_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use profiling_m
  use states_calc_m
  use states_m
  use states_dim_m

  implicit none

  private
  public ::       &
    td_calc_tacc, &
    td_calc_tvel, &
    td_calc_ionch

contains

! ---------------------------------------------------------
!> Electronic acceleration (to calculate harmonic spectrum...)
!! It is calculated as:
!!
!! \f[
!! d2<x>/dt2 = d<p>/dt + i<[H,[V_nl,x]]> =
!!           = i<[V_l,p]> + i<[V_nl,p]> - E(t)N + i<[H,[V_nl,x]]>
!! \f]
!! \warning This subroutine only works if ions are not
!!          allowed to move
! ---------------------------------------------------------
subroutine td_calc_tacc(gr, geo, st, hm, acc, time)
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(inout) :: geo
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  FLOAT,               intent(in)    :: time
  FLOAT,               intent(out)   :: acc(MAX_DIM)

  FLOAT :: field(MAX_DIM), x(MAX_DIM)
  CMPLX, allocatable :: hzpsi(:,:), hhzpsi(:,:), xzpsi(:,:,:), vnl_xzpsi(:,:)
  integer  :: j, k, ik, ist, idim

#if defined(HAVE_MPI)
  FLOAT   :: y(MAX_DIM)
#endif

  PUSH_SUB(td_calc_tacc)

  ! The term i<[V_l,p]> + i<[V_nl,p]> may be considered as equal but opposite to the
  ! force exerted by the electrons on the ions. COMMENT: This has to be thought about.
  ! Maybe we are forgetting something....
  call total_force_calculate(gr, geo, hm%ep, st, acc)

  ! Adds the laser contribution : i<[V_laser, p]>
  ! WARNING: this ignores the possibility of non-electric td external fields.
  field = M_ZERO
  do j = 1, hm%ep%no_lasers
    call laser_electric_field(hm%ep%lasers(j), field(1:gr%sb%dim), time, CNST(0.001))
    acc(1:gr%mesh%sb%dim) = acc(1:gr%mesh%sb%dim) - st%qtot*field(1:gr%mesh%sb%dim)
  end do

  if(.not. hm%ep%non_local) then
    POP_SUB(td_calc_tacc)
    return
  end if

  ! And now, i<[H,[V_nl,x]]>
  x = M_ZERO
  SAFE_ALLOCATE(hzpsi (1:gr%mesh%np_part, 1:st%d%dim))
  SAFE_ALLOCATE(hhzpsi(1:3, 1:gr%mesh%np_part))

  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end

      call zhamiltonian_apply(hm, gr%der, st%zpsi(:, :, ist, ik), hzpsi(:,:), ist, ik, time)

      SAFE_ALLOCATE(xzpsi    (1:gr%mesh%np_part, 1:st%d%dim, 1:3))
      SAFE_ALLOCATE(vnl_xzpsi(1:gr%mesh%np_part, 1:st%d%dim))
      xzpsi = M_z0
      do k = 1, gr%mesh%np
        do j = 1, gr%mesh%sb%dim
          xzpsi(k, 1:st%d%dim, j) = gr%mesh%x(k, j)*st%zpsi(k, 1:st%d%dim, ist, ik)
        end do
      end do

      do j = 1, gr%mesh%sb%dim
        call zhamiltonian_apply(hm, gr%der, xzpsi(:, :, j), vnl_xzpsi, ist, ik, time, terms = TERM_NON_LOCAL_POTENTIAL)

        do idim = 1, st%d%dim
          x(j) = x(j) - 2*st%occ(ist, ik)*real(zmf_dotp(gr%mesh, hzpsi(1:gr%mesh%np, idim), vnl_xzpsi(:, idim)), REAL_PRECISION)
        end do
      end do

      xzpsi = M_z0
      do k = 1, gr%mesh%np
        do j = 1, gr%mesh%sb%dim
          xzpsi(k, 1:st%d%dim, j) = gr%mesh%x(k, j)*hzpsi(k, 1:st%d%dim)
        end do
      end do

      do j = 1, gr%mesh%sb%dim
        call zhamiltonian_apply(hm, gr%der, xzpsi(:, :, j), vnl_xzpsi, ist, ik, time, terms = TERM_NON_LOCAL_POTENTIAL)

        do idim = 1, st%d%dim
          x(j) = x(j) + 2*st%occ(ist, ik)* &
            real(zmf_dotp(gr%mesh, st%zpsi(1:gr%mesh%np, idim, ist, ik), vnl_xzpsi(:, idim)), REAL_PRECISION)
        end do
      end do
      SAFE_DEALLOCATE_A(xzpsi)
      SAFE_DEALLOCATE_A(vnl_xzpsi)

    end do
  end do
  SAFE_DEALLOCATE_A(hzpsi)
  SAFE_DEALLOCATE_A(hhzpsi)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(x(1), y(1), gr%mesh%sb%dim, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    x = y
  end if
#endif
  acc = acc + x

  POP_SUB(td_calc_tacc)
end subroutine td_calc_tacc

! ---------------------------------------------------------
!> Electronic velocity (to calculate harmonic spectrum...)
!! It is calculated as:
!! \f[
!! d<x>/dt = <p>
!! \f]
! ---------------------------------------------------------
subroutine td_calc_tvel(gr, st, vel)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(inout) :: st
  FLOAT,               intent(out)   :: vel(MAX_DIM)

  FLOAT, allocatable :: momentum(:,:,:)
  
  PUSH_SUB(td_calc_tvel)

  SAFE_ALLOCATE(momentum(1:gr%mesh%sb%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
  call states_calc_momentum(st, gr%der, momentum)

  momentum(1:gr%mesh%sb%dim, st%st_start:st%st_end, 1) = & 
    sum(momentum(1:gr%mesh%sb%dim, st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end), 3)
  momentum(1:gr%mesh%sb%dim, 1, 1) = & 
    sum(momentum(1:gr%mesh%sb%dim, st%st_start:st%st_end, 1), 2)
  vel = momentum(1:gr%mesh%sb%dim, 1, 1)

  SAFE_DEALLOCATE_A(momentum)
  POP_SUB(td_calc_tvel)
end subroutine td_calc_tvel

! ---------------------------------------------------------
!> Multiple ionization probabilities calculated form the KS orbital densities 
!! C. Ullrich, Journal of Molecular Structure: THEOCHEM 501, 315 (2000).
! ---------------------------------------------------------
subroutine td_calc_ionch(gr, st, ch, Nch)
  type(grid_t),        intent(in)    :: gr
  type(states_t),      intent(in)    :: st
  FLOAT,               intent(out)   :: ch(0:Nch)
  integer,             intent(in)    :: Nch

  integer :: ik, ist, ii, jj, idim, Nid
  FLOAT   :: prod, prod0
  FLOAT, allocatable :: N(:), Nnot(:)
  CMPLX,  allocatable :: zpsi(:)
  !combinations   
  integer :: next
  type(c_ptr) :: c
  integer, allocatable :: idx0(:), idx(:), idxref(:)
  
#if defined(HAVE_MPI) 
  FLOAT, allocatable :: Nbuf(:)
#endif
  
  
  PUSH_SUB(td_calc_ionch)
  
  SAFE_ALLOCATE(   N(1: Nch)) 
  SAFE_ALLOCATE(Nnot(1: Nch)) 
  SAFE_ALLOCATE(zpsi(1:gr%mesh%np))
  
  N(:)   = M_ZERO
  Nnot(:)= M_ZERO

  
  ii = 1
  do ik = 1, st%d%nik
    do ist = 1, st%nst
      do idim = 1, st%d%dim        

        if (st%st_start <= ist .and. ist <= st%st_end .and. &
              st%d%kpt%start <= ik .and. ik <= st%d%kpt%end) then
          call states_get_state(st, gr%mesh, idim, ist, ik, zpsi)
          N(ii) = zmf_integrate(gr%mesh, zpsi(:) * conjg(zpsi(:)) ) 
          Nnot(ii) = M_ONE - N(ii)
        end if

        ii = ii + 1
        
      end do
    end do
  end do
  
#if defined(HAVE_MPI) 


  if(st%parallel_in_states) then
    SAFE_ALLOCATE(Nbuf(1: Nch)) 
    Nbuf(:) = M_ZERO
    call MPI_Allreduce(N(1), Nbuf(1), Nch, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    N(:) = Nbuf(:)

    Nbuf(:) = M_ZERO
    call MPI_Allreduce(Nnot(1), Nbuf(1), Nch, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    Nnot(:) = Nbuf(:)

    SAFE_DEALLOCATE_A(Nbuf) 
  end if
#endif

! print * ,mpi_world%rank, "N    =", N(:)
! print * ,mpi_world%rank, "Nnot =", Nnot(:)

  ch(:) = M_ZERO

  Nid = Nch - 1 
  SAFE_ALLOCATE(idx(0:Nid))
  SAFE_ALLOCATE(idxref(0:Nid))
  idxref = (/ (ii, ii = 0, Nid) /)
  do ii = 0, Nch  
!     print *, "Nch= ", Nch, "ii", ii
    call loct_combination_init(c, Nch, ii)
    if (ii == 0 ) then
      SAFE_ALLOCATE(idx0(0:0))
    else 
      SAFE_ALLOCATE(idx0(0:ii-1))
    end if 
!     print *,"size(idx0,1)=",size(idx0,1)
    if(in_debug_mode) then
      call messages_write("P(")
      call messages_write(ii)
      call messages_write(") = 0")
      call messages_new_line()
    end if
    ! Loop over combinations   
    do
      prod  = M_ONE
      prod0 = M_ONE
      if(in_debug_mode) then
        call messages_write("P(")
        call messages_write(ii)
        call messages_write(") += ")
      end if               
      call loct_get_combination(c, idx0)
      idx(:) = idxref(:)
      do jj = 0, ii-1
        idx(idx0(jj))= -1
        if(in_debug_mode) then
          call messages_write(" No(")
          call messages_write(idx0(jj)+1)
          call messages_write(") ")
        end if
        prod0 = prod0 * Nnot(idx0(jj)+1)
      end do

      do jj = 0, Nid
        if(idx(jj)>=0) then
          if(in_debug_mode) then
            call messages_write(" N(")
            call messages_write(idx(jj)+1)
            call messages_write(") ")
          end if
          prod = prod * N(idx(jj)+1)
        end if
      end do
      
      if(in_debug_mode) call messages_new_line()

      ch(ii) = ch(ii) + prod*prod0

      call loct_combination_next(c, next)
      if( next /=  0) exit
    end do
    SAFE_DEALLOCATE_A(idx0)
    call loct_combination_end(c)

    if(in_debug_mode) call messages_info()
  end do 

  SAFE_DEALLOCATE_A(idx)
  SAFE_DEALLOCATE_A(idxref)

  
  SAFE_DEALLOCATE_A(N)
  SAFE_DEALLOCATE_A(Nnot)
  SAFE_DEALLOCATE_A(zpsi)

  POP_SUB(td_calc_ionch)
end subroutine td_calc_ionch


end module td_calc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
