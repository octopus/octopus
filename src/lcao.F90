!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

module lcao
  use global
  use messages
  use lib_oct
  use lib_oct_gsl_spline
  use lib_basic_alg
  use lib_adv_alg
  use functions
  use mesh
  use simul_box
  use specie
  use geometry
  use states
  use system
  use hamiltonian
  use grid

  use output

  implicit none

  private
  public ::          &
    lcao_type,       &
    lcao_init,       &
    lcao_wf,         &
    lcao_initial_wf, &
    lcao_end

  type lcao_type
    integer           :: state ! 0 => non-initialized;
                               ! 1 => initialized (k, s and v1 matrices filled)
    type(states_type) :: st

    R_TYPE,  pointer  :: hamilt(:, :, :) ! hamilt stores the Hamiltonian in the LCAO subspace;
    R_TYPE,  pointer  :: s     (:, :, :) ! s is the overlap matrix;
    R_TYPE,  pointer  :: k     (:, :, :) ! k is the kinetic + spin orbit operator matrix;
    R_TYPE,  pointer  :: v     (:, :, :) ! v is the potential.
  end type lcao_type

contains

  ! ---------------------------------------------------------
  subroutine lcao_initial_wf(n, gr, psi, ispin, ik, err)
    integer,         intent(in)  :: n
    type(grid_type), intent(in)  :: gr
    R_TYPE,          intent(out) :: psi(:, :)
    integer,         intent(in)  :: ispin
    integer,         intent(in)  :: ik
    integer,         intent(out) :: err

    integer :: norbs, ia, i, j, idim, k, wf_dim
    type(specie_type), pointer :: s
    FLOAT :: x(calc_dim), r

    err = 0
    psi = R_TOTYPE(M_ZERO)

    norbs = 0
    do ia = 1, gr%geo%natoms
      norbs = norbs + gr%geo%atom(ia)%spec%niwfs
    end do

    wf_dim = 1
    if( ispin.eq.SPINORS ) then
      wf_dim = 2
      norbs = norbs * 2
    end if

    if( (n>norbs) .or. (n<1) ) then
      err = 1; return
    end if

    i = 1
    do ia = 1, gr%geo%natoms
      s => gr%geo%atom(ia)%spec
      do idim = 1, wf_dim
        do j = 1, s%niwfs
          if(n == i) then
            do k = 1, gr%m%np
              x(:) = gr%m%x(k, :) - gr%geo%atom(ia)%x(:)
              psi(k, idim) =  specie_get_iwf(s, j, calc_dim, states_spin_channel(ispin, ik, idim), x)
            end do
            r = X(states_nrm2)(gr%m, wf_dim, psi)
            psi = psi/r
            return
          end if
          i = i + 1
        end do
      end do
    end do

  end subroutine lcao_initial_wf


  ! ---------------------------------------------------------
  subroutine lcao_init(gr, lcao_data, st, h)
    type(lcao_type),         intent(out)   :: lcao_data
    type(grid_type), target, intent(inout) :: gr
    type(states_type),       intent(in)    :: st
    type(hamiltonian_type),  intent(in)    :: h

    type(geometry_type), pointer :: geo
    integer :: norbs, ik, n1, n2, ia, n, ierr
    R_TYPE, allocatable :: hpsi(:,:)

    if(lcao_data%state == 1) return

    call push_sub('lcao.lcao_init')

    geo => gr%geo
    call states_null(lcao_data%st)

    ! Fix the dimension of the LCAO problem (lcao_data%dim)
    norbs = 0
    do ia = 1, geo%natoms
      norbs = norbs + geo%atom(ia)%spec%niwfs
    end do
    if( (st%d%ispin.eq.SPINORS) ) norbs = norbs * 2
    lcao_data%st%nst = norbs
    if(lcao_data%st%nst < st%nst) then
      lcao_data%state = 0
      write(message(1),'(a)') 'Cannot do LCAO initial calculation because there are not enough atomic orbitals.'
      call write_warning(1)
      nullify(geo)
      call pop_sub()
      return
    end if

    ALLOCATE(lcao_data%st%X(psi)(NP_PART, st%d%dim, norbs, st%d%nik), NP_PART*st%d%dim*norbs*st%d%nik)
    lcao_data%st%X(psi) = R_TOTYPE(M_ZERO)

    do ik = 1, st%d%nik
      do n = 1, lcao_data%st%nst
        call lcao_initial_wf(n, gr, lcao_data%st%X(psi)(:, :, n, ik), st%d%ispin, ik, ierr)
        if(ierr.ne.0) then
          write(message(1),'(a)') 'Internal error in lcao_wf.'
          call write_fatal(1)
        end if
      end do
    end do

    ! Allocation of variables
    ALLOCATE(lcao_data%hamilt(norbs, norbs, st%d%nik), norbs*norbs*st%d%nik)
    ALLOCATE(lcao_data%s     (norbs, norbs, st%d%nik), norbs*norbs*st%d%nik)
    ALLOCATE(lcao_data%k     (norbs, norbs, st%d%nik), norbs*norbs*st%d%nik)
    ALLOCATE(lcao_data%v     (norbs, norbs, st%d%nik), norbs*norbs*st%d%nik)

    ! Overlap and kinetic+so matrices.
    ALLOCATE(hpsi(NP, st%d%dim), NP*st%d%dim)
    do ik = 1, st%d%nik
      do n1 = 1, lcao_data%st%nst
        call X(kinetic) (h, gr, lcao_data%st%X(psi)(:, :, n1, ik), hpsi(:, :), ik)
        ! Relativistic corrections...
        select case(h%reltype)
        case(NOREL)
#if defined(COMPLEX_WFNS) && defined(R_TCOMPLEX)
        case(SPIN_ORBIT)
          call zso (h, gr, lcao_data%st%X(psi)(:, :, n1, ik), hpsi(:, :), ik)
#endif
        end select

        do n2 = n1, lcao_data%st%nst
          lcao_data%k(n1, n2, ik) = X(states_dotp)(gr%m, st%d%dim, hpsi, lcao_data%st%X(psi)(:, : ,n2, ik))
          lcao_data%s(n1, n2, ik) = X(states_dotp)(gr%m, st%d%dim, lcao_data%st%X(psi)(:, :, n1, ik), &
            lcao_data%st%X(psi)(:, : ,n2, ik))
        end do

      end do
    end do
    deallocate(hpsi)
    lcao_data%state = 1

    call pop_sub()
  end subroutine lcao_init


  ! ---------------------------------------------------------
  subroutine lcao_end(lcao_data, nst)
    type(lcao_type), intent(inout) :: lcao_data
    integer,            intent(in) :: nst

    call push_sub('lcao.lcao_end')

    if(lcao_data%st%nst >= nst) then
      if(associated(lcao_data%hamilt)) deallocate(lcao_data%hamilt)
      if(associated(lcao_data%s     )) deallocate(lcao_data%s)
      if(associated(lcao_data%k     )) deallocate(lcao_data%k)
      if(associated(lcao_data%v     )) deallocate(lcao_data%v)
    endif

    if(associated(lcao_data%st%X(psi))) then
      deallocate(lcao_data%st%X(psi)); nullify(lcao_data%st%X(psi))
    end if

    lcao_data%state = 0
    call pop_sub()
  end subroutine lcao_end


  ! ---------------------------------------------------------
  subroutine lcao_wf(lcao_data, st, m, sb, h, start)
    type(lcao_type),        intent(inout) :: lcao_data
    type(states_type),      intent(inout) :: st
    type(mesh_type),        intent(in)    :: m
    type(simul_box_type),   intent(in)    :: sb
    type(hamiltonian_type), intent(in)    :: h
    integer, optional,      intent(in)    :: start

    integer :: np, dim, nst, ik, n1, n2, idim, norbs, start_
    R_TYPE, allocatable :: hpsi(:,:)
    FLOAT, allocatable :: ev(:)

    ASSERT(lcao_data%state == 1)

    call push_sub('lcao.lcao_wf')

    norbs = lcao_data%st%nst
    np = m%np
    dim = st%d%dim
    nst = st%nst

    ! Hamiltonian and overlap matrices.
    ALLOCATE(hpsi(np, dim), np*dim)
    do ik = 1, st%d%nik
      do n1 = 1, lcao_data%st%nst
        hpsi = M_ZERO
        call X(vlpsi) (h, m, lcao_data%st%X(psi)(:,:, n1, ik), hpsi(:,:), ik)
        if (h%ep%nvnl > 0) call X(vnlpsi) (h, m, sb, lcao_data%st%X(psi)(:,:, n1, ik), hpsi(:,:), ik)
        do n2 = n1, lcao_data%st%nst
          lcao_data%v(n1, n2, ik) = X(states_dotp)(m, dim, hpsi, lcao_data%st%X(psi)(1:, : ,n2, ik))
          lcao_data%hamilt(n1, n2, ik) = lcao_data%k(n1, n2, ik) + lcao_data%v(n1 , n2, ik)
          lcao_data%hamilt(n2, n1, ik) = R_CONJ(lcao_data%hamilt(n1, n2, ik))
        end do
      end do
    end do

    do ik = 1, st%d%nik
      ALLOCATE(ev(norbs), norbs)
      call lalg_geneigensolve(norbs, lcao_data%hamilt(1:norbs, 1:norbs, ik), &
        lcao_data%s(1:norbs, 1:norbs, ik), ev)

      start_ = 1
      if(present(start)) start_ = start

      do n1 = start_, nst
        st%eigenval(n1, ik) = ev(n1)
        st%X(psi)(:, :, n1, ik) = R_TOTYPE(M_ZERO)
      end do
      deallocate(ev)

      ! Change of base
      do n1 = start_, nst
        do idim = 1, dim
          do n2 = 1, norbs
            call lalg_axpy(np, lcao_data%hamilt(n2, n1, ik), lcao_data%st%X(psi)(:, idim, n2, ik), &
              st%X(psi)(:, idim, n1, ik))
          end do
        end do
      end do

    end do

    deallocate(hpsi)
    call pop_sub()
  end subroutine lcao_wf

end module lcao
