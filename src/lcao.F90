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

module lcao_m
  use global_m
  use messages_m
  use datasets_m
  use lib_oct_parser_m
  use lib_oct_m
  use lib_oct_gsl_spline_m
  use lib_basic_alg_m
  use lib_adv_alg_m
  use functions_m
  use mesh_m
  use simul_box_m
  use specie_m
  use geometry_m
  use states_m
  use system_m
  use hamiltonian_m
  use grid_m

  use output_m

  implicit none

  private
  public ::          &
    lcao_t,          &
    lcao_init,       &
    lcao_wf,         &
    lcao_initial_wf, &
    lcao_end

  type lcao_t
    integer           :: state ! 0 => non-initialized;
                               ! 1 => initialized (k, s and v1 matrices filled)
    type(states_t) :: st

    R_TYPE,  pointer  :: hamilt(:, :, :) ! hamilt stores the Hamiltonian in the LCAO subspace;
    R_TYPE,  pointer  :: s     (:, :, :) ! s is the overlap matrix;
    R_TYPE,  pointer  :: k     (:, :, :) ! k is the kinetic + spin orbit operator matrix;
    R_TYPE,  pointer  :: v     (:, :, :) ! v is the potential.
  end type lcao_t

contains

  ! ---------------------------------------------------------
  subroutine lcao_initial_wf(n, gr, psi, ispin, ik, err)
    integer,         intent(in)  :: n
    type(grid_t),    intent(in)  :: gr
    R_TYPE,          intent(out) :: psi(:, :)
    integer,         intent(in)  :: ispin
    integer,         intent(in)  :: ik
    integer,         intent(out) :: err

    integer :: norbs, ia, i, j, idim, k, wf_dim
    type(specie_t), pointer :: s
    FLOAT :: x(MAX_DIM), r

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

    idim = 1
    i = 1; j = 0
    do
      j = j + 1
      do ia = 1, gr%geo%natoms
        s => gr%geo%atom(ia)%spec
        do idim = 1, wf_dim
          if(j > s%niwfs) cycle
          if(n == i) then
            do k = 1, gr%m%np
              x(1:calc_dim) = gr%m%x(k, 1:calc_dim) - gr%geo%atom(ia)%x(1:calc_dim)
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
    type(lcao_t),         intent(out)   :: lcao_data
    type(grid_t), target, intent(inout) :: gr
    type(states_t),       intent(in)    :: st
    type(hamiltonian_t),  intent(in)    :: h

    type(geometry_t), pointer :: geo
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

    if(norbs < st%nst) then
      lcao_data%state = 0
      write(message(1),'(a)') 'Cannot do LCAO initial calculation because there are not enough atomic orbitals.'
      call write_warning(1)
      nullify(geo)
      call pop_sub()
      return
    end if

    !%Variable LCAODimension
    !%Type integer
    !%Default 0
    !%Section SCF
    !%Description
    !% Before starting the SCF cycle, an initial LCAO calculation can be performed
    !% in order to obtain reasonable initial guesses for spin-orbitals and densities.
    !% For this purpose, the code calculates a number of atomic orbitals -- this
    !% number depends on the given species. The default dimension for the LCAO basis
    !% set will be the sum of all these numbers, unless this dimension is larger than
    !% twice the number of required orbitlas for the full calculation. 
    !%
    !% This dimension however can be reduced (never increased) by making use of the 
    !% variable LCAODimension. Note that LCAODimension cannot be smaller than the 
    !% number of orbitals needed in the full calculation -- if LCAODimension is smaller, 
    !% it will be changed silently increased to meet this requirement. In the same way, 
    !% if LCAODimension is larger than the available number of atomic orbitals, 
    !% it will be reduced. If you want to use the largest possible number, set
    !% LCAODimension to a negative number.
    !%End
    call loct_parse_int(check_inp('LCAODimension'), 0, n)
    if((n > 0) .and. (n <= st%nst)) then
      norbs = st%nst
    elseif( (n > st%nst) .and. (n < norbs) ) then
      norbs = n
    elseif( n.eq.0) then
      norbs = min(norbs, 2*st%nst)
    end if

    lcao_data%st%nst = norbs
    lcao_data%st%st_start = 1
    lcao_data%st%st_end = norbs
    lcao_data%st%d%dim = st%d%dim
    lcao_data%st%d%nik = st%d%nik
    call X(states_allocate_wfns)(lcao_data%st, gr%m)

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
    type(lcao_t), intent(inout) :: lcao_data
    integer,         intent(in) :: nst

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
    type(lcao_t),        intent(inout) :: lcao_data
    type(states_t),      intent(inout) :: st
    type(mesh_t),        intent(in)    :: m
    type(simul_box_t),   intent(in)    :: sb
    type(hamiltonian_t), intent(in)    :: h
    integer, optional,   intent(in)    :: start

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

end module lcao_m
