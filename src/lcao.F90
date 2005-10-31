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
  public :: lcao_type, &
            lcao_init, &
            lcao_wf, &
            lcao_end

  type lcao_type
    integer           :: state ! 0 => non-initialized;
                               ! 1 => initialized (k, s and v1 matrices filled)
    type(states_type) :: st

    R_TYPE,  pointer :: hamilt(:, :, :) ! hamilt stores the Hamiltonian in the LCAO subspace;
    R_TYPE,  pointer :: s     (:, :, :) ! s is the overlap matrix;
    R_TYPE,  pointer :: k     (:, :, :) ! k is the kinetic + spin orbit operator matrix;
    R_TYPE,  pointer :: v     (:, :, :) ! v is the potential.
  end type lcao_type

contains

subroutine lcao_init(gr, lcao_data, st, h)
  type(lcao_type),         intent(out)   :: lcao_data
  type(grid_type), target, intent(inout) :: gr
  type(states_type),       intent(in)    :: st
  type(hamiltonian_type),  intent(in)    :: h

  type(geometry_type), pointer :: geo
  type(specie_type), pointer :: s
  integer :: ierr, norbs, ispin, ik, n1, i1, l, l1, lm1, d1, n2, i, j, ia, n, idim, is, k
  FLOAT :: x(gr%sb%dim), r
  R_TYPE, allocatable :: hpsi(:,:)

  if(lcao_data%state == 1) return

  call push_sub('lcao.lcao_init')

  geo => gr%geo
  call states_null(lcao_data%st)

  ! Fix the dimension of the LCAO problem (lcao_data%dim)
  norbs = 0
  do ia = 1, geo%natoms
     norbs = norbs + geo%atom(ia)%spec%niwfs
  enddo
  if( (st%d%ispin.eq.SPINORS) ) norbs = norbs * 2
  lcao_data%st%nst = norbs
  if(lcao_data%st%nst < st%nst) then
    lcao_data%state = 0
    write(message(1),'(a)') 'Cannot do LCAO initial calculation because there are not enough atomic orbitals.'
    call write_warning(1)
    nullify(geo)
    call pop_sub()
    return
  endif

  allocate(lcao_data%st%X(psi)(gr%m%np, st%d%dim, norbs, st%d%nik))
  lcao_data%st%X(psi) = R_TOTYPE(M_ZERO)

  do ik = 1, st%d%nik
     n = 1
     do ia = 1, geo%natoms
        s => geo%atom(ia)%spec
        do idim = 1, st%d%dim
           do j = 1, s%niwfs
              do k = 1, gr%m%np
                 x(:) = gr%m%x(k, :) - geo%atom(ia)%x(:)
                 lcao_data%st%X(psi)(k, idim, n, ik) = specie_get_iwf(s, j, gr%sb%dim, &
                                                                 states_spin_channel(st%d%ispin, ik, idim), x)
              enddo
              n = n + 1
           enddo
        enddo
     enddo
  enddo

  !Normalize.
  do ik = 1, st%d%nik
     do n = 1, lcao_data%st%nst
        r = X(states_nrm2)(gr%m, st%d%dim, lcao_data%st%X(psi)(:, :, n, ik))
        lcao_data%st%X(psi)(:, :, n, ik) = lcao_data%st%X(psi)(:, :, n, ik)/r
     enddo
  enddo

  ! Allocation of variables
  allocate(lcao_data%hamilt (norbs, norbs, st%d%nik), &
           lcao_data%s      (norbs, norbs, st%d%nik), &
           lcao_data%k      (norbs, norbs, st%d%nik), &
           lcao_data%v      (norbs, norbs, st%d%nik))

  ! Overlap and kinetic+so matrices.
  allocate(hpsi(gr%m%np, st%d%dim))
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

subroutine lcao_end(lcao_data)
  type(lcao_type), intent(inout) :: lcao_data
  call push_sub('lcao.lcao_end')

  ASSERT(lcao_data%state == 1)

  if(associated(lcao_data%hamilt)) then
    deallocate(lcao_data%hamilt, lcao_data%s, lcao_data%k, lcao_data%v)
  endif

  call states_end(lcao_data%st)

  lcao_data%state = 0
  call pop_sub()
end subroutine lcao_end

subroutine lcao_wf(lcao_data, m, sb, st, h)
  type(lcao_type),        intent(inout) :: lcao_data
  type(mesh_type),        intent(in)    :: m
  type(simul_box_type),   intent(in)    :: sb
  type(states_type),      intent(inout) :: st
  type(hamiltonian_type), intent(in)    :: h

  integer :: np, dim, nst, ik, n1, n2, idim, norbs
  R_TYPE, allocatable :: hpsi(:,:)
  FLOAT, allocatable :: ev(:)

  ASSERT(lcao_data%state == 1)

  call push_sub('lcao.lcao_wf')

  norbs = lcao_data%st%nst
  np = m%np
  dim = st%d%dim
  nst = st%nst

  ! Hamiltonian and overlap matrices.
  allocate(hpsi(np, dim))
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
    allocate(ev(norbs))
    call lalg_geneigensolve(norbs, lcao_data%hamilt(1:norbs, 1:norbs, ik), &
         lcao_data%s(1:norbs, 1:norbs, ik), ev)

    st%eigenval(1:nst, ik) = ev(1:nst)
    deallocate(ev)

    st%X(psi)(:,:,:, ik) = R_TOTYPE(M_ZERO)

    ! Change of base
    do n1 = 1, nst
       do idim = 1, dim
          do n2 = 1, norbs
             call lalg_axpy(np, lcao_data%hamilt(n2, n1, ik), lcao_data%st%X(psi)(:, idim, n2, ik), &
                            st%X(psi)(:, idim, n1, ik))
          enddo
       enddo
    enddo

   end do

  deallocate(hpsi)
  call pop_sub()
end subroutine lcao_wf

end module lcao






