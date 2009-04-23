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
!! $Id$


! ---------------------------------------------------------
subroutine X(one_body) (gr, geo, st, hm)
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in) :: hm

  integer i, j, iunit, idir, iatom, np
  R_TYPE :: me, exp_r, exp_g, corr
  R_TYPE, allocatable :: gpsi(:,:), cpsi(:,:)

  np = NP

  call io_assign(iunit)
  iunit = io_open('matrix_elements/1-body', action='write')

  do i = 1, st%nst
    do j = 1, st%nst
      if(j > i) cycle
      
      me = st%eigenval(i,1) - X(mf_integrate) (gr%mesh, R_CONJ(st%X(psi) (1:np, 1, i, 1)) * &
           hm%Vhxc(1:np, 1) * st%X(psi) (1:np, 1, j, 1))

      write(iunit, *) i, j, me
    end do
  end do

  SAFE_ALLOCATE(gpsi(1:gr%mesh%np_part, 1:MAX_DIM))
  SAFE_ALLOCATE(cpsi(1:gr%mesh%np_part, 1:1))

  
  call io_assign(iunit)
  iunit = io_open('matrix_elements/gauge', action='write')

  do i = 1, st%nst
    do j = 1, st%nst
      if(st%occ(i, 1) < CNST(0.0001)) cycle
      if(st%occ(j, 1) > CNST(0.0001)) cycle

      call X(derivatives_grad)(gr%der, st%X(psi)(:, 1, j, 1), gpsi)
       
      do idir = 1, 3
         exp_r = X(mf_integrate) (gr%mesh, R_CONJ(st%X(psi) (1:np, 1, i, 1)) * &
              gr%mesh%x(1:np, idir) * st%X(psi) (1:np, 1, j, 1))
         
         exp_g = X(mf_integrate) (gr%mesh, R_CONJ(st%X(psi) (1:np, 1, i, 1)) * &
              gpsi(1:np, idir))
         
         corr = M_ZERO
         do iatom = 1, geo%natoms
           call X(projector_commute_r)(hm%ep%proj(iatom), gr, 1, idir, 1, st%X(psi)(1:np, 1, j, 1), cpsi)
           corr = corr + &
                X(mf_integrate)(gr%mesh, R_CONJ(st%X(psi)(1:np, 1, i, 1)) * cpsi(1:np, 1))
         end do

         me = (st%eigenval(j,1) - st%eigenval(i,1)) * exp_r
         
         write(iunit, *) i, j, idir, me, me - (exp_g + corr)

       end do
       
     end do
   end do
   
  SAFE_DEALLOCATE_A(gpsi)

  call io_close(iunit)
end subroutine X(one_body)


! ---------------------------------------------------------
subroutine X(two_body) (gr, st)
  type(grid_t),        intent(inout) :: gr
  type(states_t),      intent(in)    :: st

  integer i, j, k, l, iunit
  R_TYPE :: me
  R_TYPE, allocatable :: n(:), v(:)

  call io_assign(iunit)
  iunit = io_open('matrix_elements/2-body', action='write')

  SAFE_ALLOCATE(n(1:gr%mesh%np))
  SAFE_ALLOCATE(v(1:gr%mesh%np))

  do i = 1, st%nst
    do j = 1, st%nst
      if(j > i) cycle

      n(1:gr%mesh%np) = R_CONJ(st%X(psi) (1:gr%mesh%np, 1, i, 1)) * st%X(psi) (1:gr%mesh%np, 1, j, 1)
      call X(poisson_solve) (gr, v, n, all_nodes=.false.)

      do k = 1, st%nst
        if(k > i) cycle
        do l = 1, st%nst
          if(l > k) cycle
          if(l > j) cycle

          me = X(mf_integrate) (gr%mesh, v(1:gr%mesh%np) * &
                st%X(psi) (1:gr%mesh%np, 1, k, 1) *R_CONJ(st%X(psi) (1:gr%mesh%np, 1, l, 1)))

          write(iunit, *) i, j, k, l, me
        end do
      end do
    end do
  end do

  call io_close(iunit)
end subroutine X(two_body)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
