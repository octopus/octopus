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

! Calculates the new density out the wavefunctions and occupations...
subroutine R_FUNC(calcdens)(st, np, rho, reduce)
  type(states_type), intent(inout) :: st
  integer, intent(in) :: np
  real(r8), intent(out) :: rho(np, st%nspin)
  logical, intent(in), optional :: reduce

  integer :: i, ik, p, sp
#ifdef HAVE_MPI
  real(r8), allocatable :: reduce_rho(:,:)
  R_TYPE, allocatable :: reduce_rho_off(:)
  integer :: ierr
#endif

  sub_name = 'calc_dens'; call push_sub()

  if(st%ispin == 2) then
    sp = 2
  else
    sp = 1
  end if

  rho = 0._r8
  if(st%ispin == 3) st%R_FUNC(rho_off) = 0._r8
  do ik = 1, st%nik, sp
    do p  = st%st_start, st%st_end
      do i = 1, np
        ! spin-up density
        rho(i, 1) = rho(i, 1) + st%kweights(ik)*st%occ(p, ik)&
             * R_ABS(st%R_FUNC(psi)(i, 1, p, ik))**2

        ! spin-down density
        if(st%ispin == 2) then
          rho(i, 2) = rho(i, 2) + st%kweights(ik+1)*st%occ(p, ik+1) &
               * R_ABS(st%R_FUNC(psi)(i, 1, p, ik+1))**2
        end if

        ! off-diagonal densities
        if(st%ispin == 3) then
          rho(i, 2) = rho(i, 2) + st%kweights(ik)*st%occ(p, ik) &
               * R_ABS(st%R_FUNC(psi)(i, 2, p, ik))**2
          st%R_FUNC(rho_off)(i) = st%R_FUNC(rho_off)(i) + st%kweights(ik)*st%occ(p, ik) &
               * st%R_FUNC(psi)(i, 1, p, ik)*R_CONJ(st%R_FUNC(psi)(i, 2, p, ik))
        end if

      end do
    end do
  end do

#if defined(HAVE_MPI) && defined(MPI_TD)
  ! reduce density (assumes memory is contiguous)
  if(present(reduce) .and. reduce) then
    allocate(reduce_rho(1:np, st%nspin))
    call MPI_ALLREDUCE(rho(1, 1), reduce_rho(1, 1), np*st%nspin, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    rho = reduce_rho
    deallocate(reduce_rho)

    ! update also off diagonal...
    if(st%ispin == 3) then
      allocate(reduce_rho_off(1:np))
      call MPI_ALLREDUCE(st%R_FUNC(rho_off) (1), reduce_rho_off(1), np, &
           R_MPITYPE, MPI_SUM, MPI_COMM_WORLD, ierr)
      st%R_FUNC(rho_off) = reduce_rho_off
      deallocate(reduce_rho_off)
    end if
  end if
#endif

  call pop_sub()
  return
end subroutine R_FUNC(calcdens)
