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

!!! This routine calculates the SIC exchange functional.
subroutine X(oep_sic) (xcs, gr, st, is, oep, ex, ec)
  type(xc_type),     intent(in)    :: xcs
  type(grid_type),   intent(inout) :: gr
  type(states_type), intent(inout) :: st
  integer,           intent(in)    :: is
  type(xc_oep_type), intent(inout) :: oep
  FLOAT,             intent(inout) :: ex, ec

  integer  :: i
  FLOAT :: ex2, ec2, ex_, ec_, edummy
  FLOAT, allocatable :: vxc(:, :), rho(:,:)
#if defined(HAVE_MPI)
  integer :: ierr
#endif


  call push_sub('oep_sic')
  
  allocate(rho(NP, 2), Vxc(NP, 2))
  rho(1:NP, 2) = M_ZERO

  ! loop over states
  ex_ = M_ZERO
  ec_ = M_ZERO
  do i = st%st_start, st%st_end
    if(st%occ(i, is) .gt. small) then ! we only need the occupied states
      ! get orbital density
      rho(1:NP, 1) = oep%socc*st%occ(i, is)*R_ABS(st%X(psi)(1:NP, 1, i, is))**2
     
      ! initialize before calling get_vxc
      vxc = M_ZERO
      ex2  = M_ZERO
      ec2  = M_ZERO
      
      ! calculate LDA/GGA contribution to the SIC (does not work for LB94)
      edummy = M_ZERO
      call xc_get_vxc(gr, xcs, rho, SPIN_POLARIZED, vxc, ex2, ec2, edummy, edummy) 

      ex_ = ex_ - oep%sfact*ex2
      ec_ = ec_ - oep%sfact*ec2

      oep%X(lxc)(1:NP, i) = oep%X(lxc)(1:NP, i) - &
         vxc(1:NP, 1)*R_CONJ(st%X(psi) (1:NP, 1, i, is))      
      
      ! calculate the Hartree contribution using poissons equation
      vxc(1:NP, 1) = M_ZERO
      call dpoisson_solve(gr, vxc(:, 1), rho(:, 1))

      ex_ = ex_ - M_HALF*oep%sfact*oep%socc*st%occ(i, is)* &
         sum(vxc(1:NP, 1) * R_ABS(st%X(psi)(1:NP, 1, i, is))**2 * gr%m%vol_pp(1:NP))

      oep%X(lxc)(1:NP, i) = oep%X(lxc)(1:NP, i) - vxc(1:NP, 1)*R_CONJ(st%X(psi) (1:NP, 1, i, is))
    end if
  end do

#if defined(HAVE_MPI)
  if(st%st_end - st%st_start + 1 .ne. st%nst) then ! This holds only in the td part.
    call MPI_ALLREDUCE(ec_, edummy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr); ec_ = edummy
    call MPI_ALLREDUCE(ex_, edummy, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr); ex_ = edummy
  endif
#endif

  ec = ec + ec_
  ex = ex + ex_

  deallocate(rho, Vxc)
  call pop_sub()
end subroutine X(oep_sic)
