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

!!! This routine calculates the SIC exchange functional. note that the LDA
!!! part of the functional is already included by the LDA routines
subroutine X(oep_sic) (xcs, m, f_der, st, is, oep, vxc, ex, ec)
  type(xc_type),     intent(in)    :: xcs
  type(mesh_type),   intent(in)    :: m
  type(f_der_type),  intent(inout) :: f_der
  type(states_type), intent(inout) :: st
  integer,           intent(in)    :: is
  type(xc_oep_type), intent(inout) :: oep
  FLOAT,             intent(inout) :: vxc(:,:), ex, ec

  integer  :: i
  FLOAT :: ex2, ec2, edummy
  FLOAT, allocatable :: vxc2(:, :)
  FLOAT, pointer :: rho(:,:), rho_save(:,:)
  
  allocate(rho(m%np, 2), Vxc2(m%np, 2))
  rho(:, 2) = M_ZERO

  ! save real density
  rho_save => st%rho;    st%rho => rho

  ! loop over states
  do i = 1, st%nst
    if(st%occ(i, is) .gt. small) then ! we only need the occupied states
      ! get orbital density
      st%rho(:, 1) = oep%socc*st%occ(i, is)*R_ABS(st%X(psi)(:, 1, i, is))**2
      
      ! initialize before calling get_vxc
      vxc2 = M_ZERO
      ex2  = M_ZERO
      ec2  = M_ZERO
      
      ! calculate LDA/GGA contribution to the SIC (does not work for LB94)
      edummy = M_ZERO
      call xc_get_vxc(2, xcs%sic_aux, xcs%nlcc, m, f_der, st, vxc2, ex2, &
         ec2, edummy, edummy)
      
      ex = ex - oep%sfact*ex2
      ec = ec - oep%sfact*ec2

      oep%lxc(:, i) = oep%lxc(:, i) - &
         vxc2(:, 1)*R_CONJ(st%X(psi) (:, 1, i, is))      
      
      ! calculate the Hartree contribution using poissons equation
      vxc2(:, 1) = M_ZERO
      call poisson_solve(m, f_der, vxc2(:, 1), rho(:, 1))

      ex = ex - M_HALF*oep%sfact*oep%socc*st%occ(i, is)* &
         sum(vxc2(:, 2) * R_ABS(st%X(psi)(:, 1, i, is))**2 * m%vol_pp(:))

      oep%lxc(:, i) = oep%lxc(:, i) - &
         vxc2(:, 1)*R_CONJ(st%X(psi) (:, 1, i, is))
    end if
  end do

  st%rho => rho_save
  deallocate(rho, Vxc2)

end subroutine X(oep_sic)
