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

! Please use this functional only in spin polarized calculations
! and for finite systems
subroutine R_FUNC(kli_hju) (m, st, hartr, type, Vx, ex)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(inout) :: st
  type(hartree_type), intent(inout) :: hartr
  integer, intent(in) :: type
  real(r8), intent(out) :: Vx(m%np, st%nspin), ex

  integer :: ia, is, i, k, i1, i2
  real(r8) :: ex2, N_alpha

  real(r8), allocatable :: rho2(:,:), Vx2(:, :)
  real(r8), allocatable, target :: rho(:,:)
  real(r8), pointer :: p(:,:)

  allocate(rho(m%np, 2), rho2(m%np, 2), Vx2(m%np, 2))

  ! first we get the true density
  do is = 1, st%nspin
     do k = 1, m%np
        !rho(k, is) = sum(st%occ(:, is)*R_ABS(st%R_FUNC(psi) (k, 1, :, is))**2)
       rho(k, is) = st%rho(k, is)
     end do
  end do

  call R_FUNC(xc_lda) (type, m, st, Vx, ex)

  ! now store some stuff
  i1 = st%ispin; st%ispin = 2
  i2 = st%nspin; st%nspin = 2
  p  =>st%rho;   st%rho   =>rho

  rho2(:, 2) = 0._r8
  do is = 1, st%nspin
    do ia = 1, 2 ! loop over atoms (only dimers! ;)

      rho2(:, 1) = 0._r8
      do k = 1, m%np
        if((ia == 1.and.m%Lz(k)>0).or.(ia == 2.and.m%Lz(k)<0)) then
          rho2(k, 1) = rho(k, is)
        end if
        if(m%Lz(k) == 0) rho2(k, 1) = 0.5_r8*rho(k, is)
      end do
      N_alpha = dmesh_integrate(m, rho2(:, 1))

      ! first the lda term
      Vx2 = 0.0_r8; Ex2 = 0.0_r8
      call R_FUNC(xc_lda) (type, m, st, Vx2, ex2)
      
      ! The LDA is local, so we can just add the whole array
      Vx(1:m%np, is) = Vx(1:m%np, is) - N_alpha*Vx2(1:m%np, 1)
      Ex = Ex - N_alpha*ex2
      
      ! should we add the SI Hartree correction?
      if(type == X_FUNC_LDA_NREL) then
        call hartree_solve(hartr, m, Vx2(:, 1), rho(:, 1:1))
        do k = 1, m%np
          if((ia == 1.and.m%Lz(k)>=0).or.(ia == 2.and.m%Lz(k)<=0)) then
            Vx(k, is) = Vx(k, is) - N_alpha*Vx2(k, 1)
          end if
        end do
        
        Ex = Ex - 0.5_r8*N_alpha*sum(Vx2(1:m%np, 1)*rho2(1:m%np, 1))*m%vol_pp
      end if

    end do
  end do

  ! restore state variables
  st%ispin = i1; st%nspin = i2
  st%rho => p

  ! deallocate the rest of the vars..
  deallocate(rho, rho2, Vx2)

  return
end subroutine R_FUNC(kli_hju)
