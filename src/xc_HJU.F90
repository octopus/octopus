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
subroutine R_FUNC(kli_hju) (m, st, hartr, nlcc, type, Vx, ex)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(inout) :: st
  type(hartree_type), intent(inout) :: hartr
  logical, intent(in) :: nlcc
  integer, intent(in) :: type
  real(r8), intent(out) :: Vx(m%np, st%nspin), ex

  integer :: ia, is, i, k
  real(r8) :: ex2, N_alpha

  real(r8), allocatable :: rho(:,:), Vx2(:, :)

  allocate(rho(m%np, 2), Vx2(m%np, 2))

  ! save rho
  rho = st%rho

  ex = 0._r8
  call R_FUNC(xc_lda) (type, nlcc, m, st, Vx, ex)

  st%rho(:, 2) = 0._r8
  do is = 1, st%nspin
    do ia = 1, 2 ! loop over atoms (only dimers! ;)

      st%rho(:, 1) = 0._r8
      do k = 1, m%np
        if((ia == 1.and.m%Lxyz(3, k)>0).or.(ia == 2.and.m%Lxyz(3, k)<0)) then
          st%rho(k, 1) = rho(k, is)
        end if
        if(m%Lxyz(3, k) == 0) st%rho(k, 1) = 0.5_r8*rho(k, is)
      end do
      N_alpha = dmesh_integrate(m, st%rho(:, 1))

      ! first the lda term
      Vx2 = 0.0_r8; Ex2 = 0.0_r8
      call R_FUNC(xc_lda) (type, nlcc, m, st, Vx2, Ex2)
      
      ! The LDA is local, so we can just add the whole array
      Vx(1:m%np, is) = Vx(1:m%np, is) - N_alpha*Vx2(1:m%np, 1)
      Ex = Ex - N_alpha*Ex2
      
      ! should we add the SI Hartree correction?
      if(type == X_FUNC_LDA_NREL) then
        write(message(1), '(a,i2,a,f14.6)') 'N(', is, ') = ', N_alpha
        call write_info(1)

        call hartree_solve(hartr, m, Vx2(:, 1), st%rho(:, 1:1))
        Ex2 = 0._r8
        do k = 1, m%np
          if((ia == 1.and.m%Lxyz(3, k)>=0).or.(ia == 2.and.m%Lxyz(3, k)<=0)) then
            Vx(k, is) = Vx(k, is) - N_alpha*Vx2(k, 1)
            Ex2 = Ex2 + Vx2(k, 1)*st%rho(k, 1)
          end if
        end do
        
        Ex = Ex - 0.5_r8*N_alpha*Ex2*m%vol_pp
      end if

    end do
  end do

  ! restore rho
  st%rho = rho

  ! deallocate the rest of the vars..
  deallocate(rho, Vx2)

  return
end subroutine R_FUNC(kli_hju)
