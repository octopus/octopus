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
!! $Id: cpmd_inc.F90 3694 2008-02-15 13:37:54Z marques $

! This integration is based on Tuckermann and Parrinello, JCP 101
! 1302 (1994), using the Verlet and velocity Verlet algorithms
! described on page 1306.
subroutine X(cpmd_propagate)(this, gr, hm, st, iter, dt)
  type(cpmd_t), target, intent(inout) :: this
  type(grid_t),         intent(inout) :: gr
  type(hamiltonian_t),  intent(inout) :: hm
  type(states_t),       intent(inout) :: st
  integer,              intent(in)    :: iter
  FLOAT,                intent(in)    :: dt

  integer :: ik, ist1, ddim, idim, np, ip
  R_TYPE, allocatable :: hpsi(:, :), psi(:, :), xx(:, :)
  R_TYPE, pointer     :: oldpsi(:, :, :)
  R_TYPE :: one 
  type(batch_t) :: oldpsib

  PUSH_SUB(X(cpmd_propagate))

  one = R_TOTYPE(M_ONE)

  call profiling_in(cpmd_prop, "CP_PROPAGATION")

  np = gr%mesh%np
  ddim = st%d%dim

  SAFE_ALLOCATE(xx(1:st%nst, 1:st%nst))
  SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psi(1:gr%mesh%np, 1:st%d%dim))

  this%ecorr = M_ZERO

  do ik = 1, st%d%nik

    select case(this%method)

    case(VERLET)

      oldpsi => this%X(psi2)(:, :, :, ik)

      do ist1 = st%st_start, st%st_end

        ! give the initial conditions
        if(iter == 1) this%X(psi2)(1:np, 1:ddim, ist1, ik) = st%X(psi)(1:np, 1:ddim, ist1, ik)

        ! calculate the "force"
        call X(hamiltonian_apply)(hm, gr%der, st%X(psi)(:, :, ist1, ik), hpsi, ist1, ik)

        ! propagate the wavefunctions
        psi(1:np, 1:ddim) = M_TWO*st%X(psi)(1:np, 1:ddim, ist1, ik) - this%X(psi2)(1:np, 1:ddim, ist1, ik) &
          + dt**2/this%emass*(-st%occ(ist1, ik)*hpsi(1:np, 1:ddim)) !(4.2)

        ! calculate the velocity and the fictitious electronic energy
        hpsi(1:np, 1:ddim) = abs(psi(1:np, 1:ddim) - this%X(psi2)(1:np, 1:ddim, ist1, ik))/(M_TWO*dt) !(4.7)
        this%ecorr = this%ecorr + this%emass*X(mf_nrm2)(gr%mesh, ddim, hpsi)**2 !(2.11)

        do idim = 1, ddim
          ! store the old wavefunctions
          call lalg_copy(np, st%X(psi)(:, idim, ist1, ik), this%X(psi2)(:, idim, ist1, ik))
          call lalg_copy(np, psi(:, idim), st%X(psi)(:, idim, ist1, ik))
        end do

      end do

      call profiling_in(cpmd_orth, "CP_ORTHOGONALIZATION")

      call calc_xx()

      ! psi <= psi + X * psi2
      call states_block_matr_mul_add(gr%mesh, st, one, st%st_start, st%st_end, st%st_start, st%st_end, &
        this%X(psi2)(:, :, :, ik), xx, one, st%X(psi)(:, :, :, ik)) !(4.3)

      call profiling_out(cpmd_orth)

    case(VEL_VERLET)

      SAFE_ALLOCATE(oldpsi(1:gr%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end))

      do ist1 = st%st_start, st%st_end
        
        ! calculate the "force"
        call X(hamiltonian_apply)(hm, gr%der, st%X(psi)(:, :, ist1, ik), hpsi, ist1, ik)
        
        if(iter == 1) then 
          ! give the initial conditions
          this%X(psi2)(1:np, 1:ddim, ist1, ik) = M_ZERO
        end if

        oldpsi(1:np, 1:ddim, ist1) = st%X(psi)(1:np, 1:ddim, ist1, ik)

        !(4.8)
        this%X(psi2)(1:np, 1:ddim, ist1, ik) = &
          this%X(psi2)(1:np, 1:ddim, ist1, ik) + dt*M_HALF/this%emass*(-st%occ(ist1, ik)*hpsi(1:np, 1:ddim))

        st%X(psi)(1:np, 1:ddim, ist1, ik) = st%X(psi)(1:np, 1:ddim, ist1, ik) + dt*this%X(psi2)(1:np, 1:ddim, ist1, ik)

      end do

      call profiling_in(cpmd_orth, "CP_ORTHOGONALIZATION")
      
      call calc_xx()


      ! psi <= psi + X * oldpsi
      ! psi2 <= psi2 + 1/dt * X * oldpsi

      call batch_init(oldpsib, st%d%dim, st%st_start, st%st_end, oldpsi(:, :, :))
      call X(mesh_batch_rotate)(gr%mesh, oldpsib, xx)
      call batch_end(oldpsib)
      
      do ist1 = st%st_start, st%st_end
        do idim = 1, st%d%dim
          forall(ip = 1:gr%mesh%np)
            st%X(psi)(ip, idim, ist1, ik) = st%X(psi)(ip, idim, ist1, ik) + oldpsi(ip, idim, ist1)
            this%X(psi2)(ip, idim, ist1, ik) = this%X(psi2)(ip, idim, ist1, ik) + oldpsi(ip, idim, ist1)/dt
          end forall
        end do
      end do

      call profiling_out(cpmd_orth)
      
      SAFE_DEALLOCATE_P(oldpsi)
      
    end select

  end do
  
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(xx)

  call profiling_out(cpmd_prop)
  POP_SUB(X(cpmd_propagate))

contains

  subroutine calc_xx()
    integer :: ist1, ist2, it
    FLOAT   :: res
    FLOAT,  allocatable :: ii(:, :)
    R_TYPE, allocatable :: aa(:, :), bb(:, :), xxi(:, :)

    PUSH_SUB(X(cpmd_propagate).calc_xx)

    SAFE_ALLOCATE( aa(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE( bb(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE( ii(1:st%nst, 1:st%nst))
    SAFE_ALLOCATE(xxi(1:st%nst, 1:st%nst))

    do ist1 = 1, st%nst
      do ist2 = 1, st%nst
        ii(ist1, ist2) = ddelta(ist1, ist2)
      end do
    end do

    call X(states_calc_overlap)(st, gr%mesh, ik, aa)
    call X(states_calc_overlap)(st, gr%mesh, ik, bb, psi2 = oldpsi)

    xx = M_HALF*(ii - aa) !(4.6)

    do it = 1, 10
      xxi = M_HALF*(ii - aa + matmul(xx, ii - bb) + matmul(ii - transpose(bb), xx) - matmul(xx, xx)) !(4.5)
      res = maxval(abs(xxi - xx))
      xx = xxi
      if (res < CNST(1e-5)) exit
    end do

    SAFE_DEALLOCATE_A(aa)
    SAFE_DEALLOCATE_A(bb)
    SAFE_DEALLOCATE_A(ii)
    SAFE_DEALLOCATE_A(xxi)

    POP_SUB(X(cpmd_propagate).calc_xx)
  end subroutine calc_xx

end subroutine X(cpmd_propagate)

subroutine X(cpmd_propagate_vel)(this, gr, hm, st, dt)
  type(cpmd_t),         intent(inout) :: this
  type(grid_t),         intent(inout) :: gr
  type(hamiltonian_t),  intent(inout) :: hm
  type(states_t),       intent(inout) :: st
  FLOAT,                intent(in)    :: dt

  integer :: ik, ist1, ddim, np

  R_TYPE, allocatable :: hpsi(:, :), yy(:, :)
  R_TYPE :: one

  one = R_TOTYPE(M_ONE)

  if ( this%method == VERLET ) return

  PUSH_SUB(X(cpmd_propagate_vel))
  call profiling_in(cpmd_prop, "CP_PROPAGATION")

  np = gr%mesh%np
  ddim = st%d%dim

  SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(yy(1:st%nst, 1:st%nst))

  this%ecorr = M_ZERO

  do ik = 1, st%d%nik

    do ist1 = st%st_start, st%st_end

      ! calculate the "force"
      call X(hamiltonian_apply)(hm, gr%der, st%X(psi)(:, :, ist1, ik), hpsi, ist1, ik)

      ! we have to complete the propagation of psi2 from the previous step
      this%X(psi2)(1:np, 1:ddim, ist1, ik) = &
        this%X(psi2)(1:np, 1:ddim, ist1, ik) + dt*M_HALF/this%emass*(-st%occ(ist1, ik)*hpsi(1:np, 1:ddim)) !(4.9) 2nd part

      this%ecorr = this%ecorr + this%emass*X(mf_nrm2)(gr%mesh, ddim, this%X(psi2)(:, :, ist1, ik))**2 !(2.11)

    end do

    call profiling_in(cpmd_orth, "CP_ORTHOGONALIZATION")          

    call calc_yy()

    ! psi2 <= psi2 + Y * psi
    call states_block_matr_mul_add(gr%mesh, st, one, st%st_start, st%st_end, st%st_start, st%st_end, &
      st%X(psi)(:, :, :, ik), yy, one, this%X(psi2)(:, :, :, ik)) !(4.11)

    call calc_yy()

    call profiling_out(cpmd_orth)


  end do

  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(yy)

  call profiling_out(cpmd_prop)
  POP_SUB(X(cpmd_propagate_vel))

contains

  subroutine calc_yy()
    R_TYPE, allocatable :: cc(:, :)

    PUSH_SUB(X(cpmd_propagate_vel).calc_yy)

    SAFE_ALLOCATE(cc(1:st%nst, 1:st%nst))

    call X(states_calc_overlap)(st, gr%mesh, ik, cc, psi2 = this%X(psi2)(:, :, :, ik))

    yy = -M_HALF*(cc + transpose(cc))

    SAFE_DEALLOCATE_A(cc)
    POP_SUB(X(cpmd_propagate_vel).calc_yy)
  end subroutine calc_yy

end subroutine X(cpmd_propagate_vel)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
