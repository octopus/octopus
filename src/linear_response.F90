!! Copyright (C) 2004 E.S. Kadantsev, M. Marques
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

#include "global.h"

module linear_response
  use global
  use mesh
  use states
  use mix
  use system
  use hamiltonian

  implicit none

  type lr_type
    FLOAT   :: conv_abs_dens  ! convergence required
    FLOAT   :: abs_dens       ! convergence reached
    integer :: max_iter       ! maximum number of iterations
    
    type(mix_type) :: mixer   ! can not live without it

    FLOAT, pointer :: dl_rho(:,:)      ! response of the density
    FLOAT, pointer :: ddl_psi(:,:,:,:) ! linear change of the real KS orbitals
    CMPLX, pointer :: zdl_psi(:,:,:,:) ! linear change of the complex KS orbitals
    
    FLOAT, pointer :: dl_Vhar(:)      ! linear change of the Hartree potential
    FLOAT, pointer :: dl_Vxc(:,:,:)   ! linear change of the xc potential (fxc)
  end type lr_type

contains

  ! ---------------------------------------------------------
  subroutine lr_init(st, m, lr)
    type(states_type), intent(in)  :: st
    type(mesh_type),   intent(in)  :: m
    type(lr_type),     intent(out) :: lr
    
    ! read some parameters from the input file
    call loct_parse_float("ConvAbsDens", CNST(1e-5), lr%conv_abs_dens)
    call loct_parse_int("MaximumIter", 200, lr%max_iter)
    
    ! allocate variables
    allocate(lr%dl_rho(m%np, st%d%nspin), lr%X(dl_psi)(m%np, st%d%dim, st%nst, st%nspin))
    allocate(lr%dl_Vhar(m%np), lr%dl_Vxc(m%np, st%d%nspin, st%d%nspin))

  end subroutine lr_init


  ! ---------------------------------------------------------
  subroutine lr_end(lr)
    type(lr_type), intent(inout) :: lr

    ASSERT(associated(lr%dl_rho))

    deallocate(lr%dl_rho, lr%X(dl_psi), lr%dl_Vhar, lr%dl_Vxc)
    nullify   (lr%dl_rho, lr%X(dl_psi), lr%dl_Vhar, lr%dl_Vxc)

  end subroutine lr_end

  
  ! ---------------------------------------------------------
  subroutine lr_build_dl_rho(m, st, lr)
    type(mesh_type),   intent(in)    :: m
    type(states_type), intent(in)    :: st
    type(lr_type),     intent(inout) :: lr
    
    integer :: i, p, ik, sp
    CMPLX   :: c

    call push_sub('lr_build_dl_rho')
    
    if(st%d%ispin == SPIN_POLARIZED) then
      sp = 2
    else
      sp = 1
    end if

    lr%dl_rho = M_ZERO
    do ik = 1, st%d%nik, sp
      do p  = st%st_start, st%st_end
        do i = 1, m%np
          lr%dl_rho(i, 1) = lr%dl_rho(i, 1) + M_TWO*st%d%kweights(ik)*st%occ(p, ik) * &
             R_REAL(st%X(psi)(i, 1, p, ik)*lr%X(dl_psi)(i, 1, p, ik))

          select case(st%d%ispin)
          case(SPIN_POLARIZED)
            lr%dl_rho(i, 2) = lr%dl_rho(i, 2) + M_TWO*st%d%kweights(ik+1)*st%occ(p, ik+1) * &
               R_REAL(st%X(psi)(i, 1, p, ik+1)*lr%X(dl_psi)(i, 1, p, ik+1))
          case(SPINORS)
            lr%dl_rho(i, 2) = lr%dl_rho(i, 2) + M_TWO*st%d%kweights(ik)*st%occ(p, ik) * &
               R_REAL(st%X(psi)(i, 2, p, ik)*lr%X(dl_psi)(i, 2, p, ik))


            c = st%X(psi)(i, 1, p, ik) * R_CONJ(lr%X(dl_psi)(i, 2, p, ik)) + &
                lr%X(dl_psi)(i, 1, p, ik) * R_CONJ(st%X(psi)(i, 2, p, ik))

            lr%dl_rho(i, 3) = lr%dl_rho(i, 3) + st%d%kweights(ik)*st%occ(p, ik) * R_REAL(c)
            lr%dl_rho(i, 4) = lr%dl_rho(i, 4) + st%d%kweights(ik)*st%occ(p, ik) * R_AIMAG(c)
          end select
        end do
      end do
    end do     

    call pop_sub()
  end subroutine lr_build_dl_rho


  ! ---------------------------------------------------------
  subroutine lr_solve_AXeY(sys, h, lr, ist, ik, Y, MAXITER, tol)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(lr_type),          intent(inout) :: lr
    R_TYPE,                 intent(in)    :: Y(:,:)   ! Y(m%np, st%dim) and st%dim=1
    integer,                intent(in)    :: ist, ik
    integer,                intent(in)    :: MAXITER
    FLOAT,                  intent(in)    :: tol
    
    integer :: iter, stdim  
    R_TYPE, allocatable :: z(:,:), g(:,:), pp(:,:), t(:,:)
    FLOAT  :: sm1, sm2, sm3, alpha, beta

    call push_sub('lr_solve_AXeY')

    
    allocate(z(sys%m%np, sys%st%dim), g(sys%m%np, sys%st%dim), pp(sys%m%np,sys%st%dim), t(sys%m%np,sys%st%dim))

    pp(:,:) = lr%X(dl_psi)(:,:, ist, ik)
    call X(Hpsi)(h, sys%m, sys%f_der, pp, z, ik)

    z (:,:) = -z(:,:) + sys%st%eigenval(ist, ik)*pp(:,:) + Y(:,:)
    pp(:,:) =  z(:,:)
    iter_loop: do iter = 1, MAXITER     
      call X(Hpsi)(h, sys%m, sys%f_der, pp, t, ik)
      t(:,:) = t(:,:) - sys%st%eigenval(ist, ik)*pp(:,:)

      sm1 = X(states_dotp) (sys%m, sys%st%dim, z, z)
      sm2 = X(states_dotp) (sys%m, sys%st%dim, z, t)
      alpha = sm1/sm2

      lr%X(dl_psi)(:,:, ist, ik) = lr%X(dl_psi)(:,:, ist, ik) + alpha*pp(:,:)
      g(:,:) = z(:,:)
      g(:,:) = g(:,:) - alpha*t(:,:)
      sm3 = X(states_dotp) (sys%m, sys%st%dim, g, g) 
      beta = sm3/sm1

      if(sm3 <= tol) then 
        exit iter_loop
      else 
        pp(:,:) = beta*pp(:,:) + g(:,:)
        z(:,:) = g(:,:) 
      end if
    end do iter_loop

    if(sm3 > tol) then 
      message(1) = "Response using CG: Not converged!"
      write(message(2),*) "iter = ", iter, " sm3 = ", sm3
      call write_warning(2)
    end if
    deallocate(z, g, pp, t)
    
    call pop_sub()
  end subroutine lr_solve_AXeY
  

  ! ---------------------------------------------------------
  ! orthogonalizes response of \alpha KS orbital to all occupied 
  ! \alpha KS orbitals  
  ! ---------------------------------------------------------
  subroutine lr_orth_response(m, st, lr)
    type(mesh_type),   intent(in)    :: m
    type(states_type), intent(in)    :: st
    type(lr_type),     intent(inout) :: lr
    
    integer :: ist, ik
    call push_sub('lr_orth_response')

    do ik = 1, st%d%nspin
      do ist = 1, st%nst 
        if(st%occ(ist, ik) > M_ZERO) then  
          call lr_orth_vector(m, st, lr%X(dl_psi)(:,:, ist, ik), ik) 
        endif
      end do
    end do
    
    call pop_sub()
  end subroutine lr_orth_response

  
  ! ---------------------------------------------------------
  subroutine lr_orth_vector(m, st, v, ik)
    type(mesh_type),        intent(in)    :: m
    type(states_type),      intent(in)    :: st
    R_TYPE,                 intent(inout) :: v(:,:)
    integer,                intent(in)    :: ik
    
    R_TYPE   :: scalp
    integer :: ist

    call push_sub('lr_orth_vector')
    
    do ist = 1, st%nst
      if(st%occ(ist, ik) > M_ZERO) then
        scalp = X(states_dotp)(m, st%dim, v, st%X(psi)(:,:, ist, ik))
        v(:,:) = v(:,:) - scalp*st%X(psi)(:,:, ist, ik)
      end if
    end do

    call pop_sub()
  end subroutine lr_orth_vector

end module linear_response
