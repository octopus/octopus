!! Copyright (C) 2004 Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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

module static_pol_lr
  use global
  use messages
  use units
  use mesh
  use mesh_function
  use system
  use restart
  use hamiltonian 
  use mix
  use poisson
  use linear_response
  use io

  implicit none

  private
  public :: static_pol_lr_run

contains

  ! ---------------------------------------------------------
  integer function static_pol_lr_run(sys, h, fromScratch) result(ierr)
    type(system_type), intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    logical, intent(inout) :: fromScratch

    type(lr_type) :: lr
    FLOAT :: pol(conf%dim, conf%dim)
    integer :: err
    
    ierr = 0
    call init_()

    ! load wave-functions
    call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, sys%m, err)
    if(err.ne.0) then
      message(1) = "Could not load wave-functions in pol_lr_run: Starting from scratch"
      call write_warning(1)
      
      ierr = 1
      call end_()
      return
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)
    call X(system_h_setup) (sys, h)
    
    !if(.not.fromScratch) then ! try to load delta_psi
    !  if(X(restart_read) (trim(tmpdir)//'restart_lr_static_pol', sys%st, sys%m).ne.0) then
    fromScratch = .true.

    call lr_init(lr, "SP")
    call X(lr_alloc_fHxc) (sys%st, sys%m, lr)
    err = X(lr_alloc_psi) (sys%st, sys%m, lr)
    call lr_build_fxc(sys%m, sys%st, sys%ks%xc, lr%dl_Vxc)
    call pol_tensor(sys, h, lr, pol)
    call output()
    call lr_dealloc(lr)

    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('static_pol_lr_run')

      ! allocate wfs
      allocate(sys%st%X(psi)(sys%m%np, sys%st%d%dim, sys%st%nst, sys%st%d%nik))

    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      deallocate(sys%st%X(psi))
      
      call pop_sub()
    end subroutine end_
    

    ! ---------------------------------------------------------
    subroutine output()
      integer :: j, iunit
      FLOAT :: msp
      call io_mkdir('linear')

      iunit = io_open('linear/polarizability_lr', action='write')
      write(iunit, '(2a)', advance='no') '# Static polarizability tensor [', &
         trim(units_out%length%abbrev)
      if(conf%dim.ne.1) write(iunit, '(a,i1)', advance='no') '^', conf%dim
      write(iunit, '(a)') ']'

      msp = M_ZERO
      do j = 1, conf%dim
        write(iunit, '(3f12.6)') pol(j, 1:conf%dim) &
           / units_out%length%factor**conf%dim
        msp = msp + pol(j,j)
      end do
      msp = msp / M_THREE

      write(iunit, '(a, f12.6)')  'Mean static polarizability', msp &
         / units_out%length%factor**conf%dim

      call io_close(iunit)
      
    end subroutine output

  end function static_pol_lr_run


  ! ---------------------------------------------------------
  subroutine pol_tensor(sys, h, lr, pol)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(lr_type),          intent(inout) :: lr
    FLOAT,                  intent(out)   :: pol(:,:)

    integer :: i, j
    FLOAT :: rhov
    call push_sub('pol_tensor')

    pol = M_ZERO
    do i = 1, conf%dim
      write(message(1), '(a,i1)') 'Info: Calculating polarizability for direction ', i
      call write_info(1)

      call mix_init(lr%mixer, 1, sys%m%np, sys%st%d%nspin)
      call get_response_e(sys, h, lr, i, R_TOTYPE(M_ZERO))
      call mix_end(lr%mixer)
      
      do j = 1, sys%m%np
        rhov = sum(lr%X(dl_rho)(j,:))*sys%m%vol_pp(j)
        pol(i, :) = pol(i, :) - sys%m%x(j, :)*rhov
      end do
      
    end do
    
    call pop_sub()
  end subroutine pol_tensor


  ! ---------------------------------------------------------
  subroutine get_response_e(sys, h, lr, alpha, omega)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(lr_type),          intent(inout) :: lr
    integer,                intent(in)    :: alpha
    R_TYPE,                 intent(in)    :: omega
    
    integer :: iter, iter_max, sigma, ik, ik2, ist, i 
    FLOAT, allocatable :: tmp(:), dl_rhoin(:,:,:), dl_rhonew(:,:,:), dl_rhotmp(:,:,:)
    R_TYPE, allocatable :: Y(:, :)
    logical :: finish
    
    type(mesh_type), pointer :: m
    type(states_type), pointer :: st
    
    call push_sub('get_response_e')

    m => sys%m
    st => sys%st

    allocate(tmp(m%np), Y(m%np,1))
    allocate(dl_rhoin(1,m%np,st%d%nspin), dl_rhonew(1,m%np,st%d%nspin), dl_rhotmp(1,m%np,st%d%nspin))

    call init_response_e()

    iter_loop: do iter=1, lr%max_iter

      dl_rhoin(1,:,:) = lr%X(dl_rho)(:,:)

      if(.not.h%ip_app) then
        do i = 1, m%np
          tmp(i) = sum(lr%X(dl_rho)(i,:))
        end do
        call dpoisson_solve(m, sys%f_der, lr%ddl_Vhar, tmp) 
      end if 

      lr%X(dl_rho) = M_ZERO
      do sigma = -1, 1, 2
        if(omega==M_ZERO .and. sigma==1) cycle

        do ik = 1, st%d%nspin
          do ist = 1, st%nst
            if(st%occ(ist, ik) <= M_ZERO) cycle

            Y(:,1) = (lr%ddl_Vhar(:) + m%x(:,alpha))
            do ik2 = 1, st%d%nspin
              Y(:,1) = Y(:,1) + lr%dl_Vxc(:, ik, ik2)*lr%X(dl_rho)(:,ik2)  
            end do
            Y(:,1) = -Y(:,1)*st%X(psi)(:, 1, ist, ik)
            
            call X(lr_orth_vector)(m, st, Y, ik)
            
            iter_max = 50
            call X(lr_solve_HXeY) (lr, h, sys%m, sys%f_der, sys%st%d, ik, lr%X(dl_psi)(:,:, ist, ik), Y, &
               -sys%st%eigenval(ist, ik) + real(sigma, PRECISION)*omega)

            print *, lr%iter, sum(lr%X(dl_psi)(:,1, ist, ik)**2*sys%m%vol_pp(:))
          end do
        end do


        !call lr_orth_response(m, st, lr)
        call X(lr_build_dl_rho)(m, st, lr, 1)
        lr%X(dl_rho) = M_TWO*lr%X(dl_rho)
      end do

      ! if static perturbation, then psi^+ and psi^- are the same
      !if(omega==M_ZERO) lr%X(dl_rho) = M_TWO*lr%X(dl_rho)

      ! mix to get new density
      dl_rhonew(1,:,:) = M_ZERO
      dl_rhotmp(1,:,:) = lr%X(dl_rho)(:,:)
      call mixing(lr%mixer, iter, 1, m%np, st%d%nspin, &
         dl_rhoin, dl_rhotmp, dl_rhonew)  

      ! check for convergence
      lr%abs_dens = M_ZERO
      do ik = 1, st%d%nspin
        tmp(:) = (dl_rhoin(1,:,ik) - lr%X(dl_rho)(:,ik))**2
        lr%abs_dens = lr%abs_dens + dmf_integrate(m, tmp)
      end do
      lr%abs_dens = sqrt(lr%abs_dens)

      ! are we finished?
      finish = (lr%abs_dens <= lr%conv_abs_dens) 
      if(finish) then 
        write(message(1), '(a, i4, a)')        &
           'Info: SCF for response converged in ', &
           iter, ' iterations'  
        call write_info(1)
        exit
      else  
        lr%X(dl_rho)(:,:) = dl_rhonew(1,:,:)
      end if

    end do iter_loop

    deallocate(dl_rhoin, dl_rhonew, dl_rhotmp)
    deallocate(tmp, Y)
    call pop_sub()

  contains
    subroutine init_response_e()
      integer :: ik, ist, i
      FLOAT :: rd
      
      do ik = 1, st%d%nspin
        do ist = 1, st%nst
          if (st%occ(ist, ik) > M_ZERO) then
            do i = 1, m%np
              call mesh_r(m, i, rd)
              lr%X(dl_psi)(i, 1, ist, ik) = st%X(psi)(i, 1, ist, ik)*rd*exp(-rd)
            end do
          endif
        end do
      end do
      lr%ddl_Vhar(:) = M_ZERO
      
      call lr_orth_response(m, st, lr)
      call X(lr_build_dl_rho)(m, st, lr, 1)
      lr%X(dl_rho) = M_TWO*lr%X(dl_rho)
      
    end subroutine init_response_e

  end subroutine get_response_e
  
end module static_pol_lr
