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

#include "global.h"

module static_pol_lr
  use geometry
  use system
  use states
  use restart
  use hamiltonian 
  use fxc
  use mix
  use mesh
  use poisson
  use linear_response

  implicit none

contains

  integer function static_pol_lr_run(sys, h, fromScratch) result(ierr)
    type(system_type), intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    logical, intent(inout) :: fromScratch

    type(lr_type) :: lr
    FLOAT :: pol(conf%dim, conf%dim)

    ierr = 0
    call init_()

    ! load wave-functions
    if(X(restart_read) ("tmp/restart_gs", sys%st, sys%m).ne.sys%st%nst) then
      message(1) = "Could not load wave-functions: Starting from scratch"
      call write_warning(1)
      
      ierr = 1
      call end_()
      return
    end if

    ! setup Hamiltonian
    message(1) = 'Info: Setting up Hamiltonian.'
    call write_info(1)
    call X(system_h_setup) (sys, h)

    !if(.not.fromScratch) then ! try to load delta_psi
    !  if(X(restart_read) ("tmp/restart_lr_static_pol", sys%st, sys%m).ne.sys%st%nst) then
    fromScratch = .true.

    call lr_init(sys%st, sys%m, lr)
    call pol_tensor(sys, h, lr, pol)
    call output()
    call lr_end(lr)

    call end_()

  contains

    subroutine init_()
      call push_sub('static_pol_lr_run')

      ! allocate wfs
      allocate(sys%st%X(psi)(sys%m%np, sys%st%dim, sys%st%nst, sys%st%nik))

    end subroutine init_

    subroutine end_()
      deallocate(sys%st%X(psi))
      
      call pop_sub()
    end subroutine end_
    
    subroutine output()
      integer :: j, iunit

      call loct_mkdir("linear")

      call io_assign(iunit)
      open(iunit, file='linear/polarizability_lr', status='unknown')
      write(iunit, '(2a)', advance='no') '# Static polarizability tensor [', &
         trim(units_out%length%abbrev)
      if(conf%dim.ne.1) write(iunit, '(a,i1)', advance='no') '^', conf%dim
      write(iunit, '(a)') ']'

      do j = 1, conf%dim
        write(iunit, '(3f12.6)') pol(j, 1:conf%dim) &
           / units_out%length%factor**conf%dim
      end do
      call io_close(iunit)
      
    end subroutine output

  end function static_pol_lr_run


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
      call get_response_e(sys, h, lr, i)
      call mix_end(lr%mixer)

      do j = 1, sys%m%np
        rhov = sum(lr%dl_rho(j,:))*sys%m%vol_pp(j)
        pol(i, :) = pol(i, :) - sys%m%x(j, :)*rhov
      end do
    end do
    
    call pop_sub()
  end subroutine pol_tensor


  subroutine get_response_e(sys, h, lr, alpha)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(lr_type),          intent(inout) :: lr
    integer,                intent(in)    :: alpha
    
    integer :: iter, ik, ik2, ist, i 
    FLOAT, allocatable :: tmp(:), Y(:,:), dl_rhoin(:,:,:), dl_rhonew(:,:,:), dl_rhotmp(:,:,:)
    logical :: finish
    
    type(mesh_type), pointer :: m
    type(states_type), pointer :: st

    call push_sub('get_response_e')

    m => sys%m
    st => sys%st

    allocate(tmp(m%np), Y(m%np,st%d%dim))
    allocate(dl_rhoin(1,m%np,st%d%dim), dl_rhonew(1,m%np,st%d%dim), dl_rhotmp(1,m%np,st%d%dim))

    call init_response_e()

    do iter=1, lr%max_iter
      dl_rhoin(1,:,:) = lr%dl_rho(:,:)

      do i = 1, m%np
        tmp(i) = sum(lr%dl_rho(i,:))
      end do
      if(.not.h%ip_app) call poisson_solve(m, sys%f_der, lr%dl_Vhar, tmp) 
      
      do ik = 1, st%d%nspin
        do ist = 1, st%nst 
          if (st%occ(ist, ik) > M_ZERO) then  

            Y(:,1) = (lr%dl_Vhar(:) + m%x(:,alpha))*tmp(:) 
            do ik2 = 1, st%d%nspin
              Y(:,1) = Y(:,1) + lr%dl_Vxc(:, ik, ik2)*lr%dl_rho(:,ik2)  
            end do
            do ik2 = 1, st%d%dim
              Y(:,ik2) = -Y(:,1)*lr%X(dl_psi)(:, ik2, ist, ik)
            end do

            call lr_orth_vector(m, st, Y, ik)
            call lr_solve_AXeY(sys, h, lr, ist, ik, Y, 500, CNST(1.0e-5))
          endif
        end do
      end do

      call lr_orth_response(m, st, lr)
      call lr_build_dl_rho(m, st, lr)

      dl_rhonew(:,:,:) = M_ZERO
      dl_rhotmp(1,:,:) = lr%dl_rho(:,:)
      call mixing(lr%mixer, iter, 1, m%np, st%d%nspin, &
         dl_rhoin, dl_rhotmp, dl_rhonew)  

      lr%abs_dens = M_ZERO
      do ik = 1, st%d%nspin
        tmp(:) = (dl_rhoin(1,:,ik)-lr%dl_rho(:,ik))**2
        lr%abs_dens = lr%abs_dens + dmf_integrate(m, tmp)
      end do
      lr%abs_dens = sqrt(lr%abs_dens)
      finish = (lr%abs_dens <= lr%conv_abs_dens) 

      if(finish) then 
        write(message(1), '(a, i4, a)')        &
           'Info: SCF for response converged in ', &
           iter, 'iterations'  
        exit
      else  
        ! without mixing it seems to go faster!!!!
        !lr%dl_rho(:,:) = dl_rhonew(1,:,:)
      end if
    end do

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
      lr%dl_Vhar(:) = M_ZERO
      
      call lr_orth_response(m, st, lr)
      call lr_build_dl_rho(m, st, lr)
      
    end subroutine init_response_e

  end subroutine get_response_e
  
end module static_pol_lr
