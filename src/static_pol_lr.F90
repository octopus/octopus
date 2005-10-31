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
  public :: &
    static_pol_lr_run

contains

  ! ---------------------------------------------------------
  integer function static_pol_lr_run(sys, h, fromScratch) result(ierr)
    type(system_type), target, intent(inout) :: sys
    type(hamiltonian_type),    intent(inout) :: h
    logical,                   intent(inout) :: fromScratch

    type(lr_type) :: lr
    type(grid_type), pointer :: gr
    FLOAT :: pol(1:sys%gr%sb%dim, 1:sys%gr%sb%dim)
    FLOAT :: hpol(1:sys%gr%sb%dim, 1:sys%gr%sb%dim, 1:sys%gr%sb%dim)
    integer :: err

    ierr = 0
    call init_()

    ! load wave-functions
    call X(restart_read) (trim(tmpdir)//'restart_gs', sys%st, gr%m, err)
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
    !  if(X(restart_read) (trim(tmpdir)//'restart_lr_static_pol', sys%st, m).ne.0) then
    fromScratch = .true.

    call lr_init(lr, "SP")


    call X(lr_alloc_fHxc) (sys%st, gr%m, lr)
    err = X(lr_alloc_psi) (sys%st, gr%m, lr)

    call lr_build_fxc(gr%m, sys%st, sys%ks%xc, lr%dl_Vxc)

    !   call pol_tensor(sys, h, lr, pol)
    call hyperpol_tensor(sys, h, lr, pol, hpol)
    call output()
    call lr_dealloc(lr)

    call end_()

  contains

    ! ---------------------------------------------------------
    subroutine init_()
      call push_sub('static_pol_lr.static_pol_lr_run')

      gr => sys%gr

      ! allocate wfs
      allocate(sys%st%X(psi)(NP, sys%st%d%dim, sys%st%nst, sys%st%d%nik))

    end subroutine init_


    ! ---------------------------------------------------------
    subroutine end_()
      deallocate(sys%st%X(psi))

      call pop_sub()
    end subroutine end_


    ! ---------------------------------------------------------
    subroutine output()
      integer :: j, iunit,i,k
      FLOAT :: msp
      call io_mkdir('linear')

      !! Output polarizabilty

      iunit = io_open('linear/polarizability_lr', action='write')
      write(iunit, '(2a)', advance='no') '# Static polarizability tensor [', &
        trim(units_out%length%abbrev)
      if(NDIM.ne.1) write(iunit, '(a,i1)', advance='no') '^', NDIM
      write(iunit, '(a)') ']'

      msp = M_ZERO
      do j = 1, NDIM
        write(iunit, '(3f12.6)') pol(j, 1:NDIM) &
          / units_out%length%factor**NDIM
        msp = msp + pol(j,j)
      end do
      msp = msp / M_THREE

      write(iunit, '(a, f12.6)')  'Mean static polarizability', msp &
        / units_out%length%factor**NDIM

      call io_close(iunit)

      !! Output first hyperpolarizabilty (beta)
      iunit = io_open('linear/beta_lr', action='write')
      write(iunit, '(2a)', advance='no') '# Static hyperpolarizability tensor [', &
        trim(units_out%length%abbrev)
      if(NDIM.ne.1) write(iunit, '(a,i1)', advance='no') '^', NDIM+2
      write(iunit, '(a)') ']'

      write(iunit, '(a)') '# WARNING, this values are not reliable until further testing'

      do i=1,NDIM
        do j=1,NDIM
          do k=1,NDIM
            write(iunit,'(3i1,3f12.6)') i,j,k, hpol(i,j,k)
          end do
        end do
      end do
      !      write(iunit, '(a, f12.6)')  'Mean static polarizability', msp &
      !         / units_out%length%factor**NDIM

      call io_close(iunit)

    end subroutine output

  end function static_pol_lr_run

#if 0
  ! ---------------------------------------------------------
  ! This is not used anymore, pol is calculated at the same time that hpol
  subroutine pol_tensor(sys, h, lr, pol)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(lr_type),          intent(inout) :: lr
    FLOAT,                  intent(out)   :: pol(:,:)

    integer :: i, j
    FLOAT :: rhov
    call push_sub('static_pol_lr.pol_tensor')

    pol = M_ZERO
    do i = 1, sys%gr%sb%dim
      write(message(1), '(a,i1)') 'Info: Calculating polarizability for direction ', i
      call write_info(1)

      call mix_init(lr%mixer, sys%gr%m, 1, sys%st%d%nspin)
      call get_response_e(sys, h, lr, i, R_TOTYPE(M_ZERO))
      call mix_end(lr%mixer)

      do j = 1, sys%gr%m%np
        rhov = sum(lr%X(dl_rho)(j,1:sys%st%d%nspin))*sys%gr%m%vol_pp(j)
        pol(i, :) = pol(i, :) - sys%gr%m%x(j,:)*rhov
      end do

    end do

    call pop_sub()
  end subroutine pol_tensor
#endif

  ! ---------------------------------------------------------
  subroutine hyperpol_tensor(sys, h, lr, pol, hpol)
    type(system_type),      intent(inout) :: sys
    type(hamiltonian_type), intent(inout) :: h
    type(lr_type),          intent(inout) :: lr
    FLOAT,                  intent(out)   :: pol(:,:)
    FLOAT,                  intent(out)   :: hpol(:,:,:)

    integer :: i, j, k
    integer :: ispin, ist, ispin2, ist2,n,np, dim

    R_TYPE ::prod

    R_TYPE, allocatable :: tmp(:,:), dphide(:,:,:,:), dVde(:,:,:), drhode(:,:)

    FLOAT, allocatable :: kxc(:,:,:,:)
    FLOAT :: vol

    np=sys%gr%m%np
    dim=sys%gr%sb%dim

    call push_sub('static_pol_lr.hyperpol_tensor')

    print*, "Calculating static hyperpolarizability"

    allocate(tmp(1:np,1))

    ! here we store the derivatives of the orbitals in each direction
    allocate(dphide(1:np,1:sys%st%nst,1:sys%st%d%nspin,1:dim))

    ! the derivative of the xc potential
    allocate(dVde(1:np,1:sys%st%d%nspin,1:dim))

    allocate(drhode(1:np,1:dim))

    allocate(kxc(1:np, 1:sys%st%d%nspin, 1:sys%st%d%nspin, 1:sys%st%d%nspin ))

    call xc_get_kxc(sys%ks%xc, sys%gr%m, sys%st%rho, sys%st%d%ispin, kxc)

    ! first the derivatives in all directions are calculated and stored
    write(message(1), '(a,i1)') 'Info: Calculating derivatives of the orbitals:'
    call write_info(1)

    do i = 1, dim
      write(message(1), '(a,i1)') '   direction: ', i
      call write_info(1)

      call mix_init(lr%mixer, sys%gr%m, 1, sys%st%d%nspin)

      call get_response_e(sys, h, lr, i, R_TOTYPE(M_ZERO))
      call mix_end(lr%mixer)


      do ispin = 1, sys%st%d%nspin

        ! the potential derivatives

        ! Hartree and the potential and the derivative of the
        !  potential associated at the electric field

        dVde(1:np,ispin,i)=lr%X(dl_Vhar)(1:np) + sys%gr%m%x(1:np,i)

        ! XC
        do ispin2 = 1, sys%st%d%nspin
          dVde(1:np,ispin,i)=dVde(1:np,ispin,i)+lr%dl_Vxc(1:np, ispin, ispin2)*lr%X(dl_rho)(1:np,ispin2)
        end do


        ! the derivatives of the orbitals
        do ist = 1, sys%st%nst
          dphide(1:np,ist,ispin,i)=lr%X(dl_psi)(1:np,1,ist,ispin)
        end do

      end do

      ! the density
      do n=1,np
        drhode(n,i)=sum(lr%X(dl_rho)(n,1:sys%st%d%nspin))
      end do

    end do

    print*, "Calculating polarizability tensor"

    pol = M_ZERO
    do i = 1, dim
      !      write(message(1), '(a,i1)') 'Info: Calculating polarizability for direction ', i
      call write_info(1)

      do j = 1,np
        pol(i,1:dim)=pol(i,1:dim)-sys%gr%m%x(j,1:dim)*drhode(j,i)*sys%gr%m%vol_pp(j)
      end do

    end do


    print*, "Calculating hyperpolarizability tensor"
    hpol = M_ZERO

    do i = 1, dim
      do j = 1, dim
        do k = 1, dim

          !             print*,"Component",i,j,k

          do ispin = 1, sys%st%d%nspin
            do ist = 1, sys%st%nst

              ! <D\psi_n | P_c DV_scf P_c | D\psi_n >

              tmp(1:np,1)=dphide(1:np,ist,ispin,k)

              tmp(1:np,1)=dVde(1:np, ispin, j)*tmp(1:np,1)

              hpol(i,j,k)=M_TWO*sum(R_CONJ(dphide(1:np,ist,ispin,i))*tmp(1:np,1)*sys%gr%m%vol_pp(1:np))



              ! -<D\psi_n| P_c | D\psi_m > < \psi_m| DV_scf | \psi_n >

              do ispin2 = 1, sys%st%d%nspin
                do ist2 = 1, sys%st%nst

                  tmp(1:np,1)=dphide(1:np,ist2,ispin2,j)

                  prod=sum(R_CONJ(dphide(1:np,ist,ispin,i))*tmp(1:np,1)*sys%gr%m%vol_pp(1:np))

                  !WARNING check spin dependency of the potential

                  tmp(1:np,1)=dVde(1:np,ispin2,k)*sys%st%X(psi)(1:np,1,ist2,ispin2)
                  prod=prod*sum(R_CONJ(sys%st%X(psi)(1:np,1,ist,ispin))*tmp(1:np,1)*sys%gr%m%vol_pp(1:np))

                  hpol(i,j,k)=hpol(i,j,k)-M_TWO*prod

                end do
              end do

              !                   print*,"K3"
              hpol(i,j,k)=hpol(i,j,k)+&
                sum(kxc(1:np,1,1,1)*drhode(1:np,i)*drhode(1:np,j)*drhode(1:np,k)*sys%gr%m%vol_pp(1:np))/CNST(6.0)
              !                   print*,"done.."
            end do
          end do




        end do
      end do
    end do

    vol=sum(sys%gr%m%vol_pp(1:np))

    !    print*, hpol

    hpol=M_SIX*hpol

    !    print*, "POL", pol(1:dim,1:dim)
    !    print*, "HPOL", hpol(1:dim,1:dim,1:dim)

    deallocate(tmp)
    deallocate(dphide)
    deallocate(dVde)
    deallocate(drhode)
    deallocate(kxc)

    call pop_sub()
  end subroutine hyperpol_tensor


  ! ---------------------------------------------------------
  subroutine get_response_e(sys, h, lr, alpha, omega)
    type(system_type), target, intent(inout) :: sys
    type(hamiltonian_type),    intent(inout) :: h
    type(lr_type),             intent(inout) :: lr
    integer,                   intent(in)    :: alpha
    R_TYPE,                    intent(in)    :: omega

    integer :: iter, iter_max, sigma, ik, ik2, ist, i
    FLOAT, allocatable :: tmp(:), dl_rhoin(:,:,:), dl_rhonew(:,:,:), dl_rhotmp(:,:,:)
    R_TYPE, allocatable :: Y(:, :)
    logical :: finish

    type(mesh_type), pointer :: m
    type(states_type), pointer :: st

    call push_sub('static_pol_lr.get_response_e')

    m => sys%gr%m
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
        call dpoisson_solve(sys%gr, lr%ddl_Vhar, tmp)
      end if

      lr%X(dl_rho) = M_ZERO
      do sigma = -1, 1, 2
        if(omega==M_ZERO .and. sigma==1) cycle

        do ik = 1, st%d%nspin
          do ist = 1, st%nst
            if(st%occ(ist, ik) <= M_ZERO) cycle

            Y(:,1) = (lr%ddl_Vhar(:) + m%x(:,alpha))
            do ik2 = 1, st%d%nspin
              Y(1:m%np,1) = Y(1:m%np,1) + lr%dl_Vxc(1:m%np, ik, ik2)*lr%X(dl_rho)(1:m%np,ik2)
            end do
            Y(1:m%np,1) = -Y(1:m%np,1)*st%X(psi)(1:m%np, 1, ist, ik)

            call X(lr_orth_vector)(m, st, Y, ik)

            iter_max = 200
            call X(lr_solve_HXeY) (lr, h, sys%gr, sys%st%d, ik, lr%X(dl_psi)(:,:, ist, ik), Y, &
              -sys%st%eigenval(ist, ik) + real(sigma, PRECISION)*omega)

            print *, iter, ik, ist, sum(lr%X(dl_psi)(1:m%np,1, ist, ik)**2*sys%gr%m%vol_pp(1:m%np))
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
      call mixing(lr%mixer, m, iter, 1, st%d%nspin, &
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
          end if
        end do
      end do
      lr%ddl_Vhar(:) = M_ZERO

      call lr_orth_response(m, st, lr)
      call X(lr_build_dl_rho)(m, st, lr, 1)
      lr%X(dl_rho) = M_TWO*lr%X(dl_rho)

    end subroutine init_response_e

  end subroutine get_response_e

end module static_pol_lr
