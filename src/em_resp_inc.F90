!! Copyright (C) 2005 M. Marques, X. Andrade
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


! ---------------------------------------------------------
subroutine X(dynamic_response)(sys, h, lr, props, pol, w, status)
  type(system_t),      intent(inout) :: sys
  type(hamiltonian_t), intent(inout) :: h
  type(lr_t),          intent(inout) :: lr(:,:,:) ! dim, nsigma(=2), nfreq(=1)
  type(pol_props_t),   intent(in)    :: props
  CMPLX,               intent(inout) :: pol(:,:)
  R_TYPE,              intent(in)    :: w 
  type(status_t),      intent(out)   :: status
  
  integer :: dir, j, freq, sigma

  R_TYPE :: rhov

  call push_sub('static_pol_lr.dynamic')

  freq = 1
  pol = M_ZERO

  !iterate for each direction
  do dir = 1, sys%gr%sb%dim

    write(message(1), '(a,i1)') 'Info: Calculating polarizability for direction ', dir
    call write_info(1)
    
    do sigma = 1, 2
      call mix_init(lr(dir,sigma,freq)%mixer, sys%gr%m, 1, sys%st%d%nspin)
    end do

    call X(get_response_e)(sys, h, lr(:,:,freq), dir, 2 , w, props, status)

    do sigma = 1, 2
      call mix_end(lr(dir,sigma,freq)%mixer)
    end do

    !calculate the polarizability
    do j = 1, sys%gr%m%np
      rhov = sum(lr(dir,2,freq)%X(dl_rho)(j,1:sys%st%d%nspin))*sys%gr%m%vol_pp(j)
      pol(dir, :) = pol(dir, :) - sys%gr%m%x(j,:)*rhov
    end do

  end do


  call pop_sub()

end subroutine X(dynamic_response)


! -------------------------------------------------------------
! this is the central routine of the electromagnetic response
! it calculates the first order variations of the wavefunctions 
! for an electric field
!--------------------------------------------------------------

subroutine X(get_response_e)(sys, h, lr, dir, nsigma, omega, props, status)
  type(system_t), target, intent(inout) :: sys
  type(hamiltonian_t),    intent(inout) :: h
  type(lr_t),             intent(inout) :: lr(:,:) !ndim,nsigma
  integer,                intent(in)    :: dir 
  integer,                intent(in)    :: nsigma 
  R_TYPE,                 intent(in)    :: omega
  type(pol_props_t),      intent(in)    :: props
  type(status_t),optional,intent(out)   :: status

  FLOAT, allocatable :: diff(:,:,:)
  FLOAT :: dpsimod,freq_sign
  integer :: iter, sigma, ik, ik2, ist, i, ist2, err
  FLOAT, allocatable :: dl_rhoin(:, :, :, :), dl_rhonew(:, :, :, :), dl_rhotmp(:, :, :, :), dtmp(:)
  R_TYPE, allocatable :: Y(:, :),dV(:, :), tmp(:)

  logical :: finish(2)

  type(mesh_t), pointer :: m
  type(states_t), pointer :: st
  
  character(len=30) :: dirname

  call push_sub('static_pol_lr.get_response_e')

  ASSERT( nsigma==1 .or. nsigma ==2 )
  
  m => sys%gr%m
  st => sys%st
  
  ALLOCATE(tmp(m%np),m%np)
  ALLOCATE(dtmp(m%np),m%np)
  ALLOCATE(Y(m%np,1),m%np*1)
  ALLOCATE(dV(m%np,st%d%nspin),m%np*st%d%nspin)
  ALLOCATE(dl_rhoin(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
  ALLOCATE(dl_rhonew(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
  ALLOCATE(dl_rhotmp(1,m%np,st%d%nspin,nsigma),1*m%np*st%d%nspin*nsigma)
  ALLOCATE(diff(st%nst,st%d%nspin,nsigma),st%nst*st%d%nspin*nsigma)

  diff=M_ZERO

  call init_response_e()
  finish = .false.

  message(1)="--------------------------------------------"
  call write_info(1)

  !self consistency iteration for response
  iter_loop: do iter=1, lr(dir,1)%max_iter
    write(message(1), '(a, i3)') "LR SCF Iteration: ", iter
    write(message(2), '(a, f20.6, a, f20.6, a, i1)') "Frequency: ", R_REAL(omega), " Delta : ", R_AIMAG(omega), " Dir: ", dir
    call write_info(2)
    
    !we will use the two mixers, one for the real part and one for the complex 
    !(this is ugly, we should have complex mixers)
    dl_rhoin(1, :, :, 1) = R_REAL(lr(dir, 1)%X(dl_rho)(:, :))
    if(nsigma==2) dl_rhoin(1, :, :, 2) = R_AIMAG(lr(dir, 1)%X(dl_rho)(:, :))

    !calculate the variation of hartree term
    if (props%add_hartree) then
      do sigma=1,nsigma
        do i = 1, m%np
          tmp(i) = sum(lr(dir, sigma)%X(dl_rho)(i, :))
        end do
        call X(poisson_solve)(sys%gr, lr(dir,sigma)%X(dl_Vhar), tmp)
      end do
    end if
    
    do ik = 1, st%d%nspin
      do sigma = 1, nsigma

        if(sigma==1) then 
          freq_sign =  M_ONE
        else
          freq_sign = -M_ONE
        end if
        
        !Calculate H^(1):

        !* Vext
        dV(1:m%np, ik) = m%x(1:m%np, dir)

        !* hartree
        if (props%add_hartree) dV(1:m%np, ik) = dV(1:m%np, ik) + lr(dir, sigma)%X(dl_Vhar)(1:m%np) 

        !* fxc
        if(props%add_fxc) then 
          do ik2 = 1, st%d%nspin
            dV(1:m%np, ik) = dV(1:m%np, ik) + &
                 lr(dir, sigma)%dl_Vxc(1:m%np, ik, ik2)*lr(dir, sigma)%X(dl_rho)(1:m%np, ik2)
          end do
        end if

        !now calculate response for each state
        do ist = 1, st%nst
          if(st%occ(ist, ik) <= lr_min_occ) cycle

          !calculate the RHS of the Sternheimer eq
          Y(1:m%np, 1) = -dV(1:m%np, ik)*st%X(psi)(1:m%np, 1, ist, ik)

          !and project it into the unoccupied states
          call X(lr_orth_vector)(m, st, Y, ik)

          !solve the Sternheimer equation
          if(props%ort_each_step) then 
            call X(lr_solve_HXeY) (lr(dir, 1), h, sys%gr, sys%st%d, ik, lr(dir, sigma)%X(dl_psi)(:,:, ist, ik), Y, &
                 -sys%st%eigenval(ist, ik) + freq_sign*omega, sys%st)
          else
            call X(lr_solve_HXeY) (lr(dir, 1), h, sys%gr, sys%st%d, ik, lr(dir, sigma)%X(dl_psi)(:,:, ist, ik), Y, &
                 -sys%st%eigenval(ist, ik) + freq_sign*omega)
          end if

          !altough the dl_psi we get should be orthogonal to psi
          !a re-orthogonalization is sometimes necessary 
          if(.not. props%ort_each_step) then 
            call X(lr_orth_vector)(m, st, lr(dir,sigma)%X(dl_psi)(:,:, ist, ik), ik)
          end if

          !calculate and print the norm of the variations and how much
          !they change, this is only to have an idea of the converge process
          dpsimod = sum(R_ABS(lr(dir,sigma)%X(dl_psi)(1:m%np, 1, ist, ik))**2 * sys%gr%m%vol_pp(1:m%np))
          write(message(1), '(i4, f20.6, e20.6)') (3-2*sigma)*ist, dpsimod, (dpsimod-diff(ist,ik,sigma))
          call write_info(1)
          diff(ist,ik,sigma)=dpsimod

        end do !ist
      end do !sigma
    end do !ik

    ! calculate dl_rho
    if(nsigma == 2 ) then
      call build_rho_dynamic()
    else 
      lr(dir, 1)%X(dl_rho) = M_ZERO
      call X(lr_build_dl_rho)(m, st, lr(dir,1), type=3)
    end if

    dl_rhonew(1, 1:m%np, 1:st%d%nspin, 1:nsigma) = M_ZERO

    !write restart info
    do sigma=1,nsigma 
      write(dirname,'(a,i1,a,i1)') RESTART_DIR, dir, "_", sigma
      call restart_write(trim(tmpdir)//dirname, st, sys%gr, err, iter=iter, lr=lr(dir, sigma))
    end do
    
    !all the rest is the mixing and checking for convergency

    if( lr(dir,1)%max_iter == iter  ) then 
      message(1) = "Self-consistent iteration for response did not converge"
      if(present(status)) then 
        status%ok = .false. 
        call write_warning(1);
      else
        call write_fatal(1);
      end if 
    end if

    finish=.true.

    dl_rhotmp(1,:,:,1) = R_REAL(lr(dir,1)%X(dl_rho)(:,:))
    if(nsigma==2) dl_rhotmp(1,:,:,2) = R_AIMAG(lr(dir,1)%X(dl_rho)(:,:))

    do sigma = 1, nsigma

      call mixing(lr(dir,sigma)%mixer, m, iter, 1, st%d%nspin, &
           dl_rhoin(:,:,:,sigma), dl_rhotmp(:,:,:,sigma), dl_rhonew(:,:,:,sigma))

      ! check for convergence
      lr(dir, sigma)%abs_dens = M_ZERO

      do ik = 1, st%d%nspin
        dtmp(:) = (dl_rhoin(1,:,ik,sigma) - dl_rhotmp(1,:,ik,sigma))**2
        lr(dir,sigma)%abs_dens = lr(dir,sigma)%abs_dens + dmf_integrate(m, dtmp)
      end do

      lr(dir,sigma)%abs_dens = sqrt(lr(dir,sigma)%abs_dens)

      ! are we finished?
      finish(sigma) = (lr(dir,sigma)%abs_dens <= lr(dir,sigma)%conv_abs_dens)
    end do


    if(nsigma == 1) then 
      write(message(1), '(a, e20.6)') "Res ", lr(dir,1)%abs_dens
      lr(dir,1)%abs_dens=M_ZERO 
    else 
      write(message(1), '(a, 2e20.6)') "Res ", lr(dir,1)%abs_dens, lr(dir,2)%abs_dens
      lr(dir,1)%abs_dens=M_ZERO 
      lr(dir,2)%abs_dens=M_ZERO 
    endif
    message(2)="--------------------------------------------"
    call write_info(2)

    
      
    if( finish(1) .and. finish(2) ) then 

      if(present(status)) status%ok = .true.
      write(message(1), '(a, i4, a)')        &
           'Info: SCF for response converged in ', &
           iter, ' iterations'
      call write_info(1)
      exit

    else
      
      if(props%complex_response) then 
        lr(dir,1)%X(dl_rho)(:,:) = cmplx(dl_rhonew(1,:,:,1),dl_rhonew(1,:,:,2),PRECISION)
        lr(dir,2)%X(dl_rho)(:,:) = cmplx(dl_rhonew(1,:,:,1),-dl_rhonew(1,:,:,2),PRECISION)
      else
        do sigma=1,nsigma
          lr(dir,sigma)%X(dl_rho)(:,:) = dl_rhonew(1,:,:,1)
        end do
      end if

    end if

    
  end do iter_loop

  deallocate(tmp)
  deallocate(dtmp)
  deallocate(Y)
  deallocate(dV)
  deallocate(dl_rhoin)
  deallocate(dl_rhonew)
  deallocate(dl_rhotmp)
  deallocate(diff)

  call pop_sub()

contains

  !------------------------------------------------------------
  subroutine init_response_e()
    integer :: ik, ist, i
    FLOAT :: rd

    call push_sub('static_pol_lr.init_response_e')

    do sigma=1,nsigma
      lr(dir,sigma)%X(dl_Vhar)(:) = M_ZERO
      call X(lr_orth_response)(m, st, lr(dir,sigma))
    end do

    if(nsigma == 2 ) then 
      call build_rho_dynamic() 
    else
      call X(lr_build_dl_rho)(m, st, lr(dir,1), 3)
    end if

    call pop_sub()
  end subroutine init_response_e


  !------------------------------------------------------------
  subroutine build_rho_dynamic()
    R_TYPE  :: a


    do sigma=1,nsigma 
      lr(dir,sigma)%X(dl_rho)=M_ZERO
    end do

    do ik = 1, st%d%nspin
      do ist = 1, st%nst
        do i=1,m%np
          
          a=st%d%kweights(ik)*st%occ(ist, ik)*(&
               R_CONJ(st%X(psi)(i, 1, ist, ik))*lr(dir,1)%X(dl_psi)(i,1,ist,ik) +&
               R_CONJ(lr(dir,2)%X(dl_psi)(i,1,ist,ik))*st%X(psi)(i, 1, ist, ik))

          lr(dir,1)%X(dl_rho)(i,ik)=lr(dir,1)%X(dl_rho)(i,ik)+ a

          lr(dir,2)%X(dl_rho)(i,ik)=lr(dir,2)%X(dl_rho)(i,ik)+R_CONJ(a)

        end do
      end do
    end do
    
  end  subroutine build_rho_dynamic

end subroutine X(get_response_e)


! ---------------------------------------------------------
subroutine X(static_response) (sys, h, lr, props, pol, hpol, hpol_density)
  type(system_t),      intent(inout) :: sys
  type(hamiltonian_t), intent(inout) :: h
  type(lr_t),          intent(inout) :: lr(:,:) ! lr(NDIM,1,1)
  type(pol_props_t),     intent(in)    :: props
  FLOAT,               intent(out)   :: pol(:,:)
  FLOAT,               intent(out)   :: hpol(:,:,:)
  FLOAT,               intent(inout) :: hpol_density(:,:,:,:) !(1:np,MAX_DIM,MAX_DIM,MAX_DIM)

  integer :: i, j, k
  integer :: ispin, ist, ispin2, ist2, n, np, dim, ik

  R_TYPE :: prod

  R_TYPE, allocatable :: dVde(:,:,:), drhode(:,:)
  FLOAT, allocatable :: tmp(:,:)
  FLOAT,  allocatable :: kxc(:,:,:,:), hp_tmp(:,:,:,:)

  np  = sys%gr%m%np
  dim = sys%gr%sb%dim

  call push_sub('em_resp.static_response')

  message(1) = "Info: Calculating static properties"
  call write_info(1)

  ALLOCATE(tmp(1:np, 1), np)
  ALLOCATE(dVde(1:np, 1:sys%st%d%nspin, 1:dim), np*sys%st%d%nspin*dim)
  ALLOCATE(drhode(1:np, 1:dim), np*dim)
  ALLOCATE(kxc(1:np, 1:sys%st%d%nspin, 1:sys%st%d%nspin, 1:sys%st%d%nspin), np*sys%st%d%nspin**3)

  call xc_get_kxc(sys%ks%xc, sys%gr%m, sys%st%rho, sys%st%d%ispin, kxc)

  ! first the derivatives in all directions are calculated and stored
  write(message(1), '(a)') 'Info: Calculating derivatives of the orbitals:'
  call write_info(1)

  do i = 1, dim
    write(message(1), '(a,i1)') 'Info: Derivative direction: ', i
    call write_info(1)

    call mix_init(lr(i, 1)%mixer, sys%gr%m, 1, sys%st%d%nspin)

    call X(get_response_e)(sys, h, lr(:,:), i, 1, R_TOTYPE(M_ZERO), props)
    call mix_end(lr(i, 1)%mixer)

    do ispin = 1, sys%st%d%nspin
      ! the potential derivatives

      ! Hartree and the potential and the derivative of the
      !  potential associated at the electric field

      dVde(1:np,ispin, i) = sys%gr%m%x(1:np, i)

      if(props%add_hartree) dVde(1:np, ispin, i) = dVde(1:np, ispin, i) + lr(i,1)%X(dl_Vhar)(1:np) 

      if(props%add_fxc) then 
        ! xc
        do ispin2 = 1, sys%st%d%nspin
          dVde(1:np, ispin,i) = dVde(1:np, ispin,i) + &
               lr(i, 1)%dl_Vxc(1:np, ispin, ispin2)*lr(i,1)%X(dl_rho)(1:np, ispin2)
        end do
      end if

    end do

    ! the density
    do n = 1, np
      drhode(n, i) = sum(lr(i, 1)%X(dl_rho)(n, 1:sys%st%d%nspin))
    end do

  end do !dim

  write(message(1), '(a)') 'Info: Calculating polarizability tensor'
  call write_info(1)

  pol = M_ZERO
  do i = 1, dim
    do j = 1, np
      pol(i, 1:dim) = pol(i, 1:dim) - &
           sys%gr%m%x(j, 1:dim) * drhode(j, i) * sys%gr%m%vol_pp(j)
    end do
  end do

  write(message(1), '(a)') 'Info: Calculating hyperpolarizability tensor'
  call write_info(1)

  if (sys%st%d%nspin /= UNPOLARIZED ) then 
    write(message(1), '(a)') 'WARNING: Hyperpolarizability has not been tested for spin polarized systems'
    call write_warning(1)
  end if

  hpol_density = M_ZERO

  do i = 1, dim
    do j = 1, dim
      do k = 1, dim

        do ik=1, sys%st%d%nik
          do ispin = 1, sys%st%d%nspin
            do ist = 1, sys%st%nst

              if( sys%st%occ(ist, ik) > lr_min_occ ) then 
                ! <D\psi_n | P_c DV_scf P_c | D\psi_n >

                tmp(1:np, 1)  = dVde(1:np, ispin, j) * lr(k,1)%X(dl_psi)(1:np, 1, ist, ispin)
                hpol_density(1:np,i, j, k) = hpol_density(1:np, i, j, k) + &
                     sys%st%d%kweights(ik)*sys%st%occ(ist, ik)* &
                     R_CONJ(lr(i, 1)%X(dl_psi)(1:np, 1, ist, ispin)) * tmp(1:np, 1)

                do ispin2 = 1, sys%st%d%nspin
                  do ist2 = 1, sys%st%nst
                    if( sys%st%occ(ist2, ik) > lr_min_occ ) then 

                      tmp(1:np, 1)=R_CONJ(sys%st%X(psi)(1:np, 1, ist2, ispin2)) * &
                           dVde(1:np, ispin, j) * sys%st%X(psi)(1:np, 1, ist, ispin)
                      prod = dmf_integrate(sys%gr%m,tmp(1:np,1))

                      hpol_density(1:np,i, j, k) = hpol_density(1:np,i, j, k) - & 
                           sys%st%d%kweights(ik)*sys%st%occ(ist, ik)* & 
                           R_CONJ(lr(i, 1)%X(dl_psi)(1:np, 1, ist, ispin)) * &
                           lr(k, 1)%X(dl_psi)(1:np, 1, ist2, ispin2)*prod

                    end if
                  end do ! ist2
                end do ! ispin2

              end if

            end do ! ist
          end do ! ispin
        end do !ik

        if(props%add_fxc) then 
          hpol_density(1:np,i, j, k) = hpol_density(1:np,i, j, k) + &
               kxc(1:np, 1, 1, 1) * drhode(1:np, i) * drhode(1:np, j)*drhode(1:np, k)/CNST(6.0)
        end if

      end do ! k
    end do ! j
  end do ! i

  ALLOCATE(hp_tmp(1:np,MAX_DIM,MAX_DIM,MAX_DIM),np*MAX_DIM**3)
  hp_tmp(1:np,1:dim,1:dim,1:dim)=hpol_density(1:np,1:dim,1:dim,1:dim)
  do k = 1, dim
    do j = 1, dim
      do i = 1, dim
        hpol_density(1:np,i,j,k) = -( &
             + hp_tmp(1:np,i,j,k) + hp_tmp(1:np,j,k,i) &
             + hp_tmp(1:np,k,i,j) + hp_tmp(1:np,k,j,i) &
             + hp_tmp(1:np,j,i,k) + hp_tmp(1:np,i,k,j))
        hpol(i,j,k)=dmf_integrate(sys%gr%m,hpol_density(1:np,i,j,k))
      end do ! k
    end do ! j
  end do ! i
  deallocate(hp_tmp)


  deallocate(tmp)
  deallocate(dVde)
  deallocate(drhode)
  deallocate(kxc)

  call pop_sub()
end subroutine X(static_response)

subroutine X(lr_calc_elf)(st, gr, lr, lr_m)
  type(states_t),   intent(inout) :: st
  type(grid_t),     intent(inout) :: gr
  type(lr_t),       intent(inout) :: lr
  type(lr_t), optional, intent(inout) :: lr_m !when this argument is present, we are doing dynamical response

  integer :: i, is, ist, idim, ik

  R_TYPE, allocatable :: gpsi(:,:), gdl_psi(:,:), gdl_psi_m(:,:)
  FLOAT,  allocatable :: rho(:), grho(:,:)
  R_TYPE, allocatable :: dl_rho(:), gdl_rho(:,:)
  FLOAT,  allocatable :: elf(:,:), de(:,:)
  FLOAT :: dl_d0, d0
  FLOAT :: f, s

  FLOAT, parameter :: dmin = CNST(1e-10)
  FLOAT :: u, ik_weight

  call push_sub('em_resp_inc.Xcalc_lr_elf')

  ALLOCATE(   gpsi(NP, NDIM), NP*NDIM)
  ALLOCATE(gdl_psi(NP, NDIM), NP*NDIM)

  if(present(lr_m)) ALLOCATE(gdl_psi_m(NP, NDIM), NP*NDIM)

  ALLOCATE(   grho(NP, NDIM), NP*NDIM)
  ALLOCATE(gdl_rho(NP, NDIM), NP*NDIM)

  ALLOCATE(   rho(NP_PART), NP)
  ALLOCATE(dl_rho(NP_PART), NP)

  ALLOCATE(   elf(NP, st%d%nspin), NP*st%d%nspin)
  ALLOCATE(    de(NP, st%d%nspin), NP*st%d%nspin)

  ALLOCATE(lr%X(dl_de)(NP, st%d%nspin), NP*st%d%nspin)  
  ALLOCATE(lr%X(dl_elf)(NP, st%d%nspin), NP*st%d%nspin)

  !calculate the gs elf
  call states_calc_elf(st, gr, elf, de)

  !calculate current and its variation
  if(st%d%wfs_type == M_CMPLX) then 
    call states_calc_physical_current(gr, st, st%j)
    if(present(lr_m)) then 
      call lr_calc_current(st, gr, lr, lr_m)
    else
      call lr_calc_current(st, gr, lr)
    end if
  end if

  ! single or double occupancy
  if(st%d%nspin == 1) then
    s = M_TWO
  else
    s = M_ONE
  end if

  lr%X(dl_de) = M_ZERO

  do is = 1, st%d%nspin
    rho  = M_ZERO
    grho = M_ZERO

    dl_rho  = M_ZERO
    gdl_rho = M_ZERO

    
    !first we calculate the denisities and its gradients, this could
    !be done directly, but it is less precise numerically
    do ik = is, st%d%nik, st%d%nspin
      do ist = 1, st%nst
        ik_weight = st%d%kweights(ik)*st%occ(ist, ik)/s

        do idim = 1, st%d%dim

          call X(f_gradient)(gr%sb, gr%f_der, st%X(psi)   (:, idim, ist, is), gpsi)
          call X(f_gradient)(gr%sb, gr%f_der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

          ! sum over states to obtain the spin-density
          rho(1:NP)    = rho(1:NP)    + ik_weight * abs(st%X(psi)(1:NP, idim, ist, is))**2

          !the gradient of the density
          do i = 1, NDIM
            grho(1:NP,i)    = grho(1:NP, i)   + ik_weight *  &
                 M_TWO * R_REAL(R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * gpsi(1:NP,i))
          end do

          !the variation of the density and its gradient

          if(present(lr_m)) then 

            call X(f_gradient)(gr%sb, gr%f_der, lr_m%X(dl_psi)(:, idim, ist, is), gdl_psi_m)

            dl_rho(1:NP) = dl_rho(1:NP) + ik_weight * ( &
                 R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * lr%X(dl_psi)(1:NP, idim, ist, is)+ & 
                 st%X(psi)(1:NP, idim, ist, is) * R_CONJ(lr_m%X(dl_psi)(1:NP, idim, ist, is)) )

            do i=1,NDIM

              gdl_rho(1:NP,i) = gdl_rho(1:NP,i) + ik_weight * ( &
                   R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * gdl_psi(1:NP,i) +      &
                   R_CONJ(gpsi(1:NP,i))* lr%X(dl_psi)(1:NP, idim, ist, is)  +      &
                   st%X(psi)(1:NP, idim, ist, is) * R_CONJ(gdl_psi_m(1:NP,i)) +      &
                   gpsi(1:NP,i) * R_CONJ(lr_m%X(dl_psi)(1:NP, idim, ist, is))  )

            end do

          else
            
            dl_rho(1:NP) = dl_rho(1:NP) + ik_weight * M_TWO *  &
                 R_REAL(R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * lr%X(dl_psi)(1:NP, idim, ist, is))

            do i=1,NDIM
              gdl_rho(1:NP,i) = gdl_rho(1:NP,i) + ik_weight * M_TWO * ( &
                   R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * gdl_psi(1:NP,i) +      &
                   gpsi(1:NP,i) * R_CONJ(lr%X(dl_psi)(1:NP, idim, ist, is))  )
            end do

          end if

        end do !idim
      end do !ist
    end do !ik

    !now we start to calculate the elf

    !first the term that depends on the orbitals
    !this is the only term that is different for the dynamical case
    do ik = is, st%d%nik, st%d%nspin
      do ist = 1, st%nst
        ik_weight = st%d%kweights(ik)*st%occ(ist, ik)/s
        do idim = 1, st%d%dim

          call X(f_gradient)(gr%sb, gr%f_der, st%X(psi)   (:, idim, ist, is), gpsi)
          call X(f_gradient)(gr%sb, gr%f_der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

          if(present(lr_m)) then 

            call X(f_gradient)(gr%sb, gr%f_der, lr_m%X(dl_psi)(:, idim, ist, is), gdl_psi_m)
            do i = 1, NP
              lr%X(dl_de)(i, is) = lr%X(dl_de)(i, is) +                             &
                   dl_rho(i) * ik_weight * sum(R_ABS(gpsi(i, 1:NDIM))**2) + &
                   rho(i)    * ik_weight * sum( & 
                   R_CONJ(gpsi(i,1:NDIM))*gdl_psi(i,1:NDIM) &
                   + gpsi(i,1:NDIM)*R_CONJ(gdl_psi_m(i,1:NDIM)) )
            end do

          else 

            do i = 1, NP
              lr%X(dl_de)(i, is) = lr%X(dl_de)(i, is) +                             &
                   dl_rho(i) * ik_weight * sum(R_ABS(gpsi(i, 1:NDIM))**2) + &
                   rho(i)    * ik_weight * M_TWO*(sum(R_CONJ(gpsi(i,1:NDIM))*gdl_psi(i,1:NDIM)))
            end do

          end if

        end do
      end do
    end do

    !the density term
    do i= 1, NP
      if(abs(st%rho(i, is)) >= dmin) then
        lr%X(dl_de)(i, is) = lr%X(dl_de)(i, is)                       &
             - M_HALF * sum(grho(i, 1:NDIM)*gdl_rho(i, 1:NDIM))
      end if
    end do

    !the current term
    if(st%d%wfs_type == M_CMPLX) then       
      do i= 1, NP
        if(abs(st%rho(i, is)) >= dmin) then
          lr%X(dl_de)(i, is) = lr%X(dl_de)(i, is)                       &
               +M_TWO*sum(st%j(i, 1:NDIM, is)*lr%dl_j(i, 1:NDIM, is))
        end if
      end do
    end if

    !now the normalization 
    f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
    do i = 1, NP

      if(abs(st%rho(i, is)) >= dmin) then

        d0    = f * rho(i)**(M_EIGHT/M_THREE)
        dl_d0 = M_EIGHT/M_THREE * f * dl_rho(i) * rho(i)**(M_FIVE/M_THREE)

        lr%X(dl_elf)(i,is) = M_TWO*d0*dl_d0/(d0**2+de(i,is)**2)*(1-elf(i,is))& 
             -M_TWO*de(i,is)*lr%X(dl_de)(i,is)/(d0**2+de(i,is)**2)*elf(i,is)

      else

        lr%X(dl_elf)(i, is) = M_ZERO

      end if

    end do

  end do

  deallocate(gpsi)
  deallocate(gdl_psi)
  if(present(lr_m)) deallocate(gdl_psi_m)

  deallocate(grho)
  deallocate(gdl_rho)

  deallocate(elf)
  deallocate(de)

  call pop_sub()

end subroutine X(lr_calc_elf)

