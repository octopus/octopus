!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id: em_resp.F90 2548 2006-11-06 21:42:27Z xavier $


! ---------------------------------------------------------
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
  FLOAT :: ik_weight

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

  if( .not. associated(lr%X(dl_de)) ) ALLOCATE(lr%X(dl_de)(NP, st%d%nspin), NP*st%d%nspin)  
  if( .not. associated(lr%X(dl_elf))) ALLOCATE(lr%X(dl_elf)(NP, st%d%nspin), NP*st%d%nspin)

  !calculate the gs elf
  call elf_calc(st, gr, elf, de)

  !calculate current and its variation
  if(st%d%wfs_type == M_CMPLX) then 
    call calc_physical_current(gr, st, st%j)
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

  deallocate(rho)
  deallocate(dl_rho)

  deallocate(grho)
  deallocate(gdl_rho)

  deallocate(elf)
  deallocate(de)

  call pop_sub()

end subroutine X(lr_calc_elf)


! ---------------------------------------------------------
subroutine X(lr_calc_polarizability)(sys, lr, zpol)
  type(system_t), target, intent(inout) :: sys
  type(lr_t),             intent(inout) :: lr(:,:)
  CMPLX,                  intent(out)   :: zpol(1:MAX_DIM, 1:MAX_DIM)

  integer :: dir1, dir2, j
  CMPLX :: zpol_tmp(1:MAX_DIM, 1:MAX_DIM)
  CMPLX, allocatable  :: dl_rho_tot(:), tmp(:)

  ALLOCATE(dl_rho_tot(1:sys%gr%m%np), sys%gr%m%np)
  ALLOCATE(tmp(1:sys%gr%m%np), sys%gr%m%np)

  zpol_tmp = M_ZERO
  
  do dir1 = 1, sys%gr%sb%dim

    ! sum dl_rho for all spin components
    do j = 1, sys%gr%m%np
      dl_rho_tot(j) = sum(lr(dir1, 1)%X(dl_rho)(j, 1:sys%st%d%nspin))
    end do
    
    do dir2 = 1, sys%gr%sb%dim
      tmp(1:sys%gr%m%np) = sys%gr%m%x(1:sys%gr%m%np, dir2) * dl_rho_tot(1:sys%gr%m%np)
      zpol_tmp(dir1, dir2) = zpol_tmp(dir1, dir2) - zmf_integrate(sys%gr%m, tmp)
    end do
    
  end do
  
  deallocate(dl_rho_tot, tmp)

  ! symmetrize
  do dir1 = 1, sys%gr%sb%dim
    do dir2 = 1, sys%gr%sb%dim
      zpol(dir1, dir2) = M_HALF*(zpol_tmp(dir1, dir2) + zpol_tmp(dir2, dir1))
    end do
  end do

end subroutine X(lr_calc_polarizability)


! ---------------------------------------------------------
subroutine X(lr_calc_beta) (sh, sys, h, lr, beta)
  type(sternheimer_t), intent(inout) :: sh
  type(system_t),      intent(inout) :: sys
  type(hamiltonian_t), intent(inout) :: h
  type(lr_t),          intent(inout) :: lr(:,:,:)
  CMPLX,               intent(out)   :: beta(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM)

  integer :: ifreq, isigma, idim, ispin, ispin2, np, ndim, idir, ist
  integer :: ii, jj, kk, iperm, op_sigma, idim2, ist2, n
  integer :: perm(1:3), u(1:3), w(1:3), ijk(1:3)

  integer, parameter :: ik = 1 !kpoints not supported

  R_TYPE :: prod

  R_TYPE, allocatable :: hvar(:, :, :, :, :, :)
  R_TYPE, allocatable :: tmp(:)
  FLOAT,  allocatable :: kxc(:,:,:,:)
  R_TYPE, allocatable :: hpol_density(:)

  np  = sys%gr%m%np
  ndim = sys%gr%sb%dim

  call push_sub('em_resp_inc.Xlr_calc_beta')

  write(message(1), '(a)') 'Info: Calculating hyperpolarizability tensor'
  call write_info(1)

  !calculate kxc, the derivative of fxc
  ALLOCATE(kxc(1:np, 1:sys%st%d%nspin, 1:sys%st%d%nspin, 1:sys%st%d%nspin), np*sys%st%d%nspin**3)
  kxc(:,:,:,:) = M_ZERO
  call xc_get_kxc(sys%ks%xc, sys%gr%m, sys%st%rho, sys%st%d%ispin, kxc)


  ALLOCATE(tmp(1:np), np)
  ALLOCATE(hvar(1:np, 1:sys%st%d%nspin, 1:2, 1:sys%st%d%dim, 1:ndim, 1:3), np*sys%st%d%nspin*2*sys%st%d%dim*ndim*3)
  ALLOCATE(hpol_density(1:np), np)


  do ifreq = 1, 3
    do idir = 1, ndim
      do idim = 1, sys%st%d%dim
        
        call X(sternheimer_calc_hvar)(sh, sys, h, lr(idir, :, ifreq), 2, sys%gr%m%x(:, idir), &
             hvar(:, :, :, idim, idir, ifreq))
        
      end do !idim
    end do !idir
  end do !ifreq

  beta(1:MAX_DIM, 1:MAX_DIM, 1:MAX_DIM) = M_ZERO


  do ii = 1, ndim
    do jj = 1, ndim
      do kk = 1, ndim

        ijk = (/ ii, jj, kk /)

        do isigma = 1, 2

          op_sigma = 2 
          if(isigma==2) op_sigma = 1

          do iperm = 1, 6

            call get_permutation(iperm, perm)

            u(1:3) = ijk(perm(1:3))
            w(1:3) = perm(1:3)

            hpol_density(1:np)=M_ZERO


            do ispin = 1, sys%st%d%nspin
              do ist = 1, sys%st%nst
                do idim = 1, sys%st%d%dim

                  if( sys%st%occ(ist, ik) > lr_min_occ ) then 

                    hpol_density(1:np) = hpol_density(1:np) &
                         + sys%st%d%kweights(ik)*sys%st%occ(ist, ik) &
                         * R_CONJ(lr(u(1), op_sigma, w(1))%X(dl_psi)(1:np, 1, ist, ispin)) &
                         * R_REAL(hvar(1:np, ispin, isigma, idim, u(2), w(2) )) &
                         * lr(u(3), isigma, w(3))%X(dl_psi)(1:np, idim, ist, ispin)

                    do ispin2 = 1, sys%st%d%nspin
                      do ist2 = 1, sys%st%nst
                        do idim2 = 1, sys%st%d%dim
                          if( sys%st%occ(ist2, ik) > lr_min_occ ) then 

                            tmp(1:np)= &
                                 R_CONJ(sys%st%X(psi)(1:np, idim2, ist2, ispin2)) * &
                                 R_REAL(hvar(1:np, ispin2, isigma, idim2, u(2), w(2))) * &
                                 sys%st%X(psi)(1:np, idim, ist, ispin)

                            prod = X(mf_integrate)(sys%gr%m, tmp(1:np))

                            hpol_density(1:np) = hpol_density(1:np) - & 
                                 sys%st%d%kweights(ik)*sys%st%occ(ist, ik)* & 
                                 R_CONJ(lr(u(1), op_sigma, w(1))%X(dl_psi)(1:np, 1, ist, ispin)) * &
                                 lr(u(3), isigma, w(3))%X(dl_psi)(1:np, idim2, ist2, ispin2)*prod

                          end if
                        end do !idim2
                      end do ! ist2
                    end do ! ispin2

                  end if

                end do !idim
              end do !ist
            end do !ispin

            if(sternheimer_add_fxc(sh)) then 
              do n = 1, np
                hpol_density(n) = hpol_density(n) &
                     + kxc(n, 1, 1, 1) &
                     * sum(lr(u(1), isigma, w(1))%X(dl_rho)(n, 1:sys%st%d%nspin)) &
                     * sum(lr(u(2), isigma, w(2))%X(dl_rho)(n, 1:sys%st%d%nspin)) &
                     * sum(lr(u(3), isigma, w(3))%X(dl_rho)(n, 1:sys%st%d%nspin)) &
                     / CNST(6.0)
              end do
            end if

            beta(ii, jj, kk)= beta(ii, jj, kk) - M_HALF * X(mf_integrate)(sys%gr%m, hpol_density(1:np))


          end do ! iperm

        end do !isigma

      end do !kk
    end do !jj
  end do !ii

  deallocate(hpol_density)
  deallocate(tmp)
  deallocate(hvar)

  call pop_sub()

contains

  subroutine get_permutation(i, p)
    integer, intent(in)  :: i
    integer, intent(out) :: p(1:3)

    ASSERT( i>=1 .and. i <= 6)

    select case(i)

    case(1) ; p(1)=1 ; p(2)=2 ; p(3)=3
    case(2) ; p(1)=2 ; p(2)=3 ; p(3)=1
    case(3) ; p(1)=3 ; p(2)=1 ; p(3)=2
    case(4) ; p(1)=3 ; p(2)=2 ; p(3)=1
    case(5) ; p(1)=1 ; p(2)=3 ; p(3)=2
    case(6) ; p(1)=2 ; p(2)=1 ; p(3)=3

    end select

  end subroutine get_permutation

end subroutine X(lr_calc_beta)
