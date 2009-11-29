

subroutine X(xc_oep_kernel_calc)(sys, hm, lr, nsigma, fxcn)
  type(system_t),     intent(inout):: sys
  type(hamiltonian_t),intent(in)   :: hm
  type(lr_t),         intent(in)   :: lr(:)
  integer,            intent(in)   :: nsigma
  R_TYPE,             intent(out)  :: fxcn(:)

  R_TYPE, allocatable :: r(:), s(:), t(:)
  
  FLOAT, parameter :: tol=1e-4
  integer :: np, ii, kk, ip, iis, kks
  
  np = sys%gr%mesh%np

  SAFE_ALLOCATE(r(1:np))
  SAFE_ALLOCATE(s(1:np))
  SAFE_ALLOCATE(t(1:np))

  do ip = 1, np
    if ( abs(sys%st%rho(ip, 1)) > tol) then 
      fxcn(ip) = -lr(1)%X(dl_rho)(ip, 1)/sys%st%rho(ip, 1)
    end if
  end do

  fxcn(1:np) = fxcn(1:np) * hm%vxc(1:np, 1)

  do ii = 1, sys%st%nst
    do iis = 1, nsigma
      do kk = 1, sys%st%nst
        do kks = 1, nsigma

          r(1:np) = R_CONJ(sys%st%X(psi)(1:np, 1, ii, 1))*&
               (sys%st%X(psi)(1:np, 1, kk, 1)+lr(kks)%X(dl_psi)(1:np, 1, kk, 1))

          call X(poisson_solve)(sys%gr%der, s, r, all_nodes=.false.)

          do ip = 1, np
            if ( abs(sys%st%rho(ip, 1)) > tol) then 
              fxcn(ip) = fxcn(ip) - M_HALF*&
                   (R_CONJ(sys%st%X(psi)(ip, 1, kk, 1))*&
                   (sys%st%X(psi)(ip, 1, ii, 1)+lr(iis)%X(dl_psi)(ip, 1, ii, 1))*s(ip))&
                   /sys%st%rho(ip, 1)
            end if
          end do
        end do
      end do
    end do
  end do
  
  SAFE_DEALLOCATE_A(r)
  SAFE_DEALLOCATE_A(s)
  SAFE_DEALLOCATE_A(t)

end subroutine X(xc_oep_kernel_calc)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
