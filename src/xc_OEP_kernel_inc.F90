

subroutine X(xc_oep_kernel_calc)(sys, h, lr, nsigma, fxcn)
  type(system_t),     intent(inout):: sys
  type(hamiltonian_t),intent(in)   :: h
  type(lr_t),         intent(in)   :: lr(:)
  integer,            intent(in)   :: nsigma
  R_TYPE,             intent(out)  :: fxcn(:)

  R_TYPE, allocatable :: r(:), s(:), t(:)
  
  FLOAT, parameter :: tol=1e-4
  integer :: np, ii, kk, ip, iis, kks
  
  np = sys%gr%m%np

  ALLOCATE(r(1:np), np)
  ALLOCATE(s(1:np), np)
  ALLOCATE(t(1:np), np)

  do ip = 1, np
    if ( abs(sys%st%rho(ip, 1)) > tol) then 
      fxcn(ip) = -lr(1)%X(dl_rho)(ip, 1)/sys%st%rho(ip, 1)
    end if
  end do

  fxcn(1:np) = fxcn(1:np) * h%vxc(1:np, 1)

  do ii = 1, sys%st%nst
    do iis = 1, nsigma
      do kk = 1, sys%st%nst
        do kks = 1, nsigma

          r(1:np) = R_CONJ(sys%st%X(psi)(1:np, 1, ii, 1))*&
               (sys%st%X(psi)(1:np, 1, kk, 1)+lr(kks)%X(dl_psi)(1:np, 1, kk, 1))

          call X(poisson_solve)(sys%gr, s, r, all_nodes=.false.)

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
  
  deallocate(r, s, t)

end subroutine X(xc_oep_kernel_calc)
