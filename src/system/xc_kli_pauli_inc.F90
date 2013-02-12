!! Copyright (C) 2012-2013 M. Gruning, P. Melo, M. Oliveira
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
! ---------------------------------------------------------
subroutine xc_kli_pauli_solve(mesh, st, oep)
  type(mesh_t),   intent(in)    :: mesh
  type(states_t), intent(in)    :: st
  type(xc_oep_t), intent(inout) :: oep
  !
  integer :: is, ip, ii, jj, ist, eigen_n, it, kssi 
  FLOAT, allocatable :: rho(:,:), lambda(:), n(:), t_rho(:,:)
  FLOAT, allocatable :: t_v(:,:), vloc(:,:), rhov(:), v_m1(:,:), p_i(:,:,:), t_vi(:,:,:), delta_v(:), vs(:,:)
  CMPLX, allocatable :: weighted_hf(:,:,:), rho_i(:,:,:,:)
  FLOAT :: reached_threshold(4)

  call profiling_in(C_PROFILING_XC_KLI)
  PUSH_SUB(xc_kli_pauli_solve)

  ! Density related quantities
  SAFE_ALLOCATE(rho(1:mesh%np, 4))
  SAFE_ALLOCATE(n(1:mesh%np))
  SAFE_ALLOCATE(lambda(1:mesh%np))
  rho(1:mesh%np, 1:4) = st%rho(1:mesh%np, 1:4)
  do ip = 1,mesh%np
    do is = 1, 2
      if (rho(ip,is) .lt. CNST(1e-20)) rho(ip,is) = CNST(1e-20)
    end do
  end do
  n(1:mesh%np) = rho(1:mesh%np,1) + rho(1:mesh%np,2)
  lambda(1:mesh%np) = rho(1:mesh%np,1)*rho(1:mesh%np,2) - (rho(1:mesh%np,3)**2 + rho(1:mesh%np,4)**2)

  ! Potential related quantities 
  ! (Built from HF potentials weighted with orbital densities)
  SAFE_ALLOCATE(weighted_hf(mesh%np,st%d%dim,st%d%dim))
  weighted_hf = M_Z0

  ! w_{up,down} = \sum_i \phi_{i,down}^* u_x^{i,up}^* \phi_{i,up}
  do ii = 1,st%d%dim
    do jj = 1,st%d%dim
      do ist = st%st_start,st%st_end
        weighted_hf(1:mesh%np,ii,jj) = weighted_hf(1:mesh%np,ii,jj) + &
             &oep%socc*st%occ(ist,1)*conjg(st%zpsi(1:mesh%np,jj,ist,1)*oep%zlxc(1:mesh%np,ist,ii)) ! oep%zlxc => (\phi_j)^*u_x^j 
      end do
    end do
  end do

  SAFE_ALLOCATE(t_v(mesh%np,4))
  t_v = M_ZERO
  t_v(1:mesh%np,1) = real(weighted_hf(1:mesh%np,2,2), REAL_PRECISION)
  t_v(1:mesh%np,2) = real(weighted_hf(1:mesh%np,1,1), REAL_PRECISION)
  t_v(1:mesh%np,3) = -real(weighted_hf(1:mesh%np,1,2) + weighted_hf(1:mesh%np,2,1), REAL_PRECISION)
  t_v(1:mesh%np,4) = -aimag(weighted_hf(1:mesh%np,1,2) - weighted_hf(1:mesh%np,2,1))
  SAFE_DEALLOCATE_A(weighted_hf)

  SAFE_ALLOCATE(vloc(mesh%np,4))
  vloc = M_ZERO
  vloc(1:mesh%np,1) = t_v(1:mesh%np,2) - t_v(1:mesh%np,1)
  vloc(1:mesh%np,2) = -vloc(1:mesh%np,1)
  vloc(1:mesh%np,3) = -t_v(1:mesh%np,3)
  vloc(1:mesh%np,4) = -t_v(1:mesh%np,4)

  ! Combine them to obtain Slater part
  SAFE_ALLOCATE(rhov(mesh%np))
  rhov = M_ZERO
  forall (ip = 1:mesh%np) rhov(ip) = sum(rho(ip,1:4)*t_v(ip,1:4))
  rhov(1:mesh%np) = rhov(1:mesh%np)/lambda(1:mesh%np)
  SAFE_DEALLOCATE_A(t_v)

  SAFE_ALLOCATE(t_rho(1:mesh%np, 4))
  t_rho(1:mesh%np,1) = rho(1:mesh%np,2)
  t_rho(1:mesh%np,2) = rho(1:mesh%np,1)
  t_rho(1:mesh%np,3:4) = -rho(1:mesh%np,3:4)
  forall (ip = 1:mesh%np) vloc(ip,1:4) = (vloc(ip,1:4) + t_rho(ip,1:4)*rhov(ip))/n(ip)


  select case (oep%level)
  case (XC_OEP_SLATER)
    
    oep%vxc = vloc

  case (XC_OEP_KLI)

    SAFE_ALLOCATE(vs(mesh%np,4))
    vs = vloc ! Slater part

    ! iteration criteria
    call scf_tol_init(oep%scftol,"KLI",def_maximumiter=50)

    ! get the HOMO state
    call xc_oep_AnalyzeEigen(oep, st, 1)
    eigen_n = oep%eigen_n
    if (eigen_n == 0) then 
      oep%vxc = vs

    else

      ! orbital densities
      SAFE_ALLOCATE(rho_i(mesh%np,st%d%dim,st%d%dim,eigen_n))
      rho_i = M_Z0
      do ii = 1, st%d%dim
        do jj = ii, st%d%dim
          do ist = 1, eigen_n
            kssi = oep%eigen_index(ist)
            rho_i(1:mesh%np,ii,jj,ist) = oep%socc*st%occ(kssi,1)*conjg(st%zpsi(1:mesh%np,jj,kssi,1))* &
                 st%zpsi(1:mesh%np,ii,kssi,1)
            rho_i(1:mesh%np,jj,ii,ist) = conjg(rho_i(1:mesh%np,ii,jj,ist))
          end do
        end do
      end do

      ! arrange them in a 4-vector
      SAFE_ALLOCATE(p_i(mesh%np,4,eigen_n))
      p_i = M_ZERO
      p_i(1:mesh%np,1,:) = real(rho_i(1:mesh%np,1,1,:)) 
      p_i(1:mesh%np,2,:) = real(rho_i(1:mesh%np,2,2,:))  
      p_i(1:mesh%np,3,:) = M_TWO*real(rho_i(1:mesh%np,1,2,:),REAL_PRECISION)  
      p_i(1:mesh%np,4,:) = M_TWO*aimag(rho_i(1:mesh%np,1,2,:))  
      
      SAFE_DEALLOCATE_A(rho_i)

      ! Calculate iteratively response part
      SAFE_ALLOCATE(v_m1(mesh%np,4))
      SAFE_ALLOCATE(delta_v(eigen_n)) 
      SAFE_ALLOCATE(t_vi(mesh%np,4,eigen_n)) 

      vloc = M_ZERO
      KLI_iteration: do it = 1,oep%scftol%max_iter
        v_m1 = vs + vloc

        ! delta_v^KLI
        delta_v = M_ZERO
        do ist=1,eigen_n
          kssi = oep%eigen_index(ist)
          do is = 1,st%d%nspin
            delta_v(ist) = delta_v(ist)+ dmf_dotp(mesh,p_i(1:mesh%np,is,ist),v_m1(1:mesh%np,is))
          end do
          delta_v(ist) = delta_v(ist) - real(sum(oep%uxc_bar(kssi,:)))
        end do

        !
        t_vi(1:mesh%np,1,:) = p_i(1:mesh%np,2,:) 
        t_vi(1:mesh%np,2,:) = p_i(1:mesh%np,1,:)
        t_vi(1:mesh%np,3,:) =-p_i(1:mesh%np,3,:) 
        t_vi(1:mesh%np,4,:) =-p_i(1:mesh%np,4,:)
        forall (ip=1:mesh%np,is=1:st%d%nspin) t_vi(ip,is,:) = t_vi(ip,is,:)*delta_v(:) 

        vloc = M_ZERO
        do ip = 1, mesh%np
          vloc(ip,1) = sum(t_vi(ip,2,:) - t_vi(ip,1,:)) 
          vloc(ip,2) = -vloc(ip,1)
          vloc(ip,3) = -sum(t_vi(ip,3,:))
          vloc(ip,4) = -sum(t_vi(ip,4,:))
        end do

        !
        rhov = M_ZERO
        do ip = 1, mesh%np
          do ist = 1, eigen_n
            rhov(ip) = rhov(ip) + sum(rho(ip,:)*t_vi(ip,:,ist))
          end do
        end do
        rhov = rhov/lambda

        forall (ip = 1:mesh%np) vloc(ip,:) = (vloc(ip,:) + t_rho(ip,:)*rhov(ip))/n(ip)
        !
        do is = 1, 4 
          reached_threshold(is) = dmf_nrm2(mesh,(vs(1:mesh%np,is) + vloc(1:mesh%np,is) - v_m1(1:mesh%np,is))) 
        end do
        if (all(reached_threshold(:) .le. oep%scftol%conv_abs_dens)) exit

      end do KLI_iteration

      write(message(1), '(a,i4,a,es14.6)') &
           &"Info: After ", it, " iterations, KLI converged to ", maxval(reached_threshold(:))
      message(2) = ''
      call messages_info(2)
      !
      oep%vxc = v_m1

      SAFE_DEALLOCATE_A(vs)
      SAFE_DEALLOCATE_A(v_m1)
      SAFE_DEALLOCATE_A(delta_v)
      SAFE_DEALLOCATE_A(t_vi)
      SAFE_DEALLOCATE_A(p_i)

    end if

  end select

  SAFE_DEALLOCATE_A(rho)
  SAFE_DEALLOCATE_A(lambda)
  SAFE_DEALLOCATE_A(n)
  SAFE_DEALLOCATE_A(t_rho)
  SAFE_DEALLOCATE_A(rhov)

  call profiling_out(C_PROFILING_XC_KLI)

  POP_SUB(xc_kli_pauli_solve)
end subroutine xc_KLI_Pauli_solve
