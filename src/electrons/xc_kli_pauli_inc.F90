!! Copyright (C) 2012-2013 M. Gruning, P. Melo, M. Oliveira
!! Copyright (C) 2021 N. Tancogne-Dejean
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! ---------------------------------------------------------
subroutine xc_kli_pauli_solve(mesh, st, oep)
  type(mesh_t),        intent(in)    :: mesh
  type(states_elec_t), intent(in) :: st
  type(xc_oep_t),      intent(inout) :: oep
  !
  integer :: ip, ist, jst, kssi, kssj, ik, proc, eigen_n
  CMPLX, allocatable :: psi(:,:), bij(:,:)
  FLOAT, allocatable :: sqphi(:, :, :), dd(:, :, :), v_bar_S(:), xx(:,:), yy(:,:), Ma(:,:)
  CMPLX :: tmp
  FLOAT :: nup, ndn, sqmod_updn, weight
  FLOAT :: rr, nn, mm, alpha, betar, betai, alpha2, beta2
  FLOAT :: global_b(4), local_b(4), local_v(4), global_v(4)


  call profiling_in(C_PROFILING_XC_KLI, TOSTRING(X(XC_KLI)))
  PUSH_SUB(xc_kli_pauli_solve)

  ASSERT(st%d%ispin == SPINORS)

  oep%vxc(:,:) = M_ZERO
  
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(bij(1:mesh%np, 1:3))
  bij(:,:) = M_ZERO
  SAFE_ALLOCATE(sqphi(1:mesh%np, 1:4, 1:st%nst))
  sqphi(:,:,:) = M_ZERO

  !We now construct the right-hand side of the equation defining the Slater potential
  do ik = st%d%kpt%start, st%d%kpt%end
    do ist = st%st_start, st%st_end
      weight = st%occ(ist, ik) * st%d%kweights(ik)  
      if(abs(weight) < M_EPSILON) cycle

      call states_elec_get_state(st, mesh, ist, ik, psi)
      !Here we accumulate the result for the potential
      do ip = 1, mesh%np
        bij(ip, 1) = bij(ip, 1) + weight * TOFLOAT(conjg(oep%zlxc(ip, ist, 1))*conjg(psi(ip, 1)))
        bij(ip, 2) = bij(ip, 2) + weight * TOFLOAT(conjg(oep%zlxc(ip, ist, 2))*conjg(psi(ip, 2)))
        bij(ip, 3) = bij(ip, 3) + weight * conjg(oep%zlxc(ip, ist, 1)) * conjg(psi(ip, 2))
        !The last component is simply the complex conjuguate of bij(3), so we do not compute it.

        ! We store \phi_{i,\alpha}(r)\phi^*_{i,\beta}(r). Needed for the KLI part
        sqphi(ip, 1, ist) = TOFLOAT(conjg(psi(ip, 1))*psi(ip, 1))
        sqphi(ip, 2, ist) = TOFLOAT(conjg(psi(ip, 2))*psi(ip, 2))
        tmp = conjg(psi(ip, 2))*psi(ip, 1)
        sqphi(ip, 3, ist) = TOFLOAT(tmp)
        sqphi(ip, 4, ist) = aimag(tmp)
      end do
    end do
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states .or. st%d%kpt%parallel) then
    call comm_allreduce(st%st_kpt_mpi_grp, bij)
    call comm_allreduce(st%st_kpt_mpi_grp, sqphi)
  end if
#endif

  !We now construct the Slater potential from the density and the quantity b_ij
  !Note that the final potential is real and stored as v_upup, v_downdown, Re(v_updown), and Im(updown)
  !See Eq. 22 in SI of PRB 98, 035140 (2018)
  do ip = 1, mesh%np

    nn = st%rho(ip, 1) + st%rho(ip, 2) + M_EPSILON
    sqmod_updn = st%rho(ip, 3)**2 + st%rho(ip, 4)**2

    ! 1/(2nD), where n is the charge density and D = n_uu*n_dd - n_ud*n_du
    ! where D is the determinant of the spin-density matrix
    rr = (st%rho(ip, 1) * st%rho(ip, 2)  - sqmod_updn)
    if(abs(rr) < CNST(1.0e-13) ) then !The matrix is singular
      !For determining the potential, we go to the local frame given by the local magnetization
      !In this frame, the spin density matrix has a single non-zero element on the diagonal (because the matrix is singular)
      !This allows us to determine the potential in this frame, and to rotate it back in the original frame
      !We start by finding the non-zero component in the local rotated frame
      mm = sqrt((st%rho(ip, 1) - st%rho(ip, 2))**2 + M_FOUR*sqmod_updn)
      if((st%rho(ip, 1) - st%rho(ip, 2)) > 0) then
        nup = max(M_HALF*(nn + mm), M_ZERO)
        ndn = max(M_HALF*(nn - mm), M_ZERO)
      else
        nup = max(M_HALF*(nn - mm), M_ZERO)
        ndn = max(M_HALF*(nn + mm), M_ZERO)
      end if

      !Coefficients of the rotation matrix
      call get_rotation_matrix(st%rho(ip, :), alpha, betar, betai)
      alpha2 = alpha**2
      beta2 = betar**2 + betai**2

      !We rotate to the local frame
      !We add here the factor of two that is originally in the Sylvester equation
      !In the non singular case, this is cancelled by a factor one half
      global_b(1) = M_TWO*TOFLOAT(bij(ip, 1))
      global_b(2) = M_TWO*TOFLOAT(bij(ip, 2))
      global_b(3) = M_TWO*TOFLOAT(bij(ip, 3))
      global_b(4) = M_TWO*aimag(bij(ip, 3))
      call rotate_to_local(global_b, alpha, betar, betai, alpha2, beta2, local_b)
     
      local_v(1:4) = M_ZERO
      if(nup > ndn) then !We are up
        local_v(1) = local_b(1) / (M_TWO * nup)
        local_v(3) = local_b(3) / nup
        local_v(4) = local_b(4) / nup
      else !We are down
        local_v(2) = local_b(2) / (M_TWO * ndn)
        local_v(3) = local_b(3) / ndn
        local_v(4) = local_b(4) / ndn
      end if

      !We rotate to the original frame and we accumulate the result
      call rotate_to_global(local_v, alpha, betar, betai, alpha2, beta2, global_v)
      oep%vxc(ip, 1:4) = oep%vxc(ip, 1:4) + global_v(1:4)   

    else

      rr = M_ONE/(nn * rr)

      oep%vxc(ip, 1) = oep%vxc(ip, 1) + rr * ( &
         (nn * st%rho(ip, 2) - sqmod_updn) * TOFLOAT(bij(ip, 1)) + sqmod_updn * TOFLOAT(bij(ip, 2)) &
        - M_TWO * st%rho(ip,2) * ( st%rho(ip,3) * TOFLOAT(bij(ip,3)) + st%rho(ip,4) * aimag(bij(ip,3))))

      oep%vxc(ip, 2) = oep%vxc(ip, 2) + rr * ( &
         (nn * st%rho(ip, 1) - sqmod_updn) * TOFLOAT(bij(ip, 2)) + sqmod_updn * TOFLOAT(bij(ip, 1)) &
        - M_TWO * st%rho(ip,1) * ( st%rho(ip,3) * TOFLOAT(bij(ip,3)) + st%rho(ip,4) * aimag(bij(ip,3)))) 

      tmp = -TOCMPLX(st%rho(ip, 3), st%rho(ip,4)) * (st%rho(ip, 2) * bij(ip, 1) + st%rho(ip, 1) * bij(ip,2)) &
            + (M_TWO *st%rho(ip, 1) * st%rho(ip, 2)  - sqmod_updn) * bij(ip, 3) &
            + (TOCMPLX(st%rho(ip, 3),st%rho(ip,4)))**2 * conjg(bij(ip,3))

      oep%vxc(ip, 3) = oep%vxc(ip, 3) + rr * TOFLOAT(tmp)
      oep%vxc(ip, 4) = oep%vxc(ip, 4) + rr * aimag(tmp)
    end if
  end do 
  
  SAFE_DEALLOCATE_A(psi)


  
  ! If there is more than one state, then solve linear equation.
  eigen_n = oep%eigen_n
  if(eigen_n > 0) then

    SAFE_ALLOCATE(v_bar_S(1:st%nst))
    v_bar_S = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
      weight = st%occ(ist, ik) * st%d%kweights(ik)
  
      v_bar_S(ist) = dmf_dotp(mesh, 2, sqphi(:, 1:2, ist), oep%vxc(:,1:2), reduce = .false.)
      v_bar_S(ist) = v_bar_S(ist) + M_TWO * dmf_dotp(mesh, 2, sqphi(:, 3:4, ist), oep%vxc(:,3:4), reduce = .false.)
      end do
    end do
    if(mesh%parallel_in_domains) call mesh%allreduce(v_bar_S, dim = st%st_end)

#if defined(HAVE_MPI)
    ASSERT(.not. st%d%kpt%parallel) ! Not yet implemented here
    if(st%parallel_in_states) then
      ! Broadcast the vector v_bar_S  and sqphi to all processors
      do ist = 1, st%nst
        call MPI_Bcast(v_bar_S(ist), 1, MPI_FLOAT, st%node(ist), st%mpi_grp%comm, mpi_err)
      end do
      do ist = 1, eigen_n
        kssi = oep%eigen_index(ist)
        call MPI_Bcast(sqphi(1, 1, kssi), 4*mesh%np, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
      end do
    end if
#endif

    SAFE_ALLOCATE(dd(1:mesh%np, 1:4, 1:eigen_n))
    SAFE_ALLOCATE(xx(1:eigen_n, 1:1))
    SAFE_ALLOCATE(Ma(1:eigen_n, 1:eigen_n))
    SAFE_ALLOCATE(yy(1:eigen_n, 1:1))
    xx = M_ZERO
    yy = M_ZERO
    Ma = M_ZERO
    dd = M_ZERO
    proc = st%mpi_grp%rank

    do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      if(proc /= st%node(kssi)) cycle

      do ip = 1, mesh%np
        nn = st%rho(ip, 1) + st%rho(ip, 2) + M_EPSILON
        sqmod_updn = st%rho(ip, 3)**2 + st%rho(ip, 4)**2

        ! 1/(2nD), where n is the charge density and D = n_uu*n_dd - n_ud*n_du
        ! where D is the determinant of the spin-density matrix
        rr = (st%rho(ip, 1) * st%rho(ip, 2)  - sqmod_updn)
        if(abs(rr) < CNST(1.0e-13) ) then !The matrix is singular
          ! For determining the potential, we go to the local frame given by the local magnetization
          ! In this frame, the spin density matrix has a single non-zero element on the diagonal (because the matrix is singular)
          ! This allows us to determine the potential in this frame, and to rotate it back in the original frame
          ! We start by finding the non-zero component in the local rotated frame
          mm = sqrt((st%rho(ip, 1) - st%rho(ip, 2))**2 + M_FOUR*sqmod_updn)
          if((st%rho(ip, 1) - st%rho(ip, 2)) > 0) then
            nup = max(M_HALF*(nn + mm), M_ZERO)
            ndn = max(M_HALF*(nn - mm), M_ZERO)
          else
            nup = max(M_HALF*(nn - mm), M_ZERO)
            ndn = max(M_HALF*(nn + mm), M_ZERO)
          end if

          !Coefficients of the rotation matrix
          call get_rotation_matrix(st%rho(ip, :), alpha, betar, betai)
          alpha2 = alpha**2
          beta2 = betar**2 + betai**2

          ! We rotate to the local frame
          call rotate_to_local(sqphi(ip, 1:4, kssi), alpha, betar, betai, alpha2, beta2, local_b)
   
          local_v(1:4) = M_ZERO
          if(nup > ndn) then !We are up
            local_v(1) = local_b(1) / (M_TWO * nup)
            local_v(3) = local_b(3) / nup
            local_v(4) = local_b(4) / nup
          else !We are down
            local_v(2) = local_b(2) / (M_TWO * ndn)
            local_v(3) = local_b(3) / ndn
            local_v(4) = local_b(4) / ndn
          end if

          !We rotate to the original frame and we accumulate the result
          call rotate_to_global(local_v, alpha, betar, betai, alpha2, beta2, global_v)
          dd(ip, 1:4, ist) = global_v(1:4)   

        else
          rr = M_HALF/(nn * rr)

          dd(ip, 1, ist) = rr * ((nn * st%rho(ip, 2) - sqmod_updn) * sqphi(ip, 1, kssi) + sqmod_updn * sqphi(ip, 2, kssi) &
        - M_TWO * st%rho(ip,2) * ( st%rho(ip,3) * sqphi(ip, 3, kssi) + st%rho(ip,4) * sqphi(ip, 4, kssi)))

          dd(ip, 2, ist) = rr * ((nn * st%rho(ip, 1) - sqmod_updn) * sqphi(ip, 2, kssi) + sqmod_updn * sqphi(ip, 1, kssi) &
        - M_TWO * st%rho(ip,1) * ( st%rho(ip,3) * sqphi(ip, 3, kssi) + st%rho(ip,4) * sqphi(ip, 4, kssi))) 

          tmp = -TOCMPLX(st%rho(ip, 3), st%rho(ip,4)) * (st%rho(ip, 2) * sqphi(ip, 1, kssi) + st%rho(ip, 1) * sqphi(ip, 2, kssi)) &
            + (M_TWO *st%rho(ip, 1) * st%rho(ip, 2)  - sqmod_updn) * TOCMPLX(sqphi(ip, 3, kssi), sqphi(ip, 4, kssi)) &
            + (TOCMPLX(st%rho(ip, 3),st%rho(ip,4)))**2 * TOCMPLX(sqphi(ip, 3, kssi), -sqphi(ip, 4, kssi))

          dd(ip, 3, ist) = rr * TOFLOAT(tmp)
          dd(ip, 4, ist) = rr * aimag(tmp)
        end if 
      end do
    end do 
   
    i_loop: do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      if(proc /= st%node(kssi)) cycle

      j_loop: do jst = ist, eigen_n
        kssj = oep%eigen_index(jst)
        Ma(ist, jst) = - M_TWO * dmf_dotp(mesh, 2, dd(:,1:2, jst), sqphi(:, 1:2, kssi) )
        Ma(ist, jst) = Ma(ist, jst) - M_FOUR * dmf_dotp(mesh, 2, dd(:,3:4, jst), sqphi(:, 3:4, kssi))
      end do j_loop
      Ma(ist, ist) = M_ONE + Ma(ist, ist)
      yy(ist, 1) = M_TWO * v_bar_S(kssi) - M_TWO * (oep%uxc_bar(kssi, 1) + oep%uxc_bar(kssi, 2)) 

    end do i_loop

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      do ist = 1, eigen_n
        kssi = oep%eigen_index(ist)
        call MPI_Bcast(yy(ist, 1), 1, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
        do jst = 1, eigen_n
           call MPI_Bcast(Ma(ist, jst), 1, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
        end do
      end do
    end if
#endif

    do ist = 1, eigen_n
      do jst = ist+1, eigen_n
        Ma(jst, ist) = Ma(ist, jst)
      end do
    end do

    call lalg_linsyssolve(eigen_n, 1, Ma, yy, xx)

    do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      call lalg_axpy(mesh%np, 4, st%occ(kssi, 1) * xx(ist, 1), dd(:, 1:4, ist), oep%vxc(:, 1:4))
    end do

    SAFE_DEALLOCATE_A(dd)
    SAFE_DEALLOCATE_A(xx)
    SAFE_DEALLOCATE_A(Ma)
    SAFE_DEALLOCATE_A(yy)
    SAFE_DEALLOCATE_A(v_bar_S)

  end if

  SAFE_DEALLOCATE_A(bij)

  call profiling_out(C_PROFILING_XC_KLI)

  POP_SUB(xc_kli_pauli_solve)
end subroutine xc_kli_pauli_solve

  subroutine get_rotation_matrix(dens, alpha, betar, betai)
    FLOAT,  intent(in)  :: dens(:)
    FLOAT,  intent(out) :: alpha, betar, betai

    FLOAT :: mz, mm

    mz = dens(1) - dens(2) 

    mm = sqrt(mz**2 + M_FOUR*(dens(3)**2 + dens(4)**2))

    !Fully spin unpolarized system
    if(mm < CNST(1.0e-12)) then
      alpha = M_ONE
      betar = M_ZERO
      betai = M_ZERO
      return
    end if

    alpha = sqrt( (mm + abs(mz))/(M_TWO * mm) )
    !We find the absolute values of real and imaginary parts of beta
    betar = M_TWO * dens(3) / sqrt( M_TWO * mm * (mm + abs(mz)))
    betai = M_TWO * dens(4) / sqrt( M_TWO * mm * (mm + abs(mz)))

    if(mz < M_ZERO) then
      betar = -betar
      betai = -betai
    end if

  end subroutine get_rotation_matrix

  !Given a matrix in spin space, this routine rotates is according to the rotation
  !matrix R defined by the alpha and beta coefficients
  !rotmat = R mat R^T
  subroutine rotate_to_local(mat, alpha, betar, betai, alpha2, beta2, rot_mat)
    FLOAT,  intent(in)  :: mat(:)
    FLOAT,  intent(in)  :: alpha, betar, betai, alpha2, beta2
    FLOAT,  intent(out) :: rot_mat(:)

    CMPLX :: cross

    rot_mat(1) = alpha2 * mat(1) + beta2 * mat(2) + M_TWO * alpha * (betar * mat(3) + betai * mat(4))
    rot_mat(2) = alpha2 * mat(2) + beta2 * mat(1) - M_TWO * alpha * (betar * mat(3) + betai * mat(4))
    cross = (TOCMPLX(betar, betai))**2 * TOCMPLX(mat(3), -mat(4))
    rot_mat(3) = alpha2 * mat(3) + alpha * betar * (mat(2)-mat(1)) - real(cross)
    rot_mat(4) = alpha2 * mat(4) + alpha * betai * (mat(2)-mat(1)) - aimag(cross)

  end subroutine rotate_to_local

  !Given a matrix in spin space, this routine rotates is according to the rotation
  !matrix R defined by the alpha and beta coefficients
  !rotmat = R^T mat R
  subroutine rotate_to_global(mat, alpha, betar, betai, alpha2, beta2, rot_mat)
    FLOAT,  intent(in)  :: mat(:)
    FLOAT,  intent(in)  :: alpha, betar, betai, alpha2, beta2
    FLOAT,  intent(out) :: rot_mat(:)

    CMPLX :: cross

    rot_mat(1) = alpha2 * mat(1) + beta2 * mat(2) - M_TWO * alpha * (betar * mat(3) + betai * mat(4))
    rot_mat(2) = alpha2 * mat(2) + beta2 * mat(1) + M_TWO * alpha * (betar * mat(3) + betai * mat(4))
    cross = (TOCMPLX(betar, betai))**2 * TOCMPLX(mat(3), -mat(4))
    rot_mat(3) = alpha2 * mat(3) - alpha * betar * (mat(2)-mat(1)) - real(cross)
    rot_mat(4) = alpha2 * mat(4) - alpha * betai * (mat(2)-mat(1)) - aimag(cross)

  end subroutine rotate_to_global

