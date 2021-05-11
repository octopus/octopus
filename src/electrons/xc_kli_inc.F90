!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

! ---------------------------------------------------------
subroutine X(xc_KLI_solve) (mesh, hm, st, is, oep, first)
  type(mesh_t),             intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  integer,                  intent(in)    :: is
  type(xc_oep_t),           intent(inout) :: oep
  logical,                  intent(in)    :: first

  integer :: ist, ip, jst, eigen_n, kssi, kssj, proc
  FLOAT, allocatable :: rho_sigma(:), v_bar_S(:), sqphi(:, :, :), dd(:)
  FLOAT, allocatable :: Ma(:,:), xx(:,:), yy(:,:)
  R_TYPE, allocatable :: psi(:, :), bb(:,:)
  R_TYPE, allocatable :: phi1(:,:,:)

  ASSERT(.not. oep%has_photons) 
 
  call profiling_in(C_PROFILING_XC_KLI, TOSTRING(X(XC_KLI)))

  PUSH_SUB(X(xc_KLI_solve))
  ! some intermediate quantities
  ! vxc contains the Slater part!
  SAFE_ALLOCATE(rho_sigma(1:mesh%np))
  SAFE_ALLOCATE(sqphi(1:mesh%np, 1:st%d%dim, 1:st%nst))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))

  do ist = st%st_start, st%st_end
    call states_elec_get_state(st, mesh, ist, is, psi)
    sqphi(1:mesh%np, 1:st%d%dim, ist) = abs(psi(1:mesh%np, 1:st%d%dim))**2
  end do

  do ip = 1, mesh%np
    rho_sigma(ip) = max(sum(oep%socc * st%occ(st%st_start:st%st_end, is) * sqphi(ip, 1, st%st_start:st%st_end)), CNST(1e-20))
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    SAFE_ALLOCATE(dd(1:mesh%np))
    call MPI_Allreduce(rho_sigma(1), dd(1), mesh%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    rho_sigma(1:mesh%np) = dd(1:mesh%np)
  end if
#endif

  ! Comparing to KLI paper 1990, oep%vxc corresponds to V_{x \sigma}^S in Eq. 8
  ! The n_{i \sigma} in Eq. 8 is partitioned in this code into \psi^{*} (included in lxc) and \psi (explicitly below)

  oep%vxc(1:mesh%np, is) = CNST(0.0)

  do ist = st%st_start, st%st_end
    call states_elec_get_state(st, mesh, ist, is, psi)
    do ip = 1, mesh%np
      oep%vxc(ip, is) = oep%vxc(ip, is) + oep%socc*st%occ(ist, is)*R_REAL(oep%X(lxc)(ip, ist, is)*psi(ip, 1))
    end do
  end do

  do ip = 1, mesh%np
    oep%vxc(ip, is) = oep%vxc(ip, is)/rho_sigma(ip)
  end do

  SAFE_DEALLOCATE_A(psi)
  
#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(oep%vxc(1,is), dd(1), mesh%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    oep%vxc(1:mesh%np,is) = dd(1:mesh%np)
    SAFE_DEALLOCATE_A(dd)
  end if
#endif
  
  eigen_n = oep%eigen_n

  SAFE_ALLOCATE(v_bar_S(1:st%nst))
  v_bar_S = M_ZERO
  do ist = st%st_start, st%st_end
    if(st%occ(ist, is) > M_EPSILON) then
      v_bar_S(ist) = dmf_dotp(mesh, sqphi(:, 1, ist) , oep%vxc(:,is), reduce = .false.)
    end if
  end do
  if(mesh%parallel_in_domains) call mesh%allreduce(v_bar_S, dim = st%st_end)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! Broadcast the vector v_bar_S  and sqphi to all processors
    do ist = 1, st%nst
      call MPI_Bcast(v_bar_S(ist), 1, MPI_FLOAT, st%node(ist), st%mpi_grp%comm, mpi_err)
    end do
    do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      call MPI_Bcast(sqphi(1, 1, kssi), mesh%np, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
    end do
  end if
#endif
  ! If there is more than one state, then solve linear equation.
  linear_equation: if(eigen_n > 0) then
    SAFE_ALLOCATE(dd(1:mesh%np))
    SAFE_ALLOCATE(xx(1:eigen_n, 1:1))
    SAFE_ALLOCATE(Ma(1:eigen_n, 1:eigen_n))
    SAFE_ALLOCATE(yy(1:eigen_n, 1:1))
    xx = M_ZERO
    yy = M_ZERO
    Ma = M_ZERO
    dd = M_ZERO
    proc = st%mpi_grp%rank

    i_loop: do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      if(proc  ==  st%node(kssi)) then
        dd(1:mesh%np) = sqphi(1:mesh%np, 1, kssi) / rho_sigma(1:mesh%np)
        j_loop: do jst = ist, eigen_n
          kssj = oep%eigen_index(jst)
          Ma(ist, jst) = - dmf_dotp(mesh, dd, sqphi(:, 1, kssj) )
        end do j_loop
        Ma(ist, ist) = M_ONE + Ma(ist, ist)
        yy(ist, 1) = v_bar_S(kssi) - oep%uxc_bar(kssi, is)
      end if
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
      do jst = ist, eigen_n
        Ma(jst, ist) = Ma(ist, jst)
      end do
    end do

    call lalg_linsyssolve(eigen_n, 1, Ma, yy, xx)

    do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      oep%vxc(1:mesh%np,is) = oep%vxc(1:mesh%np,is) + &
        oep%socc * st%occ(kssi, is) * xx(ist, 1) * sqphi(1:mesh%np, 1, kssi) / rho_sigma(1:mesh%np)
    end do

    SAFE_DEALLOCATE_A(dd)
    SAFE_DEALLOCATE_A(xx)
    SAFE_DEALLOCATE_A(Ma)
    SAFE_DEALLOCATE_A(yy)
    if (oep%has_photons) then
      SAFE_DEALLOCATE_A(phi1)
      SAFE_DEALLOCATE_A(bb)
    end if

  end if linear_equation
  ! The previous stuff is only needed if eigen_n>0.

  SAFE_DEALLOCATE_A(v_bar_S)
  SAFE_DEALLOCATE_A(rho_sigma)
  SAFE_DEALLOCATE_A(sqphi)
  POP_SUB(X(xc_KLI_solve))
  call profiling_out(C_PROFILING_XC_KLI)
end subroutine X(xc_KLI_solve)

! ---------------------------------------------------------
subroutine X(xc_KLI_solve_photon) (namespace, mesh, hm, st, is, oep, first)
  type(namespace_t),        intent(in)    :: namespace
  type(mesh_t),             intent(in)    :: mesh
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  integer,                  intent(in)    :: is
  type(xc_oep_t),           intent(inout) :: oep
  logical,                  intent(in)    :: first

  integer :: ist, ip, jst, eigen_n, kssi, kssj, proc
  FLOAT, allocatable :: rho_sigma(:), v_bar_S(:), sqphi(:, :, :), dd(:), coctranslation(:)
  FLOAT, allocatable :: Ma(:,:), xx(:,:), yy(:,:)
  R_TYPE, allocatable :: psi(:, :), bb(:,:)
  R_TYPE, allocatable :: phi1(:,:,:)
  
  call profiling_in(C_PROFILING_XC_KLI, TOSTRING(X(XC_KLI_PHOTON)))

  ASSERT(oep%has_photons)
  if(st%parallel_in_states) call messages_not_implemented("Photonic KLI not parallel in states")

  PUSH_SUB(X(xc_KLI_solve_photon))
  ! some intermediate quantities
  ! vxc contains the Slater part!
  SAFE_ALLOCATE(rho_sigma(1:mesh%np))
  SAFE_ALLOCATE(sqphi(1:mesh%np, 1:st%d%dim, 1:st%nst))
  SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))

#ifdef R_TREAL
  if (oep%coc_translation) then
    SAFE_ALLOCATE(coctranslation(1:mesh%np))
    coctranslation(1:mesh%np) = oep%pt%pol_dipole(1:mesh%np, 1)
    oep%pt%pol_dipole(1:mesh%np,1) = oep%pt%pol_dipole(1:mesh%np, 1) - &
      dmf_dotp(mesh, SUM(st%rho(1:mesh%np, :), dim=2),oep%pt%pol_dipole(1:mesh%np, 1))/abs(st%qtot)
  end if

  SAFE_ALLOCATE(phi1(1:mesh%np,1:st%d%dim,1:st%nst))
  SAFE_ALLOCATE(bb(1:mesh%np, 1:1))
  if (is == 1) oep%pt%ex = M_ZERO

  if(.not. lr_is_allocated(oep%photon_lr)) then
    call lr_allocate(oep%photon_lr, st, mesh)
    oep%photon_lr%ddl_psi(:,:, :, :) = M_ZERO
  end if

  if (.not. first) call xc_oep_pt_phi(namespace, mesh, hm, st, is, oep, phi1)
#else
  ! Photons with OEP are only implemented for real states
  ASSERT(.false.)
#endif

  do ist = st%st_start, st%st_end
    call states_elec_get_state(st, mesh, ist, is, psi)
    sqphi(1:mesh%np, 1:st%d%dim, ist) = abs(psi(1:mesh%np, 1:st%d%dim))**2
  end do

  do ip = 1, mesh%np
    rho_sigma(ip) = max(sum(oep%socc * st%occ(st%st_start:st%st_end, is) * sqphi(ip, 1, st%st_start:st%st_end)), CNST(1e-20))
  end do

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    SAFE_ALLOCATE(dd(1:mesh%np))
    call MPI_Allreduce(rho_sigma(1), dd(1), mesh%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    rho_sigma(1:mesh%np) = dd(1:mesh%np)
  end if
#endif

  ! Comparing to KLI paper 1990, oep%vxc corresponds to V_{x \sigma}^S in Eq. 8
  ! The n_{i \sigma} in Eq. 8 is partitioned in this code into \psi^{*} (included in lxc) and \psi (explicitly below)

  oep%vxc(1:mesh%np, 1) = CNST(0.0)

  do ist = st%st_start, st%st_end
    call states_elec_get_state(st, mesh, ist, is, psi)
#ifdef R_TREAL
    if (ist>(oep%eigen_n + 1)) exit ! included to guarantee that the photonic KLI finishes correctly but the parallel in states feature of the normal KLI works still
    bb(:,1) = oep%X(lxc)(1:mesh%np, ist, is)
    if (ist /= oep%eigen_n + 1) bb(:,1) = bb(:,1) - oep%uxc_bar(ist, is)*R_CONJ(psi(:, 1))
    if (.not.first) then
      call xc_oep_pt_rhs(mesh, st, is, oep, phi1, ist, bb)
    end if
    oep%vxc(:, is) = oep%vxc(:, is) + oep%socc*st%occ(ist, is)*bb(:, 1)*psi(:, 1)
#else
    ! Photons with OEP are only implemented for real states
    ASSERT(.false.)
#endif
  end do

  if (oep%coc_translation) then
    oep%pt%pol_dipole(1:mesh%np, 1) = oep%pt%pol_dipole(1:mesh%np,1 ) + &
        dmf_dotp(mesh, sum(st%rho(1:mesh%np, :), dim=2), coctranslation(1:mesh%np))/abs(st%qtot)
    SAFE_DEALLOCATE_A(coctranslation)
  end if

  do ip = 1, mesh%np
    oep%vxc(ip, is) = oep%vxc(ip, is)/rho_sigma(ip)
  end do

  SAFE_DEALLOCATE_A(psi)
  
#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    call MPI_Allreduce(oep%vxc(1,is), dd(1), mesh%np, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    oep%vxc(1:mesh%np,is) = dd(1:mesh%np)
    SAFE_DEALLOCATE_A(dd)
  end if
#endif
  
  eigen_n = oep%eigen_n

  SAFE_ALLOCATE(v_bar_S(1:st%nst))
  v_bar_S = M_ZERO
  do ist = st%st_start, st%st_end
    if(st%occ(ist, is) > M_EPSILON) then
      v_bar_S(ist) = dmf_dotp(mesh, sqphi(:, 1, ist) , oep%vxc(:,is), reduce = .false.)
    end if
  end do
  if(mesh%parallel_in_domains) call mesh%allreduce(v_bar_S, dim = st%st_end)

#if defined(HAVE_MPI)
  if(st%parallel_in_states) then
    ! Broadcast the vector v_bar_S  and sqphi to all processors
    do ist = 1, st%nst
      call MPI_Bcast(v_bar_S(ist), 1, MPI_FLOAT, st%node(ist), st%mpi_grp%comm, mpi_err)
    end do
    do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      call MPI_Bcast(sqphi(1, 1, kssi), mesh%np, MPI_FLOAT, st%node(kssi), st%mpi_grp%comm, mpi_err)
    end do
  end if
#endif
  ! If there is more than one state, then solve linear equation.
  linear_equation: if(eigen_n > 0) then
    SAFE_ALLOCATE(dd(1:mesh%np))
    SAFE_ALLOCATE(xx(1:eigen_n, 1:1))
    SAFE_ALLOCATE(Ma(1:eigen_n, 1:eigen_n))
    SAFE_ALLOCATE(yy(1:eigen_n, 1:1))
    xx = M_ZERO
    yy = M_ZERO
    Ma = M_ZERO
    dd = M_ZERO
    proc = st%mpi_grp%rank

    i_loop: do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      if(proc  ==  st%node(kssi)) then
        dd(1:mesh%np) = sqphi(1:mesh%np, 1, kssi) / rho_sigma(1:mesh%np)
        j_loop: do jst = ist, eigen_n
          kssj = oep%eigen_index(jst)
          Ma(ist, jst) = - dmf_dotp(mesh, dd, sqphi(:, 1, kssj) )
        end do j_loop
        Ma(ist, ist) = M_ONE + Ma(ist, ist)
        yy(ist, 1) = v_bar_S(kssi)
      end if
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
      do jst = ist, eigen_n
        Ma(jst, ist) = Ma(ist, jst)
      end do
    end do

    call lalg_linsyssolve(eigen_n, 1, Ma, yy, xx)

    do ist = 1, eigen_n
      kssi = oep%eigen_index(ist)
      oep%vxc(1:mesh%np,is) = oep%vxc(1:mesh%np,is) + &
        oep%socc * st%occ(kssi, is) * xx(ist, 1) * sqphi(1:mesh%np, 1, kssi) / rho_sigma(1:mesh%np)
    end do

    SAFE_DEALLOCATE_A(dd)
    SAFE_DEALLOCATE_A(xx)
    SAFE_DEALLOCATE_A(Ma)
    SAFE_DEALLOCATE_A(yy)
    if (oep%has_photons) then
      SAFE_DEALLOCATE_A(phi1)
      SAFE_DEALLOCATE_A(bb)
    end if

  end if linear_equation
  ! The previous stuff is only needed if eigen_n>0.

  SAFE_DEALLOCATE_A(v_bar_S)
  SAFE_DEALLOCATE_A(rho_sigma)
  SAFE_DEALLOCATE_A(sqphi)
  POP_SUB(X(xc_KLI_solve_photon))
  call profiling_out(C_PROFILING_XC_KLI)
end subroutine X(xc_KLI_solve_photon)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

