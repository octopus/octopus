!! Copyright (C) 2002-2016 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module stress_oct_m
  use boundaries_oct_m
  use comm_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use density_oct_m
  use derivatives_oct_m
  use epot_oct_m
  use fft_oct_m
  use fourier_space_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use lattice_vectors_oct_m
  use loct_math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use poisson_fft_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use projector_oct_m
  use ps_oct_m
  use simul_box_oct_m
  use species_oct_m
  use splines_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use submesh_oct_m
  use v_ks_oct_m
  use xc_f03_lib_m

  implicit none

  private
  public ::                    &
    stress_calculate

  integer, parameter ::             &
    CMD_FINISH = 1,                 &
    CMD_POISSON_SOLVE = 2


  FLOAT, allocatable, target :: rho(:, :)
  logical            :: total_density_alloc
  FLOAT, pointer     :: rho_total(:)
  CMPLX, allocatable :: rho_total_fs(:,:,:)
  FLOAT, allocatable :: FourPi_G2(:,:,:), Gvec(:,:,:,:), Gvec_G(:,:,:,:)

contains

  ! ---------------------------------------------------------
  !> This computes the total stress on the lattice
  subroutine stress_calculate(namespace, gr, hm, st, ions, ks)
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(inout) :: gr !< grid
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(states_elec_t),      intent(inout) :: st
    type(ions_t),             intent(inout) :: ions !< geometry
    type(v_ks_t),             intent(in)    :: ks

    type(profile_t), save :: stress_prof
    FLOAT :: stress(3,3) ! stress tensor in Cartecian coordinate
    FLOAT :: stress_KE(3,3), stress_Hartree(3,3), stress_xc(3,3) ! temporal
    FLOAT :: stress_ps(3,3), stress_Ewald(3,3)

    call profiling_in(stress_prof, "STRESS")
    PUSH_SUB(stress_calculate)

    SAFE_ALLOCATE(rho(1:gr%fine%mesh%np, 1:st%d%nspin))

    if(hm%kpoints%use_symmetries) then
      call messages_not_implemented("Stress tensors with k-point symmetries", namespace=namespace)
    end if

    if (geo%space%periodic_dim /= 3) then
      call messages_not_implemented("Stress tensors for periodicity different from 3D", namespace=namespace)
    end if

    if(.not.(ks%theory_level == KOHN_SHAM_DFT .and. bitand(hm%xc%family, XC_FAMILY_LDA) /= 0)) then
      write(message(1),'(a)') 'The stress tensor is currently only properly computed at the Kohn-Sham DFT at the LDA level'
      call messages_fatal(1, namespace=namespace)
    end if
    if(ks%vdw_correction /= OPTION__VDWCORRECTION__NONE) then
      write(message(1),'(a)') 'The stress tensor is currently not properly computed with vdW corrections'
      call messages_fatal(1, namespace=namespace)
    end if
  
    stress(:,:) = M_ZERO

    call calculate_density()
    call fourier_space_init(hm%psolver_fine)
    call density_rs2fs(hm%psolver_fine)
    
    ! Stress from kinetic energy of electrons    
    call stress_from_kinetic_energy_electron(gr%der, hm, st, stress, stress_KE)

    ! Stress from Hartree energy
    call stress_from_Hartree(hm, stress, stress_Hartree)

    ! Stress from exchange-correlation energy
    call stress_from_xc(gr%der, hm, stress, stress_xc)

    ! Stress from pseudopotentials
    call stress_from_pseudo(gr, hm, st, ions, stress, stress_ps)
    
    ! Stress from Ewald summation
    call stress_from_Ewald_sum(ions, stress, stress_Ewald)
    

    ! Stress from kinetic energy of ion
    ! Stress from ion-field interaction

    ! Sign changed to fit conventional definition    
    stress = - stress
    
    gr%sb%stress_tensor(1:3,1:3) = stress(1:3,1:3)

    SAFE_DEALLOCATE_A(FourPi_G2)
    SAFE_DEALLOCATE_A(Gvec)
    SAFE_DEALLOCATE_A(Gvec_G)
    SAFE_DEALLOCATE_A(rho)
    if (total_density_alloc) then
      SAFE_DEALLOCATE_P(rho_total)
    end if
    SAFE_DEALLOCATE_A(rho_total_fs)

    POP_SUB(stress_calculate)
    call profiling_out(stress_prof)

  contains
    subroutine calculate_density()
      integer :: ip
!      FLOAT                         :: amaldi_factor 

      PUSH_SUB(stress.calculate_density)

      ! get density taking into account non-linear core corrections
      call states_elec_total_density(st, gr%fine%mesh, rho)

      nullify(rho_total)

      if (allocated(st%rho_core) .or. hm%d%spin_channels > 1) then
         total_density_alloc = .true.
         
         SAFE_ALLOCATE(rho_total(1:gr%fine%mesh%np))
         
         do ip = 1, gr%fine%mesh%np
            rho_total(ip) = sum(rho(ip, 1:hm%d%spin_channels))
         end do
         
         ! remove non-local core corrections
         if (allocated(st%rho_core)) then
            do ip = 1, gr%fine%mesh%np
               rho_total(ip) = rho_total(ip) - st%rho_core(ip)
            end do
         end if
      else
         total_density_alloc = .false.
         rho_total => rho(:, 1)
      end if

      POP_SUB(stress.calculate_density)
    end subroutine calculate_density
    ! -------------------------------------------------------  

    ! ---------------------------------------------------------
    subroutine density_rs2fs(this)
      type(poisson_t), target,   intent(in)    :: this
      type(cube_t),    pointer             :: cube
      type(fourier_space_op_t), pointer    :: coulb
      type(cube_function_t) :: cf
      integer :: ii, jj, kk, iit, jjt, kkt
      FLOAT :: gx, xx(3)
      CMPLX :: zphase
      
      cube => this%cube
      coulb => this%fft_solver%coulb

      call dcube_function_alloc_RS(cube, cf, in_device = (this%fft_solver%kernel /= POISSON_FFT_KERNEL_CORRECTED))
      
      
      ! put the density in the cube
      if (cube%parallel_in_domains) then
         call dmesh_to_cube_parallel(this%der%mesh, rho_total, cube, cf, this%mesh_cube_map)
      else
         if(this%der%mesh%parallel_in_domains) then
            call dmesh_to_cube(this%der%mesh, rho_total, cube, cf, local = .true.)
         else
            call dmesh_to_cube(this%der%mesh, rho_total, cube, cf)
         end if
      end if
      
      ASSERT(allocated(cube%fft))
      ASSERT(cube%fft%library /= FFTLIB_NONE)
      
      SAFE_ALLOCATE(rho_total_fs(1:cube%rs_n_global(1),1:cube%rs_n_global(2),1:cube%rs_n_global(3)))
      
      call cube_function_alloc_fs(cube, cf)
      
      call dcube_function_rs2fs(cube, cf)
      cf%fs = cf%fs/TOFLOAT(cube%rs_n(1)*cube%rs_n(2)*cube%rs_n(3)) !Normalize

      select case(cube%fft%library)
      case(FFTLIB_PFFT)
! Not implemented yet
         write(message(1),'(a)') 'Internal error: PFFT library is not applicable for stress calculation.'
         call messages_fatal(1, namespace=namespace)
      case(FFTLIB_FFTW)
         if (allocated(cube%Lrs))then
            xx(1:3) = gr%sb%latt%red_to_cart(cube%Lrs(1,1:3))
         else
            xx(1:3) = -TOFLOAT(cube%rs_n_global(1:3)/2 )/TOFLOAT(cube%rs_n_global(1:3))
            xx(1:3) = gr%sb%latt%red_to_cart(xx(1:3))
         end if
         do kk = 1, cube%fs_n(3)
            kkt = - pad_feq(kk, cube%rs_n_global(3), .true.)
            kkt = mod(kkt+cube%rs_n_global(3),cube%rs_n_global(3)) +1
            do jj = 1, cube%fs_n(2)
               jjt = - pad_feq(jj , cube%rs_n_global(2), .true.)
               jjt = mod(jjt+cube%rs_n_global(2),cube%rs_n_global(2)) + 1
               do ii = 1, cube%fs_n(1)
                  iit = - pad_feq(ii , cube%rs_n_global(1), .true.)
                  iit = mod(iit+cube%rs_n_global(1),cube%rs_n_global(1)) + 1
                  
                  gx =  sum(xx(1:3)*Gvec(ii, jj, kk, 1:3) )
                  zphase = TOCMPLX(cos(gx), sin(gx)) 
                  rho_total_fs(ii,jj,kk) = conjg(cf%fs(ii, jj, kk))*zphase
                  rho_total_fs(iit,jjt,kkt) = cf%fs(ii, jj, kk)*conjg(zphase)
               end do
            end do
         end do
      case(FFTLIB_ACCEL)
! Not implemented yet
         write(message(1),'(a)') 'Internal error: ACCEL library is not applicable for stress calculation.'
         call messages_fatal(1, namespace=namespace)
      case default
         ASSERT(.false.)
       end select

       call cube_function_free_fs(cube, cf)
       call dcube_function_free_rs(cube, cf)
    
    end subroutine density_rs2fs

    ! -------------------------------------------------------
    subroutine fourier_space_init(this)
      type(poisson_t), target,   intent(in)    :: this
      type(cube_t),    pointer             :: cube
      type(fourier_space_op_t), pointer    :: coulb
      integer :: db(3)
      integer :: ix,iy,iz, ixx(3)
      FLOAT :: gg(3), modg2, temp(3)

      cube => this%cube
      coulb => this%fft_solver%coulb

      db(1:3) = cube%rs_n_global(1:3)

      SAFE_ALLOCATE(FourPi_G2(1:db(1),1:db(2),1:db(3)))
      SAFE_ALLOCATE(Gvec(1:db(1),1:db(2),1:db(3),3))
      SAFE_ALLOCATE(Gvec_G(1:db(1),1:db(2),1:db(3),3))
      
      
      
      if(cube%fft%library == FFTLIB_PFFT) then
!Not implemented yet
         
      else if(cube%fft%library == FFTLIB_FFTW) then

         temp(1:3) = M_TWO*M_PI/(db(1:3)*gr%fine%der%mesh%spacing(1:3))

         do ix = 1, cube%rs_n_global(1)
            ixx(1) = pad_feq(ix, db(1), .true.)
            do iy = 1, cube%rs_n_global(2)
               ixx(2) = pad_feq(iy, db(2), .true.)
               do iz = 1, cube%rs_n_global(3)
                  ixx(3) = pad_feq(iz, db(3), .true.)

                  call poisson_fft_gg_transform_l(ixx, temp, gr%fine%der%mesh%sb, gg, modg2)

                  !HH not very elegant
                  if(cube%fft%library.eq.FFTLIB_NFFT) modg2=cube%Lfs(ix,1)**2+cube%Lfs(iy,2)**2+cube%Lfs(iz,3)**2

                  if(abs(modg2) > M_EPSILON) then
                     FourPi_G2(ix,iy,iz) = 4d0*M_PI/modg2
                     Gvec_G(ix, iy, iz, 1) = gg(1)/sqrt(modg2)
                     Gvec_G(ix, iy, iz, 2) = gg(2)/sqrt(modg2)
                     Gvec_G(ix, iy, iz, 3) = gg(3)/sqrt(modg2)
                  else
                     FourPi_G2(ix,iy,iz) = M_ZERO
                     Gvec_G(ix, iy, iz, 1) = M_ZERO
                     Gvec_G(ix, iy, iz, 2) = M_ZERO
                     Gvec_G(ix, iy, iz, 3) = M_ZERO
                  end if

                  Gvec(ix, iy, iz, 1) = gg(1)
                  Gvec(ix, iy, iz, 2) = gg(2)
                  Gvec(ix, iy, iz, 3) = gg(3)

               end do
            end do
         end do

      else if(cube%fft%library == FFTLIB_ACCEL) then
!Not implemented yet
      end if
      
    end subroutine fourier_space_init
  end subroutine stress_calculate

  ! -------------------------------------------------------
  subroutine stress_from_kinetic_energy_electron(der, hm, st, stress, stress_KE)
    type(derivatives_t),  intent(in)    :: der
    type(hamiltonian_elec_t),  intent(in)    :: hm
    type(states_elec_t),  intent(inout) :: st
    FLOAT,                intent(inout) :: stress(:, :)
    FLOAT,                intent(out)   :: stress_KE(3, 3) ! temporal
    FLOAT                               :: stress_l(3, 3)
    integer :: ik, ist, idir, jdir, idim, ispin
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :)
    type(profile_t), save :: prof

    call profiling_in(prof, "STRESS_FROM_KEE")    
    PUSH_SUB(stress_from_kinetic_energy_electron)

    stress_l(:,:) = M_ZERO    

    SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:der%mesh%np, 1:der%dim, 1:st%d%dim))
    

    do ik = st%d%kpt%start, st%d%kpt%end
       ispin = st%d%get_spin_index(ik)
       do ist = st%st_start, st%st_end
          
          call states_elec_get_state(st, der%mesh, ist, ik, psi)
          
          do idim = 1, st%d%dim
             call boundaries_set(der%boundaries, psi(:, idim))
          end do
          
          if (allocated(hm%hm_base%phase)) then 
            call states_elec_set_phase(st%d, psi, hm%hm_base%phase(1:der%mesh%np_part, ik), der%mesh%np_part,.false.)  
          end if
          
          do idim = 1, st%d%dim
             call zderivatives_grad(der, psi(:, idim), gpsi(:, :, idim), set_bc = .false.)
          end do
          
          
          do idir = 1, der%dim
             do jdir = 1, der%dim
                
                do idim = 1, st%d%dim
                   stress_l(idir,jdir) = stress_l(idir,jdir) + &
                        st%d%kweights(ik)*st%occ(ist, ik) &
                        *dmf_integrate(der%mesh, TOFLOAT(conjg(gpsi(:, idir, idim))*gpsi(:, jdir, idim)))
                end do
             end do
          end do
       end do
       
    end do
    

    if(st%parallel_in_states .or. st%d%kpt%parallel) then
       ! TODO: this could take dim = (/der%mesh%np, der%dim, st%d%nspin/)) to reduce the amount of data copied
       call comm_allreduce(st%st_kpt_mpi_grp, stress_l) 
    end if

    stress_l = stress_l/der%mesh%sb%latt%rcell_volume
    stress_KE = stress_l
    stress = stress + stress_l
    

    SAFE_DEALLOCATE_A(psi)    
    SAFE_DEALLOCATE_A(gpsi)
    
    call profiling_out(prof)
    POP_SUB(stress_from_kinetic_energy_electron)
  end subroutine stress_from_kinetic_energy_electron
  
! -------------------------------------------------------
  subroutine stress_from_Hartree(hm, stress, stress_Hartree)
    type(hamiltonian_elec_t), intent(in)    :: hm
    FLOAT,                    intent(inout) :: stress(:, :)
    FLOAT,                    intent(out)   :: stress_Hartree(3, 3) ! temporal

    FLOAT :: stress_l(3, 3)
    type(cube_t), pointer :: cube
    integer :: idir, jdir, ii, jj, kk
    FLOAT :: ss
    type(profile_t), save :: prof

    cube => hm%psolver_fine%cube
    
    call profiling_in(prof, "STRESS_FROM_HARTREE")    
    PUSH_SUB(stress_from_Hartree)

    stress_l(:,:) = M_ZERO

    do idir = 1,3
       do jdir = 1,3
          ss=M_ZERO
          !$omp parallel do private(ii, jj, kk) reduction(+:ss)
          do kk = 1, cube%rs_n_global(3)
             do jj = 1, cube%rs_n_global(2)
                do ii = 1, cube%rs_n_global(1)
                   ss = ss + abs(rho_total_fs(ii,jj,kk))**2 &
                        *M_TWO*Gvec_G(ii, jj, kk,idir)*Gvec_G(ii, jj, kk,jdir) &
                        *FourPi_G2(ii, jj, kk)
                end do
             end do
          end do
          !$omp end parallel do
          stress_l(idir,jdir) = - ss 
       end do
    end do

    ss=M_ZERO
    !$omp parallel do private(ii, jj, kk) reduction(+:ss)
    do kk = 1, cube%rs_n_global(3)
       do jj = 1, cube%rs_n_global(2)
          do ii = 1, cube%rs_n_global(1)
             ss = ss + abs(rho_total_fs(ii, jj, kk))**2*FourPi_G2 (ii, jj, kk)
          end do
       end do
    end do
    !$omp end parallel do

    do idir = 1,3
       stress_l(idir,idir) = stress_l(idir,idir) + ss
    end do
    stress_l = CNST(0.5) * stress_l

    
    stress_Hartree =  stress_l
    stress = stress + stress_l

    call profiling_out(prof)
    POP_SUB(stress_from_Hartree)

  end subroutine stress_from_Hartree

  ! -------------------------------------------------------
  ! We assume hm%energy%echange, correlation, and intnvxc
  ! have already been calculated somewhere else.
  subroutine stress_from_xc(der, hm, stress, stress_xc)
    type(derivatives_t),  intent(in) :: der
    type(hamiltonian_elec_t),  intent(in)    :: hm
    FLOAT,                         intent(inout) :: stress(:, :)
    FLOAT,                         intent(out) :: stress_xc(3, 3) ! temporal
    FLOAT :: stress_l(3, 3)
    integer :: ii
    type(profile_t), save :: prof

    call profiling_in(prof, "STRESS_FROM_XC")    
    PUSH_SUB(stress_from_xc)

    stress_l = M_ZERO
    do ii = 1,3
       stress_l(ii,ii) = - hm%energy%exchange - hm%energy%correlation &
            + hm%energy%intnvxc
    end do
    stress_l = stress_l/der%mesh%sb%latt%rcell_volume
    
    stress_xc(:,:) = stress_l(:,:)
    stress(:,:) = stress(:,:) + stress_l(:,:)
    
    call profiling_out(prof)
    POP_SUB(stress_from_xc)
  end subroutine stress_from_xc

  ! -------------------------------------------------------
  subroutine stress_from_pseudo(gr, hm, st, ions, stress, stress_ps)
    type(grid_t),      target,        intent(in) :: gr !< grid
    type(hamiltonian_elec_t),  intent(inout)    :: hm
    type(states_elec_t),    intent(inout) :: st
    type(ions_t),              intent(in) :: ions !< geometry
    type(cube_t),    pointer             :: cube
    type(derivatives_t),  pointer :: der
    FLOAT,                         intent(inout) :: stress(:, :)
    FLOAT,                         intent(out) :: stress_ps(3, 3) ! temporal

    FLOAT :: stress_l(3, 3)
    FLOAT :: stress_t_SR(3, 3), stress_t_LR(3, 3), stress_t_NL(3, 3)
    CMPLX, allocatable :: gpsi(:, :, :), psi(:, :), rppsi(:, :, :)
    FLOAT :: energy_ps_SR
    FLOAT,  allocatable :: vloc(:),rvloc(:,:)
    FLOAT,  allocatable :: rho_t(:),grho_t(:,:)
    FLOAT :: sigma_erf, alpha, gx, g2
    CMPLX :: zphase, zfact
    integer :: ik, ispin, ist, idim, idir, jdir, iatom
    integer :: ii,jj,kk
    type(profile_t), save :: prof

    call profiling_in(prof, "STRESS_FROM_PSEUDO")    
    PUSH_SUB(stress_from_pseudo)

    stress_l = M_ZERO
    
    cube => hm%psolver_fine%cube
    der => gr%der


    SAFE_ALLOCATE(psi(1:der%mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(gpsi(1:der%mesh%np, 1:der%dim, 1:st%d%dim))
    SAFE_ALLOCATE(rppsi(1:der%mesh%np, 1:der%dim+1, 1:st%d%dim))
    SAFE_ALLOCATE(vloc(1:gr%mesh%np))
    SAFE_ALLOCATE(rvloc(1:gr%mesh%np, 1:der%dim))
    SAFE_ALLOCATE(rho_t(1:der%mesh%np_part))
    SAFE_ALLOCATE(grho_t(1:der%mesh%np,1:der%dim))


! calculate stress from non-local pseudopotentials
    stress_t_NL = M_ZERO
    do ik = st%d%kpt%start, st%d%kpt%end
       ispin = st%d%get_spin_index(ik)
       do ist = st%st_start, st%st_end
          
          call states_elec_get_state(st, der%mesh, ist, ik, psi)
          
          do idim = 1, st%d%dim
             call boundaries_set(der%boundaries, psi(:, idim))
          end do

          if (allocated(hm%hm_base%phase)) then 
            call states_elec_set_phase(st%d, psi, hm%hm_base%phase(1:der%mesh%np_part, ik), der%mesh%np_part, .false.)
          end if
          
          do idim = 1, st%d%dim
             call zderivatives_grad(der, psi(:, idim), gpsi(:, :, idim), set_bc = .false.)
          end do



          rppsi = M_ZERO
          do iatom = 1, ions%natoms
             if(species_is_ps(ions%atom(iatom)%species)) then
                call zr_project_psi(hm%ep%proj(iatom), der%mesh, st%d%dim, ik, psi, rppsi)
             end if
          end do

          do idir = 1,3
             do jdir = 1,3
                do idim = 1, st%d%dim

                   stress_t_NL(idir, jdir) = stress_t_NL(idir, jdir) &
                        +2d0*st%d%kweights(ik)*st%occ(ist, ik) &
                        *dmf_integrate(der%mesh, TOFLOAT(conjg(gpsi(1:der%mesh%np,idir,idim))*rppsi(1:der%mesh%np,jdir,idim)))

                   if(idir /= jdir)cycle

                   stress_t_NL(idir, jdir) = stress_t_NL(idir, jdir) &
                        +st%d%kweights(ik)*st%occ(ist, ik) &
                        *dmf_integrate(der%mesh, TOFLOAT(conjg(psi(1:der%mesh%np,idim))*rppsi(1:der%mesh%np,4,idim))) 

                end do
             end do
          end do
       end do
    end do

    
    if(st%parallel_in_states .or. st%d%kpt%parallel) then
      ! TODO: this could take dim = (/der%mesh%np, der%dim, st%d%nspin/)) to reduce the amount of data copied
       call comm_allreduce(st%st_kpt_mpi_grp, stress_t_NL)
    end if


    stress_t_NL = stress_t_NL/der%mesh%sb%latt%rcell_volume

! calculate stress from short-range local pseudopotentials
    stress_t_SR = M_ZERO

    rho_t(1:gr%mesh%np) = rho_total(1:gr%mesh%np)
    call boundaries_set(der%boundaries, rho_t(:))
    call dderivatives_grad(der, rho_t, grho_t, set_bc = .false.)

    vloc = M_ZERO
    rvloc = M_ZERO    
    do iatom = 1, ions%natoms
       call epot_local_pseudopotential_SR(gr%der, hm%ions, iatom, vloc, rvloc)
    end do

    energy_ps_SR = dmf_dotp(gr%mesh, vloc, rho_total(:))
    do idir = 1,3
       stress_t_SR(idir,idir) = energy_ps_SR
    end do
    
    do idir = 1,3
       do jdir =1,3
          
          stress_t_SR(idir, jdir) = stress_t_SR(idir, jdir) &
               +dmf_dotp(gr%mesh, rvloc(:, jdir), grho_t(:, idir))
       end do
    end do

    stress_t_SR = stress_t_SR/der%mesh%sb%latt%rcell_volume

! calculate stress from long-range local pseudopotentials
    stress_t_LR = M_ZERO


 ! We assume this value is applied for range-separation for all the species
    sigma_erf = CNST(0.625) 
    alpha = M_ONE/(sqrt(M_TWO)*sigma_erf)

    do kk = 1, cube%rs_n_global(3)
       do jj = 1, cube%rs_n_global(2)
          do ii = 1, cube%rs_n_global(1)
             
             zphase = M_Z0
             do iatom = 1, ions%natoms
                gx = sum(Gvec(ii, jj, kk, 1:3)*ions%atom(iatom)%x(1:3))
                zphase = zphase + species_zval(ions%atom(iatom)%species) &
                     *TOCMPLX(cos(gx), sin(gx)) 
             end do
             g2 = sum(Gvec(ii, jj, kk, 1:3)**2)
             zfact =  sigma_erf**2*FourPi_G2(ii,jj,kk) &
                  *exp(-M_HALF*sigma_erf**2*g2)/(M_HALF*sigma_erf**2) &
                  *(M_HALF*sigma_erf**2*g2 + M_ONE) &
                  *rho_total_fs(ii,jj,kk)*conjg(zphase)
             
             do idir =1,3
                do jdir =1,3
                   stress_t_LR(idir, jdir) = stress_t_LR(idir, jdir) &
                        + TOFLOAT(zfact)*Gvec_G(ii, jj, kk,idir)*Gvec_G(ii, jj, kk,jdir)
                end do
             end do
             
             do idir =1,3
                stress_t_LR(idir, idir) = stress_t_LR(idir, idir) &
                     - FourPi_G2(ii,jj,kk)*exp(-M_HALF*sigma_erf**2*g2) &
                     * TOFLOAT(rho_total_fs(ii,jj,kk)*conjg(zphase))
             end do

          end do
       end do
    end do

    stress_t_LR = stress_t_LR/der%mesh%sb%latt%rcell_volume

    
    stress_ps = stress_t_SR + stress_t_LR + stress_t_NL

!!! NOTE!! This part is moved to Ewald contoributoin
!! Contribition from G=0 component of the long-range part        
!    charge = M_ZERO
!    do iatom = 1, ions%natoms
!       zi = species_zval(ions%atom(iatom)%species)
!       charge = charge + zi
!    end do
!
!    do idir = 1,3
!       stress_ps(idir,idir) = stress_ps(idir,idir) &
!            + 2d0*M_PI*sigma_erf**2*charge**2 /der%mesh%sb%latt%rcell_volume**2
!    end do
    
    stress = stress + stress_ps

    
    call profiling_out(prof)
    POP_SUB(stress_from_pseudo)      
      
  end subroutine stress_from_pseudo
      
  ! -------------------------------------------------------
  subroutine epot_local_pseudopotential_SR(der, ions, iatom, vpsl, rvpsl)
    type(derivatives_t),      intent(in)    :: der
    type(ions_t),             intent(in)    :: ions
    integer,                  intent(in)    :: iatom
    FLOAT,                    intent(inout) :: vpsl(:)
    FLOAT,                    intent(inout) :: rvpsl(:,:)
      
    integer :: ip
    FLOAT :: radius, r
    FLOAT, allocatable :: vl(:)
    type(submesh_t)  :: sphere
    type(profile_t), save :: prof
    type(ps_t), pointer :: ps
    
    PUSH_SUB(epot_local_pseudopotential_sr)
    call profiling_in(prof, "EPOT_LOCAL_PSEUDOPOTENTIAL_SR")

    !the localized part
      
    if(species_is_ps(ions%atom(iatom)%species)) then
        
       ps => species_ps(ions%atom(iatom)%species)

       radius = spline_cutoff_radius(ps%vl, ps%projectors_sphere_threshold) + der%mesh%spacing(1)

       call submesh_init(sphere, ions%space, der%mesh%sb, der%mesh, ions%atom(iatom)%x, radius)
       SAFE_ALLOCATE(vl(1:sphere%np))

       do ip = 1, sphere%np
         r = sphere%x(ip, 0)
         vl(ip) = spline_eval(ps%vl, r)
       end do

       nullify(ps) 
         
       ! Cannot be written (correctly) as a vector expression since for periodic systems,
       ! there can be values ip, jp such that sphere%map(ip) == sphere%map(jp).
       do ip = 1, sphere%np
          vpsl(sphere%map(ip)) = vpsl(sphere%map(ip)) + vl(ip)
          rvpsl(sphere%map(ip) ,1) = rvpsl(sphere%map(ip), 1) + sphere%x(ip, 1) * vl(ip)
          rvpsl(sphere%map(ip) ,2) = rvpsl(sphere%map(ip), 2) + sphere%x(ip, 2) * vl(ip)
          rvpsl(sphere%map(ip) ,3) = rvpsl(sphere%map(ip), 3) + sphere%x(ip, 3) * vl(ip)
       end do
       
       SAFE_DEALLOCATE_A(vl)
       call submesh_end(sphere)

    end if
    
    call profiling_out(prof)
    POP_SUB(epot_local_pseudopotential_sr)
  end subroutine epot_local_pseudopotential_SR

  ! -------------------------------------------------------
  subroutine poisson_fft_gg_transform_l(gg_in, temp, sb, gg, modg2)
    integer,           intent(in)    :: gg_in(:)
    FLOAT,             intent(in)    :: temp(:)
    type(simul_box_t), intent(in)    :: sb
    FLOAT,             intent(inout) :: gg(:)
    FLOAT,             intent(out)   :: modg2

    ! no PUSH_SUB, called too frequently

    gg(1:3) = gg_in(1:3) * temp(1:3)
    gg(1:3) = matmul(sb%latt%klattice_primitive(1:3,1:3),gg(1:3))

    modg2 = sum(gg(1:3)**2)

  end subroutine poisson_fft_gg_transform_l

  ! ---------------------------------------------------------
  subroutine stress_from_Ewald_sum(ions, stress, stress_Ewald)
    type(ions_t),             intent(in)    :: ions
    FLOAT,                    intent(inout) :: stress(:, :)
    FLOAT,                    intent(out)   :: stress_Ewald(3, 3) ! temporal

    FLOAT :: stress_l(3, 3)

    FLOAT :: rr, xi(ions%space%dim), zi, zj, erfc, rcut
    integer :: iatom, jatom, icopy
    type(lattice_iterator_t) :: latt_iter
    integer :: ix, iy, iz, isph, ss, idim, idir, jdir
    FLOAT   :: gg(ions%space%dim), gg2, gx
    FLOAT   :: factor, charge, Hp, charge_sq
    FLOAT   :: alpha
    CMPLX   :: sumatoms, aa
    FLOAT :: sigma_erf
    type(profile_t), save :: prof

    call profiling_in(prof, "STRESS_FROM_EWALD")    
    PUSH_SUB(stress_from_Ewald_sum)

    ! Currently this is only implemented for 3D
    ASSERT(ions%space%dim == 3)

    alpha = ions%ion_interaction%alpha

    rcut = CNST(6.0)/alpha
    stress_l = M_ZERO
    latt_iter = lattice_iterator_t(ions%latt, rcut)
    ! the short-range part is calculated directly
    do iatom = ions%atoms_dist%start, ions%atoms_dist%end
      if (.not. species_represents_real_atom(ions%atom(iatom)%species)) cycle
      zi = species_zval(ions%atom(iatom)%species)

      do icopy = 1, latt_iter%n_cells
        xi = ions%atom(iatom)%x(1:ions%space%dim) + latt_iter%get(icopy)
        
        do jatom = 1,  ions%natoms
          zj = species_zval(ions%atom(jatom)%species)
          rr = norm2(xi - ions%atom(jatom)%x(1:ions%space%dim))
          
          if(rr < CNST(1e-5)) cycle
          
          erfc = M_ONE - loct_erf(alpha*rr)
          Hp = -M_TWO/sqrt(M_PI)*exp(-(alpha*rr)**2) - erfc/(alpha*rr)
          factor = M_HALF*zj*zi*alpha*Hp
          do idir = 1,3
            do jdir =1,3
              stress_l(idir, jdir) = stress_l(idir, jdir) &
                - factor*(xi(idir) - ions%atom(jatom)%x(idir))*(xi(jdir) - ions%atom(jatom)%x(jdir))/(rr**2)
            end do
          end do

        end do

      end do
    end do

    if(ions%atoms_dist%parallel) then
       call comm_allreduce(ions%atoms_dist%mpi_grp, stress_l)
    end if

    ! And the long-range part, using an Ewald sum
    charge = M_ZERO
    charge_sq = M_ZERO
    do iatom = 1, ions%natoms
       zi = species_zval(ions%atom(iatom)%species)
       charge = charge + zi
       charge_sq = charge_sq + zi**2
    end do
    
    ! get a converged value for the cutoff in g
    rcut = huge(rcut)
    do idim = 1, ions%space%dim
      rcut = min(rcut, sum(ions%latt%klattice(1:ions%space%dim, idim)**2))
    end do

    rcut = sqrt(rcut)
    
    isph = ceiling(CNST(9.5)*alpha/rcut)
    
    do ix = -isph, isph
       do iy = -isph, isph
          do iz = -isph, isph
          
             ss = ix**2 + iy**2 + iz**2
          
             if(ss == 0 .or. ss > isph**2) cycle

             gg = ix*ions%latt%klattice(:, 1) + iy*ions%latt%klattice(:, 2) + iz*ions%latt%klattice(:, 3)
             gg2 = sum(gg**2)

             ! g=0 must be removed from the sum
             if(gg2 < M_EPSILON) cycle
          
             gx = -CNST(0.25)*gg2/alpha**2
             
             if(gx < CNST(-36.0)) cycle

             factor = M_TWO*M_PI*exp(gx)/(ions%latt%rcell_volume*gg2)

             if(factor < epsilon(factor)) cycle

             sumatoms = M_Z0

             do iatom = 1, ions%natoms
                gx = sum(gg*ions%atom(iatom)%x(1:ions%space%dim))
                aa = species_zval(ions%atom(iatom)%species)*TOCMPLX(cos(gx), sin(gx))
                sumatoms = sumatoms + aa
             end do

             factor = factor*abs(sumatoms)**2
             
             do idir = 1, 3
                do jdir = 1, 3
                  stress_l(idir, jdir) = stress_l(idir, jdir) &
                        - M_TWO*factor*gg(idir)*gg(jdir)/gg2*(CNST(0.25)*gg2/alpha**2+M_ONE)
                   
                end do
                stress_l(idir, idir) = stress_l(idir, idir) + factor
             end do

          end do
       end do
    end do


    factor = M_HALF*M_PI*charge**2/(ions%latt%rcell_volume*alpha**2)
    stress_l(1, 1) = stress_l(1, 1) - factor
    stress_l(2, 2) = stress_l(2, 2) - factor
    stress_l(3, 3) = stress_l(3, 3) - factor


    ! Contribition from G=0 component of the long-range part    
    sigma_erf = CNST(0.625)
    do idir = 1,3
      stress_l(idir,idir) = stress_l(idir,idir) &
            + M_TWO*M_PI*sigma_erf**2*charge**2 /ions%latt%rcell_volume
    end do

    stress_l = stress_l/ions%latt%rcell_volume

    stress_Ewald = stress_l
    stress = stress + stress_l

    
    call profiling_out(prof)
    POP_SUB(stress_from_Ewald_sum)
    
  end subroutine stress_from_Ewald_sum
  ! -------------------------------------------------------
  ! -------------------------------------------------------
  ! -------------------------------------------------------
  

  
end module stress_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
