!! Copyright (C) 2016 H. Huebener
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

program oct_floquet
  use blas_oct_m
  use calc_mode_par_oct_m
  use comm_oct_m
  use density_oct_m
  use fft_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use ions_oct_m
  use lalg_adv_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use simul_box_oct_m
  use states_elec_oct_m
  use states_elec_restart_oct_m
  use electrons_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use v_ks_oct_m
  use xc_oct_m

  implicit none

  integer :: ierr

  type(electrons_t), pointer :: sys
  type(states_elec_t) :: st
  type(grid_t)   :: gr
  CMPLX, allocatable :: hmss(:,:), psi(:,:,:), hpsi(:,:,:), temp_state1(:,:)
  CMPLX, allocatable :: HFloquet(:,:,:), HFloq_eff(:,:), temp(:,:)
  FLOAT, allocatable :: eigenval(:), bands(:,:)
  character(len=80) :: filename
  integer :: it, nT, ik, ist, in, im, file, idim, nik, ik_count
  integer :: Forder, Fdim, m0, n0, n1, nst, ii, jj, lim_nst
  FLOAT :: dt, Tcycle,omega
  logical :: downfolding = .false.
  type(mesh_t) :: mesh
  type(restart_t) :: restart
  
  ! the usual initializations
  call global_init(is_serial = .false.)
  call calc_mode_par_init()

  call parser_init()
  
  call messages_init()

  call io_init()
  call profiling_init(global_namespace)

  call print_header()
  call messages_print_stress(stdout, "Non-interacting Floquet")
  call messages_print_stress(stdout)

  call messages_experimental("oct-floquet utility")
  call fft_all_init(global_namespace)
  call unit_system_init(global_namespace)
  call restart_module_init(global_namespace)

  call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
  sys => electrons_t(global_namespace)
  call sys%init_parallelization(mpi_world)
  ! make shortcut copies
  st = sys%st
  gr = sys%gr

  ! generate the full hamiltonian following the sequence in td_init
  call hamiltonian_elec_epot_generate(sys%hm, global_namespace, sys%space, gr, sys%ions, st, time=M_ZERO)
  call hamiltonian_elec_update(sys%hm, gr%mesh, global_namespace, sys%space, time = M_ZERO)

  call states_elec_allocate_wfns(st, gr%mesh)
  ! not sure this is needed ...
  if (gauge_field_is_applied(sys%hm%ep%gfield)) then
     !if the gauge field is applied, we need to tell v_ks to calculate the current
     call v_ks_calculate_current(sys%ks, .true.)

     ! initialize the vector field and update the hamiltonian     
     call gauge_field_init_vec_pot(sys%hm%ep%gfield, sys%ions%latt%rcell_volume, st%qtot)
     call hamiltonian_elec_update(sys%hm, gr%mesh, global_namespace, sys%space, time = M_ZERO)
  end if

  call restart_init(restart, global_namespace, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh, exact=.true.)
  if(ierr == 0) then
    call states_elec_load(restart, global_namespace, sys%space, st, gr%mesh, sys%kpoints, ierr, label = ": gs")
  end if
  if (ierr /= 0) then
     message(1) = 'Unable to read ground-state wavefunctions.'
     call messages_fatal(1)
  end if

  call density_calc(st, gr, st%rho)
  call v_ks_calc(sys%ks, global_namespace, sys%space, sys%hm, st, sys%ions, calc_eigenval=.true., time = M_ZERO)
  call hamiltonian_elec_update(sys%hm, gr%mesh, global_namespace, sys%space, time = M_ZERO)

  call floquet_init()

  call floquet_solve_non_interacting()


#if defined(HAVE_MPI)
  ! wait for all processors to finish
  call MPI_Barrier(mpi_world%comm, mpi_err)
#endif

  call fft_all_end()
  SAFE_DEALLOCATE_P(sys)
  call profiling_end(global_namespace)
  call io_end()
  call print_date("Calculation ended on ")
  call messages_end()

  call parser_end()
  call global_end()

contains

    !------------------------------------------------
    subroutine floquet_init()

      PUSH_SUB(floquet_init)

      !for now no domain distribution allowed
      ASSERT(gr%mesh%np == gr%mesh%np_global)

      ! variables documented in td/td_write.F90
      call parse_variable(global_namespace, 'TDFloquetFrequency', M_ZERO, omega, units_inp%energy)
      call messages_print_var_value(stdout,'Frequency used for Floquet analysis', omega)
      if(abs(omega)<=M_EPSILON) then
         message(1) = "Please give a non-zero value for TDFloquetFrequency"
         call messages_fatal(1)
      endif

      ! get time of one cycle
      Tcycle=M_TWO*M_PI/omega

      call parse_variable(global_namespace, 'TDFloquetSample',20 ,nt)
      call messages_print_var_value(stdout,'Number of Floquet time-sampling points', nT)
      dt = Tcycle/TOFLOAT(nT)

      call parse_variable(global_namespace, 'TDFloquetDimension',-1,Forder)
      if(Forder.ge.0) then
        call messages_print_var_value(stdout,'Order of multiphoton Floquet-Hamiltonian', Forder)
        !Dimension of multiphoton Floquet-Hamiltonian
        Fdim = 2*Forder+1
      else
        message(1) = 'Floquet-Hamiltonian is downfolded'
        call messages_info(1)
        downfolding = .true.
        Forder = 1
        Fdim = 3
      endif

      dt = Tcycle/TOFLOAT(nT)

      POP_SUB(floquet_init)

  end subroutine floquet_init

  !---------------------------------------------------
  subroutine floquet_solve_non_interacting()	
    type(states_elec_t) :: hm_st

    PUSH_SUB(floquet_solve_non_interacting)

    mesh = gr%mesh
    nst = st%nst
    
    SAFE_ALLOCATE(hmss(1:nst,1:nst))
    SAFE_ALLOCATE( psi(1:nst,1:st%d%dim,1:mesh%np))
    SAFE_ALLOCATE(hpsi(1:nst,1:st%d%dim,1:mesh%np))
    SAFE_ALLOCATE(temp_state1(1:mesh%np,1:st%d%dim))

    ! this is used to initialize the local state object
    call states_elec_copy(hm_st, st)

    ! we are only interested for k-point with zero weight
    nik = sys%kpoints%nik_skip

    ! multiphoton Floquet Hamiltonian, layout:
    !     (H_{-n,-m} ...  H_{-n,0} ...  H_{-n,m}) 
    !     (    .      .      .      .      .    )
    ! H = (H_{0,-m}  ...  H_{0,0}  ...  H_{0,m} )
    !     (    .      .      .      .      .    )
    !     (H_{n,-m}  ...  H_{n,0}  ...  H_{n,m} )    
    SAFE_ALLOCATE(HFloquet(1:nik,1:nst*Fdim, 1:nst*Fdim))
    HFloquet(1:nik,1:nst*Fdim, 1:nst*Fdim) = M_ZERO

    ! perform time-integral over one cycle
    do it=1,nT
      ! get non-interacting Hamiltonian at time (offset by one cycle to allow for ramp)
      call hamiltonian_elec_update(sys%hm, gr%mesh, global_namespace, sys%space, time=Tcycle+it*dt)
      ! get hpsi
      call zhamiltonian_elec_apply_all(sys%hm, global_namespace, gr%mesh, st, hm_st)

      ! project Hamiltonian into grounstates for zero weight k-points
      ik_count = 0

      do ik = sys%kpoints%reduced%npoints-nik+1, sys%kpoints%reduced%npoints
        ik_count = ik_count + 1

        psi(1:nst, 1:st%d%dim, 1:mesh%np)= M_ZERO
        hpsi(1:nst, 1:st%d%dim, 1:mesh%np)= M_ZERO

        do ist = st%st_start, st%st_end
          if(state_kpt_is_local(st, ist, ik)) then
            call states_elec_get_state(st, mesh, ist, ik,temp_state1 )
            do idim = 1, st%d%dim
              psi(ist, idim, 1:mesh%np) = temp_state1(1:mesh%np, idim)
            end do
            call states_elec_get_state(hm_st, mesh, ist, ik,temp_state1 )
            do idim = 1, st%d%dim
              hpsi(ist, idim, 1:mesh%np) = temp_state1(1:mesh%np, idim)
            end do
          end if
        end do
        call comm_allreduce(mpi_world, psi)
        call comm_allreduce(mpi_world, hpsi)
        hmss(1:nst, 1:nst) = M_ZERO
        call zgemm( 'n',                                  &
                    'c',                                  &
                    nst,                                  &
                    nst,                                  &
                    mesh%np_global*st%d%dim,              &
                    TOCMPLX(mesh%volume_element, M_ZERO), &
                    hpsi(1, 1, 1),                        &
                    ubound(hpsi, dim = 1),                &
                    psi(1, 1, 1),                         &
                    ubound(psi, dim = 1),                 &
                    M_z0,                                 &
                    hmss(1, 1),                           &
                    ubound(hmss, dim = 1))

        hmss(1:nst,1:nst) = CONJG(hmss(1:nst,1:nst))

        ! accumulate the Floquet integrals
        do in=-Forder,Forder
           do im=-Forder,Forder
              ii=(in+Forder)*nst
              jj=(im+Forder)*nst
              HFloquet(ik_count,ii+1:ii+nst,jj+1:jj+nst) =  &
                HFloquet(ik_count,ii+1:ii+nst,jj+1:jj+nst) + hmss(1:nst,1:nst)*exp(-(in-im)*M_zI*omega*it*dt)
              ! diagonal term
              if(in==im) then
                 do ist = 1,nst
                    HFloquet(ik_count,ii+ist,ii+ist) = HFloquet(ik_count,ii+ist,ii+ist) + in*omega
                 end do
              end if
           end do
        end do
      end do !ik

    end do ! it

    HFloquet(:,:,:) = M_ONE/nT*HFloquet(:,:,:)

    ! diagonalize Floquet Hamiltonian
    if(downfolding) then
       ! here perform downfolding
       SAFE_ALLOCATE(HFloq_eff(1:nst,1:nst))
       SAFE_ALLOCATE(eigenval(1:nst))
       SAFE_ALLOCATE(bands(1:nik,1:nst))

       HFloq_eff(1:nst,1:nst) = M_ZERO
       do ik=1,nik
          ! the HFloquet blocks are copied directly out of the super matrix
          m0 = nst ! the m=0 start position
          n0 = nst ! the n=0 start postion
          n1 = 2*nst ! the n=+1 start postion
          HFloq_eff(1:nst,1:nst) = HFloquet(ik,n0+1:n0+nst,m0+1:m0+nst) + &
               M_ONE/omega*(matmul(HFloquet(ik,1:nst,m0+1:m0+nst), HFloquet(ik,n1+1:n1+nst,m0+1:m0+nst))- &
                            matmul(HFloquet(ik,n1+1:n1+nst,m0+1:m0+nst), HFloquet(ik,1:nst,m0+1:m0+nst)))

           eigenval(1:nst) = M_ZERO
          call lalg_eigensolve(nst, HFloq_eff, eigenval)
          bands(ik,1:nst) = eigenval(1:nst)
       end do
       SAFE_DEALLOCATE_A(HFloq_eff)
    else
      ! the full Floquet 
      SAFE_ALLOCATE(eigenval(1:nst*Fdim))
      SAFE_ALLOCATE(bands(1:nik,1:nst*Fdim))
      SAFE_ALLOCATE(temp(1:nst*Fdim, 1:nst*Fdim))

      do ik=1,nik
         temp(1:nst*Fdim,1:nst*Fdim) = HFloquet(ik,1:nst*Fdim,1:nst*Fdim)
         call lalg_eigensolve(nst*Fdim, temp, eigenval)
         bands(ik,1:nst*Fdim) = eigenval(1:nst*Fdim)
      end do
    end if

    !write bandstructure to file
    if(downfolding) then
      lim_nst = nst
      filename="downfolded_floquet_bands"
    else
       lim_nst = nst*Fdim
       filename="floquet_bands"
    end if
    ! write bands (full or downfolded)
    if(mpi_world%rank==0) then
      file=987254
      file = io_open(filename, sys%namespace, action = 'write')
      do ik=1,nik
        do ist = 1,lim_nst
          write(file,'(e12.6, 1x)',advance='no') bands(ik,ist)
        end do
        write(file,'(1x)')
      end do
      call io_close(file)
    endif
    
    if(.not.downfolding) then
      ! for the full Floquet case compute also the trivially shifted
      ! Floquet bands for reference (i.e. setting H_{nm}=0 for n!=m)
      bands(1:nik,1:nst*Fdim) = M_ZERO
      do ik=1,nik
        temp(1:nst*Fdim,1:nst*Fdim) = M_ZERO
        do jj=0,Fdim-1
          ii=jj*nst
          temp(ii+1:ii+nst,ii+1:ii+nst) = HFloquet(ik,ii+1:ii+nst,ii+1:ii+nst)
        end do
        call lalg_eigensolve(nst*Fdim, temp, eigenval)
        bands(ik,1:nst*Fdim) = eigenval(1:nst*Fdim)
      end do
    
      if(mpi_world%rank==0) then
        filename='trivial_floquet_bands'
        file = io_open(filename, sys%namespace, action = 'write')
        do ik=1,nik
          do ist = 1,lim_nst
            write(file,'(e12.6, 1x)',advance='no') bands(ik,ist)
          end do
          write(file,'(1x)')
        end do
        call io_close(file)
      endif
     end if
  
    ! reset time in Hamiltonian
    call hamiltonian_elec_update(sys%hm, gr%mesh, global_namespace, sys%space, time=M_ZERO)

    SAFE_DEALLOCATE_A(hmss)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(temp_state1)
    SAFE_DEALLOCATE_A(HFloquet)
    SAFE_DEALLOCATE_A(eigenval)
    SAFE_DEALLOCATE_A(bands)
    SAFE_DEALLOCATE_A(temp)

    POP_SUB(solve_non_interacting)

  end subroutine floquet_solve_non_interacting
  
end program oct_floquet
  
!! Local Variables:
!! mode: f90				  
!! coding: utf-8 
!! End:
