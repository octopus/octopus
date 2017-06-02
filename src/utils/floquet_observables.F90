!! Copyright (C) 2017 U. De Giovannini
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

program floquet_observables
  use calc_mode_par_oct_m
  use comm_oct_m
  use command_line_oct_m
  use current_oct_m
  use geometry_oct_m
  use fft_oct_m
  use floquet_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use kpoints_oct_m
  use io_binary_oct_m
  use io_function_oct_m
  use io_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use system_oct_m
  use sort_oct_m
  use space_oct_m
  use string_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_io_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use xc_oct_m
  
  implicit none

  type fobs_t
    FLOAT   :: emax  ! max energy         
    FLOAT   :: de    ! energy step      
    integer :: ne    ! number of energy steps
    FLOAT   :: gamma ! lifetime
    integer :: nst (1:2) ! states output limits
    integer :: nkpt(1:2) ! kpoints output limits
    integer :: gauge  ! the gauge used to calculate observables 
  end type fobs_t


  integer              :: out_what  

  integer              :: ierr
  integer              :: dim, dir, how, idim, pdim
  integer              :: ii, i1,i2,i3
  type(block_t)        :: blk  
  
  type(restart_t)      :: restart
  type(system_t)      :: sys
  type(hamiltonian_t) :: hm
  character(len=512)   :: filename, str, str2
  integer              :: ist, ispin  
  type(states_t)          :: dressed_st 
  type(states_t), pointer ::    gs_st

  type(fobs_t)         :: obs


  call getopt_init(ierr)
  if(ierr /= 0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-photoelectron-spectrum command is not available."
    call messages_fatal(2)
  end if


  call global_init()
  
  call messages_init()  
  call io_init()
  
  call calc_mode_par_init()

  call fft_all_init()
  call unit_system_init()
  
  call calc_mode_par_set_parallelization(P_STRATEGY_OTHER,   default = .false.)
!   call calc_mode_par_set_parallelization(P_STRATEGY_KPOINTS, default = .true. )
  call calc_mode_par_set_parallelization(P_STRATEGY_STATES,  default = .false.)
!   call calc_mode_par_set_parallelization(P_STRATEGY_DOMAINS, default = .true. )
  
  
  
  call system_init(sys)
  
  call hamiltonian_init(hm, sys%gr, sys%geo, sys%st, sys%ks%theory_level, sys%ks%xc_family, &
                        sys%ks%xc_flags, family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))
  
  
  
  

  call floquet_init(sys,hm%F,sys%st%d%dim)
  gs_st => sys%st

  call states_init(dressed_st, sys%gr, sys%geo,floquet_dim=hm%F%floquet_dim)
  call kpoints_distribute(dressed_st%d,sys%mc)
  call states_distribute_nodes(dressed_st,sys%mc)
  call states_exec_init(dressed_st, sys%mc)
  call restart_module_init()

  !%Variable FloquetObservableCalc
  !%Type flag
  !%Default none
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Specifies what observables to calculate form the Floquet solution.
  !% Example: <tt>f_norms + f_spin</tt>
  !%Option f_norms bit(1)
  !% The norms of the Floquet states in each tensorial subspace.
  !%Option f_arpes bit(2)
  !% Calculate ARPES matrix elements for Floquet states.
  !%Option f_spin bit(3)
  !% Calculate the spin polarization of each state. (Not implemented)
  !%Option f_td_spin bit(4)
  !% Calculate the time-dependent spin projections of the Floquet eigenstates
  !% given by FloquetObservableTDspinKpoints and FloquetObservableTDspinState
  !%Option f_hhg_w bit(5)
  !% Calculate the HHG spectral weights.
  !%Option f_hhg bit(6)
  !% Calculate the HHG spectrum.
  !%Option f_wfs bit(7)
  !% Output Floquet states in the format defined by OutputFormat.
  !%Option f_conductivity bit(8)
  !% Calculate Floquet optical conductivity.
  !%End
  call parse_variable('FloquetObservableCalc', out_what, out_what)
  

  ! load floquet states only for certain tasks
  ! Note: this way it makes impossible to combine options that needs all the Floquet 
  ! states with the ones that don't UDG
  if(.not. (iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_TD_SPIN) /= 0) .and. &
     .not. (iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_WFS) /= 0) ) then
       call states_allocate_wfns(dressed_st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)
       call floquet_restart_dressed_st(hm, sys, dressed_st, ierr)
       call messages_write('Read Floquet restart files.')
       call messages_info()
  end if
  
  !%Variable FloquetObservableGauge
  !%Type flag
  !%Default guess
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Selects the gauge to be used for the calculation of the dipole matrix elements.
  !% By default is chosen according to the laser field used.
  !%Option f_velocity 1
  !% The velocity gauge (use the current operator): <i|J|j> ~ <i|p.A|j> .
  !%Option f_lenght 2
  !% The length gauge: <i|r.E|j>.
  !%End
  obs%gauge = OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGHT
  if (gauge_field_is_applied(hm%ep%gfield) .or. simul_box_is_periodic(sys%gr%sb)) & 
      obs%gauge = OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY

  call parse_variable('FloquetObservableGauge', obs%gauge, obs%gauge)
  call messages_print_var_option(stdout,'FloquetObservableGauge', obs%gauge)
  
  
  !%Variable FloquetObservableEnergyMax
  !%Type float
  !%Default 20 eV
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Observable maximum energy.
  !%End
  call parse_variable('FloquetObservableEnergyMax', CNST(0.7349968762), obs%emax, units_inp%energy)
  call messages_print_var_value(stdout,'FloquetObservableEnergyMax', obs%emax)


  !%Variable FloquetObservableEnergyStep
  !%Type float
  !%Default 0.1 eV
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Observable energy step.
  !%End
  call parse_variable('FloquetObservableEnergyStep', CNST(0.0036749843813163), obs%de, units_inp%energy)
  call messages_print_var_value(stdout,'FloquetObservableEnergyStep', obs%de)
  
  obs%ne = obs%emax/obs%de

  !%Variable FloquetObservableLifetimeBroadening
  !%Type float
  !%Default 0.001 eV
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Observable inverse lifetime used to artificially broaden spectral features.
  !%End
  call parse_variable('FloquetObservableLifetimeBroadening', CNST(0.001), obs%gamma, units_inp%energy)
  call messages_print_var_value(stdout,'FloquetObservableLifetimeBroadening', obs%gamma)
  
  

  !%Variable FloquetObservableSelectStates
  !%Type integer
  !%Default all
  !%Section Floquet
  !%Description
  !% Select a range of states and kpoints to be used for the calculation of the 
  !% wanted observable
  !%
  !% <tt>%FloquetObservableSelectStates
  !% <br>&nbsp;&nbsp; st_start| st_end 
  !% <br>&nbsp;&nbsp; kpt_start| kpt_end 
  !% <br>%</tt>
  !%
  !%End
  obs%nst(1) = 1
  obs%nst(2) = dressed_st%nst
  obs%nkpt(1) = 1 
  obs%nkpt(2) = dressed_st%d%kpt%nglobal
  
  if(parse_block('FloquetObservableSelectStates', blk) == 0) then
    if(parse_block_cols(blk,0) /= 2) call messages_input_error('FloquetObservableSelectStates')
    do idim = 1, 2
      call parse_block_integer(blk, 0, idim - 1, obs%nst(idim))
    end do
    if (parse_block_n(blk) > 1) then
      do idim = 1, 2
        call parse_block_integer(blk, 1, idim - 1, obs%nkpt(idim))
      end do
    end if
  end if

  write(message(1),'(a,i5,a,i5,a)') 'Info: FloquetObservableSelectStates [st_min, st_max]: [',&
                                     obs%nst(1),',', obs%nst(2),']'
  write(message(2),'(a,i5,a,i5,a)') '                                    [kp_min, kp_max]: [',&
                                      obs%nkpt(1),',', obs%nkpt(2),']'
  call messages_info(2)



  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_NORMS) /= 0) then
    call messages_write('Calculate norms of Floquet subspaces.')
    call messages_info()

    call calc_floquet_norms(sys%gr%der%mesh,sys%gr%sb%kpoints,gs_st,dressed_st, hm%F%iter,hm%F%floquet_dim)
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_ARPES) /= 0) then
    call messages_write('Calculate Floquet ARPES.')
    call messages_info()

    call calc_floquet_arpes()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_TD_SPIN) /= 0) then
    call messages_write('Calculate td-spin of Floquet eigenstates.')
    call messages_info()

    if(gs_st%d%dim /=2 ) then
      call messages_write('Need spin resolved calculation')
      call messages_fatal()
    end if

    call calc_floquet_td_spin()
  end if
  
  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_HHG_W) /= 0) then
    call messages_write('Calculate Floquet HH spectral weights.')
    call messages_info()

    call calc_floquet_hhg_weights()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_HHG) /= 0) then
    call messages_write('Calculate Floquet HH spectrum.')
    call messages_info()

    call calc_floquet_hhg()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_WFS) /= 0) then
    call messages_write('Output Floquet wavefunctions.')
    call messages_info()

    call out_floquet_wfs()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_CONDUCTIVITY) /= 0) then
    call messages_write('Calculate Floquet optical conductivity.')
    call messages_info()

    call calc_floquet_conductivity()
  end if


  
  call hamiltonian_end(hm)
  call system_end(sys)
  call fft_all_end()
  call io_end()
  call messages_end()
  call global_end()

contains 
    
  subroutine calc_floquet_arpes()
  
  FLOAT, allocatable :: spect(:,:), me(:,:) 
  FLOAT :: pomega, pol(1:3)
  type(block_t)        :: blk  
  
  
  !%Variable FloquetObservablePesOmega
  !%Type float
  !%Default 50 eV
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% The probe energy needed to calculate  photoemission matrix elements.
  !%End
  call parse_variable('FloquetObservablePesOmega', CNST(1.83749219065), pomega, units_inp%energy)
  call messages_print_var_value(stdout,'Frequency of PES probe field', pomega)
  
  
  !%Variable FloquetObservablePesPol
  !%Type block
  !%Default 
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Probe field polarization direction.
  !%
  !% <tt>%FloquetObservablePesPol
  !% <br> px | py | pz
  !% <br>%</tt>
  !%End
  pol(:)=M_ZERO
  pol(1)=M_ONE

  if(parse_block('FloquetObservablePesPol', blk) == 0) then
    if(parse_block_cols(blk,0) /= sys%gr%sb%dim) call messages_input_error('FloquetObservablePesPol')
    do idim = 1, sys%gr%sb%dim
      call parse_block_float(blk, 0, idim - 1, pol(idim))
    end do
    call parse_block_end(blk)
  end if
  
  write(message(1),'(a,f4.2,a,f4.2,a,f4.2,a)') 'Info: ARPES probe polarization: (', pol(1),',', pol(2),',', pol(3),')'
  call messages_info(1)
    
    
  SAFE_ALLOCATE(spect(dressed_st%nst,dressed_st%d%nik))
  SAFE_ALLOCATE(me(dressed_st%nst,dressed_st%d%nik))
  call floquet_photoelectron_spectrum(hm, sys, dressed_st, pomega, pol, spect, me)

  write(str,'(I5)') hm%F%iter
  if(simul_box_is_periodic(sys%gr%sb)) then
    filename = 'floquet_arpes_me_'//trim(adjustl(str))
    call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, sys%gr%sb, filename, vec = me)

    filename = 'floquet_arpes_'//trim(adjustl(str))
    call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, sys%gr%sb, filename, vec = spect)
  end if

  SAFE_DEALLOCATE_A(spect)
  SAFE_DEALLOCATE_A(me)
    
  end subroutine calc_floquet_arpes  

  !---------------------------------------
  subroutine floquet_photoelectron_spectrum(hm, sys, st, pomega, pol, spect, me)
    type(hamiltonian_t), intent(in) :: hm
    type(system_t), intent(in)      :: sys
    type(states_t), intent(in)      :: st
    FLOAT,          intent(in)      :: pomega     ! Probe field energy 
    FLOAT,          intent(in)      :: pol(:)     ! Probe field polarization vector
    FLOAT,          intent(out)     :: spect(:,:) ! the photoelectron spectrum
    FLOAT,          intent(out)     :: me(:,:)    ! the photoeletron matrix elements


    CMPLX, allocatable :: u_ma(:,:),  phase(:), tmp(:)
    FLOAT :: omega, dt , qq(1:3), kpt(1:3), xx(1:MAX_DIM)
    integer :: idx, im, it, ist, idim, nT, ik, ia, Fdim(2), imm, dim, pdim, ip, spindim

    type(mesh_t),   pointer :: mesh

    PUSH_SUB(floquet_photoelectron_spectrum)

    mesh  => sys%gr%der%mesh
    
    dt=hm%F%dt
    nT=hm%F%nT
    omega=hm%F%omega
    Fdim(:)=hm%F%order(:)
    spindim = hm%F%spindim 
    
    dim = mesh%sb%dim
    pdim = mesh%sb%periodic_dim
    
    SAFE_ALLOCATE(phase(1:mesh%np))
    SAFE_ALLOCATE(u_ma(1:mesh%np,hm%F%floquet_dim))
    SAFE_ALLOCATE(tmp(spindim))

    kpt(:) = M_ZERO
    qq(:)  = M_ZERO
    me(:,:) = M_ZERO
    spect(:,:) = M_ZERO
    do ik=st%d%kpt%start, st%d%kpt%end
      kpt(1:dim) = kpoints_get_point(mesh%sb%kpoints, ik) 

      do ia=st%st_start, st%st_end
        qq(dim)    = pomega - st%eigenval(ia,ik) - sum(kpt(1:pdim)**2)*M_HALF
        qq(1:pdim) = kpt(1:pdim)  
        
        do ip=1, mesh%np
          xx=mesh_x_global(mesh, ip) 
          phase(ip) = exp(-M_ZI*qq(dim)*xx(dim))
        end do
        
        call states_get_state(st, mesh, ia, ik, u_ma)
        
        tmp(:) = M_ZERO
        do idx=hm%F%flat_idx%start, hm%F%flat_idx%end
          it = hm%F%idx_map(idx,1)
          im = hm%F%idx_map(idx,2)
          imm = im - Fdim(1) + 1
          do idim=1,spindim
            tmp(idim)  =  tmp(idim) + exp(-M_zI*im*omega*it*dt)/nT * &
                          zmf_integrate(mesh, phase(1:mesh%np)*u_ma(1:mesh%np,(imm-1)*spindim+idim))
          end do
           
        end do 
        if(hm%F%is_parallel) call comm_allreduce(hm%F%mpi_grp%comm, tmp(:))   
        
        me(ia,ik) = sum(abs(tmp(:))**2) * sum((pol(1:dim)*qq(dim)))**2
        
        spect(ia,ik) =  me(ia,ik) * st%occ(ia,ik)
          
      end do
    end do
    
    call comm_allreduce(st%st_kpt_mpi_grp%comm, me)
    call comm_allreduce(st%st_kpt_mpi_grp%comm, spect)
    
    
    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(phase)
    SAFE_DEALLOCATE_A(u_ma)



    POP_SUB(floquet_photoelectron_spectrum)
    
  end subroutine floquet_photoelectron_spectrum    


  !-------------------------------------------------
  subroutine calc_floquet_td_spin()

  FLOAT :: spin(1:3)
  type(block_t)        :: blk
  integer :: nk_input, ik,iik, ist,iist, nst_input, it, iunit
  integer :: nik, dim, nst, itot
  integer, allocatable :: kpoints_input(:), states_input(:)
  CMPLX, allocatable   :: psi_t(:,:), zpsi(:), F_psi(:,:)
  FLOAT :: time

  !%Variable FloquetObservableTDspinKpoints
  !%Type block
  !%Default none
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Indices of k-points for td-spin calculation  
  !%End

  if(parse_block('FloquetObservableTDspinKpoints', blk) == 0) then
    nk_input = parse_block_cols(blk, 0)
    if(nk_input < 1 ) call messages_input_error('FloquetObservableTDspinKpoints')
    SAFE_ALLOCATE(kpoints_input(nk_input))
    do ik = 1, nk_input
      call parse_block_integer(blk, 0, ik-1, kpoints_input(ik))
    end do
    call parse_block_end(blk)
  end if
  write(message(1),'(a,i3,a)') 'Info: Read ', nk_input ,' kpoints for td-spin calculation'
  call messages_info(1)
  
  !%Variable FloquetObservableTDspinStates
  !%Type block
  !%Default none
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Indices of Floquet eigenstates for td-spin calculation
  !%End

  if(parse_block('FloquetObservableTDspinStates', blk) == 0) then
    nst_input = parse_block_cols(blk, 0)
    if(nst_input < 1 ) call messages_input_error('FloquetObservableTDspinStates')
    SAFE_ALLOCATE(states_input(nst_input))
    do ist = 1, nst_input
      call parse_block_integer(blk, 0, ist-1, states_input(ist))
    end do
    call parse_block_end(blk)
  end if
  write(message(1),'(a,i3,a)') 'Info: Read ', nst_input ,' states for td-spin calculation'
  call messages_info(1)

  ! prepare restart structure
  call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_LOAD, &
                       sys%mc, ierr,sys%gr%der%mesh)
  call states_look(restart, nik, dim, nst, ierr)
  if(dim/=dressed_st%d%dim .or. nik/=sys%gr%sb%kpoints%reduced%npoints.or.nst/=dressed_st%nst) then
     write(message(1),'(a)') 'Did not find commensurate Floquet restart structure'
     call messages_fatal(1)
  end if

  SAFE_ALLOCATE(psi_t(1:sys%gr%der%mesh%np,sys%st%d%dim))
  SAFE_ALLOCATE(zpsi(1:sys%gr%mesh%np))
  SAFE_ALLOCATE(F_psi(1:sys%gr%mesh%np,dressed_st%d%dim))

  do ik=1,nk_input
    iik=kpoints_input(ik)
    write(str,'(I5)') iik
    
    do ist=1,nst_input
      iist=states_input(ist)
      write(str2,'(I5)') iist

      ! read the floquet wavefunction for states iik,iist
      do idim = 1, dressed_st%d%dim
         itot = idim + (iist-1)*dressed_st%d%dim +  (iik-1)*dressed_st%nst*dressed_st%d%dim 
         write(filename,'(i10.10)') itot
         call zrestart_read_mesh_function(restart, trim(adjustl(filename)), sys%gr%mesh, zpsi, ierr)
         print *, ierr, trim(adjustl(filename)), iik, iist
         F_psi(1:sys%gr%mesh%np,idim) = zpsi(1:sys%gr%mesh%np)
      end do
      filename = 'floquet_td_spin_ik_'//trim(adjustl(str))//'_ist_'//trim(adjustl(str2))
      iunit = io_open(FLOQUET_DIR//filename, action='write')
      write(iunit,'(a)') '# time Sx Sy Sz'

      do it=1,hm%F%nT
        time = hm%F%Tcycle/hm%F%nT*it
        call floquet_td_state(hm%F,sys%gr%mesh,F_psi,time,psi_t)
        spin(1:3) = state_spin(sys%gr%mesh, psi_t) 
        write(iunit,'(e12.6,2x,e12.6,2x,e12.6,2x,e12.6)') time, spin(1:3)
      end do

      close(iunit)
    end do
  end do

  call restart_end(restart)

  SAFE_DEALLOCATE_A(psi_t)
  SAFE_DEALLOCATE_A(zpsi)
  SAFE_DEALLOCATE_A(F_psi)
  SAFE_DEALLOCATE_A(kpoints_input)
  SAFE_DEALLOCATE_A(states_input)
  end subroutine calc_floquet_td_spin

  !--------------------------------
  subroutine floquet_td_state(F,mesh,F_psi,time,psi_t)
  
  type(floquet_t) :: F
  type(mesh_t)    :: mesh
  CMPLX   :: F_psi(:,:)
  integer :: ik, ist
  FLOAT   :: time
  CMPLX   :: psi_t(:,:)
  
  integer :: idim, im, imm
  
  psi_t = M_ZERO
  do im= F%order(1),F%order(2)
    imm = im - F%order(1) + 1
    do idim=1,F%spindim
      psi_t(1:mesh%np,idim)  =  psi_t(1:mesh%np,idim) + exp(-M_zI*im*F%omega*time)*F_psi(1:mesh%np,(imm-1)*F%spindim+idim)
    end do
  end do

  ! normalize td-state
  psi_t(1:mesh%np,1:F%spindim)  = M_ONE/zmf_nrm2(mesh,F%spindim,psi_t)*psi_t(1:mesh%np,1:F%spindim) 

  end subroutine floquet_td_state
  
  
  
  
  !--------------------------------------------
  subroutine calc_floquet_norms(mesh,kpoints,st,dressed_st,iter,floquet_dim)
    type(mesh_t), intent(in) :: mesh
    type(kpoints_t), intent(in) :: kpoints
    type(states_t), intent(in) :: st,dressed_st
    integer :: iter , floquet_dim

    integer :: maxiter, ik, in, im, ist, idim, ierr, nik, dim, nst, iunit
    CMPLX, allocatable :: temp_state1(:,:), temp_state2(:,:)
    FLOAT, allocatable :: norms(:,:,:)
    character(len=1024):: ik_name,iter_name, filename

    SAFE_ALLOCATE(temp_state1(1:mesh%np,st%d%dim))
    SAFE_ALLOCATE(temp_state2(1:mesh%np,floquet_dim*st%d%dim))
    SAFE_ALLOCATE(norms(1:kpoints%reduced%npoints,1:dressed_st%nst,floquet_dim))

    norms = M_ZERO
   

    do ik=st%d%kpt%start,st%d%kpt%end
      do in=1,floquet_dim
        do ist=st%st_start,st%st_end
          call states_get_state(dressed_st, mesh, (in-1)*st%nst+ist, ik,temp_state2)

          do im=1,floquet_dim
            do idim=1,st%d%dim
              temp_state1(1:mesh%np,idim) = temp_state2(1:mesh%np,(im-1)*st%d%dim+idim)
            end do
            norms(ik,(in-1)*st%nst+ist,im) = zmf_nrm2(mesh,st%d%dim,temp_state1)
          end do
        enddo
      enddo
    enddo
    
    call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm, norms)

    if(mpi_world%rank==0) then
       write(iter_name,'(i4)') iter 
       do ik=kpoints%reduced%npoints-kpoints%nik_skip+1, kpoints%reduced%npoints
          write(ik_name,'(i4)') ik
          filename = FLOQUET_DIR//'/floquet_norms_ik_'//trim(adjustl(ik_name))//'_iter_'//trim(adjustl(iter_name))
          iunit = io_open(filename, action='write')
          
           do ist=1,dressed_st%nst
              do in=1,floquet_dim
                 write(iunit,'(i4,a1,e12.6,a1,i4,a1,e12.6)') ist, '', dressed_st%eigenval(ist,ik), ' ',in, ' ', norms(ik,ist,in)
             end do
             write(iunit,'(a1)') ' '
          end do

          call io_close(iunit)
       end do
    end if

    SAFE_DEALLOCATE_A(temp_state1)
    SAFE_DEALLOCATE_A(temp_state2)
    SAFE_DEALLOCATE_A(norms)
    
  end subroutine calc_floquet_norms
  
  
  subroutine calc_floquet_hhg_weights()
    
    CMPLX, allocatable   :: mel(:,:), u_ma(:,:), u_nb(:,:)
    FLOAT, allocatable   :: ediff(:), weight(:), fab(:)
    integer, allocatable :: mm(:), nn(:), alpha(:), beta(:), idx(:), kk(:)
    
    integer :: idim, im, in, ista, istb, ik, itot, ii, i, imm, inn, spindim
    FLOAT   :: omega, DE, wmax
    CMPLX   :: tmp(1:4), tmp2(1:3,1:hm%F%spindim)
    integer :: iunit, wpow
    character(len=1024):: filename, iter_name
    
    PUSH_SUB(calc_floquet_hhg_weights)
    
    spindim = hm%F%spindim 
    
    itot = dressed_st%d%kpt%nglobal * dressed_st%nst**2 * hm%F%floquet_dim**2
    
    select case (obs%gauge)
    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
      wpow = 2 
    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGHT)
      wpow = 4      
    end select    
    
    SAFE_ALLOCATE(mel(1:itot, 1:3))
    SAFE_ALLOCATE(ediff(1:itot))
    SAFE_ALLOCATE(weight(1:itot))
    SAFE_ALLOCATE(fab(1:itot))
    SAFE_ALLOCATE(mm(1:itot))
    SAFE_ALLOCATE(nn(1:itot))
    SAFE_ALLOCATE(alpha(1:itot))
    SAFE_ALLOCATE(beta(1:itot))
    SAFE_ALLOCATE(idx(1:itot))
    SAFE_ALLOCATE(kk(1:itot))
    
    SAFE_ALLOCATE(u_ma(1:sys%gr%mesh%np, hm%F%floquet_dim))
    SAFE_ALLOCATE(u_nb(1:sys%gr%mesh%np, hm%F%floquet_dim))

    mel(:,:)   = M_z0
    ediff(:) = M_ZERO
    fab(:) = M_ZERO
    weight(:) = M_ZERO
    mm(:) = 0
    nn(:) = 0
    kk(:) = 0 
    alpha(:) = 0
    beta(:) = 0
    idx(:) = 0
    
    ii=1
    do ik=1, dressed_st%d%nik

        if(dressed_st%d%kpt%start > ik .or. ik > dressed_st%d%kpt%end) cycle
          
        do ista=1, dressed_st%nst
          call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_ma)
          
          do istb=1, dressed_st%nst
            call states_get_state(dressed_st, sys%gr%mesh, istb, ik, u_nb)

            
            DE = dressed_st%eigenval(istb,ik) - dressed_st%eigenval(ista,ik)
!             print *, ista, istb, ik, DE, dressed_st%eigenval(istb,ik), dressed_st%eigenval(ista,ik)
            
            do in=hm%F%order(1),hm%F%order(2)
              inn = in - hm%F%order(1) + 1
                       
              do im=hm%F%order(1),hm%F%order(2)
                imm = im - hm%F%order(1) + 1

                alpha(ii) = ista 
                beta(ii) = istb 
                nn(ii) = in  
                mm(ii) = im 
                kk(ii) = ik

                ediff(ii) = DE + (in-im)*hm%F%omega
                
                ! get the dipole matrix elements 
                  select case (obs%gauge)
                  case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
                    call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                                               u_ma(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) ,&
                                               u_nb(:,(inn-1)*spindim+1: (inn-1)*spindim +spindim) , &
                                               ik, tmp2(:,:))
                    mel(ii,:) = M_z0
                    do idim=1,spindim
                      mel(ii,1:3) = mel(ii, 1:3) + tmp2(1:3,idim)/ (hm%F%order(2)-hm%F%order(1))
                    end do
                      
                  case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGHT)
                    mel(ii,:) = M_z0
                    do idim=1,spindim
                      call zmf_multipoles(sys%gr%mesh, conjg(u_ma(:,(imm-1)*spindim+idim)) &
                                                           * u_nb(:,(inn-1)*spindim+idim) , 1, tmp(:))

                      mel(ii,1:3) = mel(ii, 1:3) + tmp(2:4)/ (hm%F%order(2)-hm%F%order(1))
                    end do
                  end select
                  
                
                fab(ii) = dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik)
                weight(ii) = sum(abs(mel(ii,1:3))**2) * dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik) *ediff(ii)**wpow 
                
                weight(ii) = weight(ii) * dressed_st%d%kweights(ik)
                
                ii = ii + 1
              end do
            end do
            
          end do    
        end do    
        
    end do
    
    if(dressed_st%parallel_in_states .or. dressed_st%d%kpt%parallel) then
      call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  weight(:))
      call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  alpha(:))
      call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  beta(:))
      call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  nn(:))
      call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  mm(:))
      call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  kk(:))
    end if
    
    
    call sort(ediff(:), idx(:))
    
    wmax = maxval(weight(:))
    
    if(mpi_grp_is_root(mpi_world)) then
      write(iter_name,'(i4)') hm%F%iter
      filename = FLOQUET_DIR//'/floquet_hhg_w_'//trim(adjustl(iter_name))
      iunit = io_open(filename, action='write')
    
      write(iunit, '(a1,6a15)', advance ='no') '#', str_center("w", 15), str_center("strenght(w)", 15),&
                                                    str_center("alpha", 15), str_center("beta", 15),&
                                                    str_center("m", 15), str_center("n", 15)
                                               
      if (dressed_st%d%nik > 1 ) then
        write(iunit, '(a1,1a15)', advance ='no')  str_center("ikpt", 15)
      end if
    
    
      select case (obs%gauge)

      case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
        write(iunit, '(a1,4a15)')  str_center("|<u_ma|ix|u_nb>|", 15), str_center("|<u_ma|iy|u_nb>|", 15),&
                                    str_center("|<u_ma|iz|u_nb>|", 15), str_center("f_a*f_b", 15)
          
      case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGHT)

        write(iunit, '(a1,4a15)')  str_center("|<u_ma|x|u_nb>|", 15), str_center("|<u_ma|y|u_nb>|", 15),&
                                    str_center("|<u_ma|z|u_nb>|", 15), str_center("f_a*f_b", 15)
      end select
    
     
      write(iunit, '(a1,2a15)') &
        '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 15), &
        str_center('['//trim(units_abbrev(units_out%length))//'/(' &
          //trim(units_abbrev(units_out%time**2))//')]', 15)
      
      do ii = 1, itot
        if (weight(idx(ii))/wmax > CNST(1E-6) .and.  ediff(ii) > M_ZERO) then
          write(iunit, '(1x,2es15.6)', advance='no') units_from_atomic(units_out%energy, ediff(ii)) , weight(idx(ii))  
          write(iunit, '(4i15)', advance='no') alpha(idx(ii)) , beta(idx(ii)), mm(idx(ii)), nn(idx(ii))
          if (dressed_st%d%nik > 1 ) then
             write(iunit, '(1i15)', advance='no') kk(idx(ii))
          end if
          write(iunit, '(4es15.6)') abs(mel(idx(ii),1)),abs(mel(idx(ii),2)),abs(mel(idx(ii),3)), fab(idx(ii))
        end if
      end do  
    
      call io_close(iunit)
    end if
    
    SAFE_DEALLOCATE_A(mel)
    SAFE_DEALLOCATE_A(ediff)
    SAFE_DEALLOCATE_A(weight)
    SAFE_DEALLOCATE_A(fab)
    SAFE_DEALLOCATE_A(mm)
    SAFE_DEALLOCATE_A(nn)
    SAFE_DEALLOCATE_A(kk)
    SAFE_DEALLOCATE_A(alpha)
    SAFE_DEALLOCATE_A(beta)
    SAFE_DEALLOCATE_A(idx)
    
    SAFE_DEALLOCATE_A(u_ma)
    SAFE_DEALLOCATE_A(u_nb)

    POP_SUB(calc_floquet_hhg_weights)    
  end subroutine calc_floquet_hhg_weights

  subroutine calc_floquet_hhg()
    
    CMPLX, allocatable   :: u_ma(:,:), u_nb(:,:)
    FLOAT, allocatable   :: spect(:,:)
    
    integer :: idim, im, in, ista, istb, ik, itot, ii, i, imm, inn, spindim, ie, dim
    FLOAT   :: omega, DE, ediff, EE
    CMPLX   :: tmp(1:4), ampl(1:3), mel(1:3), tmp2(1:3,1:hm%F%spindim)
    integer :: iunit, wpow
    character(len=1024):: filename, iter_name
    
    PUSH_SUB(calc_floquet_hhg)
    
    spindim = hm%F%spindim 
    dim = sys%gr%sb%dim

    select case (obs%gauge)
    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
      wpow = 2 
    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGHT)
      wpow = 4      
    end select    
    
    
    itot = dressed_st%d%kpt%nglobal * dressed_st%nst**2 * hm%F%floquet_dim**2
    
    SAFE_ALLOCATE(spect(1:obs%ne,1:3))
    
    SAFE_ALLOCATE(u_ma(1:sys%gr%mesh%np, hm%F%floquet_dim))
    SAFE_ALLOCATE(u_nb(1:sys%gr%mesh%np, hm%F%floquet_dim))

    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, obs%ne)
    
    spect(:,:) = M_ZERO

    do ie = 1, obs%ne
      EE= ie * obs%de

      ampl(:) = M_z0
    
      do ik=dressed_st%d%kpt%start, dressed_st%d%kpt%end

          do ista=dressed_st%st_start, dressed_st%st_end
            call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_ma)
          
            do istb=dressed_st%st_start, dressed_st%st_end
              
              !Cut out all the component suppressed by small occupations 
              if (dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik) < 1E-14) cycle

              call states_get_state(dressed_st, sys%gr%mesh, istb, ik, u_nb)
            
              DE = dressed_st%eigenval(istb,ik) - dressed_st%eigenval(ista,ik)
              
            
              do in=hm%F%order(1),hm%F%order(2)
                inn = in - hm%F%order(1) + 1
                       
                do im=hm%F%order(1),hm%F%order(2)
                  imm = im - hm%F%order(1) + 1

                  ediff = DE + (in-im)*hm%F%omega
                  
                
                  ! get the dipole matrix elements 
                  select case (obs%gauge)
                  case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
                    call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                                               u_ma(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) ,&
                                               u_nb(:,(inn-1)*spindim+1: (inn-1)*spindim +spindim) , &
                                               ik, tmp2(:,:))
                    mel(:) = M_z0
                    do idim=1,spindim
                      mel(1:3) = mel( 1:3) + tmp2(1:3,idim)/ (hm%F%order(2)-hm%F%order(1))
                    end do
                      
                  case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGHT)                  
                    mel(:) = M_z0
                    do idim=1,spindim
                      call zmf_multipoles(sys%gr%mesh, conjg(u_ma(:,(imm-1)*spindim+idim)) &
                                                           * u_nb(:,(inn-1)*spindim+idim) , 1, tmp(:))
                  
                      mel(1:dim) = mel(1:dim) + tmp(2:2+dim-1)/(hm%F%order(2)-hm%F%order(1))
                    end do
                  end select
                  
                  ampl(1:dim) = ampl(1:dim)+ mel(1:dim) * dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik) &
                                         * aimag(1/(EE - ediff  + M_zi*obs%gamma))
!                   print *, ampl(1:3), ediff
                
                end do
                
              end do
            end do    
            
          end do    
          
          ampl(1:dim) = ampl(1:dim) * (1 + dressed_st%d%kweights(ik))

      end do
      
      if(dressed_st%parallel_in_states .or. dressed_st%d%kpt%parallel) then
        call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  ampl(:))
      end if
      
      spect(ie,1:dim) = abs(ampl(1:dim))**2 *EE**wpow 
      
      if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ie, obs%ne)
      
    end do
    if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")
    
    
    if(mpi_grp_is_root(mpi_world)) then
      write(iter_name,'(i4)') hm%F%iter
      filename = FLOQUET_DIR//'/floquet_hhg_'//trim(adjustl(iter_name))
      iunit = io_open(filename, action='write')
    
      write(iunit, '(a1,5a15)') '#', str_center("w", 15), str_center("I(w)", 15), &
                                     str_center("Ix(w)", 15), str_center("Iy(w)", 15), str_center("Ix(w)", 15)
                            
      write(iunit, '(a1,2a15)') &
        '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 15), &
        str_center('['//trim(units_abbrev(units_out%length))//'/(' &
          //trim(units_abbrev(units_out%time**2))//')]', 15)
      
      
      do ie = 1, obs%ne
        EE = ie * obs%de
        write(iunit, '(1x,5es15.6)') units_from_atomic(units_out%energy, EE) , sum(spect(ie,1:3)), &
                                    spect(ie,1), spect(ie,2), spect(ie,3)
      end do  
    
      call io_close(iunit)
    end if
    
    SAFE_DEALLOCATE_A(spect)
    
    SAFE_DEALLOCATE_A(u_ma)
    SAFE_DEALLOCATE_A(u_nb)

    POP_SUB(calc_floquet_hhg)    
  end subroutine calc_floquet_hhg


  subroutine out_floquet_wfs()
    

    integer :: ik, ist, fdim, spindim , im, imm
    integer :: nik, dim, nst, itot
    CMPLX, allocatable  ::  zpsi(:)
    character(len=512)  ::  str3, str4
    type(unit_t) :: fn_unit
    
    integer :: how
    
    PUSH_SUB(out_floquet_wfs)

    call io_function_read_how(sys%gr%sb, how, ignore_error = .true.)

    fn_unit = units_out%length**(-sys%gr%mesh%sb%dim)

    spindim = hm%F%spindim 


    ! prepare restart structure
    call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_LOAD, &
                         sys%mc, ierr,sys%gr%der%mesh)
    call states_look(restart, nik, dim, nst, ierr)
    if(dim/=dressed_st%d%dim .or. nik/=sys%gr%sb%kpoints%reduced%npoints.or.nst/=dressed_st%nst) then
       write(message(1),'(a)') 'Did not find commensurate Floquet restart structure'
       call messages_fatal(1)
    end if

    SAFE_ALLOCATE(zpsi(1:sys%gr%mesh%np))

    do ik=obs%nkpt(1), obs%nkpt(2)
      write(str,'(I5)') ik
    
      do ist=obs%nst(1), obs%nst(2)
        write(str2,'(I5)') ist

        ! read the floquet wavefunction for states ik,ist
        do im=hm%F%order(1),hm%F%order(2)
           imm = im - hm%F%order(1) + 1
           do idim=1,spindim
             fdim = (imm-1)*spindim+idim 
             itot = fdim + (ist-1)*dressed_st%d%dim +  (ik-1)*dressed_st%nst*dressed_st%d%dim
             write(filename,'(i10.10)') itot
             ierr = 0
             call zrestart_read_mesh_function(restart, trim(adjustl(filename)), sys%gr%mesh, zpsi, ierr)

             if (ierr > 0) then
               write(message(1),'(2a)') "Failed to read from restart file ", trim(filename)
               call messages_warning(2)
               cycle
             end if

             write(str3,'(I5)') im
             if(spindim>1) then
               filename = 'wfn_is_'//trim(adjustl(str4))//'_ik_'//trim(adjustl(str))//'_ist_'&
                                   //trim(adjustl(str2))//'_m_'//trim(adjustl(str3))
             else
               filename = 'wfn_ik_'//trim(adjustl(str))//'_ist_'//trim(adjustl(str2))//'_m_'//trim(adjustl(str3))
             end if
             
             call zio_function_output(how, FLOQUET_DIR, filename, sys%gr%mesh, zpsi, fn_unit, ierr)
           end do
        end do


      end do
    end do

    call restart_end(restart)

    SAFE_DEALLOCATE_A(zpsi)


    
    POP_SUB(out_floquet_wfs)
  end subroutine out_floquet_wfs


  subroutine calc_floquet_conductivity()
    
    CMPLX, allocatable   :: u_ma(:,:), u_mb(:,:), sigma(:,:,:)
    FLOAT, allocatable   :: spect(:,:)
    
    integer :: idim, im, in, ista, istb, ik, itot, ii, i, imm, inn, spindim, ie, dim
    FLOAT   :: omega, DE, ediff, EE, norm, fact
    CMPLX   :: melba(1:3), melab(1:3), tmp2(1:3,1:hm%F%spindim), ampl
    integer :: iunit, idir, jdir, dir
    character(len=1024):: filename, iter_name, str
    
    PUSH_SUB(calc_floquet_conductivity)
    
    spindim = hm%F%spindim 
    dim = sys%gr%sb%dim

    
    
    itot = dressed_st%d%kpt%nglobal * dressed_st%nst**2 * hm%F%floquet_dim**2
    
    SAFE_ALLOCATE(sigma(1:obs%ne,1:3,1:3))
    
    SAFE_ALLOCATE(u_ma(1:sys%gr%mesh%np, hm%F%floquet_dim))
    SAFE_ALLOCATE(u_mb(1:sys%gr%mesh%np, hm%F%floquet_dim))

    
    sigma(:,:,:) = M_z0

    do idir = 1, 2 
      do jdir = idir, 2 


        write(message(1),'(a,i5,a,i5,a)') 'Calculate sigma(',idir,',', jdir,'):'
        call messages_info(1)
        
        
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, obs%ne)
        
        do ie = 1, obs%ne
          EE= ie * obs%de

          sigma(ie,idir,jdir) = M_z0
    
          do ik=dressed_st%d%kpt%start, dressed_st%d%kpt%end

              do ista=dressed_st%st_start, dressed_st%st_end
                call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_ma)
      
                do istb=1, dressed_st%nst
            
                  if (ista == istb .and. idir /= jdir) cycle
                  
                  if (ista >= istb .and. idir == jdir) cycle
                  
                  DE = dressed_st%eigenval(istb,ik) - dressed_st%eigenval(ista,ik)
                  
                  ! Skip divergence 
                  if (DE < M_EPSILON) cycle                  
                  
                  !Cut out all the components suppressed by small occupations 
                  if ((dressed_st%occ(istb,ik) - dressed_st%occ(ista,ik))/DE < 1E-14) cycle

                  call states_get_state(dressed_st, sys%gr%mesh, istb, ik, u_mb)
          

                   
                  melab(:) = M_z0
                  melba(:) = M_z0
                                   
                  norm = (hm%F%order(2)-hm%F%order(1))
                  if (norm == 0) norm = M_ONE  
                                  
                  ! get the dipole matrix elements <<psi_a|J|psi_b >>
                  do im=hm%F%order(1),hm%F%order(2)
                    imm = im - hm%F%order(1) + 1
          
                    !<u_ma| J | u_mb>
                    call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                                               u_ma(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) ,&
                                               u_mb(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) , &
                                               ik, tmp2(:,:))
                    do dir=1,dim
                      melab(dir) = melab(dir) + sum(tmp2(dir,1:spindim))/ norm 
                    end do
                
                    !<u_mb| J | u_ma>
                    call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                                               u_mb(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) ,&
                                               u_ma(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) , &
                                               ik, tmp2(:,:))
                    do dir=1,dim
                      melba(dir) = melba(dir) + sum(tmp2(dir,1:spindim))/ norm 
                    end do
                
                  end do ! im loop

                  
                  ampl =  M_zI * (dressed_st%occ(istb,ik) - dressed_st%occ(ista,ik))/DE * &
                                  melab(idir)*melba(jdir) /(DE + EE + M_zi*obs%gamma) 
                          
                  if (idir == jdir) ampl = M_z2 * ampl
                  
                  ! sum over a and b         
                  sigma(ie,idir,jdir) = sigma(ie,idir,jdir) + ampl


                  
                end do ! istb loop   
              end do ! ista loop   

              ! sum over kpoints 
              sigma(ie,idir,jdir) = sigma(ie,idir,jdir) * (1 + dressed_st%d%kweights(ik))
             
            end do ! ik loop
          
          

          if(dressed_st%parallel_in_states .or. dressed_st%d%kpt%parallel) then
            call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  sigma(ie,idir,jdir))
          end if
                  
          if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(ie, obs%ne)
          
        end do ! ie loop

        if(mpi_grp_is_root(mpi_world)) write(*, "(1x)")
        
      end do            
    end do
    
    
    if(mpi_grp_is_root(mpi_world)) then
      write(iter_name,'(i4)') hm%F%iter
      filename = FLOQUET_DIR//'/floquet_conductivity_'//trim(adjustl(iter_name))
      iunit = io_open(filename, action='write')

      write(iunit, '(a1,a15)', advance = 'no') '#', str_center("w", 15)
    
      do idir = 1, 2 
        do jdir = idir, 2 
          write(str, '(a,i1,a,i1 ,a)') 'sigma(',idir,',', jdir,')'
          write(iunit, '(2a15)', advance = 'no' )  str_center('Re['//trim(str)//']', 15),& 
                                                   str_center('Im['//trim(str)//']', 15)
        end do 
      end do
      write(iunit, '(1x)')

                            
      write(iunit, '(a1,1a15)') &
            '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 15)
      
      
      do ie = 1, obs%ne
        EE = ie * obs%de
        
        write(iunit, '(1x,es15.6)', advance = 'no') units_from_atomic(units_out%energy, EE)         
        do idir = 1, 2 
          do jdir = idir, 2 
        
            write(iunit, '(2es15.6)', advance ='no') real(sigma(ie,idir,jdir)), aimag(sigma(ie,idir,jdir))
          end do
        end do
        write(iunit, '(1x)')
          
      end do  
    
      call io_close(iunit)
    end if
    
    SAFE_DEALLOCATE_A(sigma)
    
    SAFE_DEALLOCATE_A(u_ma)
    SAFE_DEALLOCATE_A(u_mb)

    POP_SUB(calc_floquet_conductivity)    
  end subroutine calc_floquet_conductivity



  end program floquet_observables

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
