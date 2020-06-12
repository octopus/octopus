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
  use cube_oct_m
  use current_oct_m
  use density_oct_m
  use geometry_oct_m
  use fft_oct_m
  use floquet_oct_m
  use forces_oct_m
  use gauge_field_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use kpoints_oct_m
  use io_binary_oct_m
  use io_function_oct_m
  use io_oct_m
  use loct_oct_m
  use loct_math_oct_m
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
  use species_oct_m
  use string_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_io_oct_m
  use states_restart_oct_m
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
    integer :: fc_method ! the method used to calculate floquet conductivity
    FLOAT   :: time0  ! probe time
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
  type(states_t), pointer :: bare_st
  type(states_t)          :: FBZ_st

  type(fobs_t)         :: obs

  integer              :: iter
  FLOAT                :: time, time_step
  
  logical              :: FBZ_st_initialized

  call getopt_init(ierr)
  if(ierr /= 0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-floquet-observables command is not available."
    call messages_fatal(2)
  end if

  FBZ_st_initialized = .false.

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
                        family_is_mgga_with_exc(sys%ks%xc, sys%st%d%nspin))

  call floquet_init(sys,hm%F,sys%st%d%dim)
  bare_st => sys%st

  dressed_st%floquet_dim = hm%F%floquet_dim
  dressed_st%floquet_FBZ = hm%F%FBZ_solver
  call states_init(dressed_st, sys%gr, sys%geo)
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
  !%Option f_forces bit(9)
  !% Calculate Floquet forces.
  !%Option f_density_plot bit(10)
  !% Plot the Floquet density.
  !%Option f_FBZ_bands bit(11)
  !% Plot the Floquet FBZ bands by filtering.
  !%End
  call parse_variable('FloquetObservableCalc', out_what, out_what)
  

  ! load floquet states only for certain tasks
  ! Note: this way it makes impossible to combine options that needs all the Floquet 
  ! states with the ones that do not UDG
  if(.not. (iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_TD_SPIN) /= 0) .and. &
     .not. (iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_WFS) /= 0) ) then
       call states_allocate_wfns(dressed_st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)
       call floquet_restart_dressed_st(hm, sys, dressed_st, ierr)
       call messages_write('Read Floquet restart files.')
       call messages_info()
       
       ! we are in the Floquet Brillouin zone (FBZ)  
       if(hm%F%FBZ_solver) then
         call messages_write('Calculate FBZ coefficients.')
         call messages_info()
         
         ! get groundstate states
         call states_allocate_wfns(bare_st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)
         call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
         if(ierr == 0) call states_load(restart, bare_st, sys%gr, ierr)
         if (ierr /= 0) then
            message(1) = 'Unable to read ground-state wavefunctions.'
            call messages_fatal(1)
         end if
         call restart_end(restart)
         
         call floquet_calc_FBZ_coefficients(hm, sys, dressed_st, bare_st, M_ZERO)
         call states_write_eigenvalues(stdout, dressed_st%nst, dressed_st, sys%gr%sb, &
                                       compact = .true.)     
       end if 
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
  !%Option f_length 2
  !% The length gauge: <i|r.E|j>.
  !%End
  obs%gauge = OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH
  if (gauge_field_is_applied(hm%ep%gfield) .or. simul_box_is_periodic(sys%gr%sb)) & 
      obs%gauge = OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY

  call parse_variable('FloquetObservableGauge', obs%gauge, obs%gauge)
  call messages_print_var_option(stdout,'FloquetObservableGauge', obs%gauge)


  !%Variable FloquetConductivityMethod
  !%Type flag
  !%Default guess
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Selects the method to be used for the calculation of optical conductivity. 
  !%Option f_Oka 1
  !% Use oka's equation
  !%Option f_FBZ 2
  !% Use 19
  !%Option f_FBZ21 3
  !% Use 21
  !%End
  obs%fc_method = OPTION__FLOQUETCONDUCTIVITYMETHOD__F_OKA

  call parse_variable('FloquetConductivityMethod', obs%fc_method, obs%fc_method)
  call messages_print_var_option(stdout,'FloquetConductivityMethod', obs%fc_method)

  
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


  !%Variable FloquetObservableProbeTime
  !%Type float
  !%Default 0.0
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Probe time.
  !%End
  call parse_variable('FloquetObservableProbeTime', M_ZERO, obs%time0, units_inp%time)
  call messages_print_var_value(stdout,'FloquetObservableProbeTime', obs%time0)
  
  

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

    call calc_floquet_norms(sys%gr%der%mesh,sys%gr%sb%kpoints,dressed_st, hm%F%iter,hm%F%floquet_dim)
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_ARPES) /= 0) then
    call messages_write('Calculate Floquet ARPES.')
    call messages_info()

    call calc_floquet_arpes()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_TD_SPIN) /= 0) then
    call messages_write('Calculate td-spin of Floquet eigenstates.')
    call messages_info()

    if(bare_st%d%dim /=2 ) then
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
    ! two different gauges
  
    if (hm%F%FBZ_solver .eqv..false. .and. obs%fc_method == OPTION__FLOQUETCONDUCTIVITYMETHOD__F_FBZ) then
      call get_FBZ_st(dressed_st, FBZ_st)
     ! get groundstate states
      call states_allocate_wfns(bare_st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)
      call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
      if(ierr == 0) call states_load(restart, bare_st, sys%gr, ierr)
      if (ierr /= 0) then
         message(1) = 'Unable to read ground-state wavefunctions.'
         call messages_fatal(1)
      end if
      call restart_end(restart)
      call floquet_calc_FBZ_coefficients(hm, sys, FBZ_st, bare_st, M_ZERO)
      call calc_floquet_conductivity(FBZ_st)
    else if (hm%F%FBZ_solver .eqv. .false. .and. obs%fc_method == OPTION__FLOQUETCONDUCTIVITYMETHOD__F_FBZ21) then
      call get_FBZ_st(dressed_st, FBZ_st)
      call calc_floquet_conductivity(FBZ_st)
    else
      call calc_floquet_conductivity(dressed_st)
    end if

  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_FORCES) /= 0) then
    call messages_write('Calculate Floquet forces.')
    call messages_info()

    call get_FBZ_st(dressed_st, FBZ_st)   
    call calc_floquet_forces()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_DENSITY_PLOT) /= 0) then
    call messages_write('Output Floquet density.')
    call messages_info()

    call plot_floquet_density()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_FBZ_BANDS) /= 0) then
    call messages_write('Output FBZ bands.')
    call messages_info()
    
    call get_FBZ_st(dressed_st, FBZ_st)     
          
    if (simul_box_is_periodic(sys%gr%sb) .and. kpoints_have_zero_weight_path(sys%gr%sb%kpoints)) then
      filename = 'FBZ_bands_occ'
      call states_write_bandstructure(FLOQUET_DIR, FBZ_st%nst, FBZ_st, sys%gr%sb, filename, vec = FBZ_st%occ)
      filename = 'FBZ_bands_coeff_Re'
      call states_write_bandstructure(FLOQUET_DIR, FBZ_st%nst, FBZ_st, sys%gr%sb, filename, vec = Real(FBZ_st%coeff))
      filename = 'FBZ_bands_coeff_Im'
      call states_write_bandstructure(FLOQUET_DIR, FBZ_st%nst, FBZ_st, sys%gr%sb, filename, vec = Aimag(FBZ_st%coeff))
    else
      filename = trim(FLOQUET_DIR)//'FBZ_eigenvalues'
      call states_write_eigenvalues(filename, FBZ_st%nst, FBZ_st, sys%gr%sb)
    end if
     
    
  end if

!   if (FBZ_st_initialized) call states_end(FBZ_st)

  
  call hamiltonian_end(hm)
  call system_end(sys)
  call fft_all_end()
  call io_end()
  call messages_end()
  call global_end()

contains 
  subroutine get_FBZ_st(dressed_st, FBZ_st)
    type(states_t), intent(in)  :: dressed_st
    type(states_t), intent(out) :: FBZ_st
    PUSH_SUB(get_FBZ_st)
    
    if (FBZ_st_initialized) then
      POP_SUB(get_FBZ_st)
      return
    end if
    
    
    if(.not.hm%F%FBZ_solver .and. hm%F%boson == OPTION__FLOQUETBOSON__CLASSICAL_FLOQUET) then
      ! this triggers FBZ allocation                                                                                                                                   
      hm%F%FBZ_solver = .true.
      FBZ_st%floquet_dim = hm%F%floquet_dim
      FBZ_st%floquet_FBZ = hm%F%FBZ_solver
      call states_init(FBZ_st, sys%gr, sys%geo)
      ! print *, 'FBZ_st dim t', FBZ_st%d%dim, dressed_st%d%dim, sys%st%d%dim
      call floquet_init_dressed_wfs(hm, sys, FBZ_st, fromScratch=.false.)
      ! reset flag
      hm%F%FBZ_solver = .false.
      call floquet_FBZ_filter(sys, hm, dressed_st, FBZ_st)

      call states_allocate_wfns(bare_st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)

      time = M_ZERO
      call messages_write('Reading td-state to calculate coefficients.')
      call messages_info()
      call restart_init(restart, RESTART_TD, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
      if(ierr == 0) call states_load(restart, bare_st, sys%gr, ierr, iter=iter, label = ": td")
      if (ierr /= 0) then
        message(1) = 'Unable to read time-dependent wavefunctions.'
        call messages_warning(1)
      end if
      if(ierr==0) then
        call parse_variable('TDTimeStep', M_ZERO, time_step, unit = units_inp%time)
        if(time_step == M_ZERO) then
           message(1) = 'Did not find time-step in Floquet init, please give a value for TDTimeStep'
          call messages_fatal(1)
        end if
        time = time_step*iter     
      end if
      call restart_end(restart)
      if(ierr/=0) then
        ierr = 0
        call messages_write('Reading gs-states to calculate coefficients.')
        call messages_info()
        call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
        if(ierr == 0) call states_load(restart, bare_st, sys%gr, ierr)
        if (ierr /= 0) then
          message(1) = 'Unable to read ground-state wavefunctions.'
          call messages_fatal(1)
        end if
        call restart_end(restart)
      end if
   
 
      call floquet_calc_FBZ_coefficients(hm, sys, FBZ_st, bare_st, time)
      
      FBZ_st_initialized = .true.

    end if
    
    POP_SUB(get_FBZ_st)
  end subroutine get_FBZ_st 


  subroutine calc_floquet_arpes()
  
  FLOAT, allocatable :: spect(:,:), me(:,:) 
  FLOAT :: pomega, pol(1:3)
  type(block_t)        :: blk  
  type(states_t)       :: d_st
  
  
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
    
  
!   if(.false.) then
!     SAFE_ALLOCATE(spect(dressed_st%nst,dressed_st%d%nik))
!     SAFE_ALLOCATE(me(dressed_st%nst,dressed_st%d%nik))
!     call floquet_photoelectron_spectrum(hm, sys, dressed_st, pomega, pol, spect, me)
!   end if


  ! need the big dressed states dimension for the output because even with the FBZ we are covering sidebands
  SAFE_ALLOCATE(spect(dressed_st%nst,dressed_st%d%nik))
  SAFE_ALLOCATE(me(dressed_st%nst,dressed_st%d%nik))
  call states_copy(d_st, dressed_st)
  call floquet_photoelectron_spectrum_FBZ(hm, sys, d_st, pomega, pol, spect, me)


  if(simul_box_is_periodic(sys%gr%sb)) then
    filename = 'floquet_arpes_me'
    call states_write_bandstructure(FLOQUET_DIR, d_st%nst, d_st, sys%gr%sb, filename, vec = me)

    filename = 'floquet_arpes'
    call states_write_bandstructure(FLOQUET_DIR, d_st%nst, d_st, sys%gr%sb, filename, vec = spect)
  end if

  SAFE_DEALLOCATE_A(spect)
  SAFE_DEALLOCATE_A(me)
  call states_end(d_st)
  
  end subroutine calc_floquet_arpes  


  !---------------------------------------
  subroutine floquet_photoelectron_spectrum_FBZ(hm, sys, dressed_st, pomega, pol, spect, me)
    type(hamiltonian_t), intent(in) :: hm
    type(system_t), intent(in)      :: sys
    type(states_t), intent(inout)   :: dressed_st
    FLOAT,          intent(in)      :: pomega     ! Probe field energy 
    FLOAT,          intent(in)      :: pol(:)     ! Probe field polarization vector
    FLOAT,          intent(out)     :: spect(:,:) ! the photoelectron spectrum
    FLOAT,          intent(out)     :: me(:,:)    ! the photoelectron matrix elements


    CMPLX, allocatable :: u_a(:,:),  phase(:), tmp(:), mec(:,:,:),spectc(:,:,:)
    FLOAT :: Omega, dt , pp(1:3), kpt(1:3), xx(1:MAX_DIM), gg(1:3), AA(1:3)
    integer :: idx, im, it, ist, idim, nT, ik, ia, Fdim(2), imm, dim, pdim, ip, spindim
    integer :: il, igx,igy, ll, inl, in, ial, box(1:3)

    type(cube_t):: cube
    FLOAT, allocatable :: Lg(:,:)

    type(mesh_t),   pointer :: mesh
      

    PUSH_SUB(floquet_photoelectron_spectrum_FBZ)

    mesh  => sys%gr%der%mesh
    
    dt=hm%F%dt
    nT=hm%F%nT
    Omega=hm%F%omega
    Fdim(:)=hm%F%order(:)
    spindim = hm%F%spindim 
    
    dim = mesh%sb%dim
    pdim = mesh%sb%periodic_dim
    
    SAFE_ALLOCATE(phase(1:mesh%np))
    SAFE_ALLOCATE(u_a(1:mesh%np,hm%F%floquet_dim))
    SAFE_ALLOCATE(tmp(spindim))

    SAFE_ALLOCATE(spectc(dressed_st%nst,dressed_st%d%nik, spindim))
    SAFE_ALLOCATE(mec(dressed_st%nst,dressed_st%d%nik, spindim))

    

    kpt(:) = M_ZERO
    pp(:)  = M_ZERO
    gg(:)  = M_ZERO
    me(:,:) = M_ZERO
    spect(:,:) = M_ZERO
    mec(:,:,:) = M_z0
    spectc(:,:,:) = M_z0
    
    AA(:)=M_ZERO
    AA(1)=M_ONE
    
    ! get fourier grid for gpoints 
    box(1:pdim) = mesh%idx%ll(1:pdim)
    box(dim) = 1
    call cube_init(cube,box, mesh%sb, fft_type = FFTLIB_FFTW, verbose = .true., spacing = mesh%spacing)
    SAFE_ALLOCATE(Lg(1:maxval(cube%fs_n_global(:)),1:3))
    print *, cube%fs_n_global(:)
    print *, cube%rs_n_global(:)
    print *,  cube%Lfs(:,1)
    print *,  cube%Lrs(:,1)
    do ii=1,maxval(cube%fs_n_global(:))   
      Lg(ii,1:dim)= matmul(mesh%sb%klattice_primitive(1:dim,1:dim), cube%Lfs(ii,1:dim))
    end do
    
        
    call get_FBZ_st(dressed_st, FBZ_st) 
  
  
    do ik=FBZ_st%d%kpt%start, FBZ_st%d%kpt%end
      kpt(1:dim) = kpoints_get_point(mesh%sb%kpoints, ik) 

      do ia=FBZ_st%st_start, FBZ_st%st_end
        call states_get_state(FBZ_st, mesh, ia, ik, u_a)
        
        do il = Fdim(1), Fdim(2)
          tmp(:) = M_ZERO

          do in = Fdim(1), Fdim(2)
            
            if ((in -il) > Fdim(2) .or. (in -il) < Fdim(1)) cycle
            inl = (in -il) - Fdim(1) + 1

            do igx = 1, cube%fs_n_global(1)
              gg(1)=Lg(igx,1)
              do igy = 1, cube%fs_n_global(2)
                gg(2)=Lg(igy,2)
          
                pp(dim) = M_TWO*(pomega + FBZ_st%eigenval(ia,ik) + il*Omega) - sum((kpt(1:pdim)+gg(1:pdim))**2)
                if (pp(dim) < M_ZERO) cycle  
                ! perpendicular (non-periodic) momentum component
                pp(dim) = sqrt(pp(dim))
                ! parallel (periodic) momentum component
                pp(1:pdim) = kpt(1:pdim)+gg(1:pdim)  
             
                do ip=1, mesh%np
                  xx=mesh_x_global(mesh, ip) 
                  phase(ip) = exp(-M_ZI*sum(pp(1:dim)*xx(1:dim)))
                end do
                
                ! Normalize to Jn(A/Omega.p)/sqrt(2pi)
                phase(:)=phase(:)/sqrt(M_TWO*M_PI)*loct_bessel(in, pp(1)/Omega)!sum(AA(1:dim)*pp(1:dim))/Omega)
                
                do idim=1,spindim
                  tmp(idim) = tmp(idim) + zmf_integrate(mesh, phase(1:mesh%np)*u_a(1:mesh%np,(inl-1)*spindim+idim))
                end do
              
                tmp(:) = tmp(:)* sum(pol(1:dim)*pp(1:dim))

              end do ! igy
            end do ! igx
          end do ! in
          
          ! combine ia and il to fill the full set of floquet eingenvalues (including replicas) 
          ial = ia + (il - Fdim(1))*FBZ_st%nst
!           print *, ial, ia, il, tmp(:)
          
          mec(ial,ik,:)    = tmp(:)    
          spectc(ial,ik,:) = tmp(:)*FBZ_st%coeff(ia,ik)
          dressed_st%eigenval(ial,ik)= FBZ_st%eigenval(ia,ik) + il*Omega
            
        end do ! il
      end do ! ia
    end do ! ik
    
    call comm_allreduce(FBZ_st%st_kpt_mpi_grp%comm, mec)
    call comm_allreduce(FBZ_st%st_kpt_mpi_grp%comm, spectc)
    
    do ik=dressed_st%d%kpt%start, dressed_st%d%kpt%end
      do ia=dressed_st%st_start, dressed_st%st_end
        me(ia,ik) = sum(abs(mec(ia,ik,1:spindim))**2) 
        spect(ia,ik) = sum(abs(spectc(ia,ik,1:spindim))**2) 
      end do
    end do
    
    ! don't forget to sort the spectra and me by energy for all kpoints 
    
    
    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(phase)
    SAFE_DEALLOCATE_A(u_a)
    SAFE_DEALLOCATE_A(mec)
    SAFE_DEALLOCATE_A(spectc)

    SAFE_DEALLOCATE_A(Lg)
    call cube_end(cube)
    call states_end(FBZ_st)
    


    POP_SUB(floquet_photoelectron_spectrum_FBZ)
    
  end subroutine floquet_photoelectron_spectrum_FBZ    

  !---------------------------------------
  subroutine floquet_photoelectron_spectrum(hm, sys, st, pomega, pol, spect, me)
    type(hamiltonian_t), intent(in) :: hm
    type(system_t), intent(in)      :: sys
    type(states_t), intent(in)      :: st
    FLOAT,          intent(in)      :: pomega     ! Probe field energy 
    FLOAT,          intent(in)      :: pol(:)     ! Probe field polarization vector
    FLOAT,          intent(out)     :: spect(:,:) ! the photoelectron spectrum
    FLOAT,          intent(out)     :: me(:,:)    ! the photoelectron matrix elements


    CMPLX, allocatable :: u_a(:,:),  phase(:), tmp(:)
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
    SAFE_ALLOCATE(u_a(1:mesh%np,hm%F%floquet_dim))
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
        
        call states_get_state(st, mesh, ia, ik, u_a)
        
        tmp(:) = M_ZERO
        do idx=hm%F%flat_idx%start, hm%F%flat_idx%end
          it = hm%F%idx_map(idx,1)
          im = hm%F%idx_map(idx,2)
          imm = im - Fdim(1) + 1
          do idim=1,spindim
            tmp(idim)  =  tmp(idim) + exp(-M_zI*im*omega*it*dt)/nT * &
                          zmf_integrate(mesh, phase(1:mesh%np)*u_a(1:mesh%np,(imm-1)*spindim+idim))
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
    SAFE_DEALLOCATE_A(u_a)



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
        call floquet_td_state(hm%F,sys%gr%mesh,F_psi, M_ZERO, time,psi_t)
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


  
  
  
  
  !--------------------------------------------
  subroutine calc_floquet_norms(mesh,kpoints,st,iter,floquet_dim)
    type(mesh_t), intent(in) :: mesh
    type(kpoints_t), intent(in) :: kpoints
    type(states_t), intent(in) :: st
    integer :: iter , floquet_dim

    integer :: maxiter, ik, in, im, ist, idim, ierr, nik, dim, nst, iunit
    CMPLX, allocatable :: temp_state1(:,:), temp_state2(:,:)
    FLOAT, allocatable :: norms(:,:,:)
    character(len=1024):: ik_name,iter_name, filename

    SAFE_ALLOCATE(temp_state1(1:mesh%np,hm%F%spindim))
    SAFE_ALLOCATE(temp_state2(1:mesh%np,st%d%dim))
    SAFE_ALLOCATE(norms(1:kpoints%reduced%npoints,1:st%nst,floquet_dim))

    norms = M_ZERO
   

    do ik=st%d%kpt%start,st%d%kpt%end
      do ist=st%st_start,st%st_end
          call states_get_state(st, mesh, ist, ik,temp_state2)

          do im=1,floquet_dim
            do idim=1,hm%F%spindim
              temp_state1(1:mesh%np,idim) = temp_state2(1:mesh%np,(im-1)*hm%F%spindim+idim)
            end do
            norms(ik,ist,im) = zmf_nrm2(mesh,hm%F%spindim,temp_state1)
          end do
          
      end do
    end do
    
    call comm_allreduce(st%st_kpt_mpi_grp%comm, norms)

    if(mpi_world%rank==0) then
      if (kpoints_have_zero_weight_path(kpoints)) then
        do ik=kpoints%reduced%npoints-kpoints%nik_skip+1, kpoints%reduced%npoints
          write(ik_name,'(i4)') ik
          filename = FLOQUET_DIR//'/norms_ik_'//trim(adjustl(ik_name))
          iunit = io_open(filename, action='write')
      
          do ist=1,st%nst
            do in=1,floquet_dim
              write(iunit,'(i4,a1,e12.6,a1,i4,a1,e12.6)') ist, '', st%eigenval(ist,ik), ' ',in, ' ', norms(ik,ist,in)
            end do
            write(iunit,'(a1)') ' '
          end do

          call io_close(iunit)
        end do
      else
        do ik=1,kpoints%reduced%npoints
          write(ik_name,'(i4)') ik
          filename = FLOQUET_DIR//'/norms_ik_'//trim(adjustl(ik_name))
          iunit = io_open(filename, action='write')
      
          do ist=1,st%nst
            do in=1,floquet_dim
              write(iunit,'(i4,a1,e12.6,a1,i4,a1,e12.6)') ist, '', st%eigenval(ist,ik), ' ',in, ' ', norms(ik,ist,in)
            end do
            write(iunit,'(a1)') ' '
          end do

          call io_close(iunit)
        end do        
        
      end if
    end if

    SAFE_DEALLOCATE_A(temp_state1)
    SAFE_DEALLOCATE_A(temp_state2)
    SAFE_DEALLOCATE_A(norms)
    
  end subroutine calc_floquet_norms
  
  
  subroutine calc_floquet_hhg_weights()
    
    CMPLX, allocatable   :: mel(:,:), u_a(:,:), u_nb(:,:)
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
    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)
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
    
    SAFE_ALLOCATE(u_a(1:sys%gr%mesh%np, hm%F%floquet_dim))
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
          call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_a)
          
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
                                               u_a(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) ,&
                                               u_nb(:,(inn-1)*spindim+1: (inn-1)*spindim +spindim) , &
                                               ik, tmp2(:,:))
                    mel(ii,:) = M_z0
                    do idim=1,spindim
                      mel(ii,1:3) = mel(ii, 1:3) + tmp2(1:3,idim)/ (hm%F%order(2)-hm%F%order(1))
                    end do
                      
                  case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)
                    mel(ii,:) = M_z0
                    do idim=1,spindim
                      call zmf_multipoles(sys%gr%mesh, conjg(u_a(:,(imm-1)*spindim+idim)) &
                                                           * u_nb(:,(inn-1)*spindim+idim) , 1, tmp(:))

                      mel(ii,1:3) = mel(ii, 1:3) + tmp(2:4)/ (hm%F%order(2)-hm%F%order(1))
                    end do
                  end select
                  
                
                fab(ii) = dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik)
                weight(ii) = sum(abs(mel(ii,1:3))**2) * dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik) !*ediff(ii)**wpow 
                
                weight(ii) = weight(ii) * dressed_st%d%kweights(ik)
                
                ii = ii + 1
              end do
            end do
            
          end do    
        end do    
        
    end do
    
    if(dressed_st%parallel_in_states .or. dressed_st%d%kpt%parallel) then
      call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  ediff(:))
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
                                               
      if (dressed_st%d%kpt%nglobal > 1 ) then
        write(iunit, '(a1,1a15)', advance ='no')  str_center("ikpt", 15)
      end if
    
    
      select case (obs%gauge)

      case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
        write(iunit, '(a1,4a15)')  str_center("|<u_a|jx|u_nb>|", 15), str_center("|<u_a|jy|u_nb>|", 15),&
                                    str_center("|<u_a|jz|u_nb>|", 15), str_center("f_a*f_b", 15)
          
      case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)

        write(iunit, '(a1,4a15)')  str_center("|<u_a|x|u_nb>|", 15), str_center("|<u_a|y|u_nb>|", 15),&
                                    str_center("|<u_a|z|u_nb>|", 15), str_center("f_a*f_b", 15)
      end select
    
     
      write(iunit, '(a1,2a15)') &
        '#', str_center('['//trim(units_abbrev(units_out%energy)) // ']', 15), &
        str_center('['//trim(units_abbrev(units_out%length))//'/(' &
          //trim(units_abbrev(units_out%time**2))//')]', 15)
      
      do ii = 1, itot
        if (weight(idx(ii))/wmax > CNST(1E-6) .and.  ediff(ii) > M_ZERO) then
          write(iunit, '(1x,2es15.6)', advance='no') units_from_atomic(units_out%energy, ediff(ii)) , weight(idx(ii))  
          write(iunit, '(4i15)', advance='no') alpha(idx(ii)) , beta(idx(ii)), mm(idx(ii)), nn(idx(ii))
          if (dressed_st%d%kpt%nglobal > 1 ) then
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
    
    SAFE_DEALLOCATE_A(u_a)
    SAFE_DEALLOCATE_A(u_nb)

    POP_SUB(calc_floquet_hhg_weights)    
  end subroutine calc_floquet_hhg_weights

  subroutine calc_floquet_hhg()
    
    CMPLX, allocatable   :: u_a(:,:), u_nb(:,:)
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
    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)
      wpow = 4      
    end select    
    
    
    itot = dressed_st%d%kpt%nglobal * dressed_st%nst**2 * hm%F%floquet_dim**2
    
    SAFE_ALLOCATE(spect(1:obs%ne,1:3))
    
    SAFE_ALLOCATE(u_a(1:sys%gr%mesh%np, hm%F%floquet_dim))
    SAFE_ALLOCATE(u_nb(1:sys%gr%mesh%np, hm%F%floquet_dim))

    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, obs%ne)
    
    spect(:,:) = M_ZERO

    do ie = 1, obs%ne
      EE= ie * obs%de

      ampl(:) = M_z0
    
      do ik=dressed_st%d%kpt%start, dressed_st%d%kpt%end

          do ista=dressed_st%st_start, dressed_st%st_end
            call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_a)
          
            do istb=dressed_st%st_start, dressed_st%st_end
              
              !Cut out all the component suppressed by small occupations 
              if (dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik) < CNST(1E-7)) cycle

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
                                               u_a(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) ,&
                                               u_nb(:,(inn-1)*spindim+1: (inn-1)*spindim +spindim) , &
                                               ik, tmp2(:,:))
                    mel(:) = M_z0
                    do idim=1,spindim
                      mel(1:3) = mel( 1:3) + tmp2(1:3,idim)/ (hm%F%order(2)-hm%F%order(1))
                    end do
                      
                  case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)                  
                    mel(:) = M_z0
                    do idim=1,spindim
                      call zmf_multipoles(sys%gr%mesh, conjg(u_a(:,(imm-1)*spindim+idim)) &
                                                           * u_nb(:,(inn-1)*spindim+idim) , 1, tmp(:))
                  
                      mel(1:dim) = mel(1:dim) + tmp(2:2+dim-1)/(hm%F%order(2)-hm%F%order(1))
                    end do
                  end select
                  
                  ampl(1:dim) = ampl(1:dim)+ mel(1:dim) * dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik) &
                                         * aimag(1/(ediff - EE  + M_zi*obs%gamma))
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
      
      if (obs%gauge == OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY) then
        ! Add the contribute of the fundamental harmonic
        ampl(1:dim) = ampl(1:dim)+ mel(1:dim) * dressed_st%occ(istb,ik) * dressed_st%occ(ista,ik) &
                               * aimag(1/(EE - ediff  + M_zi*obs%gamma))
        
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
    
    SAFE_DEALLOCATE_A(u_a)
    SAFE_DEALLOCATE_A(u_nb)

    POP_SUB(calc_floquet_hhg)    
  end subroutine calc_floquet_hhg


  subroutine out_floquet_wfs()
    

    integer :: ik, ist, fdim, spindim , im, imm
    integer :: nik, dim, nst, itot
    CMPLX, allocatable  ::  zpsi(:)
    character(len=512)  ::  str3, str4
    type(unit_t) :: fn_unit
    
    integer(8) :: how
    
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
 

  subroutine calc_floquet_conductivity(dressed_st)
    type(states_t), intent(in) :: dressed_st 
    CMPLX, allocatable   :: u_a(:,:), u_b(:,:), u_c(:,:),sigma(:,:,:)
    FLOAT, allocatable   :: spect(:,:)
    
    integer :: idim, im, in,il, ista, istb,istc, ik, itot, ii, i, imm, inn, ill, spindim, ie, dim
    FLOAT   :: omega, DE, DE_ca, DE_ab, DE_bc, ediff, EE, norm, fact
    CMPLX   :: melba(1:3), melab(1:3),dm_ac(1:3), dn_cb(1:3), dm_ab(1:3), dn_ba(1:3)
    CMPLX   :: tmp2(1:3,1:hm%F%spindim), ampl, tmp(1:4)
    integer :: iunit, idir, jdir, dir
    character(len=1024):: filename, iter_name, str, gauge, methodname
    
    PUSH_SUB(calc_floquet_conductivity)
    
    spindim = hm%F%spindim 
    dim = sys%gr%sb%dim
    SAFE_ALLOCATE(sigma(1:obs%ne,1:3,1:3))
    print *, 'forder start end',hm%F%order(1), hm%F%order(2)

    select case (obs%fc_method)

    case (OPTION__FLOQUETCONDUCTIVITYMETHOD__F_OKA) 
      SAFE_ALLOCATE(u_a(1:sys%gr%mesh%np, hm%F%floquet_dim))
      SAFE_ALLOCATE(u_b(1:sys%gr%mesh%np, hm%F%floquet_dim))
      
      sigma(:,:,:) = M_z0

      do ik=dressed_st%d%kpt%start, dressed_st%d%kpt%end
        do ista=dressed_st%st_start, dressed_st%st_end
          call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_a)
          do istb=1, dressed_st%nst
            call states_get_state(dressed_st, sys%gr%mesh, istb, ik, u_b)
            DE = dressed_st%eigenval(istb,ik) - dressed_st%eigenval(ista,ik)

            melab(:) = M_z0
            melba(:) = M_z0
                             
            if (ista == istb) cycle
            !Cut out all the components suppressed by small occupations 
            if (abs((dressed_st%occ(istb,ik) - dressed_st%occ(ista,ik))/DE) < 1E-14) cycle
            ! get the dipole matrix elements <<psi_a|J|psi_b >>
            ! get the dipole matrix elements <<psi_a|r|psi_b >>
            do im=hm%F%order(1),hm%F%order(2)
              imm = im - hm%F%order(1) + 1
      
              select case (obs%gauge)

              case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
                !<u_a| J | u_b>
                call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                                           u_a(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) ,&
                                           u_b(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) , &
                                           ik, tmp2(:,:))
                do dir=1,dim
                  melab(dir) = melab(dir) + sum(tmp2(dir,1:spindim))/DE  
                end do
                !<u_b| J | u_a>
                call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                                           u_b(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) ,&
                                           u_a(:,(imm-1)*spindim+1: (imm-1)*spindim +spindim) , &
                                           ik, tmp2(:,:))
                do dir=1,dim
                  melba(dir) = melba(dir) + sum(tmp2(dir,1:spindim))/DE 
                end do
              case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)
                !<u_a| r | u_b> and !<u_b| r | u_a>
                do idim=1,spindim
                  call zmf_multipoles(sys%gr%mesh, conjg(u_a(:,(imm-1)*spindim+idim)) &
                                                       * u_b(:,(imm-1)*spindim+idim) , 1, tmp(:))

                  melab(:) = melab(:) + tmp(2:4)
                  call zmf_multipoles(sys%gr%mesh, conjg(u_b(:,(imm-1)*spindim+idim)) &
                                                       * u_a(:,(imm-1)*spindim+idim) , 1, tmp(:))

                  melba(:) = melba(:) + tmp(2:4)
                end do 
                    
              end select
            end do ! im loop

            do ie = 1, obs%ne
              EE= ie * obs%de
              do idir = 1,dim
                do jdir = idir, dim
            
                  ampl =  M_zI * &
                           (dressed_st%occ(istb,ik) - dressed_st%occ(ista,ik)) * &
                            melab(idir)*melba(jdir) * &
                           1/(DE + EE + M_zi*obs%gamma) 
                          
                  if (idir == jdir) ampl = M_z2 * ampl
                  
                  ! sum over a and b         
              
                  sigma(ie,idir,jdir) = sigma(ie,idir,jdir) + ampl * dressed_st%d%kweights(ik)
                end do ! jdir
              end do ! idir
            end do ! ie loop

          end do ! istb loop   
        end do ! ista loop   
      end do ! ik loop
      
      SAFE_DEALLOCATE_A(u_a)
      SAFE_DEALLOCATE_A(u_b)
      
      if(dressed_st%parallel_in_states .or. dressed_st%d%kpt%parallel) then
        call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  sigma(:,:,:))
      end if



    case (OPTION__FLOQUETCONDUCTIVITYMETHOD__F_FBZ)
          
      SAFE_ALLOCATE(u_a(1:sys%gr%mesh%np, hm%F%floquet_dim))
      SAFE_ALLOCATE(u_b(1:sys%gr%mesh%np, hm%F%floquet_dim))
      SAFE_ALLOCATE(u_c(1:sys%gr%mesh%np, hm%F%floquet_dim))
      
      sigma(:,:,:) = M_z0
     
      do ik=dressed_st%d%kpt%start, dressed_st%d%kpt%end

        do ista=dressed_st%st_start, dressed_st%st_end
          ! call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_a)
          call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_a)
          do istb=1, dressed_st%nst
            !istb = ista
            call states_get_state(dressed_st, sys%gr%mesh, istb, ik, u_b)
            !print *, 'occ',conjg(dressed_st%coeff(ista,ik)) *dressed_st%coeff(istb,ik)
            do istc=1, dressed_st%nst
              
              DE_bc = dressed_st%eigenval(istb,ik) - dressed_st%eigenval(istc,ik)
              DE_ca = dressed_st%eigenval(istc,ik) - dressed_st%eigenval(ista,ik)
              
              call states_get_state(dressed_st, sys%gr%mesh, istc, ik, u_c)
      
                               
              ! get the dipole matrix elements dn_cb and dm_ac
              do im=hm%F%order(1),hm%F%order(2)
                !if (abs(DE_ca-im*hm%F%omega)<1E-14) cycle
                imm = im - hm%F%order(1) + 1
                
                !\sum_l <u_(l-m)a| r | u_lc>
                dm_ac = M_z0
                do il=hm%F%order(1), hm%F%order(2)
                  ill = il - hm%F%order(1) + 1
                  if ( -im+il .ge. hm%F%order(1) .and. -im+il .le. hm%F%order(2)) then
                    select case (obs%gauge)
                    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
                      !<u_a| J | u_b>
                      call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                              u_a(:,(-im+il-hm%F%order(1))*spindim+1: (-im+il-hm%F%order(1))*spindim +spindim) ,&
                              u_c(:,(ill-1)*spindim+1: (ill-1)*spindim +spindim) , ik, tmp2(:,:))
                      do dir=1,dim
                        dm_ac(dir) = dm_ac(dir) + sum(tmp2(dir,1:spindim)) !/(DE_ca -im * hm%F%omega) 
                      end do
                    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)
                      do idim=1,spindim
                        call zmf_multipoles(sys%gr%mesh, conjg(u_a(:,(-im+il-hm%F%order(1))*spindim+idim)) &
                                                             * u_c(:,(ill-1)*spindim+idim) , 1, tmp(:))

                        dm_ac(:) = dm_ac(:) + tmp(2:4)
                      end do ! idim loop
                    end select
 
                 end if
                
                end do ! il loop

                do in=hm%F%order(1),hm%F%order(2)
                  !if (abs(DE_bc-in*hm%F%omega)<1e-14) cycle
                  inn = in - hm%F%order(1) + 1
                  !\sum_l <u_(l-n)c| r | u_lb>
                  dn_cb = M_z0 
                  do il=hm%F%order(1), hm%F%order(2)
                    ill = il - hm%F%order(1) + 1

                    if (-in+il .ge. hm%F%order(1) .and. -in+il .le. hm%F%order(2)) then
                    select case (obs%gauge)
                      case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
                        !<u_a| J | u_b>
                        call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                                u_c(:,(-in+il-hm%F%order(1))*spindim+1: (-in+il-hm%F%order(1))*spindim +spindim),&
                                u_b(:,(ill-1)*spindim+1: (ill-1)*spindim +spindim) , ik, tmp2(:,:))
                        do dir=1,dim
                          dn_cb(dir) = dn_cb(dir) + sum(tmp2(dir,1:spindim)) !/(DE_bc-in*hm%F%omega)  
                        end do
                      case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)
                        do idim=1,spindim
                          call zmf_multipoles(sys%gr%mesh, conjg(u_c(:,(-in+il-hm%F%order(1))*spindim+idim)) &
                                                               * u_b(:,(ill-1)*spindim+idim) , 1, tmp(:))
  
                          dn_cb(:) = dn_cb(:) + tmp(2:4)
                        end do ! idim loop
                      end select
                    end if
                        

                  end do ! il loop
              
                  do ie = 1, obs%ne
                    EE= ie * obs%de
                    do idir = 1, dim 
                      do jdir = idir, dim 
                        ampl =  M_zI * &
                                 conjg(dressed_st%coeff(ista,ik)) *dressed_st%coeff(istb,ik) * &
                                 dm_ac(idir)*dn_cb(jdir)* &
                                (1/(DE_bc - in * hm%F%omega  + EE + M_zI*obs%gamma)- &
                                 1/(DE_ca - im * hm%F%omega  + EE + M_zI*obs%gamma) )!- & 
                        if (obs%gauge==OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY) then
                          ampl = ampl/EE
                        else 
                          ampl = ampl * EE           
                        end if 
                        if (idir == jdir) ampl = M_z2 * ampl
                        
                        ! sum over a and b         
                        sigma(ie,idir,jdir) = sigma(ie,idir,jdir) + ampl* (dressed_st%d%kweights(ik))
                          
                        end do    !jdir      
                      end do    !idir
      
                  end do ! ie loop

                end do ! in loop
              end do ! im loop

            end do ! istc loop
           end do ! istb loop   
        end do ! ista loop   
 
      end do ! ik loop
                 
          
      if(dressed_st%parallel_in_states .or. dressed_st%d%kpt%parallel) then
        call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm,  sigma(:,:,:))
      end if
      SAFE_DEALLOCATE_A(u_a)
      SAFE_DEALLOCATE_A(u_b)
      SAFE_DEALLOCATE_A(u_c)

    case (OPTION__FLOQUETCONDUCTIVITYMETHOD__F_FBZ21)
      
      SAFE_ALLOCATE(u_a(1:sys%gr%mesh%np, hm%F%floquet_dim))
      SAFE_ALLOCATE(u_b(1:sys%gr%mesh%np, hm%F%floquet_dim))
     
      
      sigma(:,:,:) = M_z0
      print *, 'FBZ21 nst', dressed_st%nst
      do ik=dressed_st%d%kpt%start, dressed_st%d%kpt%end

        do ista=dressed_st%st_start, dressed_st%st_end
          call states_get_state(dressed_st, sys%gr%mesh, ista, ik, u_a)
          
          do istb=1, dressed_st%nst
            ! if (ista == istb) cycle
              
            DE_ab = dressed_st%eigenval(ista,ik) - dressed_st%eigenval(istb,ik)
            ! print *, ista, istb, 'energy diff', DE_ab
            call states_get_state(dressed_st, sys%gr%mesh, istb, ik, u_b)
                             
            ! get the dipole matrix elements dn_cb and dm_ac
            do im=hm%F%order(1),hm%F%order(2)
              ! if (abs(-DE_ab-im * hm%F%omega)<1E-14) cycle
              imm = im - hm%F%order(1) + 1
              dm_ab(:) = M_z0

              do il=hm%F%order(1), hm%F%order(2)
                ill = il - hm%F%order(1) + 1
                !\sum_n <u_(l-m)a| r | u_lb> AND \sum_n <u_(l-m)b| r | u_la>
                if (-im+il .ge. hm%F%order(1) .and. -im+il .le. hm%F%order(2)) then
                  
                  select case (obs%gauge)
                  case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
                    !<u_a| J | u_b>
                    call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                            u_a(:,(-im+il-hm%F%order(1))*spindim+1: (-im+il-hm%F%order(1))*spindim +spindim) ,&
                            u_b(:,(ill-1)*spindim+1: (ill-1)*spindim +spindim) , ik, tmp2(:,:))
                    !if (abs(-DE_ab-im*hm%F%omega)>1E-5) then
                    do dir=1,dim
                      dm_ab(dir) = dm_ab(dir) + sum(tmp2(dir,1:spindim)) !/(-DE_ab-im*hm%F%omega) 
                    end do
                    !end if 
                  case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)
                    do idim=1,spindim
                      call zmf_multipoles(sys%gr%mesh, conjg(u_a(:,(-im+il-hm%F%order(1))*spindim+idim)) &
                                                           * u_b(:,(ill-1)*spindim+idim) , 1, tmp(:))  
                      dm_ab(:) = dm_ab(:) + tmp(2:4)
                    end do
                  end select
                end if
              end do ! il

              do in=hm%F%order(1),hm%F%order(2)
                ! if (abs(DE_ab-in * hm%F%omega)<1E-14) cycle
                inn = in - hm%F%order(1) + 1
                dn_ba(:) = M_z0

                do il = hm%F%order(1), hm%F%order(2)
                  ill = il - hm%F%order(1) + 1
                  if (-in+il .ge. hm%F%order(1) .and. -in+il .le. hm%F%order(2)) then
                    select case (obs%gauge)
                    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY)
                      !<u_a| J | u_b>
                      call current_calculate_mel(sys%gr%der, hm, sys%geo, &
                              u_b(:,(-in+il-hm%F%order(1))*spindim+1: (-in+il-hm%F%order(1))*spindim +spindim) ,&
                              u_a(:,(ill-1)*spindim+1: (ill-1)*spindim +spindim) , ik, tmp2(:,:))
                    !if (abs(DE_ab-in*hm%F%omega)>1E-5) then
                      do dir=1,dim
                        dn_ba(dir) = dn_ba(dir) + sum(tmp2(dir,1:spindim))!/(DE_ab-in*hm%F%omega) 
                      end do
                    !end if
                    case (OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH)
                      do idim=1,spindim
                        call zmf_multipoles(sys%gr%mesh, conjg(u_b(:,(-in+il-hm%F%order(1))*spindim+idim)) &
                                                             * u_a(:,(il-hm%F%order(1))*spindim+idim) , 1, tmp(:))
                        dn_ba(:) = dn_ba(:) + tmp(2:4)
                      end do
                    end select
                  end if
                    
                end do ! il loop
           !print *, 'matrix element', ista, istb, im, in, abs(dn_ba(1)*dm_ab(1))
           
                do ie = 1, obs%ne
                  EE= ie * obs%de
                  do idir = 1, dim
                    do jdir = idir, dim   
                      ampl = M_zI *&
                                 (dressed_st%occ(ista,ik)- dressed_st%occ(istb,ik)) *&
                                 dn_ba(idir)*dm_ab(jdir)* &
                                 (1/(DE_ab - in * hm%F%omega  + EE + M_zI*obs%gamma)) ! + &
                                  !1/(DE_ab - in * hm%F%omega  - EE - M_zI*obs%gamma) )
                      if (obs%gauge==OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY) then
                         ampl = ampl / EE        
                      else
                         ampl = ampl * EE
                      end if
                         !ampl = ampl * exp(M_zI*(2*im*hm%F%omega - EE)*obs%time0)
                      if (idir == jdir) ampl = M_z2 * ampl
                   
                      ! sum over a, b, im and ik        
                      sigma(ie,idir,jdir) = sigma(ie,idir,jdir) + ampl * (dressed_st%d%kweights(ik))
                    end do    !jdir      
                  end do    !idir
                end do ! ie loop
              end do ! in loop
            end do ! im loop
          end do ! istb loop   
        end do ! ista loop   
         
      end do ! ik loop
      
      if(dressed_st%d%kpt%parallel) then
        call comm_allreduce(dressed_st%d%kpt%mpi_grp%comm,  sigma(:,:,:))
      end if

      SAFE_DEALLOCATE_A(u_a)
      SAFE_DEALLOCATE_A(u_b)
    end select
    
    print *, 'gauge variable name', obs%gauge, obs%fc_method
    if(mpi_grp_is_root(mpi_world)) then
      write(iter_name,'(i4)') hm%F%iter

      if (obs%fc_method==OPTION__FLOQUETCONDUCTIVITYMETHOD__F_OKA) then
        methodname='oka'
      else if (obs%fc_method==OPTION__FLOQUETCONDUCTIVITYMETHOD__F_FBZ) then
        methodname='fbz'
      else if (obs%fc_method==OPTION__FLOQUETCONDUCTIVITYMETHOD__F_FBZ21) then
        methodname='fbz21'
      end if

      if (obs%gauge==OPTION__FLOQUETOBSERVABLEGAUGE__F_VELOCITY) then
        gauge='v'
      else if (obs%gauge==OPTION__FLOQUETOBSERVABLEGAUGE__F_LENGTH) then
        gauge='l'
      end if

      filename = FLOQUET_DIR//'/conductivity_'//trim(gauge)//'_'//trim(methodname)
      ! filename = FLOQUET_DIR//'/floquet_conductivity_'//trim(gauge)//'_'//trim(adjustl(iter_name))
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
    

    POP_SUB(calc_floquet_conductivity)    
  end subroutine calc_floquet_conductivity
  


  subroutine calc_floquet_forces()

    type(geometry_t) :: geo
    type(grid_t) :: gr
    type(unit_t) :: fn_unit
    
    integer :: idim, iatom, ii, idir, iunit, ierr, ib
    character(len=1024):: filename


    gr=sys%gr
    geo=sys%geo

    ! re-calculate occupations in case 
    if(hm%F%occ_cut > 0) then
       call states_allocate_wfns(bare_st,gr%der%mesh, wfs_type = TYPE_CMPLX)
       call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh, exact=.true.)
       if(ierr == 0) call states_load(restart, bare_st, gr, ierr)
       if (ierr /= 0) then
          message(1) = 'Unable to read ground-state wavefunctions.'
          call messages_fatal(1)
       end if
       call restart_end(restart)
       call floquet_calc_occupations(hm, sys, dressed_st, bare_st)
    end if

    filename = FLOQUET_DIR//'/floquet_forces'
    iunit = io_open(filename, action='write')

    call forces_calculate(gr, geo, hm, dressed_st, sys%ks)

    if(mpi_grp_is_root(mpi_world)) then
       write(iunit,'(3a)') 'Floquet-forces on the ions [', trim(units_abbrev(units_out%force)), "]"
       write(iunit,'(a,10x,99(14x,a))') ' Ion', (index2axis(idir), idir = 1, gr%sb%dim)
       do iatom = 1, geo%natoms
          write(iunit,'(i4,a10,10f15.6)') iatom, trim(species_label(geo%atom(iatom)%species)), &
               (units_from_atomic(units_out%force, geo%atom(iatom)%f(idir)), idir=1, gr%sb%dim)
       end do
       write(iunit,'(1x,100a1)') ("-", ii = 1, 13 + gr%sb%dim * 15)
       write(iunit,'(a14, 10f15.6)') " Max abs force", &
            (units_from_atomic(units_out%force, maxval(abs(geo%atom(1:geo%natoms)%f(idir)))), idir=1, gr%sb%dim)
       write(iunit,'(a14, 10f15.6)') " Total force", &
            (units_from_atomic(units_out%force, sum(geo%atom(1:geo%natoms)%f(idir))), idir=1, gr%sb%dim)
    end if

  end subroutine calc_floquet_forces


  subroutine plot_floquet_density()
    type(grid_t) :: gr

    integer :: idim, iatom, ii, idir, iunit, ierr, ib
    character(len=1024):: filename

    FLOAT, allocatable :: rho(:,:), rhoF(:,:)

    gr=sys%gr

    SAFE_ALLOCATE(rho(1:gr%mesh%np,1))
    SAFE_ALLOCATE(rhoF(1:gr%mesh%np,1))

    ! get groundstate density
    call states_allocate_wfns(bare_st,gr%der%mesh, wfs_type = TYPE_CMPLX)
    call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=gr%mesh, exact=.true.)
    if(ierr == 0) call states_load(restart, bare_st, gr, ierr)
    if (ierr /= 0) then
       message(1) = 'Unable to read ground-state wavefunctions.'
       call messages_fatal(1)
    end if
    call restart_end(restart)

    rho = M_ZERO
    call density_calc(bare_st, gr,rho)

    ! get Floquet density
    rhoF = M_ZERO
    call density_calc(dressed_st, gr,rhoF)
    call dio_function_output(io_function_fill_how("PlaneX"), FLOQUET_DIR, "rho_floquet",  gr%mesh, rhoF(:,1), unit_one, ierr)
    call dio_function_output(io_function_fill_how("PlaneY"), FLOQUET_DIR, "rho_floquet",  gr%mesh, rhoF(:,1), unit_one,ierr)
    call dio_function_output(io_function_fill_how("PlaneZ"), FLOQUET_DIR, "rho_floquet",  gr%mesh, rhoF(:,1), unit_one,ierr)
    ! difference
    rho = rho -rhoF
    call dio_function_output(io_function_fill_how("PlaneX"), FLOQUET_DIR, "rho_diff",  gr%mesh, rho(:,1), unit_one, ierr)
    call dio_function_output(io_function_fill_how("PlaneY"), FLOQUET_DIR, "rho_diff",  gr%mesh, rho(:,1), unit_one, ierr)
    call dio_function_output(io_function_fill_how("PlaneZ"), FLOQUET_DIR, "rho_diff",  gr%mesh, rho(:,1), unit_one, ierr)


    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(rhoF)
  end subroutine plot_floquet_density


  end program floquet_observables

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
