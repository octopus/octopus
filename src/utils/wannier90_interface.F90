!! Copyright (C) 2017-2019 H. Huebener, N. Tancogne-Dejean
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

program wannier90_interface
  use batch_oct_m
  use calc_mode_par_oct_m
  use comm_oct_m
  use command_line_oct_m
  use cube_oct_m
  use cube_function_oct_m
  use fft_oct_m
  use global_oct_m
  use grid_oct_m
  use io_binary_oct_m
  use io_function_oct_m
  use io_oct_m
  use ions_oct_m
  use iso_fortran_env
  use kpoints_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use lattice_vectors_oct_m
  use loct_oct_m
  use loct_math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m  
  use multicomm_oct_m
  use namespace_oct_m
  use orbitalset_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use states_elec_calc_oct_m
  use electrons_oct_m
  use space_oct_m
  use string_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use states_elec_io_oct_m
  use states_elec_restart_oct_m
  use submesh_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use wfs_elec_oct_m
  use ylm_wannier_oct_m

  implicit none

  integer              :: w90_what, w90_mode, w90_what_default

  integer              :: ierr
  integer              :: dim, idim
  integer              :: ii, nik, iter, nst

  type(restart_t)      :: restart
  type(electrons_t), pointer :: sys
  logical              :: w90_spinors, scdm_proj, w90_scdm
  integer              :: w90_nntot, w90_num_bands, w90_num_kpts   ! w90 input parameters
  integer, allocatable :: w90_nnk_list(:,:)                        !
  character(len=80)    :: w90_prefix                               ! w90 input file prefix
  integer              :: w90_num_wann                             ! input paramter
  FLOAT, allocatable   :: w90_proj_centers(:,:)                    ! projections centers
  integer, allocatable :: w90_proj_lmr(:,:)                        ! definitions for real valued Y_lm*R_r
  integer :: w90_nproj                                             ! number of such projections        
  integer, allocatable :: w90_spin_proj_component(:)               ! up/down flag 
  FLOAT, allocatable   :: w90_spin_proj_axis(:,:)                  ! spin axis (not implemented)
  integer              :: w90_num_exclude
  logical, allocatable :: exclude_list(:)                          ! list of excluded bands
  integer, allocatable :: band_index(:)                            ! band index after exclusion
  logical              :: read_td_states
  integer              :: w90_spin_channel                         !> For spin-polarized cases, the selected spin channel

  ! scdm variables
  integer, allocatable :: jpvt(:)
  CMPLX, allocatable   :: uk(:,:,:)                                ! SCDM-Wannier gauge matrices U(k)
  CMPLX, allocatable   :: psi(:,:)
  CMPLX, allocatable   :: chi(:,:), chi_diag(:,:),chi2(:,:)
  FLOAT, allocatable   :: chi_eigenval(:), occ_temp(:)
  FLOAT                :: scdm_mu, scdm_sigma, smear,  kvec(MAX_DIM)
  integer :: ist, jst, ik

  call global_init()
  call parser_init()

  call messages_init()
  call io_init()

  call calc_mode_par_init()

  call profiling_init(global_namespace)

  call restart_module_init(global_namespace)

  call fft_all_init(global_namespace)
  call unit_system_init(global_namespace)

  call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
  sys => electrons_t(global_namespace)
  call sys%init_parallelization(mpi_world)

  !%Variable Wannier90Prefix
  !%Type string
  !%Default w90
  !%Section Utilities::oct-wannier90
  !%Description
  !% Prefix for wannier90 files
  !%End
  call parse_variable(global_namespace, 'Wannier90Prefix', 'w90', w90_prefix)
  if(w90_prefix=='w90') then
    message(1) = "oct-wannier90: the prefix is set by default to w90"
    call  messages_info(1)
  end if

  !%Variable Wannier90Mode
  !%Type integer
  !%Default 0
  !%Section Utilities::oct-wannier90
  !%Description
  !% Specifies which stage of the Wannier90 interface to use
  !%Option none 0
  !% Nothing is done.
  !%Option w90_setup 1
  !% Writes parts of the wannier90 input file <tt>w90_prefix.win</tt> corresponding to
  !% the octopus inp file. Importantly it generates the correct form of Monkhorst-Pack mesh
  !% written to the file w90_kpoints that has to be used in a gs calculation of Octopus by
  !% as <tt> include w90_kpoints </tt> instead of the <tt>%KpointsGrid</tt> block.
  !%Option w90_output 2
  !% Generates the relevant files for a wannier90 run, specified by the variable <tt>W90_interface_files</tt>.
  !% This needs files previously generated
  !% by <tt>wannier90.x -pp w90 </tt>
  !%Option w90_wannier 3
  !% Parse the output of wannier90 to generate the Wannier states on the real-space grid. 
  !% The states will be written in the folder wannier. By default, the states are written as
  !% binary files, similar to the Kohn-Sham states.
  !%
  !% Not implemented for spinor states.
  !%End
  call parse_variable(global_namespace, 'Wannier90Mode', 0, w90_mode)

  if(w90_mode == 0) then
    message(1) = "Wannier90Mode must be set to a value different from 0."
    call messages_fatal(1)
  end if

  !%Variable Wannier90Files
  !%Type flag
  !%Default w90_mmn + w90_amn + w90_eig
  !%Section Utilities::oct-wannier90
  !%Description
  !% Specifies which files to generate.
  !% Example: <tt>w90_mmn + w90_unk</tt>
  !%Option w90_mmn bit(1)
  !% (see Wannier90 documentation)
  !%Option w90_unk bit(2)
  !% (see Wannier90 documentation)
  !%Option w90_amn bit(3)
  !% (see Wannier90 documentation)
  !%Option w90_eig bit(4)
  !% Eigenvalues. See Wannier90 documentation for more details.
  !%End
  w90_what_default = OPTION__WANNIER90FILES__W90_MMN + OPTION__WANNIER90FILES__W90_AMN + OPTION__WANNIER90FILES__W90_EIG
  call parse_variable(global_namespace, 'Wannier90Files', w90_what_default, w90_what)

  !%Variable Wannier90UseTD
  !%Type logical
  !%Default no
  !%Section Utilities::oct-wannier90
  !%Description
  !% By default oct-wannier90 uses the ground-state states to compute the necessary information.
  !% By setting this variable to yes, oct-wannier90 will use the TD states instead. 
  !%End
  call parse_variable(global_namespace, 'Wannier90UseTD', .false., read_td_states)

  !%Variable Wannier90UseSCDM
  !%Type logical
  !%Default no
  !%Section Utilities::oct-wannier90
  !%Description
  !% By default oct-wannier90 uses the projection method to generate the .amn file.
  !% By setting this variable to yes, oct-wannier90 will use SCDM method instead. 
  !%End
  call parse_variable(global_namespace, 'Wannier90UseSCDM', .false., w90_scdm)
  if(w90_scdm) then
    !%Variable SCDMsigma
    !%Type float
    !%Default 0.2
    !%Section Utilities::oct-wannier90  
    !%Description
    !% Broadening of SCDM smearing function.
    !%End
    call parse_variable(global_namespace, 'SCDMsigma', CNST(0.2), scdm_sigma)
  
    !%Variable SCDMmu
    !%Type float
    !%Section Utilities::oct-wannier90
    !%Description
    !% Energy range up to which states are considered for SCDM. 
    !%End
    call parse_variable(global_namespace, 'SCDMmu', M_HUGE, scdm_mu)
  end if

  if(sys%kpoints%use_symmetries) then
    message(1) = 'oct-wannier90: k-points symmetries are not allowed'
    call messages_fatal(1)
  end if
  if(sys%kpoints%use_time_reversal) then
    message(1) = 'oct-wannier90: time-reversal symmetry is not allowed'
    call messages_fatal(1)
  end if
  if(sys%kpoints%reduced%nshifts > 1) then
    message(1) = 'oct-wannier90: Wannier90 does not allow for multiple shifts of the k-point grid'
    call messages_fatal(1)
  end if


  if(sys%st%d%ispin /= UNPOLARIZED) then
    call messages_experimental("oct-wannier90 with SpinComponnents /= unpolarized") 
  end if

  w90_spinors = .false.
  w90_spin_channel = 1

  ! create setup files
  select case(w90_mode) 
  case(OPTION__WANNIER90MODE__W90_SETUP)
    call wannier90_setup(sys%ions, sys%kpoints)

  ! load states and calculate interface files
  case(OPTION__WANNIER90MODE__W90_OUTPUT)
    call wannier90_output()

  case(OPTION__WANNIER90MODE__W90_WANNIER)
    call read_wannier90_files()

   ! normal interface run
    call states_elec_allocate_wfns(sys%st, sys%gr%mesh, wfs_type = TYPE_CMPLX, skip=exclude_list)
    if(read_td_states) then
      call restart_init(restart, global_namespace, RESTART_TD, RESTART_TYPE_LOAD, &
                       sys%mc, ierr, sys%gr%mesh)
    else
      call restart_init(restart, global_namespace, RESTART_GS, RESTART_TYPE_LOAD, &
                       sys%mc, ierr, sys%gr%mesh)
    end if

    if(ierr == 0) then
      call states_elec_look(restart, nik, dim, nst, ierr)
      if(dim == sys%st%d%dim .and. nik == sys%kpoints%reduced%npoints .and. nst == sys%st%nst) then
        call states_elec_load(restart, global_namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, &
                                 ierr, iter, label = ": wannier90", skip=exclude_list)
      else
         write(message(1),'(a)') 'Restart structure not commensurate.'
         call messages_fatal(1)
      end if
    end if
    call restart_end(restart)

    call generate_wannier_states(sys%space, sys%gr%mesh, sys%ions, sys%st, sys%kpoints)
  case default
    message(1) = "Wannier90Mode is set to an unsupported value."
    call messages_fatal(1)
  end select

  SAFE_DEALLOCATE_A(exclude_list)
  SAFE_DEALLOCATE_A(band_index)
  SAFE_DEALLOCATE_A(w90_nnk_list)

  SAFE_DEALLOCATE_P(sys)
  call fft_all_end()
  call io_end()
  call profiling_end(global_namespace)
  call messages_end()
  call parser_end()
  call global_end()

contains

  subroutine wannier90_setup(ions, kpoints)
    type(ions_t),      intent(in) :: ions
    type(kpoints_t),   intent(in) :: kpoints

    character(len=80) :: filename
    integer :: w90_win, ia, axis(3), npath

    PUSH_SUB(wannier90_setup)

    ! open win file
    filename = trim(adjustl(w90_prefix)) //'.win'
    w90_win = io_open(trim(filename), global_namespace, action='write')

    write(w90_win,'(a)') '# this file has been created by the Octopus wannier90 utility'
    write(w90_win,'(a)') ' '

    ! write direct lattice vectors (in angstrom)
    write(w90_win,'(a)') 'begin unit_cell_cart'
    write(w90_win,'(a)') 'Ang'
    do idim = 1,3
      write(w90_win,'(f13.8,f13.8,f13.8)') units_from_atomic(unit_angstrom, ions%latt%rlattice(1:3,idim))
    end do
    write(w90_win,'(a)') 'end unit_cell_cart'
    write(w90_win,'(a)') ' '

    write(w90_win,'(a)') 'begin atoms_frac'
    do ia = 1, ions%natoms
       write(w90_win,'(a,2x,f13.8,f13.8,f13.8)') trim(ions%atom(ia)%label), ions%latt%cart_to_red(ions%pos(:, ia))
    end do
    write(w90_win,'(a)') 'end atoms_frac'
    write(w90_win,'(a)') ' '

    ! This is a default value. In practice, one should use projections
    write(w90_win,'(a)') 'use_bloch_phases = .true.'
    write(w90_win,'(a)') ' '

    write(w90_win,'(a10,i4)') 'num_bands ', sys%st%nst
    write(w90_win,'(a9,i4)') 'num_wann ', sys%st%nst
    write(w90_win,'(a)') ' '

    if(sys%st%d%ispin == SPINORS) then
       write(w90_win,'(a)') 'spinors = .true.'
    end if
    if(sys%st%d%ispin == SPIN_POLARIZED) then
      write(w90_win, '(a)') 'spin = up'
    end if

    ! This is for convenience. This is needed for plotting the Wannier states, if requested.
    write(w90_win,'(a)')  'write_u_matrices = .true.'
    write(w90_win,'(a)')  'translate_home_cell = .true.'
    write(w90_win,'(a)')  'write_xyz = .true.'
    write(w90_win,'(a)') ' '

    if(kpoints%reduced%npoints == 1) then
      write(w90_win,'(a)')  'gamma_only = .true.'
      write(w90_win,'(a)') ' '
    else
      if(.not.parse_is_defined(global_namespace, 'KPointsGrid')) then
        message(1) = 'oct-wannier90: need Monkhorst-Pack grid. Please specify %KPointsGrid'
        call messages_fatal(1)
      end if

      ! In case the user used also a k-point path, we ignore it
      npath = kpoints_nkpt_in_path(kpoints)

      axis(1:3) = kpoints%nik_axis(1:3)
      ASSERT(product(kpoints%nik_axis(1:3)) == kpoints%reduced%npoints - npath)

      write(w90_win,'(a8,i4,i4,i4)')  'mp_grid =', axis(1:3)
      write(w90_win,'(a)') ' '
      write(w90_win,'(a)')  'begin kpoints '
      ! Put a minus sign here for the wrong convention in Octopus

      do ii = 1, kpoints%reduced%npoints-npath
        write(w90_win,'(f13.8,f13.8,f13.8)') - kpoints%reduced%red_point(1:3,ii) 
      end do
      write(w90_win,'(a)')  'end kpoints '
    end if

    call io_close(w90_win)

    POP_SUB(wannier90_setup)

  end subroutine wannier90_setup

  subroutine wannier90_output()
    integer :: ik_real

    PUSH_SUB(wannier90_output)

    call read_wannier90_files()

    ! normal interface run
    call states_elec_allocate_wfns(sys%st, sys%gr%mesh, wfs_type = TYPE_CMPLX, skip=exclude_list)
    if(read_td_states) then
      call restart_init(restart, global_namespace, RESTART_TD, RESTART_TYPE_LOAD, &
                       sys%mc, ierr, sys%gr%mesh)
    else
      call restart_init(restart, global_namespace, RESTART_GS, RESTART_TYPE_LOAD, &
                       sys%mc, ierr, sys%gr%mesh)
    end if

    if(ierr == 0) then
      call states_elec_look(restart, nik, dim, nst, ierr)
      if(sys%st%d%ispin == SPIN_POLARIZED) then
        nik = nik / 2
      end if
      if(dim == sys%st%d%dim .and. nik == sys%kpoints%reduced%npoints .and. nst == sys%st%nst) then
         call states_elec_load(restart, global_namespace, sys%space, sys%st, sys%gr%mesh, sys%kpoints, &
                    ierr, iter, label = ": wannier90", skip=exclude_list)
      else
         write(message(1),'(a)') 'Restart structure not commensurate.'
         call messages_fatal(1)
      end if
    end if
    call restart_end(restart)

    if(w90_scdm) then
      nik = w90_num_kpts
      SAFE_ALLOCATE(jpvt(1:sys%gr%mesh%np_global))
      SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))
      SAFE_ALLOCATE(occ_temp(1:w90_num_bands))

      ! smear the states at gamma
      do ist = 1, w90_num_bands
         occ_temp(ist)= sys%st%occ(ist, 1) 
         sys%st%occ(ist, 1)=M_HALF*loct_erfc((sys%st%eigenval(ist, 1)-scdm_mu) / scdm_sigma)
      end do

      call zstates_elec_rrqr_decomposition(sys%st, sys%namespace, sys%gr%mesh, w90_num_bands, .true., 1, jpvt)

      ! reset occupations at gamma
      do ist = 1, w90_num_bands
         sys%st%occ(ist, 1) = occ_temp(ist)
      end do

      SAFE_ALLOCATE(uk(1:w90_num_bands, 1:w90_num_bands, 1:nik))

      ! auxiliary arrays for scdm procedure
      SAFE_ALLOCATE(chi(1:w90_num_bands, 1:w90_num_bands))
      SAFE_ALLOCATE(chi_diag(1:w90_num_bands, 1:w90_num_bands))
      SAFE_ALLOCATE(chi2(1:w90_num_bands, 1:w90_num_bands))
      SAFE_ALLOCATE(chi_eigenval(1:w90_num_bands))

      chi(1:w90_num_bands, 1:w90_num_bands) = M_ZERO

      do ik = 1, nik
        kvec(:) = sys%kpoints%reduced%point(:, ik)

        if(sys%st%d%ispin == SPIN_POLARIZED) then
          ik_real = (ik-1)*2 + w90_spin_channel
        else
          ik_real = ik
        end if
  
        do ist = 1, w90_num_bands
          call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik_real, psi)
          smear=M_HALF * loct_erfc((sys%st%eigenval(ist, ik_real) - scdm_mu) / scdm_sigma)
          ! NOTE: here check for domain parallelization
          do jst = 1, w90_num_bands
             chi(ist, jst) = smear * conjg(psi(jpvt(jst), 1)) &
                 * exp(M_zI * dot_product(sys%gr%mesh%x(jpvt(jst), 1:3), kvec(1:3)))
          end do
        end do

        ! loewdin orhtogonalization of chi.chi
        ! this can also be done with SVD, which might be more stable!?
        chi_diag = matmul(conjg(transpose(chi)), chi)
        call lalg_eigensolve(w90_num_bands, chi_diag, chi_eigenval)
        chi2 = conjg(transpose(chi_diag))
        
        !we need the eigenvalues to be >0
        if( any(chi_eigenval(:) .lt. M_ZERO) ) then
           message(1) = 'SCDM Wannierization failed because chi matrix is'
           message(2) = 'ill conditioned. Try increasingin scdm_sigma and/or'
           message(3) = 'change scdm_mu.'
           call messages_fatal(3)
        end if

        do ist = 1, w90_num_bands
          chi_eigenval(ist) = M_ONE / sqrt(chi_eigenval(ist))
          chi2(ist, 1:w90_num_bands) = chi_eigenval(ist) * chi2(ist, 1:w90_num_bands)
        end do
        ! the loewdin result would be: matmul(chi_diag,chi2)
        ! to get the wannier gauge U(k) we multiply this with the original chi
        uk(:,:,ik) = matmul(chi, matmul(chi_diag,chi2))

      end do
      
      SAFE_DEALLOCATE_A(chi)
      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(chi_diag)
      SAFE_DEALLOCATE_A(chi2)
      SAFE_DEALLOCATE_A(chi_eigenval)
      SAFE_DEALLOCATE_A(jpvt)
      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(occ_temp)
    end if

    ! ---- actual interface work ----------
    if(bitand(w90_what, OPTION__WANNIER90FILES__W90_MMN) /= 0) then
      call create_wannier90_mmn(sys%gr%mesh, sys%st)
    end if

    if(bitand(w90_what, OPTION__WANNIER90FILES__W90_UNK) /= 0) then
      call write_unk(sys%space, sys%gr%mesh, sys%st)
    end if

    if(bitand(w90_what, OPTION__WANNIER90FILES__W90_AMN) /= 0) then
      call create_wannier90_amn(sys%space, sys%gr%mesh, sys%ions%latt, sys%st, sys%kpoints)
    end if

    if(bitand(w90_what, OPTION__WANNIER90FILES__W90_EIG) /= 0) then
      call create_wannier90_eig()
    end if

    SAFE_DEALLOCATE_A(uk)
    SAFE_DEALLOCATE_A(w90_proj_centers)
    SAFE_DEALLOCATE_A(w90_spin_proj_component)
    SAFE_DEALLOCATE_A(w90_spin_proj_axis)

    POP_SUB(wannier90_output)
  end subroutine wannier90_output

  subroutine read_wannier90_files()
    integer ::  w90_nnkp, itemp, dummyint, io
    character(len=80) :: filename, dummy, dummy1, dummy2, line
    logical :: exist, parse_is_ok
    FLOAT :: dummyr(7)

    PUSH_SUB(read_wannier90_files)

    w90_num_kpts = product(sys%kpoints%nik_axis(1:3))
    w90_num_exclude = 0 

    ! open nnkp file
    filename = trim(adjustl(w90_prefix)) //'.nnkp'

    message(1) = "oct-wannier90: Parsing "//filename
    call messages_info(1)

    inquire(file=filename,exist=exist)
    if(.not. exist) then
      message(1) = 'oct-wannier90: Cannot find specified Wannier90 nnkp file.'
      write(message(2),'(a)') 'Please run wannier90.x -pp '// trim(adjustl(w90_prefix)) // ' first.' 
      call messages_fatal(2)
    end if

    parse_is_ok = .false.

    ! check number of k-points
    w90_nnkp = io_open(trim(filename), global_namespace, action='read')
    do
      read(w90_nnkp, *, iostat=io) dummy, dummy1
      if(io == iostat_end) exit

      if(dummy =='begin' .and. dummy1 == 'kpoints' ) then
        read(w90_nnkp,*) itemp
        if(itemp /= w90_num_kpts) then
          message(1) = 'oct-wannier90: wannier90 setup seems to have been done with a different number of k-points.'
          call messages_fatal(1)
        else
          parse_is_ok = .true.
          exit
        end if
      end if
    end do
    call io_close(w90_nnkp)

    if(.not. parse_is_ok) then
      message(1) = 'oct-wannier90: Did not find the kpoints block in nnkp file'
      call messages_fatal(1)
    end if
    parse_is_ok = .false.

    ! read from nnkp file
    ! find the nnkpts block
    w90_nnkp = io_open(trim(filename), global_namespace, action='read', position='rewind')
    do
      read(w90_nnkp, *, iostat=io) dummy, dummy1
      if(io  == iostat_end) exit !End of file

      if(dummy =='begin' .and. dummy1 == 'nnkpts' ) then
        read(w90_nnkp,*) w90_nntot
        SAFE_ALLOCATE(w90_nnk_list(1:5, 1:w90_num_kpts * w90_nntot))
        do ii = 1, w90_num_kpts * w90_nntot
          read(w90_nnkp,*) w90_nnk_list(1:5, ii)
        end do
        !make sure we are at the end of the block
        read(w90_nnkp,*) dummy
        if(dummy /= 'end') then
          message(1) = 'oct-wannier90: There dont seem to be enough k-points in nnkpts file to.'
          call messages_fatal(1)
        end if
        parse_is_ok = .true.
        exit
      end if
    end do
   
    if(.not. parse_is_ok) then
      message(1) = 'oct-wannier90: Did not find nnkpts block in nnkp file'
      call messages_fatal(1)
    end if

    ! read from nnkp file
    ! find the exclude block
    SAFE_ALLOCATE(exclude_list(1:sys%st%nst))
    !By default we use all the bands
    exclude_list(1:sys%st%nst) = .false.
    rewind(w90_nnkp)
    do
      read(w90_nnkp, *, iostat=io) dummy, dummy1
      if(io  == iostat_end) exit !End of file
      if(dummy =='begin' .and. dummy1 == 'exclude_bands') then
        read(w90_nnkp, *) w90_num_exclude
        do ii = 1, w90_num_exclude
          read(w90_nnkp, *) itemp
          exclude_list(itemp) = .true.
        end do
        !make sure we are at the end of the block
        read(w90_nnkp, *) dummy
        if(dummy /= 'end') then
          message(1) = 'oct-wannier90: There dont seem to be enough bands in exclude_bands list.'
          call messages_fatal(1)
        end if
        exit
      end if
    end do
    call io_close(w90_nnkp)

    !We get the number of bands
    w90_num_bands = sys%st%nst - w90_num_exclude

    SAFE_ALLOCATE(band_index(1:sys%st%nst))
    itemp = 0
    do ii = 1, sys%st%nst
      if(exclude_list(ii)) cycle
      itemp = itemp + 1
      band_index(ii) = itemp
    end do

    if(bitand(w90_what, OPTION__WANNIER90FILES__W90_AMN) /= 0) then
      ! parse file again for definitions of projections
      w90_nnkp = io_open(trim(filename), global_namespace, action='read', position='rewind')

      do
        read(w90_nnkp, *, iostat=io) dummy, dummy1
        if(io == iostat_end) then !End of file
          message(1) = 'oct-wannier90: Did not find projections block in w90.nnkp file'
          call messages_fatal(1)
        end if

        if(dummy =='begin' .and. (dummy1 == 'projections' .or. dummy1 == 'spinor_projections')) then

          if(dummy1 == 'spinor_projections') then
            w90_spinors = .true.
            if(sys%st%d%ispin /= SPINORS) then
              message(1) = 'oct-wannier90: Spinor = .true. is only valid with spinors wavefunctions.'
              call messages_fatal(1)
            end if

            message(1) = 'oct-wannier90: Spinor interface incomplete. Note there is no quantization axis implemented'
            call messages_warning(1)
          else
            if(sys%st%d%ispin == SPINORS) then
              message(1) = 'oct-wannier90: Octopus has spinors wavefunctions but spinor_projections is not defined.'
              message(2) = 'oct-wannier90: Please check the input file for wannier 90.'
              call messages_fatal(2)
            end if
          end if

          read(w90_nnkp, *) w90_nproj
          ! num_wann is given in w90.win, not double checked here
          ! I assume that the wannier90.x -pp run has checked this
          w90_num_wann = w90_nproj

          SAFE_ALLOCATE(w90_proj_centers(w90_nproj, 3))
          SAFE_ALLOCATE(w90_proj_lmr(w90_nproj, 3))
          if(w90_spinors) then
            SAFE_ALLOCATE(w90_spin_proj_component(w90_nproj)) 
          end if
          if(w90_spinors) then
            SAFE_ALLOCATE(w90_spin_proj_axis(w90_nproj, 3))
          end if

          do ii = 1, w90_nproj
             read(w90_nnkp, *) w90_proj_centers(ii, 1:3), w90_proj_lmr(ii, 1:3)
             ! skip a line for now
             read(w90_nnkp, *) dummyr(1:7)
             if(w90_spinors) then
                read(w90_nnkp, *) w90_spin_proj_component(ii), w90_spin_proj_axis(ii, 1:3)
                ! use octopus spindim conventions
                if(w90_spin_proj_component(ii) == -1) w90_spin_proj_component(ii) = 2
             end if
          end do
          !make sure we are at the end of the block
          read(w90_nnkp, *) dummy
          if(dummy /= 'end') then
             message(1) = 'oct-wannier90: There dont seem to be enough projections in nnkpts file to.'
             call messages_fatal(1)
          end if
          exit
        end if
      end do

      ! look for auto_projection block
      scdm_proj = .false.
      do
        read(w90_nnkp, *, iostat=io) dummy, dummy1
        if(io == iostat_end) exit !End of file

        if(dummy =='begin' .and. dummy1 =='auto_projections') then
          scdm_proj = .true.
          read(w90_nnkp, *) w90_nproj
          w90_num_wann = w90_nproj

          if(.not. w90_scdm) then
            message(1) = 'oct-wannier90: Found auto_projections block. Currently the only implemeted automatic way'
            message(2) = 'oct-wannier90: to compute projections is the SCDM method.'
            message(3) = 'oct-wannier90: Please set Wannier90UseSCDM = yes in the inp file.'
            call messages_fatal(3)
          end if

          if(w90_nproj /= w90_num_bands) then
            message(1) = 'oct-wannier90: In auto_projections block first row needs to be equal to num_bands.'
            call messages_fatal(1)
          end if
          read(w90_nnkp, *) dummyint
          if(dummyint /= 0) then
            message(1) = 'oct-wannier90: The second row in auto_projections has to be 0, per Wannier90 documentation.'
            call messages_fatal(1)
          end if
        end if
      end do
      call io_close(w90_nnkp)
 
    end if

    message(1) = "oct-wannier90: Finished parsing "//filename
    call messages_info(1)

    ! Look extra variables variable
    ! open win file
    filename = trim(adjustl(w90_prefix)) //'.win'
    message(1) = "oct-wannier90: Parsing "//filename
    call messages_info(1)
    w90_nnkp = io_open(trim(filename), global_namespace, action='read', position='rewind')
    do
      read(w90_nnkp, fmt='(a)', iostat=io) line
      if(io  == iostat_end) exit !End of file
      if(index(line, '=') > 0) then
        read(line, *, iostat=io) dummy, dummy2, dummy1
      else
        read(line, *, iostat=io) dummy, dummy1
      end if

      !Spin
      if(dummy =='spin') then
        if(sys%st%d%ispin /= SPIN_POLARIZED) then
          message(1) = 'oct-wannier90: The variable spin is set for a non spin-polarized calculation.'
          call messages_fatal(1)
        end if
  
        if(dummy1 == 'up') then
          w90_spin_channel = 1
        else if (dummy1 == 'down') then
          w90_spin_channel = 2
        else
          message(1) = 'oct-wannier90: Error parsing the variable spin.'
          call messages_fatal(1)
        end if
      end if
    end do
    call io_close(w90_nnkp)

    if(sys%st%d%ispin == SPIN_POLARIZED) then
      write(message(1), '(a,i1)') 'oct-wannier90: Using spin channel ', w90_spin_channel
      call messages_info(1)
    end if

    message(1) = "oct-wannier90: Finished parsing "//filename
    call messages_info(1)

    POP_SUB(read_wannier90_files)

  end subroutine read_wannier90_files

  subroutine create_wannier90_mmn(mesh, st)
    type(mesh_t),                intent(in) :: mesh
    type(states_elec_t), target, intent(in) :: st

    integer ::  ist, jst, ik, ip, w90_mmn, iknn, G(3), idim, ibind
    FLOAT   ::  Gr(3)
    character(len=80) :: filename
    CMPLX, allocatable :: overlap(:)
    CMPLX, allocatable :: psim(:,:), psin(:,:), phase(:)
    type(profile_t), save :: prof, reduce_prof
    type(wfs_elec_t), pointer :: batch

    PUSH_SUB(create_wannier90_mmn)

    call profiling_in(prof, "W90_MMN")

    if(st%d%kpt%parallel) then
      call messages_not_implemented("w90_mmn output with k-point parallelization")
    end if

    if(st%parallel_in_states) then
      call messages_not_implemented("w90_mmn output with states parallelization")
    end if

    message(1) = "Info: Computing the overlap matrix";
    call messages_info(1)


    filename = './'// trim(adjustl(w90_prefix))//'.mmn'
    w90_mmn = io_open(trim(filename), global_namespace, action='write')

    ! write header
    if(mpi_grp_is_root(mpi_world)) then
       write(w90_mmn,*) 'Created by oct-wannier90'
       write(w90_mmn,*) w90_num_bands, w90_num_kpts, w90_nntot
    end if

    SAFE_ALLOCATE(psim(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(psin(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(phase(1:mesh%np))
    SAFE_ALLOCATE(overlap(1:w90_num_bands))

    ! loop over the pairs specified in the nnkp file (read before in init)
    do ii = 1, w90_num_kpts * w90_nntot
       ik = w90_nnk_list(1, ii)
       iknn = w90_nnk_list(2, ii)
       G(1:3) = w90_nnk_list(3:5, ii)
       if(mpi_grp_is_root(mpi_world)) write(w90_mmn, '(I10,2x,I10,2x,I3,2x,I3,2x,I3)') ik, iknn, G(1:3)

       !For spin-polarized calculations, we select the right k-point
       if(sys%st%d%ispin == SPIN_POLARIZED) then
         ik = (ik-1)*2 + w90_spin_channel
         iknn = (iknn-1)*2 + w90_spin_channel
       end if

       Gr(1:3) = matmul(G(1:3), sys%ions%latt%klattice(1:3,1:3))

       if(any(G(1:3) /= 0)) then
         do ip = 1, mesh%np
           phase(ip) = exp(-M_zI*dot_product(mesh%x(ip,1:3), Gr(1:3)))
         end do
       end if

       ! loop over bands
       do jst = 1, st%nst
         if(exclude_list(jst)) cycle
         call states_elec_get_state(st, mesh, jst, iknn, psin)

         !Do not apply the phase if the phase factor is null
         if(any(G(1:3) /= 0)) then
           ! add phase
           do idim = 1, st%d%dim
             do ip = 1, mesh%np
               psin(ip, idim) = psin(ip, idim) * phase(ip)
             end do
           end do
         end if


         !See Eq. (25) in PRB 56, 12847 (1997)
         do ist = 1, st%nst
           if(exclude_list(ist)) cycle

           batch => st%group%psib(st%group%iblock(ist, ik), ik)
           select case(batch%status())
           case(BATCH_NOT_PACKED)
             overlap(band_index(ist)) = M_z0
             do idim = 1, st%d%dim
               ibind = batch%inv_index((/ist, idim/))
               overlap(band_index(ist)) = overlap(band_index(ist)) + &
                    zmf_dotp(mesh, batch%zff_linear(:, ibind), psin(:,idim), reduce = .false.)
             end do
           !Not properly done at the moment
           case(BATCH_PACKED, BATCH_DEVICE_PACKED)
             call states_elec_get_state(st, mesh, ist, ik, psim)
             overlap(band_index(ist)) = zmf_dotp(mesh, st%d%dim, psim, psin, reduce = .false.)
           end select
         end do

         if(mesh%parallel_in_domains) then
           call profiling_in(reduce_prof, "W90_MMN_REDUCE")
           call mesh%allreduce(overlap)
           call profiling_out(reduce_prof)
         end if
 
         ! write to W90 file
         if(mpi_grp_is_root(mpi_world)) then
           do ist = 1, st%nst
             if(exclude_list(ist)) cycle
             write(w90_mmn,'(e13.6,2x,e13.6)') overlap(band_index(ist))
           end do
         end if

       end do
    end do

    call io_close(w90_mmn)

    SAFE_DEALLOCATE_A(psim)
    SAFE_DEALLOCATE_A(psin)
    SAFE_DEALLOCATE_A(phase)
    SAFE_DEALLOCATE_A(overlap)

    call profiling_out(prof)

    POP_SUB(create_wannier90_mmn)

  end subroutine create_wannier90_mmn

  subroutine create_wannier90_eig()
    integer ::  ist, ik, w90_eig
    character(len=80) :: filename

    PUSH_SUB(create_wannier90_eig)

    if(sys%st%d%kpt%parallel) then
      call messages_not_implemented("w90_eig output with k-point parallelization")
    end if

    if(sys%st%parallel_in_states) then
      call messages_not_implemented("w90_eig output with states parallelization")
    end if

    if(mpi_grp_is_root(mpi_world)) then
      filename = './'//trim(adjustl(w90_prefix))//'.eig'
      w90_eig = io_open(trim(filename), global_namespace, action='write')
      do ik = 1, w90_num_kpts
        do ist = 1, sys%st%nst
          if(exclude_list(ist)) cycle
          if(sys%st%d%ispin /= SPIN_POLARIZED) then
            write(w90_eig,'(I5,2x,I5,2x,e13.6)') band_index(ist), ik,  &
                        units_from_atomic(unit_eV, sys%st%eigenval(ist, ik))
          else
            write(w90_eig,'(I5,2x,I5,2x,e13.6)') band_index(ist), ik,  &
                        units_from_atomic(unit_eV, sys%st%eigenval(ist, (ik-1)*2+w90_spin_channel))
          end if
        end do
      end do

      call io_close(w90_eig)
    end if

    POP_SUB(create_wannier90_eig)
  end subroutine create_wannier90_eig

  subroutine write_unk(space, mesh, st)
    type(space_t),       intent(in) :: space
    type(mesh_t),        intent(in) :: mesh
    type(states_elec_t), intent(in) :: st

    integer ::  ist, ik, unk_file, ispin
    integer ::  ix, iy, iz
    character(len=80) :: filename
    CMPLX, allocatable :: psi(:)
    type(cube_t) :: cube
    type(cube_function_t) :: cf

    PUSH_SUB(write_unk)

    if(st%d%kpt%parallel) then
      call messages_not_implemented("w90_unk output with k-point parallelization")
    end if
      
    if(sys%gr%mesh%parallel_in_domains) then
      call messages_not_implemented("w90_unk output with domain parallelization")
    end if

    if(st%parallel_in_states) then
      call messages_not_implemented("w90_unk output with states parallelization")
    end if

    call messages_experimental("Wannier90Files = w90_unk")


    SAFE_ALLOCATE(psi(1:mesh%np))

    call cube_init(cube, mesh%idx%ll, global_namespace, space, need_partition=.not.mesh%parallel_in_domains)
    call zcube_function_alloc_RS(cube, cf)

    do ik = 1, w90_num_kpts
      do ispin = 1, st%d%dim
        if(mpi_grp_is_root(mpi_world)) then
          write(filename, '(a,i5.5,a1,i1)') './UNK', ik,'.', ispin
          unk_file = io_open(trim(filename), global_namespace, action='write', form='unformatted')
          ! write header
          write(unk_file) mesh%idx%ll(1:mesh%idx%dim), ik, w90_num_bands
        end if

        ! states
        do ist = 1, st%nst
          if(exclude_list(ist)) cycle

          if(sys%st%d%ispin == SPIN_POLARIZED) then
            call states_elec_get_state(st, mesh, ispin, ist, ik, psi)
          else
            call states_elec_get_state(st, mesh, ispin, ist, (ik-1)*2+w90_spin_channel, psi)
          end if

          ! put the density in the cube
          ! Note: At the moment this does not work for domain parallelization
          if (cube%parallel_in_domains) then
            ASSERT(.not. cube%parallel_in_domains)
          else
            if(mesh%parallel_in_domains) then
              call zmesh_to_cube(mesh, psi, cube, cf, local = .true.)
            else
              call zmesh_to_cube(mesh, psi, cube, cf)
            end if
          end if

          if(mpi_grp_is_root(mpi_world)) then
            write(unk_file) (((cf%zrs(ix,iy,iz), ix=1,cube%rs_n_global(1)), iy=1,cube%rs_n_global(2)), iz=1,cube%rs_n_global(3))
          end if
        end do
        if(mpi_grp_is_root(mpi_world)) call io_close(unk_file)
      end do
    end do

    call dcube_function_free_RS(cube, cf)
    call cube_end(cube)

    SAFE_DEALLOCATE_A(psi)

    POP_SUB(write_unk)

  end subroutine write_unk

  subroutine create_wannier90_amn(space, mesh, latt, st, kpoints)
    type(space_t),           intent(in) :: space
    type(mesh_t),            intent(in) :: mesh
    type(lattice_vectors_t), intent(in) :: latt
    type(states_elec_t),     intent(in) :: st
    type(kpoints_t),         intent(in) :: kpoints

    integer ::  ist, ik, w90_amn, idim, iw, ip, ik_real
    FLOAT   ::  center(3),  kpoint(1:MAX_DIM), threshold
    character(len=80) :: filename
    CMPLX, allocatable :: psi(:,:), phase(:), projection(:)
    FLOAT, allocatable ::  rr(:,:), ylm(:)
    type(orbitalset_t), allocatable :: orbitals(:)
    type(profile_t), save :: prof, reduce_prof

    PUSH_SUB(create_wannier90_amn)
    call profiling_in(prof, "W90_AMN")

    if(st%d%kpt%parallel) then
      call messages_not_implemented("w90_amn output with k-point parallelization")
    end if

    if(st%parallel_in_states) then
      call messages_not_implemented("w90_amn output with states parallelization")
    end if

    filename = './'// trim(adjustl(w90_prefix))//'.amn'
    w90_amn = io_open(trim(filename), global_namespace, action='write')

    ! write header
    if(mpi_grp_is_root(mpi_world)) then
      write(w90_amn,*) 'Created by oct-wannier90'
      write(w90_amn,*)  w90_num_bands, w90_num_kpts, w90_num_wann
    end if

    if(scdm_proj) then

      message(1) = "Info: Writing projections obtained from SCDM."
      call messages_info(1)

      do ik = 1, w90_num_kpts
        do ist = 1, st%nst
          if(exclude_list(ist)) cycle
          if(mpi_grp_is_root(mpi_world)) then
            do iw = 1, w90_nproj
              write (w90_amn,'(I5,2x,I5,2x,I5,2x,e13.6,2x,e13.6)') band_index(ist), iw, ik, uk(band_index(ist),iw,ik) 
            end do
          end if
        end do !ist
      end do! ik

    else

      message(1) = "Info: Computing the projection matrix";
      call messages_info(1)
      
      !We use the variabel AOThreshold to deterine the threshold on the radii of the atomic orbitals
      call parse_variable(global_namespace, 'AOThreshold', CNST(0.01), threshold)
      
      SAFE_ALLOCATE(orbitals(1:w90_nproj))
      ! precompute orbitals
      do iw=1, w90_nproj
        call orbitalset_init(orbitals(iw))
        call orbitalset_init(orbitals(iw))
      
        orbitals(iw)%norbs = 1
        orbitals(iw)%ndim = 1
        orbitals(iw)%radius = -log(threshold)
        orbitals(iw)%submesh = .false.
      
        ! cartesian coordinate of orbital center
        center(1:3) = latt%red_to_cart(w90_proj_centers(iw,1:3))
        call submesh_init(orbitals(iw)%sphere, space, mesh, latt, center, orbitals(iw)%radius)
      
        ! make transpose table of submesh points for use in pwscf routine
        SAFE_ALLOCATE(rr(1:3,orbitals(iw)%sphere%np))
        do ip = 1,orbitals(iw)%sphere%np
          rr(1:3,ip) = orbitals(iw)%sphere%x(ip,1:3)
        end do
      
        ! get dorb as submesh points
        SAFE_ALLOCATE(orbitals(iw)%zorb(1:orbitals(iw)%sphere%np, 1:1, 1:1))
        SAFE_ALLOCATE(ylm(1:orbitals(iw)%sphere%np))
        ! (this is a routine from pwscf)
        call ylm_wannier(ylm, w90_proj_lmr(iw,1), w90_proj_lmr(iw,2), &
                              rr, orbitals(iw)%sphere%np)
        if(w90_proj_lmr(iw,3) == 1) then
          ! apply radial function
          do ip = 1,orbitals(iw)%sphere%np
            ylm(ip) = ylm(ip)*M_TWO*exp(-orbitals(iw)%sphere%r(ip))
          end do
        else
          call messages_not_implemented("oct-wannier90: r/=1 for the radial part")
        end if
      
        orbitals(iw)%zorb(1:orbitals(iw)%sphere%np, 1, 1) = ylm(1:orbitals(iw)%sphere%np) 
        SAFE_DEALLOCATE_A(ylm)
      
        SAFE_ALLOCATE(orbitals(iw)%phase(1:orbitals(iw)%sphere%np, 1:w90_num_kpts))
        orbitals(iw)%phase(:,:) = M_Z0
        SAFE_ALLOCATE(orbitals(iw)%eorb_mesh(1:mesh%np, 1:1, 1:1, 1:w90_num_kpts))
        orbitals(iw)%eorb_mesh(:,:,:,:) = M_Z0
      
        call orbitalset_update_phase(orbitals(iw), space%dim, st%d%kpt, kpoints, st%d%ispin == SPIN_POLARIZED, &
                                        kpt_max = w90_num_kpts)
      
        SAFE_DEALLOCATE_A(rr)
      end do
      
      SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))
      SAFE_ALLOCATE(phase(1:mesh%np))
      SAFE_ALLOCATE(projection(1:w90_nproj))
      
      do ik = 1, w90_num_kpts
        kpoint(1:space%dim) = kpoints%get_point(ik)
      
        do ip = 1, mesh%np
          phase(ip) = exp(-M_zI* sum(mesh%x(ip, 1:space%dim) * kpoint(1:space%dim)))
        end do

        !For spin-polarized calculations, we select the right k-point
        if(st%d%ispin == SPIN_POLARIZED) then
          ik_real = (ik-1)*2 + w90_spin_channel
        else
          ik_real = ik
        end if
      
        do ist = 1, st%nst
          if(exclude_list(ist)) cycle

          call states_elec_get_state(st, mesh, ist, ik_real, psi)

          do idim = 1, st%d%dim
            !The minus sign is here is for the wrong convention of Octopus
            do ip = 1, mesh%np
              psi(ip, idim) = psi(ip, idim)*phase(ip)
            end do
          end do
      
          do iw = 1, w90_nproj
            idim = 1
            if(w90_spinors) idim = w90_spin_proj_component(iw)
      
            !At the moemnt the orbitals do not depend on idim
            !The idim index for eorb_mesh would be for a spin-resolved orbital like j=1/2
            projection(iw) = zmf_dotp(mesh, psi(1:mesh%np,idim), &
                                    orbitals(iw)%eorb_mesh(1:mesh%np,1,1,ik_real), reduce = .false.)
          end do
      
          if(mesh%parallel_in_domains) then
            call profiling_in(reduce_prof, "W90_AMN_REDUCE")
            call mesh%allreduce(projection)
            call profiling_out(reduce_prof)
          end if
      
          if(mpi_grp_is_root(mpi_world)) then
            do iw = 1, w90_nproj
              write (w90_amn,'(I5,2x,I5,2x,I5,2x,e13.6,2x,e13.6)') band_index(ist), iw, ik, projection(iw)
            end do
          end if
        end do !ik
      end do !ist
      
      SAFE_DEALLOCATE_A(psi)
      SAFE_DEALLOCATE_A(phase)
      SAFE_DEALLOCATE_A(projection)
      
      do iw = 1, w90_nproj
        call orbitalset_end(orbitals(iw))
      end do
      SAFE_DEALLOCATE_A(orbitals)
    end if
    
    call io_close(w90_amn)

    call profiling_out(prof)

    POP_SUB(create_wannier90_amn)

  end subroutine create_wannier90_amn

  subroutine generate_wannier_states(space, mesh, ions, st, kpoints)
    type(space_t),          intent(in) :: space
    type(mesh_t),           intent(in) :: mesh
    type(ions_t),           intent(in) :: ions
    type(states_elec_t),    intent(in) :: st
    type(kpoints_t),        intent(in) :: kpoints

    integer :: w90_u_mat, w90_xyz, nwann, nik
    integer :: ik, iw, iw2, ip, ipmax
    logical :: what(MAX_OUTPUT_TYPES)
    integer(8) :: how(0:MAX_OUTPUT_TYPES)
    integer :: output_interval(0:MAX_OUTPUT_TYPES) 
    FLOAT, allocatable :: centers(:,:), dwn(:)
    CMPLX, allocatable :: Umnk(:,:,:)
    CMPLX, allocatable :: zwn(:), psi(:,:)
    character(len=MAX_PATH_LEN) :: fname
    FLOAT :: kpoint(1:MAX_DIM), wmod, wmodmax
    character(len=2) :: dum
    logical :: exist

    PUSH_SUB(generate_wannier_states)

    ASSERT(st%d%ispin /= SPINORS)

    inquire(file=trim(trim(adjustl(w90_prefix))//'_centres.xyz'),exist=exist)
    if(.not. exist) then
      message(1) = 'oct-wannier90: Cannot find the Wannier90 file seedname_centres.xyz.'
      write(message(2),'(a)') 'Please run wannier90.x with "write_xyz=.true." in '// trim(adjustl(w90_prefix)) // '.'
      call messages_fatal(2)
    end if

    w90_xyz = io_open(trim(trim(adjustl(w90_prefix))//'_centres.xyz'), global_namespace, action='read')
    SAFE_ALLOCATE(centers(1:3, 1:w90_num_wann))
    !Skip two lines
    read(w90_xyz, *)
    read(w90_xyz, *)
    do iw = 1, w90_num_wann
      read(w90_xyz, *) dum, centers(1:3, iw)
    end do    
    call io_close(w90_xyz) 


    inquire(file=trim(trim(adjustl(w90_prefix))//'_u_dis.mat'),exist=exist)
    if(exist) then
      message(1) = 'oct-wannier90: Calculation of Wannier states with disentanglement is not yet supported.'
      call messages_fatal(1)
    end if

    inquire(file=trim(trim(adjustl(w90_prefix))//'_u.mat'),exist=exist)
    if(.not. exist) then
      message(1) = 'oct-wannier90: Cannot find the Wannier90 seedname_u.mat file.'
      write(message(2),'(a)') 'Please run wannier90.x with "write_u_matrices=.true." in '// trim(adjustl(w90_prefix)) // '.'
      call messages_fatal(2)
    end if
    w90_u_mat = io_open(trim(trim(adjustl(w90_prefix))//'_u.mat'), global_namespace, action='read')    

    !To be read later
    w90_num_wann = w90_num_bands
    
    !Skip one line
    read(w90_u_mat, *)
    !Read num_kpts, num_wann, num_wann for consistency check
    read(w90_u_mat, *) nik, nwann, nwann
    if(nik /= w90_num_kpts .or. nwann /= w90_num_wann) then
      message(1) = "The file contains U matrices is inconsistent with the .win file."
      call messages_fatal(1)
    end if
   
    SAFE_ALLOCATE(Umnk(1:w90_num_wann, 1:w90_num_wann, 1:w90_num_kpts)) 

    do ik = 1, w90_num_kpts
      !Skip one line
      read(w90_u_mat, *)
      !Skip one line (k-point coordinate)
      read(w90_u_mat, *)
      read(w90_u_mat, '(f15.10,sp,f15.10)') ((Umnk(iw, iw2, ik), iw=1, w90_num_wann), iw2=1, w90_num_wann)
    end do
    
    call io_close(w90_u_mat)

    !We read the output format for the Wannier states
    what = .false.
    call io_function_read_what_how_when(global_namespace, space, what, how, output_interval)

    call io_mkdir('wannier', global_namespace)

    !Computing the Wannier states in the primitive cell, from the U matrices
    SAFE_ALLOCATE(zwn(1:mesh%np))
    SAFE_ALLOCATE(dwn(1:mesh%np))
    SAFE_ALLOCATE(psi(1:mesh%np, 1:st%d%dim))

    do iw = 1, w90_num_wann
      write(fname, '(a,i3.3)') 'wannier-', iw

      zwn(:) = M_Z0

      do ik = 1, w90_num_kpts
        kpoint(1:space%dim) = kpoints%get_point(ik, absolute_coordinates=.true.)

        do iw2 = 1, st%nst
          if(exclude_list(iw2)) cycle
          
          if(st%d%ispin == SPIN_POLARIZED) then
            call states_elec_get_state(st, mesh, iw2, ik, psi)
          else
            call states_elec_get_state(st, mesh, iw2, (ik-1)*2+w90_spin_channel, psi)
          end if
    
          !The minus sign is here is for the wrong convention of Octopus
          do ip = 1, mesh%np
            zwn(ip) = zwn(ip) + Umnk(band_index(iw2), iw, ik)/w90_num_kpts * psi(ip, 1) * &
                      exp(-M_zI* sum((mesh%x(ip, 1:space%dim)-centers(1:space%dim, iw)) * kpoint(1:space%dim)))
          end do
        end do!ik   
      end do!iw2


      !Following what Wannier90 is doing, we fix the global phase by setting the max to be real
      ipmax = 0
      wmodmax = M_z0
      do ip = 1, mesh%np
        wmod = TOFLOAT(zwn(ip)*conjg(zwn(ip)))
        if(wmod > wmodmax) then
          ipmax = ip
          wmodmax = wmod
        end if
      end do
      call lalg_scal(mesh%np, sqrt(wmodmax)/zwn(ipmax), zwn)
     

      do ip = 1, mesh%np
        dwn(ip) = TOFLOAT(zwn(ip))
      end do
      call dio_function_output(how(0), 'wannier', trim(fname), global_namespace, space, mesh, &
          dwn, unit_one, ierr, ions = ions, grp = st%dom_st_kpt_mpi_grp)
    end do

    SAFE_DEALLOCATE_A(Umnk)
    SAFE_DEALLOCATE_A(zwn)
    SAFE_DEALLOCATE_A(dwn)
    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(centers)

    POP_SUB(generate_wannier_states)
  end subroutine generate_wannier_states

end program wannier90_interface

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
