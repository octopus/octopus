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
  use calc_mode_par_oct_m
  use comm_oct_m
  use command_line_oct_m
  use geometry_oct_m
  use fft_oct_m
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
  use mpi_oct_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2 
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
  use states_restart_oct_m
  use submesh_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m

  implicit none

  integer              :: w90_what

  integer              :: ierr
  integer              :: dim, dir, how, idim, pdim
  integer              :: ii, i1,i2,i3, nik, iter, nst
  type(block_t)        :: blk

  type(restart_t)      :: restart
  type(system_t)       :: sys
  type(hamiltonian_t)  :: hm
  character(len=512)   :: filename, str, str2
  integer              :: ist, ispin
  type(states_t)       :: st
  logical              :: w90_setup, w90_output, w90_wannier
  logical              :: w90_spinors
  integer              :: w90_nntot, w90_num_bands, w90_num_kpts   ! w90 input parameters
  integer, allocatable :: w90_nnk_list(:,:)                        !
  character(len=80)    :: w90_prefix                               ! w90 input file prefix
  integer              :: w90_num_wann                             ! input paramter
  FLOAT, allocatable   :: w90_proj_centers(:,:)                    ! projections centers
  integer, allocatable :: w90_proj_lmr(:,:)                        ! definitions for real valued Y_lm*R_r
  integer :: w90_nproj                                             ! number of such projections        
  integer, allocatable :: w90_spin_proj_component(:)               ! up/down flag 
  FLOAT, allocatable   :: w90_spin_proj_axis(:,:)                  ! spin axis (not implemented)

  call getopt_init(ierr)
  if(ierr /= 0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-wannier90 command is not available."
    call messages_fatal(2)
  end if


  call global_init()

  call messages_init()
  call io_init()
  call calc_mode_par_init()

  call fft_all_init()
  call unit_system_init()

  call system_init(sys)

  !%Variable wannier90_prefix
  !%Type string
  !%Default w90
  !%Section Utilities::oct-wannier90
  !%Description
  !% Prefix for wannier90 files
  !%End
  call parse_variable('wannier90_prefix', 'w90', w90_prefix)
  if(w90_prefix=='w90') then
    message(1) = "Did not find wannier90_prefix keyword, will use default: w90"
    call  messages_warning(1)
  end if

  !%Variable wannier90_mode
  !%Type flag
  !%Default none
  !%Section Utilities::oct-wannier90
  !%Description
  !% Specifies which stage of the Wannier90 interface to use
  !%Option w90_setup bit(1)
  !% Writes parts of the wannier90 input file <tt>w90_prefix.win</tt> corresponding to
  !% the octopus inp file. Importantly it generates the correct form of Monkhorst-Pack mesh
  !% written to the file w90_kpoints that has to be used in a gs calculation of Octopus by
  !% as <tt> include w90_kpoints </tt> instead of the <tt>%KpointsGrid</tt> block
  !%Option w90_output bit(2)
  !% Generates the relevant files for a wannier90 run, specified by the variable <tt>W90_interface_files</tt>.
  !% This needs files previously generated
  !% by <tt>wannier90.x -pp w90 </tt>
  !% Option w90_wannier bit(3)
  !% Parse the output of wannier90 to generate the Wannier states on the real-space grid. 
  !% The states will be written in the folder wannier. By default, the states are written as
  !% binary files, similar to the Kohn-Sham states.
  !%End
  call parse_variable('wannier90_mode', 0, w90_what)
  w90_setup = iand(w90_what, OPTION__WANNIER90_MODE__W90_SETUP) /= 0
  w90_output = iand(w90_what, OPTION__WANNIER90_MODE__W90_OUTPUT) /= 0

  if(w90_what == 0) then
    message(1) = "wannier90_mode must be set to a value different from 0."
    call messages_fatal(1)
  end if

  !%Variable wannier90_files
  !%Type flag
  !%Default w90_mmn + w90_amn + w90_eig
  !%Section Utilities::oct-wannier90
  !%Description
  !% Specifies which files to generate
  !% Example: <tt>w90_mmn + w90_unk< /tt>
  !%Option w90_mmn bit(1)
  !% (see Wannier90 documentation)
  !%Option w90_unk bit(2)
  !% (see Wannier90 documentation)
  !%Option w90_amn bit(3)
  !% (see Wannier90 documentation)
  !%Option w90_eig bit(4)
  !% Eigenvalues. See Wannier90 documentation for more details.
  !%End
  w90_what = OPTION__WANNIER90_FILES__W90_MMN + OPTION__WANNIER90_FILES__W90_AMN + OPTION__WANNIER90_FILES__W90_EIG
  call parse_variable('wannier90_files', w90_what, w90_what)


  ! sanity checks
  if(w90_setup .and. w90_output) then
    message(1) = 'oct-wannier90: wannier90_setup and wannier90_output are mutually exclusive'
    call messages_fatal(1)
  end if

  if(.not.sys%gr%sb%kpoints%w90_compatible) then
    message(1) = 'oct-wannier90: Only Wannier90KPointsGrid=yes is allowed for running oct-wannier90'
    call messages_fatal(1)
  end if


  ! create setup files
  if(w90_setup) then
    call wannier90_setup()

  ! load states and calculate interface files
  elseif(w90_output) then

    ! normal interface run
    call states_init(st, sys%gr, sys%geo)
    call kpoints_distribute(st%d,sys%mc)
    call states_distribute_nodes(st,sys%mc)
    call states_exec_init(st, sys%mc)
    call restart_module_init()
    call states_allocate_wfns(st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)
    call restart_init(restart, RESTART_GS, RESTART_TYPE_LOAD, &
                       sys%mc, ierr, sys%gr%der%mesh)
    if(ierr == 0) then
      call states_look(restart, nik, dim, nst, ierr)
      if(dim==st%d%dim .and. nik==sys%gr%sb%kpoints%reduced%npoints .and. nst==st%nst) then
         call states_load(restart, st, sys%gr, ierr, iter)
      else
         write(message(1),'(a)') 'Restart structure not commensurate.'
         call messages_fatal(1)
      end if
    end if
    call restart_end(restart)

    ! ---- actual interface work ----------
    call read_wannier90_files()

    if(iand(w90_what, OPTION__WANNIER90_FILES__W90_MMN) /= 0) then
      call create_wannier90_mmn()
    end if

    if(iand(w90_what, OPTION__WANNIER90_FILES__W90_UNK) /= 0) then
      call write_unk()
    end if

    if(iand(w90_what, OPTION__WANNIER90_FILES__W90_AMN) /= 0) then
      call create_wannier90_amn()
    end if

    if(iand(w90_what, OPTION__WANNIER90_FILES__W90_EIG) /= 0) then
      call create_wannier90_eig()
    end if

  else if(w90_wannier) then
    call messages_not_implemented("The option wannier90_mode = w90_wannier is not yet implemented.")
  end if

  call system_end(sys)
  call fft_all_end()
  call io_end()
  call messages_end()
  call global_end()

contains

  subroutine wannier90_setup()
    character(len=80) :: filename
    integer :: w90_win, oct_kpts, ia, axis(3), jj, kk

    PUSH_SUB(wannier90_setup)

    call states_init(st, sys%gr, sys%geo)
    ! open win file
    filename = trim(adjustl(w90_prefix)) //'.win'
    w90_win = io_open(trim(filename), action='write')

    write(w90_win,'(a)') '# this file has been created by the Octopus wannier90 utility'
    write(w90_win,'(a)') ' '

    ! write direct lattice vectors (in angstrom)
    write(w90_win,'(a)') 'begin unit_cell_cart'
    write(w90_win,'(a)') 'Ang'
    do idim=1,3
      write(w90_win,'(f12.8,f12.8,f12.8)') units_from_atomic(unit_angstrom, sys%gr%sb%rlattice(idim,1:3))
    end do
    write(w90_win,'(a)') 'end unit_cell_cart'
    write(w90_win,'(a)') ' '

    ! for the moment projections are not implemented
    write(w90_win,'(a)') 'use_bloch_phases = .true.'
    write(w90_win,'(a)') ' '

    write(w90_win,'(a)') 'begin atoms_frac'
    do ia=1,sys%geo%natoms
       write(w90_win,'(a,2x,f12.8,f12.8,f12.6)') trim(sys%geo%atom(ia)%label), & 
         matmul(sys%geo%atom(ia)%x(1:3), sys%gr%sb%klattice(1:3, 1:3))/(M_TWO*M_PI) 
    end do
    write(w90_win,'(a)') 'end atoms_frac'
    write(w90_win,'(a)') ' '

    write(w90_win,'(a10,i4)') 'num_bands ', st%nst
    write(w90_win,'(a9,i4)') 'num_wann ', st%nst
    write(w90_win,'(a)') ' '

    if(.not.parse_is_defined('KPointsGrid')) then
       message(1) = 'oct-wannier90: need Monkhorst-Pack grid. Please specify %KPointsGrid'
       call messages_fatal(1)
    end if
    if(.not.sys%gr%sb%kpoints%w90_compatible) then
      message(1) = 'oct-wannier90: Only Wannier90KPointsGrid=yes is allowed for running oct-wannier90'
      call messages_fatal(1)
    end if

    axis(1:3) = sys%gr%sb%kpoints%nik_axis(1:3)
    write(w90_win,'(a8,i4,i4,i4)')  'mp_grid ', axis(1:3)
    write(w90_win,'(a)') ' '

    ! make wannier90 compliant MonkhorstPack mesh
    ! and write simultaneously to w90_prefix.win file and w90_kpoints for octopus input
    write(w90_win,'(a)')  'begin kpoints '

    do ii=0,axis(1)-1
       do jj=0,axis(2)-1
          do kk=0,axis(3)-1
             write(w90_win,'(f12.8,f12.8,f12.8)') ii*M_ONE/(axis(1)*M_ONE), jj*M_ONE/(axis(2)*M_ONE), kk*M_ONE/(axis(3)*M_ONE)
          end do
       end do
    end do
    write(w90_win,'(a)')  'end kpoints '

    call io_close(w90_win)

    POP_SUB(wannier90_setup)

  end subroutine wannier90_setup

  subroutine read_wannier90_files()
    integer ::  w90_nnkp, itemp
    character(len=80) :: filename, dummy, dummy1
    logical :: exist
    FLOAT :: dummyr(7)

    PUSH_SUB(read_wannier90_files)

    ! assume to use all bands and number of k-points is consistent with Wannier90
    ! input files. Consistncy is checked later
    w90_num_bands = st%nst

    w90_num_kpts = sys%gr%sb%kpoints%full%npoints

    ! open nnkp file
    filename = trim(adjustl(w90_prefix)) //'.nnkp'

    inquire(file=filename,exist=exist)
    if(.not. exist) then
       message(1) = 'oct-wannier90: Cannot find specified Wannier90 nnkp file.'
       write(message(2),'(a)') 'Please run wannier90.x -pp '// trim(adjustl(w90_prefix)) // ' first.' 
       call messages_fatal(2)
    end if

    w90_nnkp = io_open(trim(filename), action='read')

    ! check number of k-points
    do while(.true.)
       read(w90_nnkp,*) dummy, dummy1
       if(dummy =='begin' .and. dummy1 == 'kpoints' ) then
          read(w90_nnkp,*) itemp
          if(itemp /= w90_num_kpts) then
             message(1) = 'oct-wannier90: wannier90 setup seems to have been done with a different number of k-points.'
             call messages_fatal(1)
          else
             exit
          end if
       end if
    end do
   close(w90_nnkp)

    w90_nnkp = io_open(trim(filename), action='read')
    ! read from nnkp file
    ! find the nnkpts block
    do while(.true.)
       read(w90_nnkp,*,end=101) dummy, dummy1
       if(dummy =='begin' .and. dummy1 == 'nnkpts' ) then
          read(w90_nnkp,*) w90_nntot
          SAFE_ALLOCATE(w90_nnk_list(w90_num_kpts*w90_nntot,5))
          do ii=1,w90_num_kpts*w90_nntot
             read(w90_nnkp,*) w90_nnk_list(ii,1:5)
          end do
          !make sure we are at the end of the block
          read(w90_nnkp,*) dummy
          if(dummy /= 'end') then
             message(1) = 'oct-wannier90: There dont seem to be enough k-points in nnkpts file to.'
             call messages_fatal(1)
          end if
          goto 102
       end if
    end do

    ! jump point when EOF found while looking for nnkpts block
101 message(1) = 'oct-wannier90: Did not find nnkpts block in nnkp file'
    call messages_fatal(1)

102 continue

    call io_close(w90_nnkp)

    if(iand(w90_what, OPTION__WANNIER90_FILES__W90_AMN) /= 0) then
       ! parse file again for definitions of projections
       w90_nnkp = io_open(trim(filename), action='read')

       do while(.true.)
          read(w90_nnkp,*,end=201) dummy, dummy1

          if(dummy =='begin' .and. (dummy1 == 'projections' .or. dummy1 == 'spinor_projections')) then

              if(dummy1 == 'spinor_projections') then
                 w90_spinors = .true.
                 message(1) = 'oct-wannier90: Spinor interface incomplete. Note there is no quantization axis implemented'
                 call messages_warning(1)
              end if

              read(w90_nnkp,*) w90_nproj
              ! num_wann is given in w90.win, not double checked here
              ! I assume that the wannier90.x -pp run has checked this
              w90_num_wann = w90_nproj

              SAFE_ALLOCATE(w90_proj_centers(w90_nproj,3))
              SAFE_ALLOCATE(w90_proj_lmr(w90_nproj,3))
              if(w90_spinors) SAFE_ALLOCATE(w90_spin_proj_component(w90_nproj))
              if(w90_spinors) SAFE_ALLOCATE(w90_spin_proj_axis(w90_nproj,3))

              do ii=1,w90_nproj
                 read(w90_nnkp,*) w90_proj_centers(ii,1:3), w90_proj_lmr(ii,1:3)
                 ! skip a line for now
                 read(w90_nnkp,*) dummyr(1:7)
                 if(w90_spinors) then
                    read(w90_nnkp,*) w90_spin_proj_component(ii), w90_spin_proj_axis(ii,1:3)
                    ! use octopus spindim conventions
                    if(w90_spin_proj_component(ii) == -1) w90_spin_proj_component(ii) = 2
                 end if
              end do
              !make sure we are at the end of the block
              read(w90_nnkp,*) dummy
              if(dummy /= 'end') then
                 message(1) = 'oct-wannier90: There dont seem to be enough projections in nnkpts file to.'
                 call messages_fatal(1)
              end if
              goto 202
          end if
       end do

       ! jump point when EOF found while looking for projections block
201    message(1) = 'oct-wannier90: Did not find projections block in w90.nnkp file'
       call messages_fatal(1)

       ! jump point when projections is found in file
202    continue

    end if

    call io_close(w90_nnkp)

    POP_SUB(read_wannier90_files)

  end subroutine read_wannier90_files

  subroutine create_wannier90_mmn()
    integer ::  ist, jst, ik, w90_mmn, iknn, G(3), ii, jj, idim, idim2
    FLOAT   ::  Gr(3), t1, t2
    character(len=80) :: filename
    CMPLX   :: overlap
    CMPLX, allocatable :: state1(:,:), state2(:,:)

    PUSH_SUB(create_wannier90_mmn)

    if(st%d%kpt%parallel) then
      call messages_not_implemented("w90_mmn output with k-point parallelization.")
    end if

    if(st%parallel_in_states) then
      call messages_not_implemented("w90_mmn output with states parallelization.")
    end if

    filename = './'// trim(adjustl(w90_prefix))//'.mmn'
    w90_mmn = io_open(trim(filename), action='write')

    ! write header
    if(mpi_grp_is_root(mpi_world)) then
       write(w90_mmn,*) ' '
       write(w90_mmn,*) w90_num_bands, w90_num_kpts, w90_nntot
    end if

    SAFE_ALLOCATE(state1(1:sys%gr%der%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(state2(1:sys%gr%der%mesh%np, 1:st%d%dim))

    ! loop over the pairs specified in the nnkp file (read before in init)
    do ii=1,w90_num_kpts*w90_nntot
       ik = w90_nnk_list(ii,1)
       iknn = w90_nnk_list(ii,2)
       G(1:3) = w90_nnk_list(ii,3:5)
       if(mpi_grp_is_root(mpi_world)) write(w90_mmn,'(I10,2x,I10,2x,I3,2x,I3,2x,I3)') ik, iknn, G(1:3)

       Gr(1:3) = matmul(G(1:3),sys%gr%sb%klattice(1:3,1:3))

       ! loop over bands
       do jst=1,w90_num_bands
          call states_get_state(st, sys%gr%der%mesh, jst, iknn, state2)

         ! add phase
         do jj=1,sys%gr%der%mesh%np
           state2(jj,1:st%d%dim) = state2(jj,1:st%d%dim)*exp(M_zI*dot_product(sys%gr%der%mesh%x(jj,1:3),Gr(1:3)))
         end do

         do ist=1,w90_num_bands
           call states_get_state(st, sys%gr%der%mesh, ist, ik, state1)

           overlap = M_ZERO
           do idim=1,st%d%dim
              overlap =overlap +  zmf_dotp(sys%gr%der%mesh, state1(:,idim), state2(:,idim))
           end do

           ! write to W90 file
           if(mpi_grp_is_root(mpi_world)) write(w90_mmn,'(e12.6,2x,e12.6)') overlap
         end do
       end do

    end do

    close(w90_mmn)

    SAFE_DEALLOCATE_A(state1)
    SAFE_DEALLOCATE_A(state2)

    POP_SUB(create_wannier90_mmn)

  end subroutine create_wannier90_mmn

  subroutine create_wannier90_eig()
    integer ::  ist, ik, w90_eig
    character(len=80) :: filename

    PUSH_SUB(create_wannier90_eig)

    if(st%d%kpt%parallel) then
      call messages_not_implemented("w90_eig output with k-point parallelization.")
    end if

    if(st%parallel_in_states) then
      call messages_not_implemented("w90_eig output with states parallelization.")
    end if

    if(mpi_grp_is_root(mpi_world)) then
      filename = './'//trim(adjustl(w90_prefix))//'.eig'
      w90_eig = io_open(trim(filename), action='write')
      do ik=1,w90_num_kpts
        do ist=1,w90_num_bands
          write(w90_eig,'(I5,2x,I5,2x,e12.6)') ist, ik, units_from_atomic(units_out%energy, st%eigenval(ist, ik))
        end do
      end do

      close(w90_eig)
    end if

    POP_SUB(create_wannier90_eig)
  end subroutine create_wannier90_eig

  subroutine write_unk()
    integer ::  ist, jst, ik, unk_file,  iknn, G(3), ii, jj, idim, idim2,  ispin
    integer :: nr(2,3), ix, iy, iz, ip
    FLOAT   ::  Gr(3), t1, t2
    character(len=80) :: filename
    CMPLX   :: overlap
    CMPLX, allocatable :: state1(:,:), state2(:)

    PUSH_SUB(write_unk)

    if(st%d%kpt%parallel) then
      call messages_not_implemented("w90_unk output with k-point parallelization.")
    end if
      
    if(sys%gr%mesh%parallel_in_domains) then
      call messages_not_implemented("w90_unk output with domain parallelization")
    end if

    if(st%parallel_in_states) then
      call messages_not_implemented("w90_unk output with states parallelization.")
    end if


    SAFE_ALLOCATE(state1(1:sys%gr%der%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(state2(1:sys%gr%der%mesh%np))

    ! set boundaries of inner box
    nr(1,1:3) = sys%gr%mesh%idx%nr(1,1:3) + sys%gr%mesh%idx%enlarge(1:3)
    nr(2,1:3) = sys%gr%mesh%idx%nr(2,1:3) - sys%gr%mesh%idx%enlarge(1:3)

    do ik=1,w90_num_kpts
       do ispin=1,st%d%dim

          if(mpi_grp_is_root(mpi_world)) then
             write(filename,'(a,i5.5,a1,i1)') './UNK', ik,'.', ispin
             unk_file = io_open(trim(filename), action='write',form='unformatted')
             ! write header
             write(unk_file) sys%gr%mesh%idx%ll(1:sys%gr%mesh%idx%dim), ik,  w90_num_bands
             ! states
             do ist=1,w90_num_bands
                call states_get_state(st, sys%gr%der%mesh, ist, ik, state1)
                ! reorder state
                ip=0
                do iz=nr(1,3),nr(2,3)
                   do iy=nr(1,2),nr(2,2)
                      do ix=nr(1,1),nr(2,1)
                         ip=ip+1
                         state2(ip) =  state1(sys%gr%mesh%idx%lxyz_inv(ix, iy, iz),ispin)
                      end do
                   end do
                end do
                write(unk_file) state2(:)
             end do
             call io_close(unk_file)
          end if
       end do
    end do

    POP_SUB(write_unk)

  end subroutine write_unk

  subroutine create_wannier90_amn()
    integer ::  ist, ik, w90_amn, idim, iw, ip
    FLOAT   ::  center(3)
    character(len=80) :: filename
    CMPLX   :: projection
    CMPLX, allocatable :: psi(:,:), orbital(:,:)
    type(submesh_t) :: submesh
    FLOAT, allocatable ::  rr(:,:), ylm(:)

    PUSH_SUB(create_wannier90_amn)

    if(st%d%kpt%parallel) then
      call messages_not_implemented("w90_amn output with k-point parallelization.")
    end if

    if(st%parallel_in_states) then
      call messages_not_implemented("w90_amn output with states parallelization.")
    end if

    ! precompute orbitals
    SAFE_ALLOCATE(orbital(1:sys%gr%mesh%np, 1:w90_nproj))
    orbital = M_ZERO
    do iw=1,w90_nproj
      ! cartesian coordinate of orbital center
      center(1:3) =  matmul(w90_proj_centers(iw,1:3), sys%gr%sb%rlattice(1:3,1:3))
      call submesh_init(submesh, sys%gr%sb, sys%gr%mesh, center, CNST(2.0)*maxval(sys%gr%sb%lsize(:)))

      ! make transpose table of submesh points for use in pwscf routine
      SAFE_ALLOCATE(rr(1:3,submesh%np))
      do ip=1,submesh%np
        rr(1:3,ip) = submesh%x(ip,1:3)
      end do

      ! get ylm as submesh points
      SAFE_ALLOCATE(ylm(1:submesh%np))
      ! (this is a routine from pwscf)
      call ylm_wannier(ylm, w90_proj_lmr(iw,1), w90_proj_lmr(iw,2), rr, submesh%np)

      ! apply radial function
      do ip=1,submesh%np
        ylm(ip) = ylm(ip)*M_TWO*exp(-submesh%x(ip,0))
      end do

      call submesh_add_to_mesh(submesh, ylm, orbital(1:sys%gr%mesh%np, iw))

      SAFE_DEALLOCATE_A(ylm)
      SAFE_DEALLOCATE_A(rr)
      call submesh_end(submesh)
    end do

    filename = './'// trim(adjustl(w90_prefix))//'.amn'
    w90_amn = io_open(trim(filename), action='write')

    ! write header
    if(mpi_grp_is_root(mpi_world)) then
      write(w90_amn,*) ' '
      write(w90_amn,*)  w90_num_bands, w90_num_kpts, w90_num_wann
    end if

    SAFE_ALLOCATE(psi(1:sys%gr%der%mesh%np, 1:st%d%dim))

    do ist=1, w90_num_bands
      do iw=1, w90_nproj
        do ik=1, w90_num_kpts
          call states_get_state(st, sys%gr%der%mesh, ist, ik, psi)
          projection = M_ZERO
          do idim=1,st%d%dim
            if(w90_spinors .and. idim /= w90_spin_proj_component(iw)) cycle
            projection = projection + zmf_dotp(sys%gr%mesh,psi(1:sys%gr%mesh%np,idim),orbital(1:sys%gr%mesh%np, iw))
          end do
          write (w90_amn,'(I5,2x,I5,2x,I5,2x,e12.6,2x,e12.6)') ist, iw, ik, projection
        end do
      end do
    end do
    call io_close(w90_amn)

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(orbital)

    POP_SUB(create_wannier90_amn)

  end subroutine create_wannier90_amn

  ! the definitions of atomic orbitals are taken from this file:
#include "wannier90_interface_defs_from_pwscf.F90"

end program wannier90_interface

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
