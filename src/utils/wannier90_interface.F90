!! Copyright (C) 2017 H. Huebener
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
  use floquet_oct_m
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
  logical              :: w90_setup, w90_output, w90_floquet, w90_unk
  integer :: w90_nntot, w90_num_bands, w90_num_kpts    ! w90 input parameters
  integer, allocatable ::  w90_nnk_list(:,:)           !
  character(len=80) :: w90_prefix                      ! w90 input file prefix
  integer :: w90_num_wann                              ! input paramter
!not used for the momment
!  logical :: w90_proj                                  ! whether to comput projections for  w90 intialization
!  FLOAT, allocatable :: w90_proj_centers(:,:)          ! projections centers
!  integer, allocatable ::  w90_proj_lmr(:,:)           ! definitions for real valued Y_lm*R_r
!  integer :: w90_nproj                                 ! number of such projections


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

  call system_init(sys)

  !%Variable W90_prefix
  !%Type string
  !%Default w90
  !%Section Utilities::oct-wannier90_interface
  !%Description
  !% Prefix for Wannier90 files
  !%End
  !    call parse_string('W90_prefix', 'w90', w90_prefix)
  !  if(w90_prefix=='w90') then
  !     message(1) = "Did not find w90_prefix keyword, will use default: w90"
  !     call  messages_warning(1)
  !  end if
  w90_prefix = 'w90'

  !%Variable W90_interface_mode
  !%Type flag
  !%Default none
  !%Section Utilities::oct-wannier90_interface
  !%Description
  !% Specifies what to do for the Wannier90 interface
  !% Example: <tt>w90_setup + w90_floquet</tt>
  !%Option w90_setup bit(1)
  !% Writes parts of the wannier90 input file <tt>w90_prefix.win</tt> corresponding to 
  !% the octopus inp file. Importantly it generates the correct form of Monkhorst-Pack mesh
  !% written to the file w90_kpoints that has to be used in a gs calculation of Octopus by
  !% as <tt> input w90_kpoints <tt> instead of the <tt>%KpointsGrid<tt> block 
  !%Option w90_output bit(2)
  !% Generates the relevant files for a wannier90 run. This needs files generated 
  !% by <tt>wannier90.x -pp w90 <tt>
  !%Option w90_floquet bit(3)
  !% Make interface for a Floquet structure. Takes care of different dimensionalities etc. 
  !%
  !%Option w90_unk bit(4)
  !% write unk files for plotting routines of wannier90 
  !%
  !%End
  call parse_variable('W90_interface_mode', w90_what, w90_what)

  w90_setup = iand(w90_what, OPTION__W90_INTERFACE_MODE__W90_SETUP) /= 0
  w90_output = iand(w90_what, OPTION__W90_INTERFACE_MODE__W90_OUTPUT) /= 0
  w90_floquet = iand(w90_what, OPTION__W90_INTERFACE_MODE__W90_FLOQUET) /= 0
  w90_unk = iand(w90_what, OPTION__W90_INTERFACE_MODE__W90_UNK) /= 0


  ! sanity checks
  if(w90_setup .and. w90_output) then
      message(1) = 'W90: w90_setup and w90_output are mutually exclusive'
      call messages_fatal(1)
  end if

  if(w90_floquet .and. .not. w90_setup .and. .not. w90_output) then
      message(1) = 'W90: w90_floquet has to be combined with either w90_setup or w90_output.'
      call messages_fatal(1)
  end if

  if(w90_unk .and. .not. w90_output) then
      message(1) = 'W90: w90_unk has to be combined with w90_output.'
      call messages_fatal(1)
  end if

  ! create setup files
  if(w90_setup) then
     call wannier90_setup()

  ! load states and calculate interface files
  elseif(w90_output) then

    ! normal interface run
    if(.not.w90_floquet ) then
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

    ! floquet interface run
    elseif(w90_floquet) then
      call floquet_init(sys,hm%F,sys%st%d%dim)

      call states_init(st, sys%gr, sys%geo,floquet_dim=hm%F%floquet_dim)
      call kpoints_distribute(st%d,sys%mc)
      call states_distribute_nodes(st,sys%mc)
      call states_exec_init(st, sys%mc)
      call restart_module_init()
      call states_allocate_wfns(st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)
      call floquet_restart_dressed_st(hm, sys, st, ierr)
      call messages_write('Read Floquet restart files.')
      call messages_info()

    end if

    ! ---- actual interface work ----------
    call read_wannier90_files()
    call create_wannier90_mmn()
    if(w90_unk) call write_unk()

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

    if(.not. w90_floquet ) then
      call states_init(st, sys%gr, sys%geo)
    else
      call floquet_init(sys,hm%F,sys%st%d%dim)
      call states_init(st, sys%gr, sys%geo,floquet_dim=hm%F%floquet_dim)
    end if

    ! open win file
    filename = trim(adjustl(w90_prefix)) //'.win'
    w90_win = io_open(trim(filename), action='write')

    write(w90_win,'(a)') '# this file has been created by the Octopus wannier90 utility'
    write(w90_win,'(a)') ' '

    ! write direct lattice vectors (in angstrom)
    write(w90_win,'(a)') 'begin unit_cell_cart'
    do idim=1,3
       write(w90_win,'(f12.8,f12.8,f12.8)') sys%gr%sb%rlattice(idim,1:3)*0.529177249
    end do
    write(w90_win,'(a)') 'end unit_cell_cart'
    write(w90_win,'(a)') ' '

    ! for the moment projections are not implemented
    write(w90_win,'(a)') 'use_bloch_phases = .true.'
    write(w90_win,'(a)') ' '

    write(w90_win,'(a)') 'begin atoms_cart'
    do ia=1,sys%geo%natoms
       write(w90_win,'(a,2x,f12.8,f12.8,f12.6)') trim(sys%geo%atom(ia)%label), sys%geo%atom(ia)%x(1:3)
    end do
    write(w90_win,'(a)') 'end atoms_cart'
    write(w90_win,'(a)') ' '

    write(w90_win,'(a10,i4)') 'num_bands ', st%nst
    write(w90_win,'(a9,i4)') 'num_wann ', st%nst
    write(w90_win,'(a)') ' '

    if(.not.parse_is_defined('KPointsGrid')) then
       message(1) = 'W90: need Monkhorst-Pack grid. Please specify %KPointsGrid'
       call messages_fatal(1)
    end if
    if(parse_is_defined('KPointsPath')) then
       message(1) = 'W90: can only run with Monkhorst-Pack grid. Please specify only %KPointsGrid'
       call messages_fatal(1)
    end if

    axis(1:3) = sys%gr%sb%kpoints%nik_axis(1:3)
    write(w90_win,'(a8,i4,i4,i4)')  'mp_grid ', axis(1:3)
    write(w90_win,'(a)') ' '

    ! make wannier90 compliant MonkhorstPack mesh
    ! and write simultaneously to w90_prefix.win file and w90_kpoints for octopus input
    filename = 'w90_kpoints'
    oct_kpts = io_open(trim(filename), action='write')
    write(oct_kpts,'(a)') '%KpointsReduced'
    write(w90_win,'(a)')  'begin kpoints '

    do ii=0,axis(1)-1
       do jj=0,axis(2)-1
          do kk=0,axis(3)-1
             write(w90_win,'(f12.8,f12.8,f12.8)')  ii/real(axis(1)), jj/real(axis(2)), kk/real(axis(3))
             write(oct_kpts,'(a6,f12.8,a3,f12.8,a3,f12.8)') ' 1. | ',  ii/real(axis(1)),' | ', jj/real(axis(2)), ' | ', kk/real(axis(3))
          end do
       end do
    end do
    write(oct_kpts,'(a)') '%'
    write(w90_win,'(a)')  'end kpoints '

  end subroutine wannier90_setup

  subroutine read_wannier90_files()
    integer ::  w90_nnkp, itemp
    character(len=80) :: filename, dummy, dummy1
    logical :: exist

    ! assume to use all bands and number of k-points is consistent with Wannier90
    ! input files. Consistncy is checked later
    w90_num_bands = st%nst

    w90_num_kpts = sys%gr%sb%kpoints%full%npoints

    ! open nnkp file
    filename = trim(adjustl(w90_prefix)) //'.nnkp'

    inquire(file=filename,exist=exist)
    if(.not. exist) then
       message(1) = 'W90: Cannot find specified Wannier90 nnkp file.'
       write(message(2),'(a)') 'W90: Please run wannier90.x -pp '// trim(adjustl(w90_prefix)) // ' first.' 
       call messages_fatal(2)
    end if

    w90_nnkp = io_open(trim(filename), action='read')

    ! check number of k-points
    do while(.true.)
       read(w90_nnkp,*) dummy, dummy1
       if(dummy =='begin' .and. dummy1 == 'kpoints' ) then
          read(w90_nnkp,*) itemp
          if(itemp /= w90_num_kpts) then
             message(1) = 'W90: wannier90 setup seems to have been done with a different number of k-points.'
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
             message(1) = 'W90: There dont seem to be enough k-points in nnkpts file to.'
             call messages_fatal(1)
          end if
          goto 102
       end if
    end do

    ! jump point when EOF found while looking for nnkpts block
101 message(1) = 'W90: Did not find nnkpts block in nnkp file'
    call messages_fatal(1)

102 continue

    close(w90_nnkp)

    ! TODO here also read projections file, if necessary

  end subroutine read_wannier90_files

  subroutine create_wannier90_mmn()
    integer ::  ist, jst, ik, w90_mmn, w90_eig, w90_amn,  iknn, G(3), ii, jj, idim, idim2
    FLOAT   ::  Gr(3), t1, t2
    character(len=80) :: filename
    CMPLX   :: overlap
    CMPLX, allocatable :: state1(:,:), state2(:,:)

    ASSERT(st%d%kpt%start==1 .and. st%d%kpt%end==sys%gr%sb%kpoints%full%npoints)

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

       Gr(1:3) = matmul(G(1:3),sys%gr%sb%kpoints%klattice(1:3,1:3))

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


    SAFE_DEALLOCATE_A(state1)
    SAFE_DEALLOCATE_A(state2)

  end subroutine create_wannier90_mmn

  subroutine write_unk()
    integer ::  ist, jst, ik, unk_file,  iknn, G(3), ii, jj, idim, idim2,  ispin
    integer :: nr(2,3), ix, iy, iz, ip
    FLOAT   ::  Gr(3), t1, t2
    character(len=80) :: filename
    CMPLX   :: overlap
    CMPLX, allocatable :: state1(:,:), state2(:)

    ASSERT(st%d%kpt%start==1 .and. st%d%kpt%end==sys%gr%sb%kpoints%full%npoints)
    ASSERT(sys%gr%mesh%np==sys%gr%mesh%np_global)

    SAFE_ALLOCATE(state1(1:sys%gr%der%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(state2(1:sys%gr%der%mesh%np))

    ! set boundaries of inner box
    nr(1,:) = sys%gr%mesh%idx%nr(1,:) + sys%gr%mesh%idx%enlarge(:)
    nr(2,:) = sys%gr%mesh%idx%nr(2,:) - sys%gr%mesh%idx%enlarge(:)

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
                         ip = sys%gr%mesh%idx%lxyz_inv(ix, iy, iz)
                         ip=ip+1
                         state2(ip) =  state1(sys%gr%mesh%idx%lxyz_inv(ix, iy, iz),ispin)
                      end do
                   end do
                end do
                write(unk_file) state2(:)
             end do
             close(unk_file)
          end if
       end do
    end do

  end subroutine write_unk

end program wannier90_interface

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
