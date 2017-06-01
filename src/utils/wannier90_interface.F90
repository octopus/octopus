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
  logical              :: w90_setup, w90_output, w90_floquet, w90_unk, w90_amn, w90_mmn
  integer              ::w90_nntot, w90_num_bands, w90_num_kpts    ! w90 input parameters
  integer, allocatable ::  w90_nnk_list(:,:)                       !
  character(len=80)    :: w90_prefix                               ! w90 input file prefix
  integer              :: w90_num_wann                             ! input paramter
  FLOAT, allocatable :: w90_proj_centers(:,:)                      ! projections centers
  integer, allocatable ::  w90_proj_lmr(:,:)                       ! definitions for real valued Y_lm*R_r
  integer :: w90_nproj                                             ! number of such projections        

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
  !% Specifies which stage of the Wannier90 interface to use
  !%Option w90_setup bit(1)
  !% Writes parts of the wannier90 input file <tt>w90_prefix.win</tt> corresponding to
  !% the octopus inp file. Importantly it generates the correct form of Monkhorst-Pack mesh
  !% written to the file w90_kpoints that has to be used in a gs calculation of Octopus by
  !% as <tt> input w90_kpoints </tt> instead of the <tt>%KpointsGrid</tt> block
  !%Option w90_output bit(2)
  !% Generates the relevant files for a wannier90 run, specified by the variable <tt>W90_interface_files</tt>.
  !% This needs files previously generated
  !% by <tt>wannier90.x -pp w90 </tt>
  !%Option w90_floquet bit(3)
  !% Make interface for a Floquet structure. Takes care of different dimensionalities etc.
  !%
  !%End
  call parse_variable('W90_interface_mode', w90_what, w90_what)

  w90_setup = iand(w90_what, OPTION__W90_INTERFACE_MODE__W90_SETUP) /= 0
  w90_output = iand(w90_what, OPTION__W90_INTERFACE_MODE__W90_OUTPUT) /= 0

  !%Variable W90_interface_files
  !%Type flag
  !%Default w90_mmn
  !%Section Utilities::oct-wannier90_interface
  !%Description
  !% Specifies which files to generate
  !% Example: <tt>w90_mmn + w90_unk< /tt>
  !%Option w90_mmn bit(1)
  !% (see Wannier90 documentation)
  !%Option w90_unk bit(2)
  !% (see Wannier90 documentation)
  !%Option w90_amn bit(3)
  !% (see Wannier90 documentation)
  !%End
  w90_what = 0
  call parse_variable('W90_interface_files', w90_what, w90_what)

  w90_unk = iand(w90_what, OPTION__W90_INTERFACE_FILES__W90_UNK) /= 0
  w90_mmn = iand(w90_what, OPTION__W90_INTERFACE_FILES__W90_MMN) /= 0
  w90_amn = iand(w90_what, OPTION__W90_INTERFACE_FILES__W90_AMN) /= 0

  ! default
  if(w90_what == 0) then
     w90_mmn = .true.
     w90_amn = .true.
  end if

  ! sanity checks
  if(w90_setup .and. w90_output) then
      message(1) = 'W90: w90_setup and w90_output are mutually exclusive'
      call messages_fatal(1)
  end if

  if(w90_floquet .and. .not. w90_setup .and. .not. w90_output) then
      message(1) = 'W90: w90_floquet has to be combined with either w90_setup or w90_output.'
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
    if(w90_mmn) call create_wannier90_mmn()
    if(w90_unk) call write_unk()
    if(w90_amn) call create_wannier90_amn()
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
             write(w90_win,'(f12.8,f12.8,f12.8)')  ii*M_ONE/(axis(1)*M_ONE), jj*M_ONE/(axis(2)*M_ONE),kk*M_ONE/(axis(3)*M_ONE)
             write(oct_kpts,'(a6,f12.8,a3,f12.8,a3,f12.8)') ' 1. | ',  ii*M_ONE/(axis(1)*M_ONE) ,' | ', jj*M_ONE/(axis(2)*M_ONE), ' | ', kk*M_ONE/(axis(3)*M_ONE)
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
    FLOAT :: dummyr(7)

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

    call io_close(w90_nnkp)

    if(w90_amn) then
       ! parse file again for definitions of projections
       w90_nnkp = io_open(trim(filename), action='read')
       do while(.true.)
          read(w90_nnkp,*,end=101) dummy, dummy1
          if(dummy =='begin' .and. dummy1 == 'projections' ) then
             goto 202
          end if
       end do
       message(1) = 'W90: Did not find projections in w90.nnkp file'
       call messages_fatal(2)

       ! jump point when projections is found in file
202    continue
       read(w90_nnkp,*) w90_nproj
       ! num_wann is given in w90.win, not double checked here
       ! I assume that the wannier90.x -pp run has checked this
       w90_num_wann = w90_nproj

       SAFE_ALLOCATE(w90_proj_centers(w90_nproj,3))
       SAFE_ALLOCATE(w90_proj_lmr(w90_nproj,3))
       do ii=1,w90_nproj
          read(w90_nnkp,*) w90_proj_centers(ii,1:3), w90_proj_lmr(ii,1:3)
          ! skip a line for now
          read(w90_nnkp,*) dummyr(1:7)
       end do
    end if

    call io_close(w90_nnkp)

  end subroutine read_wannier90_files

  subroutine create_wannier90_mmn()
    integer ::  ist, jst, ik, w90_mmn, w90_eig,  iknn, G(3), ii, jj, idim, idim2
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

  end subroutine write_unk

  subroutine create_wannier90_amn()
    integer ::  ist, jst, ik, w90_amn,  iknn, G(3), ii, jj, idim, idim2, iw, ip
    FLOAT   ::  Gr(3), t1, t2, center(3), dd
    character(len=80) :: filename
    CMPLX   :: projection
    CMPLX, allocatable :: state1(:,:), state2(:,:), orbital(:,:)
    type(submesh_t) :: submesh
    FLOAT, allocatable ::  rr(:,:), ylm(:)

    ASSERT(st%d%kpt%start==1 .and. st%d%kpt%end==sys%gr%sb%kpoints%full%npoints)

    ! precompute orbitals
    SAFE_ALLOCATE(orbital(w90_nproj,1:sys%gr%mesh%np))
    orbital = M_ZERO
    do iw=1,w90_nproj
      ! cartesian coordinate of orbital center
      center(1:3) =  matmul(w90_proj_centers(iw,1:3),sys%gr%sb%rlattice(1:3,1:3))
      call submesh_init(submesh, sys%gr%sb, sys%gr%mesh, center, CNST(2.0)*maxval(sys%gr%sb%lsize(:)))

      ! make transpose table of submesh points for use in pwscf routine
      SAFE_ALLOCATE(rr(1:3,submesh%np))
      do ip=1,submesh%np
         rr(1:3,ip) = submesh%x(ip,1:3)
      end do

      ! get ylm as submesh points
      SAFE_ALLOCATE(ylm(1:submesh%np))
      ! (this is a routine from pwscf)
      call ylm_wannier(ylm,w90_proj_lmr(iw,1),w90_proj_lmr(iw,2),rr,submesh%np)

      ! apply radial function
      do ip=1,submesh%np
         dd=sqrt(dot_product(submesh%x(ip,1:3),submesh%x(ip,1:3)))
         ylm(ip) = ylm(ip)*M_TWO*exp(-dd)
      end do

      call submesh_add_to_mesh(submesh,ylm, orbital(iw, :))

      SAFE_DEALLOCATE_A(ylm)
      SAFE_DEALLOCATE_A(rr)
    end do

    filename = './'// trim(adjustl(w90_prefix))//'.amn'
    w90_amn = io_open(trim(filename), action='write')

    ! write header
    if(mpi_grp_is_root(mpi_world)) then
       write(w90_amn,*) ' '
       write(w90_amn,*)  w90_num_bands, w90_num_kpts, w90_num_wann
    end if

    SAFE_ALLOCATE(state1(1:sys%gr%der%mesh%np, 1:st%d%dim))

    do ist=1,w90_num_bands
       do iw=1,w90_num_wann
          do ik=1, w90_num_kpts
             call states_get_state(st, sys%gr%der%mesh, ist, ik, state1)
             projection = M_ZERO
             do idim=1,st%d%dim
                projection = projection + zmf_dotp(sys%gr%mesh,state1(1:sys%gr%mesh%np,idim),orbital(iw,1:sys%gr%mesh%np))
             end do
             write (w90_amn,'(I5,2x,I5,2x,I5,2x,e12.6,2x,e12.6)') ist, iw, ik, projection
          end do
       end do
    end do
    call io_close(w90_amn)

    SAFE_DEALLOCATE_A(state1)
    SAFE_DEALLOCATE_A(orbital)

  end subroutine create_wannier90_amn


  ! routine originally from Pwscf
  subroutine ylm_wannier(ylm,l,mr,r,nr)
  ! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr)
  ! of the spherical harmonic identified  by indices (l,mr)
  ! in table 3.1 of the wannierf90 specification.
  !
  ! No reference to the particular ylm ordering internal to octopus
  ! is assumed.
  !
  ! If ordering in wannier90 code is changed or extended this should be the
  ! only place to be modified accordingly
  !
  
     IMPLICIT NONE
  ! I/O variables
  !
     INTEGER :: l, mr, nr
     FLOAT :: ylm(nr), r(3,nr)
  !
  ! local variables
  !
  !     FLOAT, EXTERNAL ::  s, p_z,px,py, dz2, dxz, dyz, dx2my2, dxy 
  !     FLOAT, EXTERNAL :: fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
     FLOAT :: rr, cost, phi
     INTEGER :: ir
     FLOAT :: bs2, bs3, bs6, bs12
     bs2 = 1.d0/sqrt(2.d0)
     bs3=1.d0/sqrt(3.d0)
     bs6 = 1.d0/sqrt(6.d0)
     bs12 = 1.d0/sqrt(12.d0)
  
     IF (l > 3 .or. l < -5 ) stop 'ylm_wannier l out of range '
     IF (l>=0) THEN
        IF (mr < 1 .or. mr > 2*l+1) stop 'ylm_wannier mr out of range'
     ELSE
        IF (mr < 1 .or. mr > abs(l)+1 ) stop 'ylm_wannier mr out of range'
     ENDIF
  
     DO ir=1, nr
        rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
        IF (rr < M_EPSILON) then !stop 'ylm_wannier rr too small'
           ylm(ir)=0
           cycle
        end IF
        cost =  r(3,ir) / rr
        !
        !  beware the arc tan, it is defined modulo M_PI
        !
        IF (r(1,ir) > M_EPSILON) THEN
           phi = atan( r(2,ir)/r(1,ir) )
        ELSEIF (r(1,ir) < -M_EPSILON ) THEN
           phi = atan( r(2,ir)/r(1,ir) ) + M_PI
        ELSE
           phi = sign( M_PI/2.d0,r(2,ir) )
        ENDIF
  
  
        IF (l==0) THEN   ! s orbital
           ylm(ir) = s(cost,phi)
        ENDIF
        IF (l==1) THEN   ! p orbitals
           IF (mr==1) ylm(ir) = p_z(cost,phi)
           IF (mr==2) ylm(ir) = px(cost,phi)
           IF (mr==3) ylm(ir) = py(cost,phi)
        ENDIF
        IF (l==2) THEN   ! d orbitals
           IF (mr==1) ylm(ir) = dz2(cost,phi)
           IF (mr==2) ylm(ir) = dxz(cost,phi)
           IF (mr==3) ylm(ir) = dyz(cost,phi)
           IF (mr==4) ylm(ir) = dx2my2(cost,phi)
           IF (mr==5) ylm(ir) = dxy(cost,phi)
        ENDIF
        IF (l==3) THEN   ! f orbitals
           IF (mr==1) ylm(ir) = fz3(cost,phi)
           IF (mr==2) ylm(ir) = fxz2(cost,phi)
           IF (mr==3) ylm(ir) = fyz2(cost,phi)
           IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
           IF (mr==5) ylm(ir) = fxyz(cost,phi)
           IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
           IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
        ENDIF
        IF (l==-1) THEN  !  sp hybrids
           IF (mr==1) ylm(ir) = bs2 * ( s(cost,phi) + px(cost,phi) )
           IF (mr==2) ylm(ir) = bs2 * ( s(cost,phi) - px(cost,phi) )
        ENDIF
        IF (l==-2) THEN  !  sp2 hybrids
           IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
           IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
           IF (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
        ENDIF
        IF (l==-3) THEN  !  sp3 hybrids
           IF (mr==1) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)+py(cost,phi)+p_z(cost,phi))
           IF (mr==2) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)-py(cost,phi)-p_z(cost,phi))
           IF (mr==3) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)+py(cost,phi)-p_z(cost,phi))
           IF (mr==4) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)-py(cost,phi)+p_z(cost,phi))
        ENDIF
        IF (l==-4) THEN  !  sp3d hybrids
           IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
           IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
           IF (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
           IF (mr==4) ylm(ir) = bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
           IF (mr==5) ylm(ir) =-bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
        ENDIF
        IF (l==-5) THEN  ! sp3d2 hybrids
           IF (mr==1) ylm(ir) = bs6*s(cost,phi)-bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
           IF (mr==2) ylm(ir) = bs6*s(cost,phi)+bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
           IF (mr==3) ylm(ir) = bs6*s(cost,phi)-bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
           IF (mr==4) ylm(ir) = bs6*s(cost,phi)+bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
           IF (mr==5) ylm(ir) = bs6*s(cost,phi)-bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
           IF (mr==6) ylm(ir) = bs6*s(cost,phi)+bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
        ENDIF
  
     ENDDO
  
     RETURN
   end subroutine ylm_wannier
  
    !======== l = 0 =====================================================================
     FUNCTION s(cost,phi)
       IMPLICIT NONE
       FLOAT :: s, cost,phi
       s = 1.d0/ sqrt((M_FOUR*M_PI))
       RETURN
    END FUNCTION s
    !======== l = 1 =====================================================================
    FUNCTION p_z(cost,phi)
       IMPLICIT NONE
       FLOAT ::p_z, cost,phi
       p_z =  sqrt(3.d0/(M_FOUR*M_PI)) * cost
       RETURN
    END FUNCTION p_z
    FUNCTION px(cost,phi)
       IMPLICIT NONE
       FLOAT ::px, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       px =  sqrt(3.d0/(M_FOUR*M_PI)) * sint * cos(phi)
       RETURN
    END FUNCTION px
    FUNCTION py(cost,phi)
       IMPLICIT NONE
       FLOAT ::py, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       py =  sqrt(3.d0/(M_FOUR*M_PI)) * sint * sin(phi)
       RETURN
    END FUNCTION py
    !======== l = 2 =====================================================================
    FUNCTION dz2(cost,phi)
       IMPLICIT NONE
       FLOAT ::dz2, cost, phi
       dz2 =  sqrt(1.25d0/(M_FOUR*M_PI)) * (3.d0* cost*cost-1.d0)
       RETURN
    END FUNCTION dz2
    FUNCTION dxz(cost,phi)
       IMPLICIT NONE
       FLOAT ::dxz, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       dxz =  sqrt(15.d0/(M_FOUR*M_PI)) * sint*cost * cos(phi)
       RETURN
    END FUNCTION dxz
    FUNCTION dyz(cost,phi)
       IMPLICIT NONE
       FLOAT ::dyz, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       dyz =  sqrt(15.d0/(M_FOUR*M_PI)) * sint*cost * sin(phi)
       RETURN
    END FUNCTION dyz
    FUNCTION dx2my2(cost,phi)
       IMPLICIT NONE
       FLOAT ::dx2my2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       dx2my2 =  sqrt(3.75d0/(M_FOUR*M_PI)) * sint*sint * cos(2.d0*phi)
       RETURN
    END FUNCTION dx2my2
    FUNCTION dxy(cost,phi)
       IMPLICIT NONE
       FLOAT ::dxy, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       dxy =  sqrt(3.75d0/(M_FOUR*M_PI)) * sint*sint * sin(2.d0*phi)
       RETURN
    END FUNCTION dxy
    !======== l = 3 =====================================================================
    FUNCTION fz3(cost,phi)
       IMPLICIT NONE
       FLOAT ::fz3, cost, phi
       fz3 =  0.25d0*sqrt(7.d0/M_PI) * ( 5.d0 * cost * cost - 3.d0 ) * cost
       RETURN
    END FUNCTION fz3
    FUNCTION fxz2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fxz2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fxz2 =  0.25d0*sqrt(10.5d0/M_PI) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
       RETURN
    END FUNCTION fxz2
    FUNCTION fyz2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fyz2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fyz2 =  0.25d0*sqrt(10.5d0/M_PI) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
       RETURN
    END FUNCTION fyz2
    FUNCTION fzx2my2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fzx2my2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fzx2my2 =  0.25d0*sqrt(105d0/M_PI) * sint * sint * cost * cos(2.d0*phi)
       RETURN
    END FUNCTION fzx2my2
    FUNCTION fxyz(cost,phi)
       IMPLICIT NONE
       FLOAT ::fxyz, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fxyz =  0.25d0*sqrt(105d0/M_PI) * sint * sint * cost * sin(2.d0*phi)
       RETURN
    END FUNCTION fxyz
    FUNCTION fxx2m3y2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fxx2m3y2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fxx2m3y2 =  0.25d0*sqrt(17.5d0/M_PI) * sint * sint * sint * cos(3.d0*phi)
       RETURN
    END FUNCTION fxx2m3y2
    FUNCTION fy3x2my2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fy3x2my2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fy3x2my2 =  0.25d0*sqrt(17.5d0/M_PI) * sint * sint * sint * sin(3.d0*phi)
       RETURN
    END FUNCTION fy3x2my2
  
end program wannier90_interface

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
