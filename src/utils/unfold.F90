!!  Copyright (C) 20182019 M. S. Mrudul, N. Tancogne-Dejean
!!  
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


program oct_unfold
  use batch_oct_m
  use batch_ops_oct_m
  use calc_mode_par_oct_m
  use cube_oct_m
  use comm_oct_m
  use command_line_oct_m
  use fft_oct_m
  use fftw_params_oct_m
  use fourier_space_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use multicomm_oct_m
  use cube_function_oct_m
  use fourier_shell_oct_m
  use io_oct_m
  use io_binary_oct_m
  use io_function_oct_m
  use loct_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use parser_oct_m
  use poisson_fft_oct_m
  use profiling_oct_m
  use restart_oct_m
  use space_oct_m
  use simul_box_oct_m
  use states_oct_m
  use states_restart_oct_m
  use states_dim_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use xc_oct_m

  implicit none

  type(system_t) :: sys
  type(simul_box_t) :: sb
  integer ::  ik, idim, nkpoints
  type(restart_t) :: restart
  type(cube_t) :: zcube
  type(cube_function_t) :: cf

  integer :: ierr, run_mode, file_Gvec, jdim
  FLOAT :: lparams(MAX_DIM), rlattice_PC(MAX_DIM, MAX_DIM), klattice_PC(MAX_DIM, MAX_DIM)
  FLOAT :: volume_element_PC
  type(block_t) :: blk
  integer :: nhighsympoints, nsegments
  integer :: icol, idir, ncols
  integer, allocatable :: resolution(:)
  FLOAT, allocatable   :: highsympoints(:,:), coord_along_path(:)
  type(kpoints_grid_t) :: path_kpoints_grid

  
  ! the usual initializations
  call global_init(is_serial = .false.)
  call calc_mode_par_init()

  call messages_init()

  call io_init()
  call profiling_init()

  call print_header()
  call messages_print_stress(stdout, "Unfolding Band-structure")
  call messages_print_stress(stdout)

  call messages_experimental("oct-unfold utility")
  call fft_all_init()
  call unit_system_init()
  call restart_module_init()

  call calc_mode_par_set_parallelization(P_STRATEGY_STATES, default = .false.)
  call system_init(sys)
  call simul_box_init(sb, sys%geo, sys%space)

  if(sb%periodic_dim == 0) then
    message(1) = "oct-unfold can only be used for periodic ystems."
    call messages_fatal(1)
  end if

  if(sys%st%parallel_in_states) then
    call messages_not_implemented("oct-unfold with states parallelization.")
  end if 

  if(sys%st%d%ispin == SPINORS) then
    call messages_not_implemented("oct-unfold for spinors")
  end if


  !%Variable UnfoldMode
  !%Type flag
  !%Default none
  !%Section Utilities::oct-unfold
  !%Description
  !% Specifies which stage of the unfolding tool to use
  !%Option unfold_setup bit(1)
  !% Writes the list of k-points corresponding to the path specified by <tt>UnfoldKPointPath</tt>.
  !% This list of k-point (unfold_kpt.dat) must be used for an unocc calculation of the supercell,
  !% adding the line "include 'unfold_kpt.dat'" to the inp file and removing the KPointGrid information.
  !%Option unfold_run bit(2)
  !% Perform the actual unfolding, based on the states obtained from the previous unocc run.
  !%End
  call parse_variable('UnfoldMode', 0, run_mode)
  if(run_mode == 0) then
    message(1) = "UnfoldMode must be set to a value different from 0."
    call messages_fatal(1)
  end if

  !%Variable UnfoldLatticeParameters
  !%Type block
  !%Default 1 | 1 | 1
  !%Section Utilities:oct-unfold
  !%Description
  !% The lattice parameters of the primitive cell, on which unfolding is performed. 
  !%End
  lparams(:) = M_ONE
  if (parse_block('UnfoldLatticeParameters', blk) == 0) then
    do idim = 1, sb%dim
      call parse_block_float(blk, 0, idim - 1, lparams(idim))
    end do
  else
    message(1) = "UnfoldLatticeParameters is not specified"
    call messages_fatal(1)
  end if

  !%Variable UnfoldLatticeVectors
  !%Type block
  !%Default simple cubic
  !%Section Utilites:oct-unfold
  !%Description
  !% Lattice vectors of the primitive cell on which the unfolding is performed. 
  !%End
  rlattice_PC = M_ZERO
  forall(idim = 1:sb%dim) rlattice_PC(idim, idim) = M_ONE

  if (parse_block('UnfoldLatticeVectors', blk) == 0) then
    do idim = 1, sb%dim
      do jdim = 1, sb%dim
        call parse_block_float(blk, idim - 1,  jdim - 1, rlattice_PC(jdim, idim))
      enddo
    end do
    call parse_block_end(blk)
  else
    message(1) = "UnfoldLatticeVectors is not specified"
    call messages_fatal(1)
  end if

  do idim = 1, sb%dim
    forall(jdim = 1:sb%dim)
      rlattice_PC(jdim, idim) = rlattice_PC(jdim, idim) * lparams(idim)
    end forall
  end do

  call reciprocal_lattice(rlattice_PC, klattice_PC, volume_element_PC, sb%dim)
  klattice_PC = klattice_PC * M_TWO*M_PI


  !%Variable UnfoldKPointsPath
  !%Type block
  !%Section Utilites:oct-unfold
  !%Description
  !% Specifies the k-point path for which the unfolding need to be done.
  !% The syntax is identical to <tt>KPointsPath</tt>.
  !%End
  if(parse_block('UnfoldKPointsPath', blk) /= 0) then
    write(message(1),'(a)') 'Internal error while reading UnfoldPointsPath.'
    call messages_fatal(1)
  end if

  ! There is one high symmetry k-point per line
  nsegments = parse_block_cols(blk, 0)
  nhighsympoints = parse_block_n(blk)-1
  if( nhighsympoints /= nsegments +1) then
    write(message(1),'(a,i3,a,i3)') 'The first row of UnfoldPointsPath is not compatible with the number of specified k-points.'
    call messages_fatal(1)
  end if

  SAFE_ALLOCATE(resolution(1:nsegments))
  do icol = 1, nsegments
    call parse_block_integer(blk, 0, icol-1, resolution(icol))
  end do
  !Total number of points in the segment
  nkpoints = sum(resolution)+1

  SAFE_ALLOCATE(highsympoints(1:sb%dim, 1:nhighsympoints))
  do ik = 1, nhighsympoints
    !Sanity check
    ncols = parse_block_cols(blk, ik)
    if(ncols /= sb%dim) then
      write(message(1),'(a,i3,a,i3)') 'UnfoldPointsPath row ', ik, ' has ', ncols, ' columns but must have ', sb%dim
      call messages_fatal(1)
    end if

    do idir = 1, sb%dim
        call parse_block_float(blk, ik, idir-1, highsympoints(idir, ik))
    end do
  end do

  call kpoints_grid_init(sb%dim, path_kpoints_grid, nkpoints, 1)
  ! For the output of band-structures
  SAFE_ALLOCATE(coord_along_path(1:nkpoints))

  call kpoints_path_generate(sb%dim, klattice_PC, nkpoints, nsegments, resolution, &
           nhighsympoints, highsympoints, path_kpoints_grid%point, coord_along_path)

  SAFE_DEALLOCATE_A(resolution)
  SAFE_DEALLOCATE_A(highsympoints)

  !We convert the k-point to the reduce coordinate of the supercell
  do ik = 1, path_kpoints_grid%npoints
    call kpoints_to_reduced(sb%rlattice, path_kpoints_grid%point(:, ik), &
                                path_kpoints_grid%red_point(:, ik), sb%dim)
  end do

  if(iand(run_mode, OPTION__UNFOLDMODE__UNFOLD_SETUP) /= 0) then

    call unfold_setup()
    
  else if(iand(run_mode, OPTION__UNFOLDMODE__UNFOLD_RUN) /= 0) then

    !Sanity check
    file_Gvec = io_open('unfold_gvec.dat', action='read')
    read(file_Gvec, *)
    read(file_Gvec, *) ik
    if(ik /= path_kpoints_grid%npoints) then
      message(1) = 'There is an inconsistency between unfold_gvec.dat and the input file'
      call messages_fatal(1)
    end if
    call io_close(file_Gvec)
 
    call states_allocate_wfns(sys%st, sys%gr%mesh)

    call restart_init(restart, RESTART_UNOCC, RESTART_TYPE_LOAD, sys%mc, ierr, mesh=sys%gr%mesh, exact=.true.)
    if(ierr == 0) call states_load(restart, sys%st, sys%gr, ierr, label = ": unfold")
    if (ierr /= 0) then
      message(1) = 'Unable to read unocc wavefunctions.'
      call messages_fatal(1)
    end if
    call restart_end(restart)  

    call cube_init(zcube, sys%gr%mesh%idx%ll, sb, fft_type = FFT_COMPLEX, dont_optimize = .true.)
    call cube_function_null(cf)
    call zcube_function_alloc_rs(zcube, cf)
    call cube_function_alloc_fs(zcube, cf)

    call wfs_extract_spec_fn(sys%st, sys%gr, zcube, cf)

    call cube_function_free_fs(zcube, cf)
    call zcube_function_free_rs(zcube, cf)

  else
    message(1) = "Unsupported or incorrect value of UnfoldMode." 
    call messages_fatal(1)
  end if

  SAFE_DEALLOCATE_A(coord_along_path)

  call kpoints_grid_end(path_kpoints_grid)

  call simul_box_end(sb)
  call fft_all_end()
  call system_end(sys)
  call profiling_output()
  call profiling_end()
  call io_end()
  call print_date("Calculation ended on ")
  call messages_end()
  call global_end()



contains
  !-----------------------------------------------------------------
  subroutine unfold_setup()
    integer :: file_Gvec, file_kpts
    integer :: Gvec(1:MAX_DIM)

    PUSH_SUB(unfold_setup)

    if(mpi_grp_is_root(mpi_world)) then
      file_Gvec = io_open('unfold_gvec.dat', action='write')
      file_kpts = io_open('unfold_kpt.dat', action='write')
      write(file_kpts,'(a)') '%KpointsReduced'
      write(file_Gvec,*) '#Created by oct-unfold'
      write(file_Gvec,*) path_kpoints_grid%npoints

      !We convert the k-point to the reduce coordinate of the supercell
      do ik = 1, path_kpoints_grid%npoints
        Gvec(1:sb%dim) = nint(path_kpoints_grid%red_point(:, ik)+M_HALF*CNST(1e-7))
        write(file_kpts,'(a6,f12.8,a3,f12.8,a3,f12.8)') ' 1. | ',  path_kpoints_grid%red_point(1, ik)-Gvec(1), &
                  ' | ', path_kpoints_grid%red_point(2, ik)-Gvec(2), ' | ', path_kpoints_grid%red_point(3, ik)-Gvec(3)
        write(file_Gvec,'(3i)') Gvec(1:3)
      end do
      write(file_kpts,'(a)') '%'
      call io_close(file_Gvec)
      call io_close(file_kpts)
    end if

    POP_SUB(unfold_setup)
  end subroutine unfold_setup

  !--------------------------------------------------------------------
  subroutine wfs_extract_spec_fn(st, gr, zcube, cf)
    type(states_t),        intent(in)    :: st
    type(grid_t),          intent(in)    :: gr
    type(cube_t),          intent(inout) :: zcube
    type(cube_function_t), intent(inout) :: cf

    FLOAT, allocatable :: PKm(:,:), AkE(:,:), eigs(:)
    CMPLX, allocatable :: zpsi(:), field_g(:)
    integer :: file_ake, iq, ist, idim, nenergy 
    integer :: ig, ix, iy, iz, ik, ie, Gmin, Gmax
    FLOAT   :: eigmin, eigmax, de, norm, tol=1e-7
    integer, parameter :: nextend = 10
    FLOAT :: vec_PC(MAX_DIM),vec_SC(MAX_DIM)
    type(fourier_shell_t) :: shell 
    character(len=100) :: filename
    FLOAT, allocatable :: gvec_abs(:,:)
    logical, allocatable :: g_select(:)

    PUSH_SUB(wfs_extract_spec_fn)
   
    SAFE_ALLOCATE(zpsi(1:gr%mesh%np))

    !%Variable UnfoldEnergyStep
    !%Type float
    !%Default 0
    !%Section Utilities::oct-unfold
    !%Description
    !% Specifies the energy resolution for the unfolded band structure.
    !% If you specify 0, the resolution will be set to be 1/1000 points between <tt>UnfoldMinEnergy</tt>
    !% and <tt>UnfoldMaxEnergy</tt> 
    !%End
    call parse_variable('UnfoldEnergyStep', M_ZERO, de)
    if(de < M_ZERO) then
      message(1) = "UnfoldEnergyStep must be positive"
      call messages_fatal(1)
    end if

    !%Variable UnfoldMinEnergy
    !%Type float
    !%Section Utilities::oct-unfold
    !%Description
    !% Specifies the start of the energy range for the unfolded band structure.
    !% The default value correspond to the samllest eigenvalue.
    !%End
    call parse_variable('UnfoldMinEnergy', minval(st%eigenval(:,:)), eigmin)

    !%Variable UnfoldMaxEnergy
    !%Type float
    !%Section Utilities::oct-unfold
    !%Description
    !% Specifies the end of the energy range for the unfolded band structure.
    !% The default value correspond to the largest eigenvalue.
    !%End
    call parse_variable('UnfoldMinEnergy', maxval(st%eigenval(:,:)), eigmax)
 
    if(de == M_ZERO) then
      de = (eigmax-eigmin)/CNST(1000)
    end if
    
    !We increase a bit the energy range
    nenergy = floor((eigmax-eigmin+2*nextend*de)/de)
    SAFE_ALLOCATE(eigs(1:nenergy))
    do ie = 1, nenergy
      eigs(ie) = eigmin - nextend*de + (ie-1)*de
    end do        

    PKm(:,:) = M_ZERO

    SAFE_ALLOCATE(gvec_abs(sb%periodic_dim, st%d%nik))
    file_Gvec = io_open('./unfold_gvec.dat', action='read')
    read(file_Gvec,*)
    read(file_Gvec,*)
    do ik = 1, st%d%nik
      read(file_Gvec,*) vec_SC
      call kpoints_to_absolute(sb%klattice, vec_SC(1:sb%periodic_dim), gvec_abs(1:sb%periodic_dim,ik), &
                 sb%periodic_dim)
    end do 
    call io_close(file_Gvec)

    if (mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, (st%d%kpt%end-st%d%kpt%start+1)*st%nst)

    SAFE_ALLOCATE(AkE(1:nenergy, 1:st%d%nik))
    AkE(:,:) = M_ZERO

    SAFE_ALLOCATE(PKm(st%d%kpt%start:st%d%kpt%end,1:st%nst))
    do ik = st%d%kpt%start, st%d%kpt%end
      iq =  states_dim_get_kpoint_index(st%d, ik) 

      call fourier_shell_init(shell, zcube, gr%mesh, kk = sb%kpoints%reduced%red_point(:, iq))  

      Gmin = minval(shell%red_gvec(:,:))
      Gmax = maxval(shell%red_gvec(:,:))

      SAFE_ALLOCATE(g_select(1:shell%ngvectors))
      g_select(:) = .false.

      select case(sb%periodic_dim)
      case(3)
        do ig = 1, shell%ngvectors
          call kpoints_to_absolute(sb%klattice, real(shell%red_gvec(1:3,ig), REAL_PRECISION), vec_SC, 3)
          do ix = Gmin,Gmax
            do iy = Gmin,Gmax
              do iz = Gmin,Gmax
                vec_PC = ix*klattice_PC(1:3,1) + iy*klattice_PC(1:3,2) + iz*klattice_PC(1:3,3)
                if(abs(vec_SC(1)-vec_PC(1)-gvec_abs(1, ik)) < tol &
                  .and. abs(vec_SC(2)-vec_PC(2)-gvec_abs(2, ik)) < tol &
                  .and. abs(vec_SC(3)-vec_PC(3)-gvec_abs(3, ik)) < tol) then
                  g_select(ig) = .true.
                end if
              end do !iz
            end do !iy
          end do !ix
        end do !ig
           
      case(2)

        do ig = 1, shell%ngvectors
          call kpoints_to_absolute(sb%klattice, real(shell%red_gvec(1:2,ig), REAL_PRECISION), vec_SC, 2)
          do ix = Gmin,Gmax
            do iy = Gmin,Gmax
              vec_PC = ix*klattice_PC(1:3,1) + iy*klattice_PC(1:3,2)
              if(abs(vec_SC(1)-vec_PC(1)-gvec_abs(1, ik)) < tol &
                  .and. abs(vec_SC(2)-vec_PC(2)-gvec_abs(2, ik)) < tol) then
                  g_select(ig) = .true.
              end if
            end do !ix
          end do !ix
        end do !ig

      case default
        call messages_not_implemented("Unfolding for dimensions other than 2 or 3") 
      end select

      if (mpi_grp_is_root(gr%mesh%mpi_grp)) then
        write(filename,"(A13,I3.3,A4)") "./static/AkE_",ik,".dat"
        file_ake = io_open(trim(filename), action='write')
        write(file_ake, *) '#Energy Ak(E)'
        write(file_ake, *) '#Number of points in energy window ',  nenergy 
      end if

      do ist = 1, st%nst 
        !loop over states
        do idim = 1, st%d%dim
          ! Getting wavefunctions 
          ! for the moment we treat all functions as complex
          call states_get_state(st, gr%mesh, idim, ist, ik, zpsi)
            
          if(gr%mesh%parallel_in_domains) then
            call zmesh_to_cube(gr%mesh, zpsi, zcube, cf, local = .true.)
          else
            call zmesh_to_cube(gr%mesh, zpsi, zcube, cf)
          end if

          !Fourier transform from real-space to fourier space
          call zcube_function_rs2fs(zcube, cf)

          ! Normalisation
          SAFE_ALLOCATE(field_g(1:shell%ngvectors))
          norm = M_ZERO
          do ig = 1, shell%ngvectors
            field_g(ig) = cf%fs(shell%coords(1, ig), shell%coords(2, ig), shell%coords(3, ig)) * &
              sqrt(gr%mesh%volume_element / product(zcube%rs_n_global(1:3)))
            norm = norm + abs(field_g(ig))**2
          end do
          field_g(:) = field_g(:) / sqrt(norm)

          !Finding sub-g and calculating the Projection Pkm
          do ig = 1, shell%ngvectors
            if(.not. g_select(ig)) cycle
            PKm(ik,ist) = PKm(ik,ist) + abs(field_g(ig))**2
          end do 

          SAFE_DEALLOCATE_A(field_g)
        end do !idim

        !ist loop end
        if (mpi_grp_is_root(mpi_world)) call loct_progress_bar((ik-st%d%kpt%start)*st%nst+ist, &
                                               (st%d%kpt%end-st%d%kpt%start+1)*st%nst)
      end do !ist

      ! Calculating the spectral function 
      !TODO: We could implement here a different broadening
      do ist = 1, st%nst
        do ie = 1, nenergy
           AkE(ie, ik)  = AkE(ie, ik) + PKm(ik,ist) * ( 3.0*de / M_PI ) / & 
                     ( (eigs(ie)-st%eigenval(ist,ik))**2 + (3.0*de)**2 ) 
        end do   
      end do 
        
      ! writing Spectral-Function
      do ie = 1, nenergy
        write(file_ake,fmt ='(1es19.12,1x,1es19.12)') eigs(ie), AkE(ie, ik)
      end do

      if (mpi_grp_is_root(gr%mesh%mpi_grp)) call io_close(file_ake) 

      call fourier_shell_end(shell)
      SAFE_DEALLOCATE_A(g_select)
       
    end do !ik

#if defined(HAVE_MPI)        
    if(st%d%kpt%parallel) then
      call comm_allreduce(st%st_kpt_mpi_grp%comm, AkE)
    end if
#endif  

    if(mpi_grp_is_root(mpi_world)) then
      file_ake = io_open("static/AkE.dat", action='write')
      write(file_ake, *) '#Energy Ak(E)'
      write(file_ake, *) '#Number of points in energy window ',  nenergy 
      do ik = 1, st%d%nik
        do ie = 1, nenergy
          write(file_ake,fmt ='(1es19.12,1x,1es19.12,1x,1es19.12)') coord_along_path(ik), &
                   eigs(ie), AkE(ie, ik)
        end do
      end do
      
      call io_close(file_ake)
    end if

    SAFE_DEALLOCATE_A(eigs)
    SAFE_DEALLOCATE_A(AkE)

    SAFE_DEALLOCATE_A(gvec_abs)
    SAFE_DEALLOCATE_A(PKm)
    SAFE_DEALLOCATE_A(zpsi)
    POP_SUB(wfs_extract_spec_fn)

  end subroutine wfs_extract_spec_fn
  !------------------------------------------------------------

end program oct_unfold
!! Local Variables:
!! mode: f90                              
!! coding: utf-8 
!! End:

