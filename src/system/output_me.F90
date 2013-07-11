!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any %later version.
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
!! $Id: unocc.F90 5320 2009-04-23 15:29:57Z marques $

#include "global.h"

module output_me_m
  use datasets_m
  use derivatives_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use kpoints_m
  use loct_math_m
  use mpi_m
  use mpi_lib_m
  use parser_m
  use poisson_m
  use profiling_m
  use projector_m
  use simul_box_m
  use states_m
  use states_calc_m
  use states_dim_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m

  implicit none

  private
  public ::           &
    output_me_t,      &
    output_me_init,   &
    output_me

  type output_me_t
    integer :: what                !< what to output 
    !> If output_ksdipole, this number sets up which matrix elements will
    !! be printed: e.g. if ksmultipoles = 3, the dipole, quadrupole and 
    !! octopole matrix elements (between Kohn-Sham or single-particle orbitals).
    integer :: ks_multipoles      
  end type output_me_t

  integer, parameter, public :: &
    OUTPUT_ME_MOMENTUM       =   1, &
    OUTPUT_ME_ANG_MOMENTUM   =   2, &
    OUTPUT_ME_ONE_BODY       =   4, &
    OUTPUT_ME_TWO_BODY       =   8, &
    OUTPUT_ME_KS_MULTIPOLES  =  16

contains
  
  ! ---------------------------------------------------------
  subroutine output_me_init(this, sb)
    type(output_me_t), intent(out) :: this
    type(simul_box_t), intent(in)  :: sb

    PUSH_SUB(output_me_init)

    !%Variable OutputMatrixElements
    !%Type flag
    !%Default no
    !%Section Output
    !%Description
    !% Specifies what matrix elements to print.
    !% Enabled only if <tt>Output</tt> includes <tt>matrix_elements</tt>.
    !% The output files go into the <tt>static</tt> directory, except when
    !% running a time-dependent simulation, when the directory <tt>td.XXXXXXX</tt> is used.
    !% Example: "momentum + ks_multipoles"
    !%Option momentum 1
    !% Momentum. Filename: <tt>ks_me_momentum</tt>.
    !%Option ang_momentum 2
    !% Dimensionless angular momentum (r x k). Filename: <tt>ks_me_angular_momentum</tt>.
    !%Option one_body 4
    !% <math>&lt;i|T + V_{ext}|j&gt;</math>. Not available with states parallelization.
    !%Option two_body 8
    !% <math>&lt;ij| 1/|r_1-r_2| |kl&gt;</math>. Not available with states parallelization.
    !%Option ks_multipoles 16
    !% See <tt>OutputMEMultipoles</tt>. Not available with states parallelization.
    !%End

    call parse_integer(datasets_check('OutputMatrixElements'), 0, this%what)
    if(.not.varinfo_valid_option('OutputMatrixElements', this%what, is_flag=.true.)) then
      call input_error('OutputMatrixElements')
    end if

    if(sb%dim /= 2 .and. sb%dim /= 3) this%what = iand(this%what, not(OUTPUT_ME_ANG_MOMENTUM))

    if(iand(this%what, OUTPUT_ME_KS_MULTIPOLES) /= 0) then
      !%Variable OutputMEMultipoles
      !%Type integer
      !%Default 1
      !%Section Output
      !%Description
      !% This variable decides which multipole moments are printed out for
      !% <tt>OutputMatrixElements = ks_multipoles</tt>: <i>e.g.</i>, if 1, then the
      !% program will print three files, <tt>ks_me_multipoles.x</tt> (<tt>x</tt>=1,2,3), containing
      !% respectively the (1,-1), (1,0) and (1,1) multipole matrix elements
      !% between Kohn-Sham states.
      !%End
      call parse_integer(datasets_check('OutputMEMultipoles'), 1, this%ks_multipoles)
    end if

    POP_SUB(output_me_init)
  end subroutine output_me_init


  ! ---------------------------------------------------------
  subroutine output_me(this, dir, st, gr, geo, hm)
    type(output_me_t),   intent(in)    :: this
    character(len=*),    intent(in)    :: dir
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(in)    :: geo
    type(hamiltonian_t), intent(in)    :: hm

    integer :: id, ll, mm, ik
    character(len=256) :: fname
    
    PUSH_SUB(output_me)

    if(iand(this%what, output_me_momentum) /= 0) then
      write(fname,'(2a)') trim(dir), '/ks_me_momentum'
      call output_me_out_momentum(fname, st, gr)
    end if

    if(iand(this%what, output_me_ang_momentum) /= 0) then
      write(fname,'(2a)') trim(dir), '/ks_me_angular_momentum'
      call output_me_out_ang_momentum(fname, st, gr)
    end if

    if(iand(this%what, output_me_ks_multipoles) /= 0) then
      ! The content of each file should be clear from the header of each file.
      id = 1
      do ik = 1, st%d%nik
        do ll = 1, this%ks_multipoles
          do mm = -ll, ll
            write(fname,'(i4)') id
            write(fname,'(a)') trim(dir)//'/ks_me_multipoles.'//trim(adjustl(fname))
            if (states_are_real(st)) then
              call doutput_me_ks_multipoles(fname, st, gr, ll, mm, ik)
            else
              call zoutput_me_ks_multipoles(fname, st, gr, ll, mm, ik)
            end if

            id = id + 1
          end do
        end do
      end do
    end if

    if(iand(this%what, output_me_one_body) /= 0) then
      message(1) = "Computing one-body matrix elements"
      call messages_info(1)
      if (states_are_real(st)) then
        call done_body(dir, gr, geo, st, hm)
      else
        call zone_body(dir, gr, geo, st, hm)
      end if
    end if

    if(iand(this%what, output_me_two_body) /= 0) then
      message(1) = "Computing two-body matrix elements"
      call messages_info(1)
      if (states_are_real(st)) then
        call dtwo_body(dir, gr, st)
      else
        call ztwo_body(dir, gr, st)
      end if
    end if

    POP_SUB(output_me)
  end subroutine output_me


  ! ---------------------------------------------------------
  subroutine output_me_out_momentum(fname, st, gr)
    character(len=*), intent(in) :: fname
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr

    integer            :: ik, ist, is, ns, iunit, idir
    character(len=80)  :: cspin, str_tmp
    FLOAT              :: kpoint(1:MAX_DIM)
    FLOAT, allocatable :: momentum(:,:,:)

    PUSH_SUB(output_me_out_momentum)

    SAFE_ALLOCATE(momentum(1:gr%sb%dim, 1:st%nst, 1:st%d%nik))

    call states_calc_momentum(st, gr%der, momentum)

    iunit = io_open(fname, action='write')

    ns = 1
    if(st%d%nspin == 2) ns = 2

    write(message(1),'(a)') 'Momentum of the KS states [a.u.]:'
    call messages_info(1, iunit)      
    if (st%d%nik > ns) then
      message(1) = 'k-points [' // trim(units_abbrev(unit_one/units_out%length)) // ']'
      call messages_info(1, iunit)
    end if

    do ik = 1, st%d%nik, ns
      kpoint = M_ZERO
      kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))

      if(st%d%nik > ns) then
        write(message(1), '(a,i4, a)') '#k =', ik, ', k = ('
        do idir = 1, gr%sb%dim
          write(str_tmp, '(f12.6, a)') units_from_atomic(unit_one/units_out%length, kpoint(idir)), ','
          message(1) = trim(message(1)) // trim(str_tmp)
          if(idir == gr%sb%dim) then
            message(1) = trim(message(1)) // ")"
          else
            message(1) = trim(message(1)) // ","
          endif
        enddo
        call messages_info(1, iunit)
      end if

      write(message(1), '(a4,1x,a5)') '#st',' Spin'
      do idir = 1, gr%sb%dim
        write(str_tmp, '(a,a1,a)') '        <p', index2axis(idir), '>'
        message(1) = trim(message(1)) // trim(str_tmp)
      enddo
      write(str_tmp, '(4x,a12,1x)') 'Occupation '
      message(1) = trim(message(1)) // trim(str_tmp)
      call messages_info(1, iunit)
      
      do ist = 1, st%nst
        do is = 0, ns-1

          if(is  ==  0) cspin = 'up'
          if(is  ==  1) cspin = 'dn'
          if(st%d%ispin  ==  UNPOLARIZED .or. st%d%ispin  ==  SPINORS) cspin = '--'
          
          write(message(1), '(i4,3x,a2,1x)') ist, trim(cspin)
          do idir = 1, gr%sb%dim
            write(str_tmp, '(f12.6)') momentum(idir, ist, ik)
            message(1) = trim(message(1)) // trim(str_tmp)
          enddo
          write(str_tmp, '(3x,f12.6)') st%occ(ist, ik+is)
          message(1) = trim(message(1)) // trim(str_tmp)
          call messages_info(1, iunit)
          
        end do
      end do
      
      write(message(1),'(a)') ''
      call messages_info(1, iunit)      
      
    end do
    
    SAFE_DEALLOCATE_A(momentum)
    call io_close(iunit)

    POP_SUB(output_me_out_momentum)
  end subroutine output_me_out_momentum


  ! ---------------------------------------------------------
  subroutine output_me_out_ang_momentum(fname, st, gr)
    character(len=*), intent(in)    :: fname
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr

    integer            :: iunit, ik, ist, is, ns, idir, kstart, kend
    character(len=80)  :: tmp_str(MAX_DIM), cspin
    FLOAT              :: angular(3), lsquare, kpoint(1:MAX_DIM)
    FLOAT, allocatable :: ang(:, :, :), ang2(:, :)
#if defined(HAVE_MPI)
    integer            :: tmp
    FLOAT, allocatable :: lang(:, :)
    integer            :: kn
#endif

    PUSH_SUB(output_me_out_ang_momentum)

    ns = 1
    if(st%d%nspin == 2) ns = 2
    ASSERT(gr%sb%dim == 3)

    iunit = io_open(fname, action='write')

    if(mpi_grp_is_root(mpi_world)) then
      write(iunit,'(a)') 'Warning: When non-local pseudopotentials are used '
      write(iunit,'(a)') '         the numbers below may be meaningless.    '
      write(iunit,'(a)') '                                                  '
      write(iunit,'(a)') 'Angular Momentum of the KS states [dimensionless]:'
      ! r x k is dimensionless. we do not include hbar.
      if (st%d%nik > ns) then
        message(1) = 'k-points [' // trim(units_abbrev(unit_one/units_out%length)) // ']'
        call messages_info(1, iunit)
      end if
    end if

    SAFE_ALLOCATE(ang (1:st%nst, 1:st%d%nik, 1:3))
    SAFE_ALLOCATE(ang2(1:st%nst, 1:st%d%nik))

    if (states_are_real(st)) then
      call dstates_angular_momentum(st, gr, ang, ang2)
    else
      call zstates_angular_momentum(st, gr, ang, ang2)
    end if

    kstart = st%d%kpt%start
    kend = st%d%kpt%end
    do idir = 1, 3
      angular(idir) = states_eigenvalues_sum(st, ang(st%st_start:st%st_end, kstart:kend, idir))
    enddo
    lsquare = states_eigenvalues_sum(st, ang2(st%st_start:st%st_end, kstart:kend))

#if defined(HAVE_MPI)
    if(st%d%kpt%parallel) then
      kn = st%d%kpt%nlocal
      
      ASSERT(.not. st%parallel_in_states)
      
      ! note: could use lmpi_gen_allgatherv here?
      SAFE_ALLOCATE(lang(1:st%lnst, 1:kn))
      do idir = 1, 3
        lang(1:st%lnst, 1:kn) = ang(st%st_start:st%st_end, kstart:kend, idir)
        call MPI_Allgatherv(lang, st%nst*kn, MPI_FLOAT, &
          ang(:, :, idir), st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
          st%d%kpt%mpi_grp%comm, mpi_err)
      end do
      lang(1:st%lnst, 1:kn) = ang2(st%st_start:st%st_end, kstart:kend)
      call MPI_Allgatherv(lang, st%nst*kn, MPI_FLOAT, &
        ang2, st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
        st%d%kpt%mpi_grp%comm, mpi_err)
      SAFE_DEALLOCATE_A(lang)
   end if

   if(st%parallel_in_states) then
      SAFE_ALLOCATE(lang(1:st%lnst, 1))
    endif
#endif

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then

        kpoint = M_ZERO
        kpoint(1:gr%sb%dim) = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))
        
        write(message(1), '(a,i4, a)') '#k =', ik, ', k = ('
        do idir = 1, gr%sb%dim
          write(tmp_str(1), '(f12.6, a)') units_from_atomic(unit_one/units_out%length, kpoint(idir)), ','
          message(1) = trim(message(1)) // trim(tmp_str(1))
          if(idir == gr%sb%dim) then
            message(1) = trim(message(1)) // ")"
          else
            message(1) = trim(message(1)) // ","
          endif
        enddo
        call messages_info(1, iunit)
      end if
      
      ! Exchange ang and ang2.
#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        ASSERT(.not. st%d%kpt%parallel)

        do is = 1, ns
          do idir = 1, 3
            lang(1:st%lnst, 1) = ang(st%st_start:st%st_end, ik+is-1, idir)
            call lmpi_gen_allgatherv(st%lnst, lang(:, 1), tmp, ang(:, ik+is-1, idir), st%mpi_grp)
          end do
          lang(1:st%lnst, 1) = ang2(st%st_start:st%st_end, ik+is-1)
          call lmpi_gen_allgatherv(st%lnst, lang(:, 1), tmp, ang2(:, ik+is-1), st%mpi_grp)
        enddo
      end if
#endif
      write(message(1), '(a4,1x,a5,4a12,4x,a12,1x)')       &
        '#st',' Spin','        <Lx>', '        <Ly>', '        <Lz>', '        <L2>', 'Occupation '
      call messages_info(1, iunit)

      if(mpi_grp_is_root(mpi_world)) then
        do ist = 1, st%nst
          do is = 0, ns-1
            
            if(is  ==  0) cspin = 'up'
            if(is  ==  1) cspin = 'dn'
            if(st%d%ispin  ==  UNPOLARIZED .or. st%d%ispin  ==  SPINORS) cspin = '--'
            
            write(tmp_str(1), '(i4,3x,a2)') ist, trim(cspin)
            write(tmp_str(2), '(1x,4f12.6,3x,f12.6)') &
              (ang(ist, ik+is, idir), idir = 1, 3), ang2(ist, ik+is), st%occ(ist, ik+is)
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))
            call messages_info(1, iunit)
          end do
        end do
      end if
      write(message(1),'(a)') ''
      call messages_info(1, iunit)
      
    end do

    write(message(1),'(a)') 'Total Angular Momentum L [dimensionless]'
    write(message(2),'(10x,4f12.6)') angular(1:3), lsquare
    call messages_info(2, iunit)

    call io_close(iunit)

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      SAFE_DEALLOCATE_A(lang)
    endif
#endif
    
    SAFE_DEALLOCATE_A(ang)
    SAFE_DEALLOCATE_A(ang2)
    
    POP_SUB(output_me_out_ang_momentum)
  end subroutine output_me_out_ang_momentum


#include "undef.F90"
#include "real.F90"
#include "output_me_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "output_me_inc.F90"


end module output_me_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
