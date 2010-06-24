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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
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
  use parser_m
  use mpi_m
  use mpi_lib_m
  use poisson_m
  use profiling_m
  use projector_m
  use simul_box_m
  use states_m
  use states_calc_m
  use states_dim_m
  use unit_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::           &
    output_me_t,      &
    output_me_init,   &
    output_me

  type output_me_t
    integer :: what                ! what to output
    integer :: ks_multipoles       ! If output_ksdipole, this number sets up which matrix elements will
                                   ! be printed: e.g. if ksmultipoles = 3, the dipole, quadrupole and 
                                   ! octopole matrix elements (between Kohn-Sham or single-particle orbitals).
  end type output_me_t

  integer, parameter, public ::               &
       output_me_momentum       =   1, &
       output_me_ang_momentum   =   2, &
       output_me_one_body       =   4, &
       output_me_two_body       =   8, &
       output_me_ks_multipoles  =  16

contains
  
  ! ---------------------------------------------------------
  subroutine output_me_init(this, sb)
    type(output_me_t), intent(out) :: this
    type(simul_box_t), intent(in)  :: sb

    call push_sub('output_me.output_me_init')

    !%Variable OutputMatrixElements
    !%Type flag
    !%Default no
    !%Section Output
    !%Description
    !% Specifies what matrix elements to print.
    !% The output files go into the <tt>static</tt> directory, except when
    !% running a time-dependent simulation, when the directory <tt>td.XXXXXXX</tt> is used.
    !% Example: "momentum + ks_multipoles"
    !%Option multipoles 1
    !% TODO
    !%Option ang_momentum 2
    !% TODO
    !%Option one_body 4
    !% <math>&lt;i|T + V_{ext}|j&gt;</math>
    !%Option two_body 8
    !% <math>&lt;ij| 1/|r_1-r_2| |kl&gt;</math>
    !%Option momentum 16
    !% TODO
    !%End

    call parse_integer(datasets_check('OutputMatrixElements'), 0, this%what)
    if(.not.varinfo_valid_option('OutputMatrixElements', this%what, is_flag=.true.)) then
      call input_error('OutputMatrixElements')
    end if

    if(sb%dim.ne.2 .and. sb%dim.ne.3) this%what = iand(this%what, not(output_me_ang_momentum))

    if(iand(this%what, output_me_ks_multipoles).ne.0) then
      !%Variable OutputMEMultipoles
      !%Type integer
      !%Default 1
      !%Section Output
      !%Description
      !% This variable decides which multipole moments
      !% are printed out: <i>e.g.</i>, if 1, then the
      !% program will print three files, <tt>ks_multipoles.x</tt> (<tt>x</tt>=1,2,3), containing
      !% respectively the (1,-1), (1,0) and (1,1) multipole matrix elements
      !% between Kohn-Sham states.
      !%End
      call parse_integer(datasets_check('OutputMatrixElementsL'), 1, this%ks_multipoles)
    end if

    call pop_sub('output_me.output_me_init')
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
    
    call push_sub('output_me.output_me')

    if(iand(this%what, output_me_momentum).ne.0) then
      write(fname,'(2a)') trim(dir), '/ks_me_momentum'
      call output_me_out_momentum(fname, st, gr)
    end if

    if(iand(this%what, output_me_ang_momentum).ne.0) then
      write(fname,'(2a)') trim(dir), '/ks_me_angular_momentum'
      call output_me_out_ang_momentum(fname, st, gr)
    end if

    if(iand(this%what, output_me_ks_multipoles).ne.0) then
      ! The files will be called matrix_elements.x. The content of each file
      ! should be clear from the header of each file.
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

    if(iand(this%what, output_me_one_body).ne.0) then
      message(1) = "Computing one-body matrix elements"
      call write_info(1)
      if (states_are_real(st)) then
        call done_body(dir, gr, geo, st, hm)
      else
        call zone_body(dir, gr, geo, st, hm)
      end if
    end if

    if(iand(this%what, output_me_two_body).ne.0) then
      message(1) = "Computing two-body matrix elements"
      call write_info(1)
      if (states_are_real(st)) then
        call dtwo_body(dir, gr, st)
      else
        call ztwo_body(dir, gr, st)
      end if
    end if

    call pop_sub('output_me.output_me')
  end subroutine output_me


  ! ---------------------------------------------------------
  subroutine output_me_out_momentum(fname, st, gr)
    character(len=*), intent(in) :: fname
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr

    integer            :: ik, j, is, ns, iunit
    character(len=80)  :: cspin
    FLOAT              :: o, kpoint(1:MAX_DIM)
    FLOAT, allocatable :: momentum(:,:,:)

    call push_sub('output_me.output_me_out_momentum')   

    SAFE_ALLOCATE(momentum(1:gr%sb%dim, 1:st%nst, 1:st%d%nik))

    if (states_are_real(st)) then
      call dstates_calc_momentum(gr, st, momentum)
    else
      call zstates_calc_momentum(gr, st, momentum)
    end if

    iunit = io_open(fname, action='write')

    ns = 1
    if(st%d%nspin == 2) ns = 2

    write(message(1),'(a)') 'Momentum of the KS states [a.u.]:'
    call write_info(1, iunit)      
    if (st%d%nik > ns) then
      message(1) = 'k-points [' // trim(units_abbrev(unit_one/units_out%length)) // ']'
      call write_info(1, iunit)
    end if

    do ik = 1, st%d%nik, ns
      kpoint = M_ZERO
      kpoint = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))

      if(st%d%nik > ns) then
        write(message(1), '(a,i4,3(a,f12.6),a)') '#k =', ik, ', k = (',  &
             units_from_atomic(unit_one/units_out%length, kpoint(1)), ',', &
             units_from_atomic(unit_one/units_out%length, kpoint(2)), ',', &
             units_from_atomic(unit_one/units_out%length, kpoint(3)), ')'
        call write_info(1, iunit)
      end if

      write(message(1), '(a4,1x,a5,3a12,4x,a12,1x)')       &
           '#st',' Spin','       <px>', '        <py>', '        <pz>', 'Occupation '
      call write_info(1, iunit)
      
      do j = 1, st%nst
        do is = 0, ns-1
          
          if(j > st%nst) then
            o = M_ZERO
          else
            o = st%occ(j, ik+is)
          end if
          
          if(is.eq.0) cspin = 'up'
          if(is.eq.1) cspin = 'dn'
          if(st%d%ispin.eq.UNPOLARIZED.or.st%d%ispin.eq.SPINORS) cspin = '--'
          
          write(message(1), '(i4,3x,a2,1x,3f12.6,3x,f12.6)')        &
               j, trim(cspin), momentum(1:min(gr%sb%dim, 3), j, ik), o
          call write_info(1, iunit)
          
        end do
      end do
      
      write(message(1),'(a)') ''
      call write_info(1, iunit)      
      
    end do
    
    SAFE_DEALLOCATE_A(momentum)
    call io_close(iunit)

    call pop_sub('output_me.output_me_out_momentum')   
  end subroutine output_me_out_momentum


  ! ---------------------------------------------------------
  subroutine output_me_out_ang_momentum(fname, st, gr)
    character(len=*), intent(in)    :: fname
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr

    integer            :: iunit, ik, ist, is, ns, j
    character(len=80)  :: tmp_str(MAX_DIM), cspin
    FLOAT              :: angular(3), lsquare, o, kpoint(1:MAX_DIM)
    FLOAT, allocatable :: ang(:, :, :), ang2(:, :)
#if defined(HAVE_MPI)
    integer            :: tmp
    FLOAT, allocatable :: lang(:, :)
    integer            :: kstart, kend, kn
#endif

    call push_sub('output_me.output_me_out_ang_momentum')

    ns = 1
    if(st%d%nspin == 2) ns = 2

    iunit = io_open(fname, action='write')

    if(mpi_grp_is_root(mpi_world)) then
      write(iunit,'(a)') 'Warning: When non-local pseudopotentials are used '
      write(iunit,'(a)') '         the numbers below may be meaningless.    '
      write(iunit,'(a)') '                                                  '
      write(iunit,'(a)') 'Angular Momentum of the KS states [dimensionless]:'
      if (st%d%nik > ns) then
        message(1) = 'k-points [' // trim(units_abbrev(unit_one/units_out%length)) // ']'
        call write_info(1, iunit)
      end if
    end if

    SAFE_ALLOCATE(ang (1:st%nst, 1:st%d%nik, 1:MAX_DIM))
    SAFE_ALLOCATE(ang2(1:st%nst, 1:st%d%nik))
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        if (states_are_real(st)) then
          call dstates_angular_momentum(gr, st%dpsi(:, :, ist, ik), ang(ist, ik, :), ang2(ist, ik))
        else
          call zstates_angular_momentum(gr, st%zpsi(:, :, ist, ik), ang(ist, ik, :), ang2(ist, ik))
        end if
      end do
    end do

    angular(1) =  states_eigenvalues_sum(st, ang (st%st_start:st%st_end, :, 1))
    angular(2) =  states_eigenvalues_sum(st, ang (st%st_start:st%st_end, :, 2))
    angular(3) =  states_eigenvalues_sum(st, ang (st%st_start:st%st_end, :, 3))
    lsquare    =  states_eigenvalues_sum(st, ang2(st%st_start:st%st_end, :))

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then

        kpoint = M_ZERO
        kpoint = kpoints_get_point(gr%sb%kpoints, states_dim_get_kpoint_index(st%d, ik))
        
        write(message(1), '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
             units_from_atomic(unit_one/units_out%length, kpoint(1)), ',', &
             units_from_atomic(unit_one/units_out%length, kpoint(2)), ',', &
             units_from_atomic(unit_one/units_out%length, kpoint(3)), ')'
        call write_info(1, iunit)
      end if
      
      ! Exchange ang and ang2.
#if defined(HAVE_MPI)
      if(st%d%kpt%parallel) then
        kstart = st%d%kpt%start
        kend = st%d%kpt%end
        kn = st%d%kpt%nlocal
        
        ASSERT(.not. st%parallel_in_states)
        
        SAFE_ALLOCATE(lang(1:st%lnst, 1:kn))
        do j = 1, 3
          lang(1:st%lnst, 1:kn) = ang(st%st_start:st%st_end, kstart:kend, j)
          call MPI_Allgatherv(lang, st%nst*kn, MPI_FLOAT, &
               ang(:, :, j), st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
               st%d%kpt%mpi_grp%comm, mpi_err)
        end do
        lang(1:st%lnst, 1:kn) = ang2(st%st_start:st%st_end, kstart:kend)
        call MPI_Allgatherv(lang, st%nst*kn, MPI_FLOAT, &
             ang2(:, :), st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
             st%d%kpt%mpi_grp%comm, mpi_err)
        SAFE_DEALLOCATE_A(lang)
      end if
      
      if(st%parallel_in_states) then
        SAFE_ALLOCATE(lang(1:st%lnst, 1:1))
        do j = 1, 3
          lang(1:st%lnst, 1) = ang(st%st_start:st%st_end, ik, j)
          call lmpi_gen_allgatherv(st%lnst, lang(:, 1), tmp, ang(:, ik, j), st%mpi_grp)
        end do
        lang(1:st%lnst, 1) = ang2(st%st_start:st%st_end, ik)
        call lmpi_gen_allgatherv(st%lnst, lang(:, 1), tmp, ang2(:, ik), st%mpi_grp)
        SAFE_DEALLOCATE_A(lang)
      end if
#endif
      write(message(1), '(a4,1x,a5,4a12,4x,a12,1x)')       &
           '#st',' Spin','        <Lx>', '        <Ly>', '        <Lz>', '        <L2>', 'Occupation '
      call write_info(1, iunit)

      if(mpi_grp_is_root(mpi_world)) then
        do j = 1, st%nst
          do is = 0, ns-1

            if(j > st%nst) then
              o = M_ZERO
            else
              o = st%occ(j, ik+is)
            end if
            
            if(is.eq.0) cspin = 'up'
            if(is.eq.1) cspin = 'dn'
            if(st%d%ispin.eq.UNPOLARIZED.or.st%d%ispin.eq.SPINORS) cspin = '--'

            write(tmp_str(1), '(i4,3x,a2)') j, trim(cspin)
            write(tmp_str(2), '(1x,4f12.6,3x,f12.6)') &
                 ang(j, ik+is, 1:3), ang2(j, ik+is), o
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))
            call write_info(1, iunit)
          end do
        end do
      end if
      write(message(1),'(a)') ''
      call write_info(1, iunit)

    end do

    write(message(1),'(a)') 'Total Angular Momentum L [dimensionless]'
    write(message(2),'(10x,4f12.6)') angular(1:3), lsquare
    call write_info(2, iunit)

    call io_close(iunit)

    SAFE_DEALLOCATE_A(ang)
    SAFE_DEALLOCATE_A(ang2)
    
    call pop_sub('output_me.output_me_out_ang_momentum')
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
