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
  use global_m
  use grid_m
  use io_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use loct_math_m
  use loct_parser_m
  use profiling_m
  use states_m
  use states_calc_m
  use states_dim_m
  use units_m
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
  subroutine output_me_init(this)
    type(output_me_t), intent(out) :: this

    !%Variable OutputMatrixElements
    !%Type flag
    !%Default no
    !%Section Output
    !%Description
    !% Specifies what matrix elements to print.
    !% The output files go into the "static" directory, except when
    !% running a time-dependent simulation, when the directory "td.XXXXXXX" is used.
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

    call loct_parse_int(datasets_check('OutputMatrixElements'), 0, this%what)
    if(.not.varinfo_valid_option('OutputMatrixElements', this%what, is_flag=.true.)) then
      call input_error('OutputMatrixElements')
    end if

    if(iand(this%what, output_me_ks_multipoles).ne.0) then
      !%Variable OutputMEMultipoles
      !%Type integer
      !%Default 1
      !%Section Output
      !%Description
      !% This variable decides which multipole moments
      !% are printed out: e.g., if 1, then the
      !% program will print three files, ks_multipoles.x (x=1,2,3), containing
      !% respectively the (1,-1), (1,0) and (1,1) multipole matrix elements
      !% between Kohn-Sham states.
      !%End
      call loct_parse_int(datasets_check('OutputMatrixElementsL'), 1, this%ks_multipoles)
    end if

  end subroutine output_me_init


  ! ---------------------------------------------------------
  subroutine output_me(this, st, gr, dir)
    type(output_me_t), intent(in) :: this
    type(states_t),    intent(inout) :: st
    type(grid_t),      intent(inout) :: gr
    character(len=*),  intent(in) :: dir

    integer :: id, ll, mm, ik
    character(len=256) :: fname
    
    if(iand(this%what, output_me_momentum).ne.0) then
      write(fname,'(2a)') trim(dir), '/ks_me_momentum'
      call output_me_out_momentum(fname, st, gr)
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
            if (st%wfs_type == M_REAL) then
              call doutput_me_ks_multipoles(fname, st, gr, ll, mm, ik)
            else
              call zoutput_me_ks_multipoles(fname, st, gr, ll, mm, ik)
            end if

            id = id + 1
          end do
        end do
      end do
    end if

  end subroutine output_me


  ! ---------------------------------------------------------
  subroutine output_me_out_momentum(fname, st, gr)
    character(len=*), intent(in) :: fname
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr

    integer            :: ik, j, is, ns, iunit
    character(len=80)  :: cspin
    FLOAT              :: o
    FLOAT, allocatable :: momentum(:,:,:)

    call push_sub('output_me.output_me_momentum')   

    SAFE_ALLOCATE(momentum(1:gr%sb%dim, 1:st%nst, 1:st%d%nik))

    if (st%wfs_type == M_REAL) then
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
      message(1) = 'Kpoints [' // trim(units_out%length%abbrev) // '^-1]'
      call write_info(1, iunit)
    end if

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then
        write(message(1), '(a,i4,3(a,f12.6),a)') '#k =', ik, ', k = (',  &
             st%d%kpoints(1, ik)*units_out%length%factor, ',',            &
             st%d%kpoints(2, ik)*units_out%length%factor, ',',            &
             st%d%kpoints(3, ik)*units_out%length%factor, ')'
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

    call pop_sub()
  end subroutine output_me_out_momentum


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
