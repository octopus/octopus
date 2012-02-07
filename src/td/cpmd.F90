!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id: td.F90 3694 2008-02-15 13:37:54Z marques $

#include "global.h"
  
module cpmd_m
  use batch_m
  use mesh_batch_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use loct_m
  use loct_math_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use states_calc_m
  use states_block_m
  use system_m
  use varinfo_m

  implicit none

  private
  public ::                 &
    cpmd_t,                 &
    cpmd_init,              &
    cpmd_end,               &
    cpmd_electronic_energy, &
    cpmd_restart_read,      &
    cpmd_restart_write,     &
    dcpmd_propagate,        & 
    zcpmd_propagate,        &
    dcpmd_propagate_vel,    & 
    zcpmd_propagate_vel

  type cpmd_t
    private
    integer        :: method
    FLOAT          :: emass
    FLOAT          :: ecorr

    !> for verlet, this stores the previous wfs
    !! for vel_verlet, the time derivative of the wfs
    FLOAT, pointer :: dpsi2(:, :, :, :)
    CMPLX, pointer :: zpsi2(:, :, :, :)
    
  end type cpmd_t

  integer, parameter ::   &
    VEL_VERLET = 1,       &
    VERLET     = 2

  type(profile_t), save :: cpmd_prop, cpmd_orth

contains

  ! ---------------------------------------------------------
  subroutine cpmd_init(this, gr, st)
    type(cpmd_t),        intent(out)   :: this
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    
    PUSH_SUB(cpmd_init)

    !%Variable CPElectronicMass
    !%Type float
    !%Default 1.0
    !%Section Time-Dependent::Propagation
    !%Description
    !% The fictitious electronic mass used to propagate the electronic
    !% wavefunctions in the Car-Parrinello formalism.
    !%End
    
    call parse_float(datasets_check('CPElectronicMass'), CNST(1.0), this%emass)

    !%Variable CPMethod
    !%Type integer
    !%Default verlet
    !%Section Time-Dependent::Propagation
    !%Description
    !% This variable defines how to integrate the Car-Parrinello
    !% equations. The default is <tt>verlet</tt>.
    !%Option verlet 2
    !% Standard Verlet.
    !%Option vel_verlet 1
    !% RATTLE/Velocity Verlet integrator.
    !%End

    call parse_integer(datasets_check('CPMethod'), VERLET, this%method)
    if(.not.varinfo_valid_option('CPMethod', this%method)) call input_error('CPMethod')
    call messages_print_var_option(stdout, 'CPMethod', this%method)

    if(st%parallel_in_states) then
      message(1) = 'CP dynamics does not work with states-parallelization.'
      call messages_fatal(1)
    end if
        
    nullify(this%dpsi2, this%zpsi2)

    if(states_are_real(st)) then
      SAFE_ALLOCATE(this%dpsi2(1:gr%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    else
      SAFE_ALLOCATE(this%zpsi2(1:gr%mesh%np_part, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik))
    end if

    POP_SUB(cpmd_init)

  end subroutine cpmd_init
  

  ! ---------------------------------------------------------
  subroutine cpmd_end(this)
    type(cpmd_t), intent(inout) :: this

    PUSH_SUB(cpmd_end)

    SAFE_DEALLOCATE_P(this%dpsi2)
    SAFE_DEALLOCATE_P(this%zpsi2)

    POP_SUB(cpmd_end)

  end subroutine cpmd_end


  ! ---------------------------------------------------------
  FLOAT function cpmd_electronic_energy(this)
    type(cpmd_t),        intent(in)    :: this

    cpmd_electronic_energy = this%ecorr
    
  end function cpmd_electronic_energy


  ! ---------------------------------------------------------
  subroutine cpmd_restart_write(this, gr, st)
    type(cpmd_t),        intent(in)    :: this
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st

    integer :: ik, ist, idim, ii, err
    character(len=80) :: filename

    PUSH_SUB(cpmd_restart_write)

    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(trim(tmpdir)//'td/cpmd', is_tmp=.true.)
    endif

    ii = 0
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do idim = 1, st%d%dim
          ii = ii + 1
          if (ist < st%st_start .or.  ist > st%st_end) cycle
          write(filename,'(i10.10)') ii
          if(states_are_real(st)) then
            call drestart_write_function(trim(tmpdir)//'td/cpmd', filename, gr%mesh, this%dpsi2(:, idim, ist, ik), err)
          else
            call zrestart_write_function(trim(tmpdir)//'td/cpmd', filename, gr%mesh, this%zpsi2(:, idim, ist, ik), err)
          end if
        end do
      end do
    end do

    POP_SUB(cpmd_restart_write)

  end subroutine cpmd_restart_write


  ! ---------------------------------------------------------
  subroutine cpmd_restart_read(this, gr, st, ierr)
    type(cpmd_t),        intent(inout) :: this
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st
    integer,             intent(out)   :: ierr

    integer :: ik, ist, idim, ii
    character(len=80) :: filename

    PUSH_SUB(cpmd_restart_read)

    ierr = 0

    ii = 0
    do ik = 1, st%d%nik
      do ist = 1, st%nst
        do idim = 1, st%d%dim
          ii = ii + 1
          if (ist < st%st_start .or.  ist > st%st_end) cycle
          write(filename,'(i10.10)') ii
          if(states_are_real(st)) then
            call drestart_read_function(trim(tmpdir)//'td/cpmd', filename, gr%mesh, this%dpsi2(:, idim, ist, ik), ierr)
          else
            call zrestart_read_function(trim(tmpdir)//'td/cpmd', filename, gr%mesh, this%zpsi2(:, idim, ist, ik), ierr)
          end if
          if(ierr /= 0) then
            POP_SUB(cpmd_restart_read)
            return
          end if
        end do
      end do
    end do

    POP_SUB(cpmd_restart_read)

  end subroutine cpmd_restart_read

#include "undef.F90"
#include "real.F90"
#include "cpmd_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "cpmd_inc.F90"


end module cpmd_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
