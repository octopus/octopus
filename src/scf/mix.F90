!! Copyright (C) 2002-2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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
!! $Id$

#include "global.h"

module mix_m
  use datasets_m
  use global_m
  use io_m
  use io_function_m
  use lalg_adv_m
  use lalg_basic_m
  use mesh_m
  use messages_m
  use mpi_m
  use parser_m
  use profiling_m
  use restart_m
  use types_m
  use unit_system_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    mix_set_mixing,             &
    mix_t,                      &
    mix_init,                   &
    mix_clear,                  &
    mix_end,                    &
    mix_dump,                   &
    mix_load,                   &
    dmixing,                    &
    zmixing

  integer, parameter, public :: &
    MIX_LINEAR  = 0,            &
    MIX_GRPULAY = 1,            &
    MIX_BROYDEN = 2

  type mix_t
    private
    integer :: scheme           !< the mixing scheme used (linear, broyden, etc)

    FLOAT :: alpha              !< the mixing coefficient (in linear mixing: vnew = (1-alpha)*vin + alpha*vout)


    type(type_t) :: func_type   !< type of the functions to be mixed
    integer :: ns               !< number of steps used to extrapolate the new vector
    integer :: ns_stored        !< number of steps that are actually stored at the moment

    integer :: d1, d2, d3, d4   !< the dimensions of the arrays that store the information from the previous iterations

    integer :: last_ipos        !< where is the information about the last iteration stored in arrays df and dv

    FLOAT, pointer :: ddf(:, :, :, :)
    FLOAT, pointer :: ddv(:, :, :, :)
    FLOAT, pointer :: df_old(:, :, :)
    FLOAT, pointer :: dvin_old(:, :, :)

    CMPLX, pointer :: zdf(:, :, :, :)
    CMPLX, pointer :: zdv(:, :, :, :)
    CMPLX, pointer :: zf_old(:, :, :)
    CMPLX, pointer :: zvin_old(:, :, :)

  end type mix_t

contains

  ! ---------------------------------------------------------
  subroutine mix_init(smix, d1, d2, d3, def_, func_type, prefix_)
    type(mix_t),                intent(out) :: smix
    integer,                    intent(in)  :: d1, d2, d3
    integer,          optional, intent(in)  :: def_
    type(type_t),     optional, intent(in)  :: func_type
    character(len=*), optional, intent(in)  :: prefix_

    integer :: def
    character(len=32) :: prefix

    PUSH_SUB(mix_init)

    def = MIX_BROYDEN
    if(present(def_)) def = def_
    if(present(func_type)) then 
      smix%func_type = func_type
    else 
      smix%func_type = TYPE_FLOAT
    end if
    prefix = ""
    if(present(prefix_)) prefix = prefix_


    !%Variable TypeOfMixing
    !%Type integer
    !%Default broyden
    !%Section SCF::Mixing
    !%Description
    !% The scheme used to produce, at each iteration in the self-consistent cycle
    !% that attempts to solve the Kohn-Sham equations, the input density from the value
    !% of the input and output densities of previous iterations.
    !%Option linear 0
    !% Simple linear mixing.
    !%Option gr_pulay 1
    !% "Guaranteed-reduction" Pulay scheme [D. R. Bowler and M. J. Gillan, <i>Chem. Phys. 
    !% Lett.</i> <b>325</b>, 473 (2000)].
    !%Option broyden 2
    !% Broyden scheme [C. G Broyden, <i>Math. Comp.</i> <b>19</b>, 577 (1965); 
    !% D. D. Johnson, <i>Phys. Rev. B</i> <b>38</b>, 12807 (1988)].
    !%End
    call parse_integer(datasets_check(trim(prefix)//'TypeOfMixing'), def, smix%scheme)
    if(.not.varinfo_valid_option('TypeOfMixing', smix%scheme)) call input_error('TypeOfMixing')
    call messages_print_var_option(stdout, "TypeOfMixing", smix%scheme)

    !%Variable Mixing
    !%Type float
    !%Default 0.3
    !%Section SCF::Mixing
    !%Description
    !% Both the linear and the Broyden scheme depend on a "mixing parameter", set by this variable. 
    !% Must be 0 < <tt>Mixing</tt> <= 1.
    !%End
    if (smix%scheme == MIX_LINEAR .or. smix%scheme == MIX_BROYDEN) then
      call parse_float(datasets_check(trim(prefix)//'Mixing'), CNST(0.3), smix%alpha)
      if(smix%alpha <= M_ZERO .or. smix%alpha > M_ONE) call input_error('Mixing')
    end if

    !%Variable MixNumberSteps
    !%Type integer
    !%Default 3
    !%Section SCF::Mixing
    !%Description
    !% In the Broyden and GR-Pulay schemes, the new input density or potential is constructed
    !% from the values of the densities/potentials of a given number of previous iterations.
    !% This number is set by this variable. Must be greater than 1.
    !%End
    if (smix%scheme == MIX_GRPULAY .or. smix%scheme == MIX_BROYDEN) then
      call parse_integer(datasets_check(trim(prefix)//'MixNumberSteps'), 3, smix%ns)
      if(smix%ns <= 1) call input_error('MixNumberSteps')
    end if


    nullify(smix%ddf)
    nullify(smix%ddv)
    nullify(smix%df_old)
    nullify(smix%dvin_old)

    nullify(smix%zdf)
    nullify(smix%zdv)
    nullify(smix%zf_old)
    nullify(smix%zvin_old)

    smix%d1 = d1
    smix%d2 = d2
    smix%d3 = d3
    select case (smix%scheme)
    case (MIX_LINEAR)
      smix%d4 = 0
    case (MIX_GRPULAY)
      smix%d4 = smix%ns + 1
    case (MIX_BROYDEN)
      smix%d4 = smix%ns
    end select

    if (smix%scheme /= MIX_LINEAR) then
      if(smix%func_type == TYPE_FLOAT) then 
        SAFE_ALLOCATE(     smix%ddf(1:d1, 1:d2, 1:d3, 1:smix%d4))
        SAFE_ALLOCATE(smix%dvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     smix%ddv(1:d1, 1:d2, 1:d3, 1:smix%d4))
        SAFE_ALLOCATE(  smix%df_old(1:d1, 1:d2, 1:d3))
      else
        SAFE_ALLOCATE(     smix%zdf(1:d1, 1:d2, 1:d3, 1:smix%d4))
        SAFE_ALLOCATE(smix%zvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     smix%zdv(1:d1, 1:d2, 1:d3, 1:smix%d4))
        SAFE_ALLOCATE(  smix%zf_old(1:d1, 1:d2, 1:d3))
      end if
    end if

    call mix_clear(smix)

    POP_SUB(mix_init)
  end subroutine mix_init


  ! ---------------------------------------------------------
  subroutine mix_clear(smix)
    type(mix_t),             intent(inout) :: smix
    
    PUSH_SUB(mix_clear)

    if (smix%scheme /= MIX_LINEAR) then
      if(smix%func_type == TYPE_FLOAT) then 
        smix%ddf = M_ZERO
        smix%ddv = M_ZERO
        smix%dvin_old = M_ZERO
        smix%df_old = M_ZERO
      else
        smix%zdf = M_z0
        smix%zdv = M_z0
        smix%zvin_old = M_z0
        smix%zf_old = M_z0
      end if
    end if

    smix%ns_stored = 0
    smix%last_ipos = 0

    POP_SUB(mix_clear)
  end subroutine mix_clear


  ! ---------------------------------------------------------
  subroutine mix_end(smix)
    type(mix_t), intent(inout) :: smix

    PUSH_SUB(mix_end)

    ! Arrays got allocated for all mixing schemes, except linear mixing
    if (smix%scheme /= MIX_LINEAR) then
      SAFE_DEALLOCATE_P(smix%ddf)
      SAFE_DEALLOCATE_P(smix%ddv)
      SAFE_DEALLOCATE_P(smix%dvin_old)
      SAFE_DEALLOCATE_P(smix%df_old)

      SAFE_DEALLOCATE_P(smix%zdf)
      SAFE_DEALLOCATE_P(smix%zdv)
      SAFE_DEALLOCATE_P(smix%zvin_old)
      SAFE_DEALLOCATE_P(smix%zf_old)
    end if

    POP_SUB(mix_end)
  end subroutine mix_end


  ! ---------------------------------------------------------
  subroutine mix_set_mixing(smix, newmixing)
    type(mix_t), intent(inout) :: smix
    FLOAT, intent(in):: newmixing

    PUSH_SUB(mix_set_mixing)
    
    if(smix%scheme == MIX_LINEAR) then
      smix%alpha = newmixing
    else
    !  message(1) = "Mixing can only be adjusted in linear mixing scheme."
    !  call messages_fatal(1)
    endif
    
    POP_SUB(mix_set_mixing)
  end subroutine mix_set_mixing


  ! ---------------------------------------------------------
  subroutine mix_dump(restart, smix, mesh, ierr)
    type(restart_t), intent(in)  :: restart
    type(mix_t),     intent(in)  :: smix
    type(mesh_t),    intent(in)  :: mesh
    integer,         intent(out) :: ierr

    integer :: iunit
    integer :: id2, id3, id4
    character(len=80) :: filename

    PUSH_SUB(mix_dump)

    if (restart_skip(restart)) then
      ierr = 0
      POP_SUB(mix_dump)
      return
    end if

    ! functions to be written need to be compatible with the mesh
    ASSERT(mesh%np == smix%d1)

    ! First we write some information about the mixing
    iunit = restart_open(restart, 'mixing')
    if (mpi_grp_is_root(mpi_world)) then
      write(iunit, '(a11,i1)')  'scheme=    ', smix%scheme
      ! Number of global mesh points have to be written, not only smix%d1
      write(iunit, '(a11,i10)') 'd1=        ', mesh%np_global 
      write(iunit, '(a11,i10)') 'd2=        ', smix%d2
      write(iunit, '(a11,i10)') 'd3=        ', smix%d3
      write(iunit, '(a11,i10)') 'd4=        ', smix%ns_stored
      write(iunit, '(a11,i10)') 'last_ipos= ', smix%last_ipos
    end if
    call restart_close(restart, iunit)

    ! Now we write the different functions. 
    ! These are not needed when using linear mixing, so we will make sure we skip this step in that case.
    if (smix%scheme /= MIX_LINEAR) then
      do id2 = 1, smix%d2
        do id3 = 1, smix%d3
          do id4 = 1, smix%ns_stored

            write(filename,'(a3,i2.2,i2.2,i2.2)') 'df_', id2, id3, id4
            if (smix%func_type == TYPE_FLOAT) then
              call drestart_write_function(restart, filename, mesh, smix%ddf(1:mesh%np, id2, id3, id4), ierr)
            else
              call zrestart_write_function(restart, filename, mesh, smix%zdf(1:mesh%np, id2, id3, id4), ierr)
            end if
            if(ierr /= 0) then
              write(message(1), '(3a)') 'Unsuccessful write of "', trim(filename), '"'
              call messages_fatal(1)
            end if

          
            write(filename,'(a3,i2.2,i2.2,i2.2)') 'dv_', id2, id3, id4
            if (smix%func_type == TYPE_FLOAT) then
              call drestart_write_function(restart, filename, mesh, smix%ddv(1:mesh%np, id2, id3, id4), ierr)
            else
              call zrestart_write_function(restart, filename, mesh, smix%zdv(1:mesh%np, id2, id3, id4), ierr)
            end if
            if(ierr /= 0) then
              write(message(1), '(3a)') 'Unsuccessful write of "', trim(filename), '"'
              call messages_fatal(1)
            end if

          end do


          write(filename,'(a6,i2.2,i2.2)') 'f_old_', id2, id3
          if (smix%func_type == TYPE_FLOAT) then
              call drestart_write_function(restart, filename, mesh, smix%df_old(1:mesh%np, id2, id3), ierr)
            else
              call zrestart_write_function(restart, filename, mesh, smix%zf_old(1:mesh%np, id2, id3), ierr)
          end if
          if(ierr /= 0) then
            write(message(1), '(3a)') 'Unsuccessful write of "', trim(filename), '"'
            call messages_fatal(1)
          end if


          write(filename,'(a8,i2.2,i2.2)') 'vin_old_', id2, id3
          if (smix%func_type == TYPE_FLOAT) then
              call drestart_write_function(restart, filename, mesh, smix%dvin_old(1:mesh%np, id2, id3), ierr)
            else
              call zrestart_write_function(restart, filename, mesh, smix%zvin_old(1:mesh%np, id2, id3), ierr)
          end if
          if(ierr /= 0) then
            write(message(1), '(3a)') 'Unsuccessful write of "', trim(filename), '"'
            call messages_fatal(1)
          end if
          
        end do
      end do

    end if

    POP_SUB(mix_dump)
  end subroutine mix_dump


  !---------------------------------------------------------
  ! Error codes:
  !  -1  - type of mixing is not the same.
  !  -2  - dimensions of the arrays are not consistent
  !  -3  - no mixing file found
  !  > 0 - unable to read ierr functions
  subroutine mix_load(restart, smix, mesh, ierr)
    type(restart_t), intent(in)    :: restart
    type(mix_t),     intent(inout) :: smix
    type(mesh_t),    intent(in)    :: mesh
    integer,         intent(out)   :: ierr

    integer :: iunit, err
    integer :: scheme, d1, d2, d3, d4, last_ipos
    integer :: id2, id3, id4, ipos
    character(len=11)  :: str
    character(len=80)  :: filename
    character(len=256) :: line

    PUSH_SUB(mix_load)

    if (restart_skip(restart)) then
      ierr = -3
      POP_SUB(mix_load)
      return
    end if

    ! First we read some information about the mixing
    iunit = restart_open(restart, 'mixing')
    if (iunit < 0) then
      write(message(1), '(3a)') 'Unsuccessful read of file "mixing".'
      call messages_warning(1)
      ierr = -3
      POP_SUB(mix_load)
      return
    end if

    call iopar_read(mesh%mpi_grp, iunit, line, err)
    read(line, *) str, scheme
    call iopar_read(mesh%mpi_grp, iunit, line, err)
    read(line, *) str, d1
    call iopar_read(mesh%mpi_grp, iunit, line, err)
    read(line, *) str, d2
    call iopar_read(mesh%mpi_grp, iunit, line, err)
    read(line, *) str, d3
    call iopar_read(mesh%mpi_grp, iunit, line, err)
    read(line, *) str, d4
    call iopar_read(mesh%mpi_grp, iunit, line, err)
    read(line, *) str, last_ipos

    call restart_close(restart, iunit)

    ! We can only use the restart information if the mixing scheme remained the same
    if (scheme /= smix%scheme) then
      message(1) = "The mixing scheme from the restart data is not the same as the one used in the current calculation."
      call messages_warning(1)
      ierr = -1
      POP_SUB(mix_load)
      return
    end if

    ! Check the dimensions of the arrays to be read
    if (mesh%np_global /= d1 .or. mesh%np /= smix%d1 .or. d2 /= smix%d2 .or. d3 /= smix%d3 ) then
      message(1) = "The dimensions of the arrays from the mixing restart data"
      message(2) = "are not the same as the ones used in this calculation."
      call messages_warning(2)
      ierr = -2
      POP_SUB(mix_load)
      return
    end if


    ! Now we read the different functions.
    ! Note that we may have more or less functions than the ones needed (d4 /= smix%d4)
    ierr = 0
    if (smix%scheme /= MIX_LINEAR) then
      smix%ns_stored = min(d4, smix%d4)
      smix%last_ipos = smix%ns_stored

      do id2 = 1, smix%d2
        do id3 = 1, smix%d3

          if (d4 > 0) ipos = mod(last_ipos, d4) + 1
          do id4 = 1, smix%ns_stored

            write(filename,'(a3,i2.2,i2.2,i2.2)') 'df_', id2, id3, ipos
            if (smix%func_type == TYPE_FLOAT) then
              call drestart_read_function(restart, trim(filename), mesh, smix%ddf(1:mesh%np, id2, id3, id4), err)
            else
              call zrestart_read_function(restart, trim(filename), mesh, smix%zdf(1:mesh%np, id2, id3, id4), err)
            end if
            if(err /= 0) then
              write(message(1), '(3a)') 'Unsuccessful read of "', trim(filename), '"'
              call messages_warning(1)
              ierr = ierr + 1
            end if
          
            write(filename,'(a3,i2.2,i2.2,i2.2)') 'dv_', id2, id3, ipos
            if (smix%func_type == TYPE_FLOAT) then
              call drestart_read_function(restart, trim(filename), mesh, smix%ddv(1:mesh%np, id2, id3, id4), err)
            else
              call zrestart_read_function(restart, trim(filename), mesh, smix%zdv(1:mesh%np, id2, id3, id4), err)
            end if
            if(err /= 0) then
              write(message(1), '(3a)') 'Unsuccessful read of "', trim(filename), '"'
              call messages_warning(1)
              ierr = ierr + 1
            end if

            ipos = mod(ipos, d4) + 1
          end do

          write(filename,'(a6,i2.2,i2.2)') 'f_old_', id2, id3
          if (smix%func_type == TYPE_FLOAT) then
            call drestart_read_function(restart, trim(filename), mesh, smix%df_old(1:mesh%np, id2, id3), err)
          else
            call zrestart_read_function(restart, trim(filename), mesh, smix%zf_old(1:mesh%np, id2, id3), err)
          end if
          if(err /= 0) then
            write(message(1), '(3a)') 'Unsuccessful read of "', trim(filename), '"'
            call messages_warning(1)
            ierr = ierr + 1
          end if

          write(filename,'(a8,i2.2,i2.2)') 'vin_old_', id2, id3
          if (smix%func_type == TYPE_FLOAT) then
            call drestart_read_function(restart, trim(filename), mesh, smix%dvin_old(1:mesh%np, id2, id3), err)
          else
            call zrestart_read_function(restart, trim(filename), mesh, smix%zvin_old(1:mesh%np, id2, id3), err)
          end if
          if(err /= 0) then
            write(message(1), '(3a)') 'Unsuccessful read of "', trim(filename), '"'
            call messages_warning(1)
            ierr = ierr + 1
          end if

        end do
      end do
    end if

    if (ierr /= 0) then
      ! Something went wront, so make sure we start from scratch
      call mix_clear(smix)
    end if

    POP_SUB(mix_load)
  end subroutine mix_load


#include "undef.F90"
#include "real.F90"

#include "mix_inc.F90"

#include "undef.F90"
#include "complex.F90"

#include "mix_inc.F90"

end module mix_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
