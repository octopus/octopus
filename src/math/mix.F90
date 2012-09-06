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
!! $Id$

#include "global.h"

module mix_m
  use datasets_m
  use global_m
  use lalg_adv_m
  use lalg_basic_m
  use messages_m
  use parser_m
  use profiling_m
  use types_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    mix_set_mixing,             &
    mix_t,                      &
    mix_init,                   &
    mix_clear,                  &
    mix_end,                    &
    dmixing,                    &
    zmixing

  integer, parameter, public :: &
    MIX_LINEAR  = 0,            &
    MIX_GRPULAY = 1,            &
    MIX_BROYDEN = 2

  type mix_t
    private
    integer  :: type_of_mixing

    FLOAT :: alpha              ! vnew = (1-alpha)*vin + alpha*vout

    integer :: ns               ! number of steps used to extrapolate the new vector

    integer :: d1, d2, d3

    FLOAT, pointer :: ddf(:, :, :, :)
    FLOAT, pointer :: ddv(:, :, :, :)
    FLOAT, pointer :: df_old(:, :, :)
    FLOAT, pointer :: dvin_old(:, :, :)

    CMPLX, pointer :: zdf(:, :, :, :)
    CMPLX, pointer :: zdv(:, :, :, :)
    CMPLX, pointer :: zf_old(:, :, :)
    CMPLX, pointer :: zvin_old(:, :, :)

    integer :: last_ipos
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
    type(type_t) :: func_type_
    character(len=32) :: prefix

    PUSH_SUB(mix_init)

    def = MIX_BROYDEN
    if(present(def_)) def = def_
    if(present(func_type)) then 
      func_type_ = func_type
    else 
      func_type_ = TYPE_FLOAT
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
    call parse_integer(datasets_check(trim(prefix)//'TypeOfMixing'), def, smix%type_of_mixing)
    if(.not.varinfo_valid_option('TypeOfMixing', smix%type_of_mixing)) call input_error('TypeOfMixing')
    call messages_print_var_option(stdout, "TypeOfMixing", smix%type_of_mixing)

    !%Variable Mixing
    !%Type float
    !%Default 0.3
    !%Section SCF::Mixing
    !%Description
    !% Both the linear and the Broyden scheme depend on a "mixing parameter", set by this variable.  Must be 0 < <tt>Mixing</tt> <= 1.
    !%End
    if (smix%type_of_mixing == MIX_LINEAR .or. smix%type_of_mixing == MIX_BROYDEN) then
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
    if (smix%type_of_mixing == MIX_GRPULAY .or. smix%type_of_mixing == MIX_BROYDEN) then
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

    select case (smix%type_of_mixing)
    case (MIX_GRPULAY)
      if(func_type_ == TYPE_FLOAT) then 
        SAFE_ALLOCATE(     smix%ddf(1:d1, 1:d2, 1:d3, 1:smix%ns + 1))
        SAFE_ALLOCATE(smix%dvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     smix%ddv(1:d1, 1:d2, 1:d3, 1:smix%ns + 1))
        SAFE_ALLOCATE(  smix%df_old(1:d1, 1:d2, 1:d3))
      else
        SAFE_ALLOCATE(     smix%zdf(1:d1, 1:d2, 1:d3, 1:smix%ns + 1))
        SAFE_ALLOCATE(smix%zvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     smix%zdv(1:d1, 1:d2, 1:d3, 1:smix%ns + 1))
        SAFE_ALLOCATE(  smix%zf_old(1:d1, 1:d2, 1:d3))
      end if

    case (MIX_BROYDEN)
      if(func_type_ == TYPE_FLOAT) then 
        SAFE_ALLOCATE(     smix%ddf(1:d1, 1:d2, 1:d3, 1:smix%ns))
        SAFE_ALLOCATE(smix%dvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     smix%ddv(1:d1, 1:d2, 1:d3, 1:smix%ns))
        SAFE_ALLOCATE(  smix%df_old(1:d1, 1:d2, 1:d3))
      else
        SAFE_ALLOCATE(     smix%zdf(1:d1, 1:d2, 1:d3, 1:smix%ns))
        SAFE_ALLOCATE(smix%zvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     smix%zdv(1:d1, 1:d2, 1:d3, 1:smix%ns))
        SAFE_ALLOCATE(  smix%zf_old(1:d1, 1:d2, 1:d3))
      end if

    end select

    call mix_clear(smix, func_type_)

    POP_SUB(mix_init)
  end subroutine mix_init


  ! ---------------------------------------------------------
  subroutine mix_clear(smix, func_type)
    type(mix_t),             intent(inout) :: smix
    type(type_t),  optional, intent(in)    :: func_type
    
    type(type_t) :: func_type_

    PUSH_SUB(mix_clear)

    if(present(func_type)) then 
      func_type_ = func_type
    else 
      func_type_ = TYPE_FLOAT
    end if

    select case (smix%type_of_mixing)
    case (MIX_GRPULAY)
      if(func_type_ == TYPE_FLOAT) then 
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

    case (MIX_BROYDEN)
      if(func_type_ == TYPE_FLOAT) then 
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

    end select

    smix%last_ipos = 0

    POP_SUB(mix_clear)
  end subroutine mix_clear


  ! ---------------------------------------------------------
  subroutine mix_end(smix)
    type(mix_t), intent(inout) :: smix

    PUSH_SUB(mix_end)

    ! Arrays got allocated for all mixing schemes, except linear mixing
    if (smix%type_of_mixing .ne. MIX_LINEAR) then
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
    
    if(smix%type_of_mixing == MIX_LINEAR) then
      smix%alpha = newmixing
    else
    !  message(1) = "Mixing can only be adjusted in linear mixing scheme."
    !  call messages_fatal(1)
    endif
    
    POP_SUB(mix_set_mixing)
  end subroutine mix_set_mixing

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
