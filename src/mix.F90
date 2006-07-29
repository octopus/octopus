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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module mix_m
  use global_m
  use messages_m
  use datasets_m
  use mesh_m
  use mesh_function_m
  use lib_oct_parser_m
  use lib_basic_alg_m
  use lib_adv_alg_m
  use varinfo_m

  implicit none

  private
  public ::                     &
    mix_t,                      &
    mix_init,                   &
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

    FLOAT :: alpha              !  vnew = (1-alpha)*vin + alpha*vout

    integer :: ns               ! number of steps used to extrapolate the new vector

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
  subroutine mix_init(smix, m, d2, d3, def_, func_type)
    type(mix_t),       intent(out) :: smix
    type(mesh_t),      intent(in)  :: m
    integer,           intent(in)  :: d2, d3
    integer, optional, intent(in)  :: def_
    integer, optional, intent(in)  :: func_type

    integer :: def, func_type_

    call push_sub('mix.mix_init')

    def = MIX_BROYDEN
    if(present(def_)) def = def_
    if(present(func_type)) then 
      func_type_=func_type
    else 
      func_type_=M_REAL
    end if

    !%Variable TypeOfMixing
    !%Type integer
    !%Default broyden
    !%Section SCF::Mixing
    !%Description
    !% The scheme scheme used to produce, at each iteration in the self consistent cycle
    !% that attempts to solve the Kohn-Sham equations, the input density from the value
    !% of the input and output densities of previous iterations.
    !%Option linear 0
    !% Simple linear mixing.
    !%Option gr_pulay 1
    !% "Guaranteed-reduction" Pulay scheme [D. R. Bowler and M. J. Gillan, Chem. Phys. 
    !% Lett. 325, 473 (2000)].
    !%Option broyden 2
    !% Broyden scheme [C. G Broyden, Math. Comp. 19, 577 (1965); 
    !% D. D. Johnson, Phys. Rev. B 38, 12807 (1988)].
    !%End
    call loct_parse_int(check_inp('TypeOfMixing'), def, smix%type_of_mixing)
    if(.not.varinfo_valid_option('TypeOfMixing', smix%type_of_mixing)) call input_error('TypeOfMixing')
    call messages_print_var_option(stdout, "TypeOfMixing", smix%type_of_mixing)

    !%Variable Mixing
    !%Type float
    !%Default 0.3
    !%Section SCF::Mixing
    !%Description
    !% Both the linear and the Broyden scheme depend on a "mixing parameter", set by this variable.
    !%End
    if (smix%type_of_mixing == MIX_LINEAR .or. smix%type_of_mixing == MIX_BROYDEN) then
      call loct_parse_float(check_inp('Mixing'), CNST(0.3), smix%alpha)
      if(smix%alpha <= M_ZERO .or. smix%alpha > M_ONE) call input_error('Mixing')
    end if

    !%Variable MixNumberSteps
    !%Type integer
    !%Default 3
    !%Section SCF::Mixing
    !%Description
    !% In the Broyden and in the GR-Pulay scheme, the new input density or potential is constructed
    !% from the values of densities/potentials of previous a given number of previous iterations.
    !% This number is set by this variable.
    !%End
    if (smix%type_of_mixing == MIX_GRPULAY .or. smix%type_of_mixing == MIX_BROYDEN) then
      call loct_parse_int(check_inp('MixNumberSteps'), 3, smix%ns)
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

    select case (smix%type_of_mixing)
    case (MIX_GRPULAY)
      if(func_type_ == M_REAL ) then 
        ALLOCATE(smix%ddf(m%np, d2, d3, smix%ns + 1), m%np*d2*d3*(smix%ns + 1))
        ALLOCATE(smix%dvin_old(m%np, d2, d3),         m%np*d2*d3)
        ALLOCATE(smix%ddv(m%np, d2, d3, smix%ns + 1), m%np*d2*d3*(smix%ns + 1))
        ALLOCATE(smix%df_old(m%np, d2, d3),           m%np*d2*d3)
        smix%ddf = M_ZERO; smix%ddv = M_ZERO; smix%dvin_old = M_ZERO; smix%df_old = M_ZERO
      else
        ALLOCATE(smix%zdf(m%np, d2, d3, smix%ns + 1), m%np*d2*d3*(smix%ns + 1))
        ALLOCATE(smix%zvin_old(m%np, d2, d3),         m%np*d2*d3)
        ALLOCATE(smix%zdv(m%np, d2, d3, smix%ns + 1), m%np*d2*d3*(smix%ns + 1))
        ALLOCATE(smix%zf_old(m%np, d2, d3),           m%np*d2*d3)
        smix%zdf = M_z0; smix%zdv = M_z0; smix%zvin_old = M_z0; smix%zf_old = M_z0
      end if

    case (MIX_BROYDEN)
      if(func_type_ == M_REAL ) then 
        ALLOCATE(smix%ddf(m%np, d2, d3, smix%ns), m%np*d2*d3*smix%ns)
        ALLOCATE(smix%dvin_old(m%np, d2, d3),     m%np*d2*d3)
        ALLOCATE(smix%ddv(m%np, d2, d3, smix%ns), m%np*d2*d3*smix%ns)
        ALLOCATE(smix%df_old(m%np, d2, d3),       m%np*d2*d3)
        smix%ddf = M_ZERO; smix%ddv = M_ZERO; smix%dvin_old = M_ZERO; smix%df_old = M_ZERO
      else
        ALLOCATE(smix%zdf(m%np, d2, d3, smix%ns), m%np*d2*d3*smix%ns)
        ALLOCATE(smix%zvin_old(m%np, d2, d3),     m%np*d2*d3)
        ALLOCATE(smix%zdv(m%np, d2, d3, smix%ns), m%np*d2*d3*smix%ns)
        ALLOCATE(smix%zf_old(m%np, d2, d3),       m%np*d2*d3)
        smix%zdf = M_z0; smix%zdv = M_z0; smix%zvin_old = M_z0; smix%zf_old = M_z0
      end if

    end select

    smix%last_ipos = 0

    call pop_sub()
  end subroutine mix_init


  ! ---------------------------------------------------------
  subroutine mix_end(smix)
    type(mix_t), intent(inout) :: smix

    call push_sub('mix.mix_end')

    ! Arrays got allocated for all mixing schemes, except linear mixing
    if (smix%type_of_mixing .ne. MIX_LINEAR) then
      if (associated(smix%ddf))      deallocate(smix%ddf)
      if (associated(smix%ddv))      deallocate(smix%ddv)
      if (associated(smix%dvin_old)) deallocate(smix%dvin_old)
      if (associated(smix%df_old))   deallocate(smix%df_old)

      if (associated(smix%zdf))      deallocate(smix%zdf)
      if (associated(smix%zdv))      deallocate(smix%zdv)
      if (associated(smix%zvin_old)) deallocate(smix%zvin_old)
      if (associated(smix%zf_old))   deallocate(smix%zf_old)
    end if

    call pop_sub()
  end subroutine mix_end

#include "undef.F90"
#include "real.F90"

#include "mix_inc.F90"

#include "undef.F90"
#include "complex.F90"

#include "mix_inc.F90"

end module mix_m
