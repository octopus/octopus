!! Copyright (C) 2019 F. Buchholz, M. Oliveira
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

module dressed_interaction_oct_m
  use global_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                      &
    dressed_interaction_t,       &
    dressed_init,                &
    dressed_write_info,          &
    dressed_add_poisson_terms

  type dressed_interaction_t
    private
    integer :: dim
    FLOAT :: lambda(1:MAX_DIM - 1)
    FLOAT, public :: omega
    FLOAT :: n_electrons
    FLOAT :: coulomb
  end type dressed_interaction_t

contains

  !-----------------------------------------------------------------
  subroutine dressed_init(this, namespace, dim, n_electrons)
    type(dressed_interaction_t), intent(out) :: this
    type(namespace_t),           intent(in)  :: namespace
    integer,                     intent(in)  :: dim
    FLOAT,                       intent(in)  :: n_electrons

    integer :: idir
    type(block_t) :: blk

    PUSH_SUB(dressed_init)

    this%dim = dim - 1
    this%n_electrons = n_electrons

    !%Variable DressedLambda
    !%Type block
    !%Section SCF::RDMFT
    !%Description
    !% Polarization vector including the interaction strength in dressed orbital formalism,
    !% in units of energy. The default is zero for all components.
    !%End
    this%lambda = M_ZERO
    if (parse_block(namespace, 'DressedLambda', blk) == 0) then
      if (parse_block_cols(blk, 0) /= this%dim) then
        call messages_input_error(' DressedLambda')
      end if
      do idir = 1, this%dim
        call parse_block_float(blk, 0, idir - 1, this%lambda(idir), units_inp%energy)
      end do
    end if
    call parse_block_end(blk)
    call messages_print_var_value(stdout, 'DressedLambda', this%lambda, unit = units_out%energy)

    !%Variable DressedOmega
    !%Type float
    !%Default 1.0 Ha
    !%Section SCF::RDMFT
    !%Description
    !% mode frequency in dressed orbital formalism.
    !%End
    call parse_variable(namespace, 'DressedOmega', CNST(1.0), this%omega, units_inp%energy)
    call messages_print_var_value(stdout, 'DressedOmega', this%omega, unit = units_out%energy)

    !%Variable DressedCoulomb
    !%Type float
    !%Default 1.0
    !%Section SCF::RDMFT
    !%Description
    !% allows to control the prefactor of the electron electron interaction
    !%End
    call parse_variable(namespace, 'DressedCoulomb', CNST(1.0), this%coulomb)
    call messages_print_var_value(stdout, 'DressedCoulomb', this%coulomb)

    POP_SUB(dressed_init)
  end subroutine dressed_init

  !-----------------------------------------------------------------
  subroutine dressed_write_info(this, iunit)
    type(dressed_interaction_t), intent(in) :: this
    integer,                     intent(in) :: iunit

    integer :: idir

    PUSH_SUB(dressed_write_info)

    write(iunit, '(a,1x)', advance='no') 'DressedLambda:   '
    write(iunit, '(f14.12)') (this%lambda(idir), idir = 1, this%dim)
    write(iunit, '(a,1x,f14.12)') 'DressedOmega:    ', this%omega
    write(iunit, '(a,1x,f14.12)') 'DressedCoulomb:  ', this%coulomb

    POP_SUB(dressed_write_info)
  end subroutine dressed_write_info

  !-----------------------------------------------------------------
  subroutine dressed_add_poisson_terms(this, mesh, rho, pot)
    type(dressed_interaction_t), intent(in)    :: this
    type(mesh_t),                intent(in)    :: mesh
    FLOAT,                       intent(in)    :: rho(:)
    FLOAT,                       intent(inout) :: pot(:)

    integer :: ip, dim_ele
    FLOAT :: lx, ld, dipole(1:MAX_DIM)

    PUSH_SUB(dressed_add_poisson_terms)

    dim_ele = mesh%sb%dim - 1

    call dmf_dipole(mesh, rho, dipole)
    ld = dot_product(this%lambda(1:dim_ele), dipole(1:dim_ele))

    do ip = 1, mesh%np
      lx = dot_product(this%lambda(1:dim_ele), mesh%x(ip, 1:dim_ele))
      pot(ip) = pot(ip)*this%coulomb - &
        this%omega/sqrt(this%n_electrons)*(mesh%x(ip, dim_ele + 1)*ld + lx*dipole(dim_ele + 1)) + lx*ld
    end do

    POP_SUB(dressed_add_poisson_terms)
  end subroutine dressed_add_poisson_terms

end module dressed_interaction_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
