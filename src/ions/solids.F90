!! Copyright (C) 2007 X. Andrade, M. Marques
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
!! $Id: simul_box.F90 3479 2007-11-09 08:36:10Z xavier $

#include "global.h"

module solids_m
  use geometry_m
  use global_m
  use io_m
  use messages_m
  use mpi_m
  use profiling_m
  use simul_box_m
  use unit_m
  use unit_system_m

  implicit none

  private

  public ::                   &
    periodic_copy_t,          &
    periodic_copy_init,       &
    periodic_copy_end,        &
    periodic_copy_position,   &
    periodic_copy_num,        &
    periodic_write_crystal

  !> parts of this module explicitly work on for 3 dimensions
  type periodic_copy_t
    private
    FLOAT :: pos(1:3)
    FLOAT :: pos_chi(1:3)
    FLOAT :: range
    integer :: nbmax(3), nbmin(3)
    integer, pointer :: icell(:, :)
  end type periodic_copy_t

contains

  subroutine periodic_copy_init(this, sb, pos, range)
    type(periodic_copy_t), intent(out) :: this
    type(simul_box_t),     intent(in)  :: sb
    FLOAT,                 intent(in)  :: pos(:)
    FLOAT,                 intent(in)  :: range

    integer :: pd, dim4syms
    integer :: icell1, icell2, icell3, jj

    PUSH_SUB(periodic_copy_init)

    ASSERT(range >= M_ZERO)

    dim4syms = min(3, sb%dim)

    this%range = range
    this%pos(1:dim4syms) = pos(1:dim4syms)
    nullify(this%icell)

    if(.not. simul_box_is_periodic(sb)) then
      this%nbmin = 0
      this%nbmax = 0
      POP_SUB(periodic_copy_init)
      return
    end if

    pd = sb%periodic_dim

    !convert the position to the orthogonal space
    this%pos_chi(1:pd) = matmul(pos(1:pd), sb%klattice_primitive(1:pd, 1:pd))

    this%nbmin(1:pd) = -int(-(this%pos_chi(1:pd) - range)/(M_TWO*sb%lsize(1:pd)) + M_HALF)
    this%nbmax(1:pd) = int((this%pos_chi(1:pd) + range)/(M_TWO*sb%lsize(1:pd)) + M_HALF)

    ! no copies in non-periodic directions
    this%nbmin(pd + 1:3) = 0
    this%nbmax(pd + 1:3) = 0

    SAFE_ALLOCATE(this%icell(1:3, 1:periodic_copy_num(this)))

    jj = 1
    do icell1 = this%nbmin(1), this%nbmax(1)
      do icell2 = this%nbmin(2), this%nbmax(2)
        do icell3 = this%nbmin(3), this%nbmax(3)
          this%icell(1:3, jj) = (/icell1, icell2, icell3/)
          jj = jj + 1
        end do
      end do
    end do

    POP_SUB(periodic_copy_init)
  end subroutine periodic_copy_init

  ! ----------------------------------------------------------------

  subroutine periodic_copy_end(this)
    type(periodic_copy_t), intent(inout) :: this

    PUSH_SUB(periodic_copy_end)

    SAFE_DEALLOCATE_P(this%icell)

    this%nbmin = 0
    this%nbmax = 0

    POP_SUB(periodic_copy_end)
  end subroutine periodic_copy_end

  ! ----------------------------------------------------------------

  integer pure function periodic_copy_num(this) result(num)
    type(periodic_copy_t), intent(in)    :: this

    ! no push_sub allowed in pure function
    num = product(this%nbmax - this%nbmin + 1)

  end function periodic_copy_num
  
  ! ----------------------------------------------------------------

  pure function periodic_copy_position(this, sb, ii) result(pcopy)
    type(periodic_copy_t),   intent(in)  :: this
    type(simul_box_t),       intent(in)  :: sb
    integer, intent(in)                  :: ii
    FLOAT                                :: pcopy(sb%dim)
    
    integer :: pd

    pd = sb%periodic_dim

    if(.not. simul_box_is_periodic(sb)) then
      pcopy = this%pos
      return
    end if

    pcopy(1:pd) = this%pos_chi(1:pd) - M_TWO*sb%lsize(1:pd)*this%icell(1:pd, ii)
    pcopy(1:pd) = matmul(sb%rlattice_primitive(1:pd, 1:pd), pcopy(1:pd))
    pcopy(pd + 1:sb%dim) = this%pos(pd+1:sb%dim)

  end function periodic_copy_position

  ! ----------------------------------------------------------------
  !> This subroutine creates a crystal by replicating the geometry and
  !! writes the result to dir//'crystal.xyz'
  subroutine periodic_write_crystal(sb, geo, dir) 
    type(simul_box_t), intent(in) :: sb
    type(geometry_t),  intent(in) :: geo 
    character(len=*),  intent(in) :: dir
    
    type(periodic_copy_t) :: pp
    FLOAT :: radius, pos(1:MAX_DIM)
    integer :: total_atoms, iatom, icopy, iunit

    PUSH_SUB(periodic_write_crystal)
    
    radius = maxval(sb%lsize)*(M_ONE + M_EPSILON)
    
    !count the number of atoms in the crystal
    total_atoms = 0
    do iatom = 1, geo%natoms
      call periodic_copy_init(pp, sb, geo%atom(iatom)%x, radius)
      total_atoms = total_atoms + periodic_copy_num(pp)
      call periodic_copy_end(pp)
    end do
    
    ! now calculate
    if(mpi_grp_is_root(mpi_world)) then
      
      iunit = io_open(trim(dir)//'/crystal.xyz', action='write')
      
      write(iunit, '(i9)') total_atoms
      write(iunit, '(a)') '#generated by Octopus'
      
      do iatom = 1, geo%natoms
        call periodic_copy_init(pp, sb, geo%atom(iatom)%x, radius)
        do icopy = 1, periodic_copy_num(pp)
          pos(1:sb%dim) = units_from_atomic(units_out%length, periodic_copy_position(pp, sb, icopy))
          write(iunit, '(a, 99f12.6)') geo%atom(iatom)%label, pos(1:sb%dim)
          
        end do
        call periodic_copy_end(pp)
      end do

      call io_close(iunit)
    end if

    POP_SUB(periodic_write_crystal)
  end subroutine periodic_write_crystal

end module solids_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
