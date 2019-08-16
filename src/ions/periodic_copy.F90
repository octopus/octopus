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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module periodic_copy_oct_m
  use geometry_oct_m
  use global_oct_m
  use io_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private

  public ::                   &
    periodic_copy_t,          &
    periodic_copy_init,       &
    periodic_copy_end,        &
    periodic_copy_position,   &
    periodic_copy_num,        &
    periodic_write_crystal

  type periodic_copy_t
    private
    integer :: num
    FLOAT :: pos(1:MAX_DIM)
    FLOAT :: pos_chi(1:MAX_DIM)
    FLOAT :: range
    integer :: nbmax(1:MAX_DIM), nbmin(1:MAX_DIM)
    integer, pointer :: icell(:, :) !< (sb%dim, num)
  end type periodic_copy_t

contains

  subroutine periodic_copy_init(this, sb, pos, range)
    type(periodic_copy_t), intent(out) :: this
    type(simul_box_t),     intent(in)  :: sb
    FLOAT,                 intent(in)  :: pos(:) !< (sb%dim)
    FLOAT,                 intent(in)  :: range

    integer :: pd, jj, kk, idir

    PUSH_SUB(periodic_copy_init)

    ASSERT(range >= M_ZERO)

    this%range = range
    this%pos(1:sb%dim) = pos(1:sb%dim)

    if(.not. simul_box_is_periodic(sb)) then
      this%num = 1
      this%nbmin = 0
      this%nbmax = 0
      nullify(this%icell)

      POP_SUB(periodic_copy_init)
      return
    end if

    pd = sb%periodic_dim

    !convert the position to the orthogonal space
    this%pos_chi(1:pd) = matmul(pos(1:pd), sb%klattice_primitive(1:pd, 1:pd))

    this%nbmin(1:pd) = -nint(-(this%pos_chi(1:pd) - range)/(M_TWO*sb%lsize(1:pd)) + M_HALF)
    this%nbmax(1:pd) = nint((this%pos_chi(1:pd) + range)/(M_TWO*sb%lsize(1:pd)) + M_HALF)
    ! no copies in non-periodic directions

    this%num = product(this%nbmax(1:sb%periodic_dim) - this%nbmin(1:sb%periodic_dim) + 1)
    SAFE_ALLOCATE(this%icell(1:sb%periodic_dim, 1:this%num))

    do jj = 1, this%num
      kk = jj - 1
      do idir = sb%periodic_dim, 1, -1
        this%icell(idir, jj) = mod(kk, this%nbmax(idir) - this%nbmin(idir) + 1) + this%nbmin(idir)
        if(idir > 1) &
          kk = kk / (this%nbmax(idir) - this%nbmin(idir) + 1)
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
    num = this%num

  end function periodic_copy_num
  
  ! ----------------------------------------------------------------

  pure function periodic_copy_position(this, sb, ii) result(pcopy)
    type(periodic_copy_t),   intent(in)  :: this
    type(simul_box_t),       intent(in)  :: sb
    integer,                 intent(in)  :: ii
    FLOAT                                :: pcopy(sb%dim)
    
    integer :: pd

    pd = sb%periodic_dim

    if(.not. simul_box_is_periodic(sb)) then
      pcopy(1:sb%dim) = this%pos(1:sb%dim)
      return
    end if

    pcopy(1:pd) = this%pos_chi(1:pd) - M_TWO*sb%lsize(1:pd)*this%icell(1:pd, ii)
    pcopy(1:pd) = matmul(sb%rlattice_primitive(1:pd, 1:pd), pcopy(1:pd))
    pcopy(pd + 1:sb%dim) = this%pos(pd+1:sb%dim)

  end function periodic_copy_position

  ! ----------------------------------------------------------------
  !> This subroutine creates a crystal by replicating the geometry and
  !! writes the result to dir//'crystal.xyz'
  subroutine periodic_write_crystal(sb, geo, dir, namespace)
    type(simul_box_t), intent(in) :: sb
    type(geometry_t),  intent(in) :: geo 
    character(len=*),  intent(in) :: dir
    type(namespace_t), intent(in) :: namespace 
    
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
      
      iunit = io_open(trim(dir)//'/crystal.xyz', namespace, action='write')
      
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

end module periodic_copy_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
