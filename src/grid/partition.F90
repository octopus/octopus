!! Copyright (C) 2010 X. Andrade
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
!! $Id: partitioner.F90 6396 2010-03-26 08:51:58Z mjv500 $

#include "global.h"

module partition_m
  use batch_m
  use c_pointer_m
  use datasets_m
  use global_m
  use iihash_m
  use index_m
  use io_m
  use loct_math_m
  use mesh_m
  use messages_m
  use parser_m
  use mpi_m
  use profiling_m
  use stencil_m

  implicit none

  private

  public ::                &
    partition_t,           &
    partition_init,        &
    partition_end,         &
    partition_build,       &
    partition_quality,     &
    partition_write_info

  type partition_t
    integer, pointer :: nghost(:)
    integer, pointer :: nbound(:)
    integer, pointer :: nlocal(:)
    integer, pointer :: nneigh(:)
    integer          :: npart
    integer          :: npoints
    integer          :: library
  end type partition_t

  integer, parameter, public :: &
       METIS     = 2,           &
       ZOLTAN    = 3,           &
       PFFT_PART = 5,           &
       PARMETIS  = 6

contains

  subroutine partition_init(this, mesh)
    type(partition_t), intent(out) :: this
    type(mesh_t),      intent(in)  :: mesh

    integer :: default
    
    PUSH_SUB(partition_init)
    
    this%npart = mesh%mpi_grp%size
    this%npoints = mesh%np_part_global

    default = ZOLTAN
#ifdef HAVE_METIS
    default = METIS
#endif
    ! this variable is documented in src/grid/mesh_partition.F90
    call parse_integer(datasets_check('MeshPartitionPackage'), default, this%library)

    SAFE_ALLOCATE(this%nghost(1:this%npart))
    SAFE_ALLOCATE(this%nbound(1:this%npart))
    SAFE_ALLOCATE(this%nlocal(1:this%npart))
    SAFE_ALLOCATE(this%nneigh(1:this%npart))

    POP_SUB(partition_init)
  end subroutine partition_init
  ! ----------------------------------------------------------------------

  subroutine partition_build(this, mesh, stencil, point_to_part)
    type(partition_t), intent(out) :: this
    type(mesh_t),      intent(in)  :: mesh
    type(stencil_t),   intent(in)  :: stencil
    integer,           intent(in)  :: point_to_part(:)

    integer :: ip, ipcoords(1:MAX_DIM)
    integer, allocatable :: jpcoords(:, :), jp(:)
    integer :: istencil, ipart, jpart
    type(profile_t), save :: prof
    logical, allocatable :: is_a_neigh(:, :), gotit(:)

    call profiling_in(prof, "PARTITION_BUILD")
    PUSH_SUB(partition_build)

    SAFE_ALLOCATE(is_a_neigh(1:this%npart, 1:this%npart))

    is_a_neigh = .false.

    SAFE_ALLOCATE(gotit(1:mesh%np_part_global))

    SAFE_ALLOCATE(jpcoords(1:MAX_DIM, 1:stencil%size))
    SAFE_ALLOCATE(jp(1:stencil%size))

    this%nghost = 0
    this%nbound = 0
    this%nlocal = 0
    this%nneigh = 0

    do ipart = 1, this%npart
      gotit = .false.
      do ip = 1, mesh%np_global
        if(ipart /= point_to_part(ip)) cycle

        INCR(this%nlocal(ipart), 1)
        call index_to_coords(mesh%idx, mesh%sb%dim, ip, ipcoords)
        
        do istencil = 1, stencil%size
          jpcoords(:, istencil) = ipcoords + stencil%points(:, istencil)
        end do
        
        call index_from_coords_vec(mesh%idx, mesh%sb%dim, stencil%size, jpcoords, jp)
        
        do istencil = 1, stencil%size
          if(stencil%center == istencil) cycle

          if(.not. gotit(jp(istencil))) then
            jpart = point_to_part(jp(istencil))
         
            if(jpart /= ipart) then
              INCR(this%nghost(ipart), 1)
              is_a_neigh(ipart, jpart) = .true.
            else if(jp(istencil) > mesh%np_global) then
              INCR(this%nbound(ipart), 1)
            end if
            
            gotit(jp(istencil)) = .true.
          end if
          
        end do
        
      end do
    end do

    forall(ipart = 1:this%npart)
      this%nneigh(ipart) = count(is_a_neigh(ipart, 1:this%npart))
    end forall

    SAFE_DEALLOCATE_A(is_a_neigh)
    SAFE_DEALLOCATE_A(gotit)
    SAFE_DEALLOCATE_A(jpcoords)
    SAFE_DEALLOCATE_A(jp)

    POP_SUB(partition_build)
    call profiling_out(prof)
  end subroutine partition_build

  ! ----------------------------------------------------------------------

  subroutine partition_write_info(this)
    type(partition_t), intent(in) :: this
    
    integer :: ipart

    PUSH_SUB(partition_write_info)

    ! Write information about partitions.
    message(1) = &
      'Info: Mesh partition:'
    message(2) = ''
    call messages_info(2)

    write(message(1),'(a,e16.6)') &
      '      Partition quality:', partition_quality(this)
    message(2) = ''
    call messages_info(2)

    write(message(1),'(a)') &
      '                 Neighbours         Ghost points'
    write(message(2),'(a,i5,a,i10)') &
      '      Average  :      ', sum(this%nneigh)/this%npart, '           ', sum(this%nghost)/this%npart
    write(message(3),'(a,i5,a,i10)') &
      '      Minimum  :      ', minval(this%nneigh),        '           ', minval(this%nghost)
    write(message(4),'(a,i5,a,i10)') &
      '      Maximum  :      ', maxval(this%nneigh),        '           ', maxval(this%nghost)
    message(5) = ''
    call messages_info(5)

    do ipart = 1, this%npart
      write(message(1),'(a,i5)')  &
        '      Nodes in domain-group  ', ipart
      write(message(2),'(a,i10,a,i10)') &
        '        Neighbours     :', this%nneigh(ipart), &
        '        Local points    :', this%nlocal(ipart)
      write(message(3),'(a,i10,a,i10)') &
        '        Ghost points   :', this%nghost(ipart), &
        '        Boundary points :', this%nbound(ipart)
      call messages_info(3)
    end do

    message(1) = ''
    call messages_info(1)

    POP_SUB(partition_write_info)
  end subroutine partition_write_info

  ! ----------------------------------------------------------------------

  subroutine partition_end(this)
    type(partition_t), intent(inout) :: this

    PUSH_SUB(partition_end)

    SAFE_DEALLOCATE_P(this%nghost)
    SAFE_DEALLOCATE_P(this%nbound)
    SAFE_DEALLOCATE_P(this%nlocal)
    SAFE_DEALLOCATE_P(this%nneigh)

    POP_SUB(partition_end)
  end subroutine partition_end

  ! ----------------------------------------------------------------------

  FLOAT pure function partition_quality(this) result(quality)
    type(partition_t), intent(in) :: this
    
    FLOAT :: scal

    scal = real(this%npart, REAL_PRECISION)/this%npoints

    quality = M_ZERO

    quality = quality + (maxval(this%nlocal) - minval(this%nlocal))**3
    quality = quality + (sum(TOFLOAT(this%nghost)**2))

    quality = M_ONE/(M_ONE + quality)

  end function partition_quality

end module partition_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
