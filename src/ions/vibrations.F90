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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module vibrations_oct_m
  use global_oct_m
  use io_oct_m
  use ions_oct_m
  use lalg_adv_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use profiling_oct_m
  use species_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public :: &
    vibrations_t,                     &
    vibrations_init,                  &
    vibrations_end,                   &
    vibrations_symmetrize_dyn_matrix, &
    vibrations_normalize_dyn_matrix,  &
    vibrations_out_dyn_matrix_row,    &
    vibrations_out_dyn_matrix_header, &
    vibrations_norm_factor,           &
    vibrations_diag_dyn_matrix,       &
    vibrations_get_index,             &
    vibrations_get_atom,              &
    vibrations_get_dir,               &
    vibrations_output,                &
    vibrations_get_suffix
  
  type vibrations_t
    private
    integer,            public :: num_modes
    integer,            public :: ndim
    integer,            public :: natoms
    FLOAT, allocatable, public :: dyn_matrix(:,:)
    FLOAT, allocatable, public :: infrared(:,:)
    FLOAT, allocatable, public :: normal_mode(:,:)
    FLOAT, allocatable, public :: freq(:)
    FLOAT,              public :: disp
    FLOAT :: total_mass
    character (len=2) :: suffix
    character (len=80) :: filename_dynmat
    type(unit_t) :: unit_dynmat
    type(namespace_t), pointer :: namespace
  end type vibrations_t

contains

  ! ---------------------------------------------------------
  subroutine vibrations_init(this, ions, suffix, namespace)
    type(vibrations_t),        intent(out) :: this
    type(ions_t),              intent(in)  :: ions
    character (len=2),         intent(in)  :: suffix
    type(namespace_t), target, intent(in)  :: namespace

    PUSH_SUB(vibrations_init)

    this%ndim = ions%space%dim
    this%natoms = ions%natoms
    this%num_modes = ions%natoms*ions%space%dim
    this%namespace => namespace
    SAFE_ALLOCATE(this%dyn_matrix(1:this%num_modes, 1:this%num_modes))
    SAFE_ALLOCATE(this%infrared(1:this%num_modes, 1:this%ndim))
    SAFE_ALLOCATE(this%normal_mode(1:this%num_modes, 1:this%num_modes))
    SAFE_ALLOCATE(this%freq(1:this%num_modes))

    this%total_mass = sum(ions%mass)

    ! Since frequencies are reported as invcm, the matrix they are derived from can be expressed in invcm**2.
    ! However, that matrix is the dynamical matrix divided by the total mass, so it has a different unit.
    this%unit_dynmat = units_out%energy / units_out%length**2

    this%suffix = suffix
    this%filename_dynmat = VIB_MODES_DIR//'dynamical_matrix_'//trim(this%suffix)
    if(mpi_grp_is_root(mpi_world)) then
      call io_mkdir(VIB_MODES_DIR, namespace)
      call io_rm(this%filename_dynmat, namespace)
      call vibrations_out_dyn_matrix_header(this)
    end if

    POP_SUB(vibrations_init)
  end subroutine vibrations_init


  ! ---------------------------------------------------------
  subroutine vibrations_end(this)
    type(vibrations_t), intent(inout) :: this

    PUSH_SUB(vibrations_end)

    SAFE_DEALLOCATE_A(this%dyn_matrix)
    SAFE_DEALLOCATE_A(this%infrared)
    SAFE_DEALLOCATE_A(this%freq)
    SAFE_DEALLOCATE_A(this%normal_mode)

    POP_SUB(vibrations_end)
  end subroutine vibrations_end


  ! ---------------------------------------------------------
  character(len=2) function vibrations_get_suffix(this)
    type(vibrations_t), intent(in) :: this

    PUSH_SUB(vibrations_get_suffix)
    vibrations_get_suffix = this%suffix

    POP_SUB(vibrations_get_suffix)
  end function vibrations_get_suffix


  ! ---------------------------------------------------------
  subroutine vibrations_symmetrize_dyn_matrix(this)
    type(vibrations_t), intent(inout) :: this

    integer :: imat, jmat
    FLOAT :: average, maxdiff

    PUSH_SUB(vibrations_symmetrize_dyn_matrix)

    ! FIXME: enforce acoustic sum rule.
    
    maxdiff = M_ZERO
    do imat = 1, this%num_modes
      do jmat = imat + 1, this%num_modes
        average = M_HALF * (this%dyn_matrix(imat, jmat) + this%dyn_matrix(jmat, imat))
        maxdiff = max(maxdiff, abs(this%dyn_matrix(imat, jmat) - this%dyn_matrix(jmat, imat)))
        this%dyn_matrix(imat, jmat) = average
        this%dyn_matrix(jmat, imat) = average
      end do
    end do

    write(message(1),'(a)') 'Info: Symmetrizing dynamical matrix.'
    write(message(2),'(a,es12.6,a,a)') 'Info: Maximum discrepancy from symmetry: ', &
      units_from_atomic(this%unit_dynmat, maxdiff), &
      " ", trim(units_abbrev(this%unit_dynmat))
    call messages_info(2)

    POP_SUB(vibrations_symmetrize_dyn_matrix)
  end subroutine vibrations_symmetrize_dyn_matrix

  ! ---------------------------------------------------------
  subroutine vibrations_normalize_dyn_matrix(this, ions)
    type(vibrations_t), intent(inout) :: this
    type(ions_t),       intent(in)    :: ions

    integer :: iatom, idir, jatom, jdir, imat, jmat

    PUSH_SUB(vibrations_normalize_dyn_matrix)

    do iatom = 1, this%natoms
      do idir = 1, this%ndim
        
        imat = vibrations_get_index(this, iatom, idir)

        do jatom = 1, this%natoms
          do jdir = 1, this%ndim
            
            jmat = vibrations_get_index(this, jatom, jdir)
            
            this%dyn_matrix(jmat, imat) = &
              this%dyn_matrix(jmat, imat) * vibrations_norm_factor(this, ions, iatom, jatom)

          end do
        end do

      end do
    end do

    POP_SUB(vibrations_normalize_dyn_matrix)
  end subroutine vibrations_normalize_dyn_matrix

  ! ---------------------------------------------------------
  FLOAT pure function vibrations_norm_factor(this, ions, iatom, jatom)
    type(vibrations_t), intent(in) :: this
    type(ions_t),       intent(in) :: ions
    integer,            intent(in) :: iatom
    integer,            intent(in) :: jatom

    vibrations_norm_factor = this%total_mass / &
      sqrt(ions%mass(iatom) * ions%mass(jatom))

  end function vibrations_norm_factor

  ! ---------------------------------------------------------
  subroutine vibrations_out_dyn_matrix_row(this, imat)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: imat

    integer :: iatom, idir, jatom, jdir, jmat, iunit

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(vibrations_out_dyn_matrix_row)

    iatom = vibrations_get_atom(this, imat)
    idir  = vibrations_get_dir (this, imat)

    iunit = io_open(this%filename_dynmat, this%namespace, action='write', position='append')

    do jmat = 1, this%num_modes
      jatom = vibrations_get_atom(this, jmat)
      jdir  = vibrations_get_dir (this, jmat)
    
      write(iunit, '(2(i8, i6), e25.12)') &
        jatom, jdir, iatom, idir, units_from_atomic(this%unit_dynmat, this%dyn_matrix(jmat, imat))
    end do
    call io_close(iunit)

    POP_SUB(vibrations_out_dyn_matrix_row)
  end subroutine vibrations_out_dyn_matrix_row

  ! ---------------------------------------------------------
  subroutine vibrations_out_dyn_matrix_header(this)
    type(vibrations_t), intent(in) :: this

    integer :: iunit

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(vibrations_out_dyn_matrix_header)

    iunit = io_open(this%filename_dynmat, namespace=this%namespace, action='write') ! start at the beginning
    write(iunit, '(2(a8, a6), a25)') 'atom', 'dir', 'atom', 'dir', &
      '[' // trim(units_abbrev(this%unit_dynmat)) // ']'
    call io_close(iunit)

    POP_SUB(vibrations_out_dyn_matrix_header)
  end subroutine vibrations_out_dyn_matrix_header

  ! ---------------------------------------------------------
  subroutine vibrations_diag_dyn_matrix(this)
    type(vibrations_t), intent(inout) :: this
    
    integer :: imode

    PUSH_SUB(vibrations_diag_dyn_matrix)

    this%normal_mode = M_ZERO
    
    this%normal_mode = this%dyn_matrix
    call lalg_eigensolve(this%num_modes, this%normal_mode, this%freq)

    ! FIXME: why not remove total mass here and in norm_factor?
    this%freq(1:this%num_modes) = -this%freq(1:this%num_modes) / this%total_mass

    if(any(this%freq(1:this%num_modes) < -M_EPSILON)) then
      message(1) = "There are imaginary vibrational frequencies (represented as negative)."
      call messages_warning(1, namespace=this%namespace)
    end if

    do imode = 1, this%num_modes
      if(this%freq(imode) > M_EPSILON) then
        this%freq(imode) =  sqrt(abs(this%freq(imode)))
      else
        this%freq(imode) = -sqrt(abs(this%freq(imode)))
      end if

      ! make the largest component positive, to specify the phase
      if( maxval(this%normal_mode(:, imode)) - abs(minval(this%normal_mode(:, imode))) < -M_EPSILON) then
        this%normal_mode(:, imode) = -this%normal_mode(:, imode)
      end if
    end do

    POP_SUB(vibrations_diag_dyn_matrix)
  end subroutine vibrations_diag_dyn_matrix


  ! ---------------------------------------------------------
  integer pure function vibrations_get_index(this, iatom, idim)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: iatom
    integer,            intent(in) :: idim

    vibrations_get_index = (iatom - 1)*this%ndim + idim
  end function vibrations_get_index


  ! ---------------------------------------------------------
  integer pure function vibrations_get_atom(this, index)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: index

    vibrations_get_atom = 1 + (index - 1)/ this%ndim 
  end function vibrations_get_atom


  ! ---------------------------------------------------------
  integer pure function vibrations_get_dir(this, index)
    type(vibrations_t), intent(in) :: this
    integer,            intent(in) :: index

    vibrations_get_dir =  1 + mod(index - 1, this%ndim)
  end function vibrations_get_dir


  ! ---------------------------------------------------------
  subroutine vibrations_output(this)
    type(vibrations_t), intent(in) :: this
    
    integer :: iunit, imat, jmat

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(vibrations_output)

    ! output frequencies and eigenvectors
    iunit = io_open(VIB_MODES_DIR//'normal_frequencies_'//trim(this%suffix), this%namespace, action='write')
    do imat = 1, this%num_modes
      write(iunit, '(i6,f17.8)') imat, units_from_atomic(unit_invcm, this%freq(imat))
    end do
    call io_close(iunit)

    ! output eigenvectors
    iunit = io_open(VIB_MODES_DIR//'normal_modes_'//trim(this%suffix), this%namespace, action='write')
    do imat = 1, this%num_modes
      write(iunit, '(i6)', advance='no') imat
      do jmat = 1, this%num_modes
        write(iunit, '(es14.5)', advance='no') this%normal_mode(jmat, imat)
      end do
      write(iunit, '(1x)')
    end do
    call io_close(iunit)

    POP_SUB(vibrations_output)
  end subroutine vibrations_output

end module vibrations_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
