!! Copyright (C) 2009 D. Strubbe
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
!! $Id: geometry.F90 5547 2009-06-03 12:27:56Z acastro $

#include "global.h"

module born_charges_m
  use datasets_m
  use geometry_m
  use global_m
  use io_m
  use io_function_m
  use loct_parser_m
  use messages_m
  use profiling_m
  use species_m
  use states_m

  implicit none

  private
  public ::                &
    Born_charges_t,        &
    Born_charges_init,     &
    Born_charges_end,      &
    out_Born_charges

  type Born_charges_t
    CMPLX, pointer :: charge(:, :, :)    ! i, j, atom: Z*(i,j) = dF(j)/dE(i) = dP(i) / dR(j)
    CMPLX :: sum_ideal(MAX_DIM, MAX_DIM) ! the sum of Born charges according to acoustic sum rule 
    CMPLX :: delta(MAX_DIM, MAX_DIM)     ! discrepancy of sum of Born charge tensors from sum rule
    logical :: correct                   ! correct according to sum rule?
  end type Born_charges_t

  contains

  ! ---------------------------------------------------------
  subroutine Born_charges_init(this, geo, st, dim)
    type(Born_charges_t), intent(out) :: this
    type(geometry_t),     intent(in)  :: geo
    type(states_t),       intent(in)  :: st
    integer,              intent(in)  :: dim

    integer idir

    call push_sub('Born_charges.Born_charges_init')

    nullify(this%charge)
    SAFE_ALLOCATE(this%charge(1:dim, 1:dim, 1:geo%natoms))
    this%charge(1:dim, 1:dim, 1:geo%natoms) = M_ZERO
    this%delta(1:dim, 1:dim) = M_ZERO

    this%sum_ideal(1:dim, 1:dim) = M_ZERO
    do idir = 1, dim
      this%sum_ideal(idir, idir) = -(st%val_charge + st%qtot) ! total charge
    enddo

    !%Variable BornChargeSumRuleCorrection
    !%Type logical
    !%Default true
    !%Section Linear Response::Polarizabilities
    !%Description
    !% Enforce the acoustic sum rule by distributing the excess sum of Born charges equally among the atoms.
    !% sum rule: sum(iatom) Z*(iatom,idir,idir2) = Z_tot delta(idir1, idir2)
    !%End

    call loct_parse_logical(datasets_check('BornChargeSumRuleCorrection'), .true., this%correct)

    call pop_sub()
  end subroutine Born_charges_init

  ! ---------------------------------------------------------
  subroutine Born_charges_end(this)
    type(Born_charges_t), intent(out) :: this

    call push_sub('Born_charges.Born_charges_end')

    SAFE_DEALLOCATE_P(this%charge)

    call pop_sub()
  end subroutine Born_charges_end

  ! ---------------------------------------------------------
  ! The sum over atoms of a given tensor component of the Born charges
  !  should be Z delta_ij to satisfy the acoustic sum rule, where Z is total charge of system
  subroutine correct_Born_charges(this, geo, dim)
    type(Born_charges_t), intent(inout) :: this
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(in)    :: dim

    CMPLX :: Born_sum(MAX_DIM, MAX_DIM)        ! the sum of Born charges from the calculation 
    integer iatom, idir

    call push_sub('Born_charges.correct_Born_charges')

    Born_sum(1:dim, 1:dim) = M_ZERO 

    do iatom = 1, geo%natoms
      Born_sum(1:dim, 1:dim) = Born_sum(1:dim, 1:dim) + this%charge(1:dim, 1:dim, iatom)
    enddo

    this%delta(1:dim, 1:dim) = Born_sum(1:dim, 1:dim) - this%sum_ideal(1:dim, 1:dim)

    if(this%correct) then
      do iatom = 1, geo%natoms
        this%charge(1:dim, 1:dim, iatom) = &
          this%charge(1:dim, 1:dim, iatom) - this%delta(1:dim, 1:dim) / geo%natoms
      enddo
    endif

    call pop_sub()
  end subroutine correct_Born_charges

  ! ---------------------------------------------------------
  subroutine out_Born_charges(this, geo, dim, dirname)
    type(Born_charges_t), intent(inout) :: this
    type(geometry_t),     intent(in)    :: geo
    integer,              intent(in)    :: dim
    character(len=*),     intent(in)    :: dirname

    integer iatom, iunit
    FLOAT :: phase(1:MAX_DIM, 1:MAX_DIM)

    call push_sub('Born_charges.out_Born_charges')

    call correct_Born_charges(this, geo, dim)

    iunit = io_open(trim(dirname)//'/Born_charges', action='write')
    write(iunit,'(a)') '# (Frequency-dependent) Born effective charge tensors'
    write(iunit,'(a)') '# Real and imaginary parts'
    do iatom = 1, geo%natoms
      write(iunit,'(a,i5,a,a5,a,f10.4)') 'Index: ', iatom, '   Label: ', trim(species_label(geo%atom(iatom)%spec)), &
        '   Ionic charge: ', species_zval(geo%atom(iatom)%spec)

      write(iunit,'(a)') 'Real:'
      call io_output_tensor(iunit, real(this%charge(:, :, iatom)), dim, M_ONE)

      write(iunit,'(a)') 'Imaginary:'
      call io_output_tensor(iunit, aimag(this%charge(:, :, iatom)), dim, M_ONE)
      write(iunit,'(a)')
    enddo

    write(iunit,'(a)') '# Magnitude and phase'
    do iatom = 1, geo%natoms
      write(iunit,'(a,i5,a,a5,a,f10.4)') 'Index: ', iatom, '   Label: ', trim(species_label(geo%atom(iatom)%spec)), &
        '   Ionic charge: ', species_zval(geo%atom(iatom)%spec)

      write(iunit,'(a)') 'Magnitude:'
      call io_output_tensor(iunit, TOFLOAT(abs(this%charge(:, :, iatom))), dim, M_ONE)

      write(iunit,'(a)') 'Phase:'
      phase(1:dim, 1:dim) = atan2(aimag(this%charge(1:dim, 1:dim, iatom)), real(this%charge(1:dim, 1:dim, iatom)))
      call io_output_tensor(iunit, phase(:, :), dim, M_ONE, write_average = .false.)
      write(iunit,'(a)')
    enddo

    write(iunit,'(a)') '# Discrepancy of Born effective charges from acoustic sum rule before correction' 
    write(iunit,'(a)') 'Real:'
    call io_output_tensor(iunit, real(this%delta(:, :)), dim, M_ONE)
    write(iunit,'(a)') 'Imaginary:'
    call io_output_tensor(iunit, aimag(this%delta(:, :)), dim, M_ONE)

    call io_close(iunit)
    call pop_sub()
  end subroutine out_Born_charges

end module Born_charges_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
