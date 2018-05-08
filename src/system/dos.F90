!! Copyright (C) 2017 N. Tancogne-Dejean
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

module dos_oct_m
  use comm_oct_m
  use geometry_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use io_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use species_oct_m
  use states_oct_m
  use states_dim_oct_m
  use submesh_oct_m
  use unit_oct_m
  use unit_system_oct_m


  implicit none

  private
  public ::                    &
    dos_t,                     &
    dos_init,                  &
    dos_end,                   & 
    dos_write_dos

  type dos_t
    FLOAT   :: emin
    FLOAT   :: emax
    integer :: epoints
    FLOAT   :: gamma
    FLOAT   :: de
  end type dos_t

contains

  subroutine dos_init(this, st)
    type(dos_t), intent(out)   :: this
    type(states_t), intent(in) :: st

    FLOAT :: evalmin, evalmax, eextend

    PUSH_SUB(dos_init)

    evalmin = minval(st%eigenval)
    evalmax = maxval(st%eigenval)
    ! we extend the energy mesh by this amount
    eextend  = (evalmax - evalmin) / M_FOUR

    !%Variable DOSEnergyMin
    !%Type float
    !%Section Output
    !%Description
    !% Lower bound for the energy mesh of the DOS.
    !% The default is the lowest eigenvalue, minus a quarter of the total range of eigenvalues.
    !%End
    call parse_variable('DOSEnergyMin', evalmin - eextend, this%emin, units_inp%energy)

    !%Variable DOSEnergyMax
    !%Type float
    !%Section Output
    !%Description
    !% Upper bound for the energy mesh of the DOS.
    !% The default is the highest eigenvalue, plus a quarter of the total range of eigenvalues.
    !%End
    call parse_variable('DOSEnergyMax', evalmax + eextend, this%emax, units_inp%energy)

    !%Variable DOSEnergyPoints
    !%Type integer
    !%Default 500
    !%Section Output
    !%Description
    !% Determines how many energy points <tt>Octopus</tt> should use for 
    !% the DOS energy grid.
    !%End
    call parse_variable('DOSEnergyPoints', 500, this%epoints)

    !%Variable DOSGamma
    !%Type float
    !%Default 0.008 Ha
    !%Section Output
    !%Description
    !% Determines the width of the Lorentzian which is used for the DOS sum.
    !%End
    call parse_variable('DOSGamma', units_from_atomic(units_inp%energy, CNST(0.008)), this%gamma)
    this%gamma = units_to_atomic(units_inp%energy, this%gamma)

    ! spacing for energy mesh
    this%de = (this%emax - this%emin) / (this%epoints - 1)

    POP_SUB(dos_init)
  end subroutine dos_init

  ! ---------------------------------------------------------
  subroutine dos_end(this)
    type(dos_t), intent(inout)   :: this

    PUSH_SUB(dos_end)


    POP_SUB(dos_end)
  end subroutine dos_end


  ! ---------------------------------------------------------
  subroutine dos_write_dos(this, dir, st, sb, geo, mesh)
    type(dos_t),               intent(in) :: this
    character(len=*),         intent(in) :: dir
    type(states_t),           intent(in) :: st
    type(simul_box_t),        intent(in) :: sb
    type(geometry_t), target, intent(in) :: geo
    type(mesh_t),             intent(in) :: mesh

    integer :: ie, ik, ist, is, ns, maxdos
    integer, allocatable :: iunit(:)
    FLOAT   :: energy
    FLOAT   :: tdos
    FLOAT, allocatable :: dos(:,:,:)
    character(len=64)  :: filename

    integer :: ii, ll, mm, nn, work, norb, work2
    integer :: ia, iorb, idim, ip
    FLOAT   :: norm
    FLOAT, allocatable :: dpsi(:,:), ddot(:,:)
    CMPLX, allocatable :: zpsi(:,:), zdot(:,:)
    FLOAT, allocatable :: weight(:,:)

    PUSH_SUB(dos_write_dos)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    if(mpi_grp_is_root(mpi_world)) then 
      ! space for state-dependent DOS
      SAFE_ALLOCATE(dos(1:this%epoints, 1:st%nst, 0:ns-1))
      SAFE_ALLOCATE(iunit(0:ns-1))    

      ! compute band/spin-resolved density of states
      do ist = 1, st%nst

        do is = 0, ns-1
          if (ns > 1) then
            write(filename, '(a,i4.4,a,i1.1,a)') 'dos-', ist, '-', is+1,'.dat'
          else
            write(filename, '(a,i4.4,a)') 'dos-', ist, '.dat'
          end if
          iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
          ! write header
          write(iunit(is), '(3a)') '# energy [', trim(units_abbrev(units_out%energy)), '], band-resolved DOS'
        end do

        do ie = 1, this%epoints
          energy = this%emin + (ie - 1) * this%de
          dos(ie, ist, :) = M_ZERO
          ! sum up Lorentzians
          do ik = 1, st%d%nik, ns
            do is = 0, ns-1
             dos(ie, ist, is) = dos(ie, ist, is) + st%d%kweights(ik+is) * M_ONE/M_Pi * &
                this%gamma / ( (energy - st%eigenval(ist, ik+is))**2 + this%gamma**2 )
            end do
          end do
          do is = 0, ns-1
            write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                          units_from_atomic(unit_one / units_out%energy, dos(ie, ist, is))
            call messages_info(1, iunit(is))
          end do
        end do

        do is = 0, ns-1
          call io_close(iunit(is))
        end do
      end do

      ! for spin-polarized calculations also output spin-resolved tDOS
      if(st%d%nspin > 1) then    
        do is = 0, ns-1
          write(filename, '(a,i1.1,a)') 'total-dos-', is+1,'.dat'
          iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
          ! write header
          write(iunit(is), '(3a)') '# energy [', trim(units_abbrev(units_out%energy)), '], total DOS (spin-resolved)'

          do ie = 1, this%epoints
            energy = this%emin + (ie - 1) * this%de
            tdos = M_ZERO
            do ist = 1, st%nst
              tdos = tdos + dos(ie, ist, is)
            end do
            write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                          units_from_atomic(unit_one / units_out%energy, tdos)
            call messages_info(1, iunit(is))
          end do

          call io_close(iunit(is))
        end do
      end if


      iunit(0) = io_open(trim(dir)//'/'//'total-dos.dat', action='write')    
      write(iunit(0), '(3a)') '# energy [', trim(units_abbrev(units_out%energy)), '], total DOS'
      
      ! compute total density of states
      do ie = 1, this%epoints
        energy = this%emin + (ie - 1) * this%de
        tdos = M_ZERO
        do ist = 1, st%nst
          do is = 0, ns-1
            tdos = tdos + dos(ie, ist, is)
          end do
        end do
        write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                      units_from_atomic(unit_one / units_out%energy, tdos)
        call messages_info(1, iunit(0))
      end do

      call io_close(iunit(0))


      ! write Fermi file
      iunit(0) = io_open(trim(dir)//'/'//'total-dos-efermi.dat', action='write')
      write(message(1), '(3a)') '# Fermi energy [', trim(units_abbrev(units_out%energy)), &
        '] in a format compatible with total-dos.dat'

      ! this is the maximum that tdos can reach
      maxdos = st%smear%el_per_state * st%nst

      write(message(2), '(2f12.6)') units_from_atomic(units_out%energy, st%smear%e_fermi), M_ZERO
      write(message(3), '(f12.6,i6)') units_from_atomic(units_out%energy, st%smear%e_fermi), maxdos

      call messages_info(3, iunit(0))
      call io_close(iunit(0))

    end if

    SAFE_DEALLOCATE_A(iunit)
    SAFE_DEALLOCATE_A(dos)

    POP_SUB(dos_write_dos)
  end subroutine dos_write_dos


end module dos_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
