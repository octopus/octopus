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
!! $Id: states.F90 2515 2006-10-24 17:13:30Z acastro $

#include "global.h"

module states_output_m
  use global_m
  use messages_m
  use states_m
  use grid_m
  use output_m
  use units_m
  use lib_oct_m
  use io_m
  use magnetic_m
  use elf_m

  private
  public ::                         &
    states_output

contains

  ! ---------------------------------------------------------
  subroutine states_output(st, gr, dir, outp)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(inout) :: gr
    character(len=*), intent(in)    :: dir
    type(output_t),   intent(in)    :: outp

    integer :: ik, ist, idim, is, id, ierr, iunit
    character(len=80) :: fname
    FLOAT :: u
    FLOAT, allocatable :: dtmp(:), elf(:,:)

    call push_sub('states.states_output')

    u = M_ONE/units_out%length%factor**NDIM

    if(iand(outp%what, output_density).ne.0) then
      do is = 1, st%d%nspin
        write(fname, '(a,i1)') 'density-', is
        call doutput_function(outp%how, dir, fname, gr%m, gr%sb, &
          st%rho(:, is), u, ierr, is_tmp = .false.)
      end do
    end if

    if(iand(outp%what, output_pol_density).ne.0) then
      ALLOCATE(dtmp(NP), NP)
      do idim=1, NDIM
        do is = 1, st%d%nspin
          dtmp(1:NP)=st%rho(1:NP,is)*gr%m%x(1:NP,idim)
          write(fname, '(a,i1,a,i1)') 'dipole_density-', is, '-',idim
          call doutput_function(outp%how, dir, fname, gr%m, gr%sb, &
            dtmp(:), u, ierr, is_tmp = .false.)
        end do
      end do
    end if

    if( (iand(outp%what, output_current).ne.0) .and. (st%d%wfs_type == M_CMPLX) ) then
      ! calculate current first
      call calc_paramagnetic_current(gr, st, st%j)
      do is = 1, st%d%nspin
        do id = 1, NDIM
          write(fname, '(a,i1,a,a)') 'current-', is, '-', index2axis(id)
          call doutput_function(outp%how, dir, fname, gr%m, gr%sb, &
            st%j(:, id, is), u, ierr, is_tmp = .false.)
        end do
      end do
    end if

    if(iand(outp%what, output_wfs).ne.0) then
      do ist = st%st_start, st%st_end
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = 1, st%d%nik
            do idim = 1, st%d%dim
              write(fname, '(a,i3.3,a,i4.4,a,i1)') 'wf-', ik, '-', ist, '-', idim
              if (st%d%wfs_type == M_REAL) then
                call doutput_function(outp%how, dir, fname, gr%m, gr%sb, &
                     st%dpsi(1:, idim, ist, ik), sqrt(u), ierr, is_tmp = .false.)
              else
                call zoutput_function(outp%how, dir, fname, gr%m, gr%sb, &
                     st%zpsi(1:, idim, ist, ik), sqrt(u), ierr, is_tmp = .false.)
              end if
            end do
          end do
        end if
      end do
    end if

    if(iand(outp%what, output_wfs_sqmod).ne.0) then
      ALLOCATE(dtmp(NP_PART), NP_PART)
      do ist = 1, st%nst
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = 1, st%d%nik
            do idim = 1, st%d%dim
              write(fname, '(a,i3.3,a,i4.4,a,i1)') 'sqm-wf-', ik, '-', ist, '-', idim
              if (st%d%wfs_type == M_REAL) then
                dtmp = abs(st%dpsi(:, idim, ist, ik))**2
              else
                dtmp = abs(st%zpsi(:, idim, ist, ik))**2
              end if
              call doutput_function (outp%how, dir, fname, gr%m, gr%sb, &
                dtmp, u, ierr, is_tmp = .false.)
            end do
          end do
        end if
      end do
      deallocate(dtmp)
    end if

    if(NDIM .eq. 3) then ! If the dimensions is not three, the ELF calculation will not work.
      if(  iand(outp%what, output_elf).ne.0  ) then ! First, ELF in real space.
        select case(st%d%ispin)
          case(UNPOLARIZED)
            ALLOCATE(elf(1:gr%m%np, 1),gr%m%np)
            call elf_calc(st, gr, elf)
            write(fname, '(a)') 'elf_rs'
            call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
              elf(:,1), M_ONE, ierr, is_tmp = .false.)
            deallocate(elf)

          case(SPIN_POLARIZED, SPINORS)
            ALLOCATE(elf(1:gr%m%np, 3), 3*gr%m%np)
            call elf_calc(st, gr, elf)
            write(fname, '(a)') 'elf_rs'
            call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
              elf(:, 3), M_ONE, ierr, is_tmp = .false.)
            do is = 1, 2
              write(fname, '(a,a,i1)') 'elf_rs', '-', is
              call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
                elf(:, is), M_ONE, ierr, is_tmp = .false.)
            end do
            deallocate(elf)
        end select
      end if

      if(  iand(outp%what, output_elf_fs).ne.0  ) then ! Second, ELF in Fourier space.
        ALLOCATE(elf(1:gr%m%np,1:st%d%nspin),gr%m%np*st%d%nspin)
        call elf_calc_fs(st, gr, elf)
        do is = 1, st%d%nspin
          write(fname, '(a,a,i1)') 'elf_fs', '-', is
          call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
            elf(:,is), M_ONE, ierr, is_tmp = .false.)
        end do
        deallocate(elf)
      end if
    end if

    if(iand(outp%what, output_ksdipole).ne.0) then
      ! The files will be called dir/dipole.k.d, where d will be the
      ! direction ("x", "y" or "z"), and k is the k-point -- this may
      ! also mean the spin subspace in spin-polarized calculations.
      do idim = 1, gr%sb%dim
        do ik = 1, st%d%nik
          select case(idim)
            case(1); write(fname,'(i3.3,a2)') ik, '.x'
            case(2); write(fname,'(i3.3,a2)') ik, '.y'
            case(3); write(fname,'(i3.3,a2)') ik, '.z'
          end select
          write(fname,'(a)') trim(dir)//'/dipole.'//trim(adjustl(fname))
          iunit = io_open(file = fname, action = 'write')
          call states_write_dipole_matrix(st, gr, idim, ik, iunit)
          call io_close(iunit)
        end do
      end do
    end if

    call pop_sub()
  end subroutine states_output


  ! ---------------------------------------------------------
  ! Prints out the dipole matrix elements between KS states.
  ! It prints the states to the file opened in iunit.
  ! It prints the dipole moment in the direction "k", for
  ! the Kohn-Sham states in the irreducible subspace ik.
  ! ---------------------------------------------------------
  subroutine states_write_dipole_matrix(st, gr, k, ik, iunit)
    type(states_t), intent(in) :: st
    type(grid_t), intent(in) :: gr
    integer, intent(in) :: k, ik, iunit

    integer :: ii, jj
    FLOAT, allocatable :: dipole(:, :)
    CMPLX :: dip_element

    call push_sub('states.states_write_dipole_matrix')

    ALLOCATE(dipole(gr%m%np_part, st%d%dim), gr%m%np_part)
    dipole = M_ZERO
    do ii = 1, st%d%dim
      dipole(1:gr%m%np, ii) = gr%m%x(1:gr%m%np, k)
    end do

    do ii = 1, st%nst
      do jj = 1, st%nst
        if (st%d%wfs_type == M_REAL) then 
          write(iunit,fmt = '(f20.10)', advance = 'no') dstates_dotp(gr%m, st%d%dim, &
            st%dpsi(:, :, ii, 1), &
            st%dpsi(:, :, jj, 1) * dipole(:, :))
        else
          dip_element = zstates_dotp(gr%m, st%d%dim, &
            st%zpsi(:, :, ii, 1), &
            st%zpsi(:, :, jj, 1) * dipole(:, :))
          write(iunit,fmt = '(f20.12,a1,f20.12,3x)', advance = 'no') &
            real(dip_element, PRECISION), ',', aimag(dip_element)
        end if
        if(jj==st%nst) write(iunit, '(a)') 
      end do
    end do

    deallocate(dipole)
    call pop_sub()
  end subroutine states_write_dipole_matrix

end module states_output_m
