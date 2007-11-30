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
!! $Id: states.F90 2515 2006-10-24 17:13:30Z acastro $

#include "global.h"

module states_output_m
  use basins_m
  use elf_m
  use global_m
  use grid_m
  use io_m
  use loct_m
  use loct_math_m
  use magnetic_m
  use messages_m
  use output_m
  use states_m
  use mesh_m
  use math_m
  use units_m

  implicit none

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

    integer :: ik, ist, idim, is, id, ierr, iunit, l, m
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
      deallocate(dtmp)
    end if

    if( (iand(outp%what, output_current).ne.0) .and. (st%wfs_type == M_CMPLX) ) then
      ! calculate current first
      call states_paramagnetic_current(gr, st, st%j)
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
              if (st%wfs_type == M_REAL) then
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
      do ist = st%st_start, st%st_end
        if(loct_isinstringlist(ist, outp%wfs_list)) then
          do ik = 1, st%d%nik
            do idim = 1, st%d%dim
              write(fname, '(a,i3.3,a,i4.4,a,i1)') 'sqm-wf-', ik, '-', ist, '-', idim
              if (st%wfs_type == M_REAL) then
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

    call out_elf()

    if(iand(outp%what, output_ked).ne.0) then
      ALLOCATE(elf(gr%m%np, st%d%nspin),gr%m%np*st%d%nspin)
      call kinetic_energy_density(st, gr, elf)
      select case(st%d%ispin)
        case(UNPOLARIZED)
          write(fname, '(a)') 'tau'
          call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
            elf(:,1), M_ONE, ierr, is_tmp = .false.)
        case(SPIN_POLARIZED, SPINORS)
          do is = 1, 2
            write(fname, '(a,a,i1)') 'tau', '-', is
            call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
              elf(:, is), M_ONE, ierr, is_tmp = .false.)
          end do
      end select
      deallocate(elf)
    end if

    if(iand(outp%what, output_ksdipole).ne.0) then
      ! The files will be called matrix_elements.x. The content of each file
      ! should be clear from the header of each file.
      id = 1
      do ik = 1, st%d%nik
        do l = 1, outp%ksmultipoles
          do m = -l, l
            write(fname,'(i4)') id
            write(fname,'(a)') trim(dir)//'/matrix_elements.'//trim(adjustl(fname))
            iunit = io_open(file = fname, action = 'write')
            call states_write_multipole_matrix(st, gr, l, m, ik, iunit)
            call io_close(iunit)
            id = id + 1
          end do
        end do
      end do
    end if

    call pop_sub()

  contains
    ! ---------------------------------------------------------
    subroutine out_elf()
      FLOAT, allocatable :: elf(:,:)
      type(basins_t) :: basins
      integer :: imax, iunit

      ! If it is a one-dimensional problem, the ELF calculation will not work.
      if(NDIM .eq. 1) return

      if(iand(outp%what, output_elf).ne.0 .or. iand(outp%what, output_elf_basins).ne.0) then
        imax = 1
        if(st%d%ispin.ne.UNPOLARIZED) imax = 3

        ALLOCATE(elf(1:NP, imax), NP*imax)
        call elf_calc(st, gr, elf)
      end if

      ! output ELF in real space
      if(iand(outp%what, output_elf).ne.0) then
        write(fname, '(a)') 'elf_rs'
        call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
           elf(:,imax), M_ONE, ierr, is_tmp = .false.)

        if(st%d%ispin.ne.UNPOLARIZED) then
          do is = 1, 2
            write(fname, '(a,a,i1)') 'elf_rs', '-', is
            call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
               elf(:, is), M_ONE, ierr, is_tmp = .false.)
          end do
        end if
      end if

      if(iand(outp%what, output_elf_basins).ne.0) then
        call basins_init(basins, gr%m)
        call basins_analyze(basins, gr%m, st%d%nspin, elf(:,1), st%rho, CNST(0.01))

        write(fname, '(a)') 'elf_rs_basins'
        call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, &
           real(basins%map, REAL_PRECISION), M_ONE, ierr, is_tmp = .false.)
        
        write(fname,'(2a)') trim(dir), '/elf_rs_basins.info'
        iunit = io_open(file = fname, action = 'write')
        call basins_write(basins, gr%m, iunit)
        call io_close(iunit)

        call basins_end(basins)
      end if

      ! clean up
      if(iand(outp%what, output_elf).ne.0 .or. iand(outp%what, output_elf_basins).ne.0) then
        deallocate(elf)
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
    end subroutine out_elf

  end subroutine states_output


  ! ---------------------------------------------------------
  ! Prints out the multipole matrix elements between KS states.
  ! It prints the states to the file opened in iunit.
  ! It prints the (l,m) multipole moment, for
  ! the Kohn-Sham states in the irreducible subspace ik.
  ! ---------------------------------------------------------
  subroutine states_write_multipole_matrix(st, gr, l, m, ik, iunit)
    type(states_t), intent(in) :: st
    type(grid_t), intent(in) :: gr
    integer, intent(in) :: l, m, ik, iunit

    integer :: ii, jj, i
    FLOAT, allocatable :: multipole(:, :)
    CMPLX :: multip_element
    FLOAT :: r, x(MAX_DIM), ylm

    call push_sub('states.states_write_multipole_matrix')

    write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | r**l * Y_{lm}(theta,phi) | Phi_j>' 
    write(iunit, fmt = '(a,i2,a,i2)') '# l =', l, '; m =', m
    write(iunit, fmt = '(a,i4)')      '# ik =', ik
    if(l>1) then
      write(iunit, fmt = '(a,i1)') '# Units = ['//trim(units_out%length%abbrev)//']^',l
    else
      write(iunit, fmt = '(a)')    '# Units = ['//trim(units_out%length%abbrev)//']'
    end if

    ALLOCATE(multipole(gr%m%np_part, st%d%dim), gr%m%np_part)
    multipole = M_ZERO
    do ii = 1, st%d%dim
      do i = 1, NP
        call mesh_r(gr%m, i, r, x = x)
        ylm = loct_ylm(x(1), x(2), x(3), l, m)
        multipole(i, ii) = r**l * ylm
      end do
    end do

    do ii = 1, st%nst
      do jj = 1, st%nst
        if (st%wfs_type == M_REAL) then 
          write(iunit,fmt = '(f20.10)', advance = 'no') dstates_dotp(gr%m, st%d%dim, &
            st%dpsi(:, :, ii, 1), &
            st%dpsi(:, :, jj, 1) * multipole(:, :)) / units_out%length%factor**l

        else
          multip_element = zstates_dotp(gr%m, st%d%dim, &
            st%zpsi(:, :, ii, 1), &
            st%zpsi(:, :, jj, 1) * multipole(:, :)) / units_out%length%factor**l

          write(iunit,fmt = '(f20.12,a1,f20.12,3x)', advance = 'no') &
            real(multip_element, REAL_PRECISION), ',', aimag(multip_element)
        end if
        if(jj==st%nst) write(iunit, '(a)') 
      end do
    end do

    deallocate(multipole)
    call pop_sub()
  end subroutine states_write_multipole_matrix

end module states_output_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
