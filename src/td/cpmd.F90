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
!! $Id: td.F90 3694 2008-02-15 13:37:54Z marques $

#include "global.h"
  
module cpmd_m
  use global_m
  use io_m
  use datasets_m
  use io_function_m
  use lalg_basic_m
  use loct_math_m
  use loct_parser_m
  use units_m
  use math_m
  use messages_m
  use mesh_m
  use external_pot_m
  use geometry_m
  use hamiltonian_m
  use loct_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use v_ks_m
  use grid_m

  implicit none

  private
  public ::                 &
    cpmd_t,                 &
    cpmd_init,              &
    cpmd_end,               &
    cpmd_electronic_energy, &
    cpmd_restart_read,      &
    cpmd_restart_write,     &
    cpmd_propagate

  type cpmd_t
    private
    FLOAT          :: emass
    FLOAT          :: ecorr
    CMPLX, pointer :: oldpsi(:, :, :, :)
  end type cpmd_t

  type(profile_t), save :: cpmd_prop, cpmd_orth

contains

  subroutine cpmd_init(this, gr, st)
    type(cpmd_t),        intent(out)   :: this
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    
    integer :: size

    !%Variable CPElectronicMass
    !%Type float
    !%Default 1.0
    !%Section Time Dependent::Propagation
    !%Description
    !% The fictious electronic mass used to propagate the electronic
    !% wavefunctions in the Carr-Parrinelo formalism.
    !%End
    
    call loct_parse_float(check_inp('CPElectronicMass'), CNST(1.0), this%emass)

    size = gr%m%np * st%d%dim * st%lnst * st%d%nik
    ALLOCATE(this%oldpsi(gr%m%np, st%d%dim, st%st_start:st%st_end, st%d%nik), size)

  end subroutine cpmd_init
  
  subroutine cpmd_end(this)
    type(cpmd_t), intent(inout) :: this

    deallocate(this%oldpsi)

  end subroutine cpmd_end
  
  subroutine cpmd_propagate(this, gr, h, st, iter, dt)
    type(cpmd_t),        intent(inout) :: this
    type(grid_t),        intent(inout) :: gr
    type(hamiltonian_t), intent(inout) :: h
    type(states_t),      intent(inout) :: st
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt

    ! this integration is based on Tuckermann and Parrinello, JCP 101
    ! 1302 (1994), using the verlet algorithm described in page 1306,
    ! eqs 4.2 to 4.7.

    integer :: ik, ist1, ist2, ddim, np

    CMPLX, allocatable :: hpsi(:, :), psi(:, :), xx(:, :)

    call profiling_in(cpmd_prop, "CP_PROPAGATION")

    np = gr%m%np
    ddim = st%d%dim

    ALLOCATE(xx(1:st%nst, 1:st%nst), st%nst**2)
    ALLOCATE(hpsi(1:gr%m%np, 1:st%d%dim), gr%m%np*st%d%dim)
    ALLOCATE(psi(1:gr%m%np, 1:st%d%dim), gr%m%np*st%d%dim)

    this%ecorr = M_ZERO

    do ik = 1, st%d%nik
      do ist1 = st%st_start, st%st_end

        ! give the initial conditions
        if(iter == 1) this%oldpsi(1:np, 1:ddim, ist1, ik) = st%zpsi(1:np, 1:ddim, ist1, ik)

        ! calculate the "force"
        call zHpsi(h, gr, st%zpsi(:, :, ist1, ik), hpsi, ist1, ik)
        
        ! propagate the wavefunctions
        psi(1:np, 1:ddim) = M_TWO*st%zpsi(1:np, 1:ddim, ist1, ik) - this%oldpsi(1:np, 1:ddim, ist1, ik) &
             + dt**2/this%emass*(-st%occ(ist1, ik)*hpsi(1:np, 1:ddim)) !(4.2)
        
        ! calculate the velocity and the fictitious electronic energy
        hpsi(1:np, 1:ddim) = abs(psi(1:np, 1:ddim) - this%oldpsi(1:np, 1:ddim, ist1, ik))/(M_TWO*dt)
        this%ecorr = this%ecorr + this%emass*zstates_nrm2(gr%m, ddim, hpsi)**2 !(2.11)

        ! store the old wavefunctions
        this%oldpsi(1:np, 1:ddim, ist1, ik) = st%zpsi(1:np, 1:ddim, ist1, ik)
        st%zpsi(1:np, 1:ddim, ist1, ik) = psi(1:np, 1:ddim)
        
      end do

      call profiling_in(cpmd_orth, "CP_ORTHOGONALIZATION")

      call calc_lambda
      
      do ist1 = st%st_start, st%st_end
        do ist2 = 1, st%nst
          st%zpsi(1:np, 1:ddim, ist1, ik) =  st%zpsi(1:np, 1:ddim, ist1, ik) &
               + xx(ist1, ist2)*this%oldpsi(1:np, 1:ddim, ist2, ik) !(4.3)
        end do
      end do
      
      call profiling_out(cpmd_orth)

    end do

    deallocate(hpsi, psi, xx)

    call profiling_out(cpmd_prop)

  contains

    subroutine calc_lambda
      !this routine should be modified for states parallelization
      integer :: ist1, ist2, it
      FLOAT   :: res
      FLOAT, allocatable :: ii(:, :)
      CMPLX, allocatable :: aa(:, :), bb(:, :), xxi(:, :)
      
      ALLOCATE(aa(1:st%nst, 1:st%nst), st%nst**2)
      ALLOCATE(bb(1:st%nst, 1:st%nst), st%nst**2)
      ALLOCATE(ii(1:st%nst, 1:st%nst), st%nst**2)
      ALLOCATE(xxi(1:st%nst, 1:st%nst), st%nst**2)
            
      do ist1 = 1, st%nst
        do ist2 = 1, st%nst
          ii(ist1, ist2) = ddelta(ist1, ist2)
          aa(ist1, ist2) = zstates_dotp(gr%m, ddim,     st%zpsi(:, :, ist1, ik), st%zpsi(:, :, ist2, ik))
          bb(ist1, ist2) = zstates_dotp(gr%m, ddim, this%oldpsi(:, :, ist1, ik), st%zpsi(:, :, ist2, ik))
        end do
      end do
      
      xx = M_HALF*(ii - aa) !(4.6)
      
      do it = 1, 10
        xxi = M_HALF*(ii - aa + matmul(xx, ii - bb) + matmul(ii - transpose(bb), xx) - matmul(xx, xx)) !(4.5)
        res = maxval(abs(xxi - xx))
        xx = xxi
        if (res < CNST(1e-5)) exit
      end do

      deallocate(aa, bb, ii, xxi)

    end subroutine calc_lambda

  end subroutine cpmd_propagate

  FLOAT function cpmd_electronic_energy(this)
    type(cpmd_t),        intent(in)    :: this

    cpmd_electronic_energy = this%ecorr
    
  end function cpmd_electronic_energy

  subroutine cpmd_restart_write(this, gr, st)
    type(cpmd_t),        intent(in)    :: this
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st

    integer :: ik, ist, idim, ii, err
    character(len=80) :: filename

    call io_mkdir(trim(tmpdir)//'td/cpmd')

    ii = 1
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          write(filename,'(i10.10)') ii
          call zrestart_write_function(trim(tmpdir)//'td/cpmd', filename, gr, this%oldpsi(:, idim, ist, ik), err, gr%m%np)
          ii = ii + 1
        end do
      end do
    end do

  end subroutine cpmd_restart_write

  subroutine cpmd_restart_read(this, gr, st, ierr)
    type(cpmd_t),        intent(inout) :: this
    type(grid_t),        intent(in)    :: gr
    type(states_t),      intent(in)    :: st
    integer,             intent(out)   :: ierr

    integer :: ik, ist, idim, ii
    character(len=80) :: filename

    ierr = 0

    ii = 1
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          write(filename,'(i10.10)') ii
          call zrestart_read_function(trim(tmpdir)//'td/cpmd', filename, gr%m, this%oldpsi(:, idim, ist, ik), ierr)
          if(ierr /= 0) return
          ii = ii + 1
        end do
      end do
    end do

  end subroutine cpmd_restart_read

end module cpmd_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
