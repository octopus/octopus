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
!! $Id: run.F90 5113 2009-03-26 02:06:16Z dstrubbe $

#include "global.h"

module invert_ks_m
  use datasets_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use loct_parser_m
  use mesh_function_m
  use messages_m
  use mix_m
  use profiling_m
  use restart_m
  use scf_m
  use states_m
  use system_m
  use v_ks_m
  
  implicit none

  private
  public :: invert_ks_run

contains

  ! ---------------------------------------------------------
  subroutine invert_ks_run(sys, hm)
    type(system_t),              intent(inout) :: sys
    type(hamiltonian_t),         intent(inout) :: hm

    integer :: ii, ierr, np, ndim, nspin
    FLOAT, allocatable :: target_rho(:,:), rhoin(:,:,:), rhoout(:,:,:), rhonew(:,:,:)
    type(scf_t) :: scfv
    type(mix_t) :: smix
    character(len=256) :: fname

    ! load wave-functions
    call states_allocate_wfns(sys%st, sys%gr%mesh)
    call restart_read(trim(tmpdir)//'gs', sys%st, sys%gr, sys%geo, ierr)
    if(ierr.ne.0) then
      message(1) = "Could not load wave-functions"
      call write_fatal(1)
    end if

    np = sys%gr%mesh%np
    ndim    = sys%gr%mesh%sb%dim ! shortcut
    nspin   = sys%st%d%nspin

    ! generate density
    call states_calc_dens(sys%st, np)

    ! read target density
    SAFE_ALLOCATE(target_rho(1:np, 1:nspin))
    call read_target_rho()

    call scf_init(scfv, sys%gr, sys%geo, sys%st, hm)
    call mix_init(smix, np, nspin, 1, prefix_="InvertKS")

    !call v_ks_calc(sys%gr, sys%ks, hm, sys%st, .true.)
    !hm%vhxc(1:np, 1) =  hm%vhxc(1:np, 1) + hm%ep%vpsl(1:np)

    hm%epot     = M_ZERO
    hm%ehartree = M_ZERO
    hm%ex       = M_ZERO
    hm%ec       = M_ZERO
    hm%vhartree = M_ZERO
    hm%vxc      = M_ZERO
    hm%ep%vpsl  = M_ZERO
    hm%ep%local_potential = M_ZERO

    sys%ks%frozen_hxc = .true.
    hm%vhxc(1:np, 1) = M_ONE

    SAFE_ALLOCATE(rhoin (1:np, 1:nspin, 1:1))
    SAFE_ALLOCATE(rhoout(1:np, 1:nspin, 1:1))
    SAFE_ALLOCATE(rhonew(1:np, 1:nspin, 1:1))

    rhoin(1:np,1:nspin,1) = hm%vhxc(1:np,1:nspin)

    do ii = 1, 400
      !write(fname,'(i3.3)') ii
      !call doutput_function(io_function_fill_how("AxisX"), &
      !     ".", "vhxc"//fname, sys%gr%mesh, sys%gr%sb, hm%vhxc(:,1), M_ONE, ierr)
      !call doutput_function(io_function_fill_how("AxisX"), &
      !     ".", "rho"//fname, sys%gr%mesh, sys%gr%sb, sys%st%rho(:,1), M_ONE, ierr)

      call hamiltonian_update_potential(hm, sys%gr%mesh)

      call scf_run(scfv, sys%gr, sys%geo, sys%st, sys%ks, hm, sys%outp, gs_run = .false.)
      call states_calc_dens(sys%st, np)

      rhoout(1:np, 1:nspin, 1) = (sys%st%rho(1:np,1:nspin))/(target_rho(1:np,1:nspin)) &
           * hm%vhxc(1:np, 1:nspin)
      call dmixing(smix, ii, rhoin, rhoout, rhonew, dmf_dotp_aux)

      hm%vhxc(1:np,1:nspin) = rhonew(1:np, 1:nspin, 1)
      rhoin(1:np, 1:nspin, 1) = hm%vhxc(1:np, 1:nspin)
    end do

    call mix_end(smix)
    call scf_end(scfv)

    SAFE_DEALLOCATE_A(target_rho)

  contains

  ! ---------------------------------------------------------
    subroutine read_target_rho()
      character(len=256) :: filename
      integer :: pass, iunit, ierr, ii, npoints
      FLOAT   :: l_xx(MAX_DIM), l_ff(4), rr
      FLOAT, allocatable :: xx(:,:), ff(:,:)

      call loct_parse_string(datasets_check('InvertKSTargetDensity'), "target_density.dat", filename)

      iunit = io_open(filename, action='read', status='old')

      npoints = 0
      do pass = 1, 2
        ii = 0
        rewind(iunit)
        do
          read(iunit, fmt=*, iostat=ierr) l_xx(1:ndim), l_ff(1:nspin)
          if(ierr.ne.0) exit
          ii = ii + 1
          if(pass == 1) npoints = npoints + 1
          if(pass == 2) then
            xx(ii, 1:ndim)  = l_xx(1:ndim)
            ff(ii, 1:nspin) = l_ff(1:nspin)
          end if
        end do
        if(pass == 1) then
          SAFE_ALLOCATE(xx(1:npoints, 1:ndim))
          SAFE_ALLOCATE(ff(1:npoints, 1:nspin))
        end if
      end do

      do ii = 1, nspin
        call dmf_interpolate_points(ndim, npoints, xx, ff(:,ii), &
             np, sys%gr%mesh%x, target_rho(:, ii))
      end do
 
      ! we now renormalize the density (necessary if we have a charged system)
      rr = M_ZERO
      do ii = 1, sys%st%d%spin_channels
        rr = rr + dmf_integrate(sys%gr%mesh, target_rho(:, ii))
      end do
      rr = sys%st%qtot/rr
      target_rho(:,:) = rr*target_rho(:,:)

      call doutput_function(io_function_fill_how("AxisX"), &
           ".", "func", sys%gr%mesh, sys%gr%sb, target_rho(:,1), M_ONE, ierr)

      SAFE_DEALLOCATE_A(xx)
      SAFE_DEALLOCATE_A(ff)

    end subroutine read_target_rho

  end subroutine invert_ks_run

end module invert_ks_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
