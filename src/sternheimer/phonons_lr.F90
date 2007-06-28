!! Copyright (C) 2004 Xavier Andrade, Eugene S. Kadantsev (ekadants@mjs1.phy.queensu.ca)
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
!! $Id: em_resp.F90 2686 2007-02-03 22:10:51Z xavier $

#include "global.h"
#define RESTART_DIR "restart_phn_lr/"

module phonons_lr_m
  use datasets_m
  use external_pot_m
  use functions_m
  use geometry_m
  use global_m
  use grid_m
  use geometry_m
  use hamiltonian_m
  use io_m
  use lib_basic_alg_m
  use lib_oct_parser_m
  use linear_response_m
  use math_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use output_m
  use phonons_m
  use projector_m
  use pert_m
  use restart_m
  use specie_m
  use states_m
  use sternheimer_m
  use string_m
  use system_m
  use units_m

  implicit none

  private
  public :: &
       phonons_lr_run

contains

  ! ---------------------------------------------------------
  subroutine phonons_lr_run(sys, h, fromScratch)
    type(system_t),         intent(inout) :: sys
    type(hamiltonian_t),    intent(inout) :: h
    logical,                intent(inout) :: fromScratch

    type(sternheimer_t) :: sh
    type(lr_t)          :: lr(1:1)
    type(phonons_t)     :: ph
    type(pert_t)   :: perturbation

    integer :: natoms, ndim, iatom, idir, jatom, jdir

    natoms = sys%geo%natoms
    ndim = sys%NDIM

    !CONSTRUCT

    call push_sub('phonons_lr.phonons_lr_run')

    call parse_input()
    call read_wfs(sys%st, sys%gr, sys%geo, .false.)

    message(1) = 'Info: Setting up Hamiltonian for linear response'
    call write_info(1)

    call system_h_setup(sys, h)
    call sternheimer_init(sh, sys, h, "Phn")
    call phonons_init(ph, sys)

    call lr_init(lr(1))
    call lr_allocate(lr(1), sys%st, sys%gr%m)

    !CALCULATE

    !the ionic contribution
    call build_ionic_dm()

    call pert_init(perturbation, PERTURBATION_IONIC, sys%gr, sys%geo)

    do iatom = 1, natoms
      do idir = 1, ndim

        call pert_setup_atom(perturbation, iatom)
        call pert_setup_dir(perturbation, idir)

        call dsternheimer_solve(sh, sys, h, lr, 1, M_ZERO, perturbation, &
             RESTART_DIR, phn_rho_tag(iatom, idir), phn_wfs_tag(iatom, idir))


        do jatom = 1, natoms
          do jdir = 1, ndim

            call pert_setup_atom(perturbation, jatom, iatom)
            call pert_setup_dir(perturbation, jdir, idir)

            ph%dm(phn_idx(ph, iatom, idir), phn_idx(ph, jatom, jdir)) &
                 = ph%dm(phn_idx(ph, iatom, idir), phn_idx(ph, jatom, jdir))&
                 -dpert_expectation_value(perturbation,sys%gr,sys%geo,h,sys%st, lr(1)%ddl_psi, sys%st%dpsi)&
                 -dpert_expectation_value(perturbation,sys%gr,sys%geo,h,sys%st, sys%st%dpsi, lr(1)%ddl_psi)&
                 -dpert_expectation_value(perturbation,sys%gr,sys%geo,h,sys%st, sys%st%dpsi, sys%st%dpsi, pert_order = 2)
                 
          end do
        end do

      end do
    end do
    
    call pert_end(perturbation)

    call phonons_normalize_dm(ph, sys%geo)
    call phonons_diagonalize_dm(ph)
    call phonons_output(ph, "_lr")

    !DESTRUCT

    call lr_dealloc(lr(1))
    
    call phonons_end(ph)

    call sternheimer_end(sh)

    call states_deallocate_wfns(sys%st)

    call pop_sub()

  contains

    ! ---------------------------------------------------------
    subroutine parse_input()
      
      call push_sub('phonons_lr.parse_input')

      call pop_sub()

    end subroutine parse_input
    
    subroutine build_ionic_dm()
      FLOAT :: ac, xi(1:MAX_DIM), xj(1:MAX_DIM), xk(1:MAX_DIM), r2
      integer :: katom

      do iatom = 1, natoms
        do idir = 1, sys%NDIM
          
          do jatom = 1, natoms
            do jdir = 1, sys%NDIM         
              
              xi(1:MAX_DIM) = sys%geo%atom(iatom)%x(1:MAX_DIM)

              !ion - ion
              if( iatom == jatom) then 
                
                ac = M_ZERO
                do katom = 1, natoms
                  if ( katom == iatom ) cycle
                  
                  xk(1:MAX_DIM) = sys%geo%atom(katom)%x(1:MAX_DIM)
                  r2 = sum((xi(1:sys%NDIM) - xk(1:sys%NDIM))**2)
                  
                  ac = ac + sys%geo%atom(iatom)%spec%Z_val * sys%geo%atom(katom)%spec%Z_val &
                       /(r2**CNST(1.5)) *(&
                       -ddelta(idir, jdir) + &
                       (M_THREE*(xi(idir)-xk(idir))*(xi(jdir)-xk(jdir)))/r2 &
                       )
                  
                end do

              else ! iatom /= jatom
                
                xj(1:MAX_DIM) = sys%geo%atom(jatom)%x(1:MAX_DIM)
                
                r2 = sum((xi(1:sys%NDIM) - xj(1:sys%NDIM))**2)
                ac = sys%geo%atom(iatom)%spec%Z_val * sys%geo%atom(jatom)%spec%Z_val &
                     /(r2**CNST(1.5))*(&
                     ddelta(idir, jdir) - (M_THREE*(xi(idir)-xj(idir))*(xi(jdir)-xj(jdir)))/r2)

              end if
                
            
              ph%dm(phn_idx(ph, iatom, idir), phn_idx(ph, jatom, jdir)) = -ac

            end do
          end do
        end do
      end do
      
    end subroutine build_ionic_dm
    

  end subroutine phonons_lr_run


  ! ---------------------------------------------------------
  subroutine read_wfs(st, gr, geo, complex_wfs)
    type(states_t),   intent(inout) :: st
    type(grid_t),     intent(in)    :: gr
    type(geometry_t), intent(in)    :: geo
    logical,          intent(in)    :: complex_wfs

    integer :: kpoints, nst, ierr, dim
      
    !check how many wfs we have

    call push_sub('em_resp.read_wfs')

    call states_look(trim(tmpdir)//'restart_gs', gr%m, kpoints, dim, nst, ierr)
    if(ierr.ne.0) then
      message(1) = 'Could not properly read wave-functions from "'//trim(tmpdir)//'restart_gs".'
      call write_fatal(1)
    end if

    st%nst    = nst
    st%st_end = nst
    deallocate(st%eigenval, st%occ)

    if ( complex_wfs ) then 
      call states_allocate_wfns(st, gr%m, M_CMPLX)
    else 
      call states_allocate_wfns(st, gr%m, M_REAL)
    end if

    ALLOCATE(st%eigenval(st%nst, st%d%nik), st%nst*st%d%nik)
    ALLOCATE(st%occ(st%nst, st%d%nik), st%nst*st%d%nik)

    if(st%d%ispin == SPINORS) then
      ALLOCATE(st%spin(3, st%nst, st%d%nik), st%nst*st%d%nik*3)
      st%spin = M_ZERO
    end if
    st%eigenval = huge(REAL_PRECISION)
    st%occ      = M_ZERO

    ! load wave-functions
    call restart_read(trim(tmpdir)//'restart_gs', st, gr, geo, ierr)  
    if(ierr.ne.0) then
      message(1) = "Could not read KS orbitals from '"//trim(tmpdir)//"restart_gs'"
      message(2) = "Please run a calculation of the ground state first!"
      call write_fatal(2)
    end if

    call pop_sub()

  end subroutine read_wfs

  character(len=100) function phn_rho_tag(iatom, dir) result(str)
    integer, intent(in) :: iatom, dir
    
    call push_sub('phonons_lr.phn_rho_tag')
    
    write(str, '(a,i1,a,i1)') 'phn_rho_', iatom, '_',  dir

    call pop_sub()

  end function phn_rho_tag
  
  character(len=100) function phn_wfs_tag(iatom, dir) result(str)
    integer, intent(in) :: iatom, dir

    call push_sub('phonons_lr.phn_wfs_tag')

    write(str, '(a,i1,a,i1)') "phn_wfs_", iatom, "_", dir

    call pop_sub()
    
  end function phn_wfs_tag

end module phonons_lr_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
