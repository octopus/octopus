!! Copyright (C) 2016 N. Tancogne-Dejean 
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
#include "global.h"

module orbitalset_utils_oct_m
  use atomic_orbital_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use distributed_oct_m
  use global_oct_m
  use io_function_oct_m
  use ions_oct_m
  use lattice_vectors_oct_m
  use loct_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use orbitalset_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use species_oct_m
  use submesh_oct_m
  use unit_oct_m
 
  implicit none

  private

  public ::                            &
       orbitalset_utils_count,         &
       dorbitalset_utils_getorbitals,  &
       zorbitalset_utils_getorbitals,  &
       orbitalset_init_intersite

  integer, public, parameter ::     &
    SM_POISSON_DIRECT          = 0, &
    SM_POISSON_ISF             = 1, &
    SM_POISSON_PSOLVER         = 2, &
    SM_POISSON_FFT             = 3


contains

  integer function orbitalset_utils_count(ions, ia, iselect) result(norb)
    type(ions_t),         intent(in) :: ions
    integer,              intent(in) :: ia
    integer, optional,    intent(in) :: iselect

    integer :: iorb, ii, ll, mm

    !We count the number of orbital sets we have for a given atom
    !If iselect is present, this routine return instead the number of orbital for a given
    !value of i
    norb = 0
    do iorb = 1, species_niwfs(ions%atom(ia)%species)
      call species_iwf_ilm(ions%atom(ia)%species, iorb, 1, ii, ll, mm)
      if(present(iselect)) then
        if(ii == iselect) norb = norb + 1
      else
        norb = max(norb, ii)
      end if
    end do
  end function orbitalset_utils_count

  subroutine orbitalset_init_intersite(this, namespace, ind, ions, der, psolver, os, nos, maxnorbs, &
                rcut, kpt, has_phase, sm_poisson)
    type(orbitalset_t),           intent(inout) :: this
    type(namespace_t),            intent(in)    :: namespace
    integer,                      intent(in)    :: ind
    type(ions_t),                 intent(in)    :: ions
    type(derivatives_t),          intent(in)    :: der
    type(poisson_t),              intent(in)    :: psolver
    type(orbitalset_t),           intent(inout) :: os(:) !> inout as this is also in orbs
    integer,                      intent(in)    :: nos, maxnorbs
    FLOAT,                        intent(in)    :: rcut
    type(distributed_t),          intent(in)    :: kpt
    logical,                      intent(in)    :: has_phase
    integer,                      intent(in)    :: sm_poisson


    type(lattice_iterator_t) :: latt_iter
    FLOAT :: xat(ions%space%dim), xi(ions%space%dim)
    FLOAT :: rr
    integer :: inn, ist, jst
    integer :: np_sphere, ip, ios
    type(submesh_t) :: sm
    FLOAT, allocatable :: tmp(:), vv(:), nn(:)
    FLOAT, allocatable :: orb(:,:,:)
    FLOAT, parameter :: TOL_INTERSITE = CNST(1.e-5)
    type(distributed_t) :: dist

    PUSH_SUB(orbitalset_init_intersite)

    call messages_print_stress(stdout, "Intersite Coulomb integrals")

    !From the formula and the coding point of view, it is interesting to split the intersite terms
    !into the neighbors (between two different atoms) and periodic copies (same atom).
    !In the last case, due to periodicity of Bloch states, projections on the atomic orbitals
    !or its periodic copy is the same, making easier to implement the formulas

    this%nneighbors = 0
    if(this%iatom /= -1)then
      xat = ions%pos(:, this%iatom)

      latt_iter = lattice_iterator_t(ions%latt, rcut)

      !We first count first the number of neighboring atoms at a distance max rcut 
      do ios = 1, nos
        do inn = 1, latt_iter%n_cells

          xi = os(ios)%sphere%center(1:ions%space%dim) + latt_iter%get(inn)
          rr = norm2(xi - xat)

          !This atom is too far
          if( rr >rcut + TOL_INTERSITE ) cycle
          !Intra atomic interaction
          if( ios == ind .and. rr < TOL_INTERSITE) cycle

          this%nneighbors = this%nneighbors +1
        end do
      end do

      !The first three values are the position of the periodic copies
      !and the zero value one is used to store the actual value of V_ij
      SAFE_ALLOCATE(this%V_ij(1:this%nneighbors, 0:ions%space%dim+1))
      this%V_ij(1:this%nneighbors, 0:ions%space%dim+1) = M_ZERO
      SAFE_ALLOCATE(this%map_os(1:this%nneighbors))
      this%map_os(1:this%nneighbors) = 0
      if(has_phase) then
        SAFE_ALLOCATE(this%phase_shift(1:this%nneighbors, kpt%start:kpt%end))
      end if
 
      this%nneighbors = 0
      do ios = 1, nos
        do inn = 1, latt_iter%n_cells
          xi = os(ios)%sphere%center(1:ions%space%dim) + latt_iter%get(inn)
          rr = norm2(xi - xat)

          if( rr > rcut + TOL_INTERSITE ) cycle
          if( ios == ind .and. rr < TOL_INTERSITE) cycle

          this%nneighbors = this%nneighbors +1

          this%V_ij(this%nneighbors, 1:ions%space%dim) = xi(1:ions%space%dim) -os(ios)%sphere%center(1:ions%space%dim)
          this%V_ij(this%nneighbors, ions%space%dim+1) = rr
          
          this%map_os(this%nneighbors) = ios
        end do
      end do

      write(message(1),'(a, i3, a)')    'Intersite interaction will be computed for ', this%nneighbors, ' neighboring atoms.'
      call messages_info(1)


      SAFE_ALLOCATE(this%coulomb_IIJJ(1:this%norbs,1:this%norbs,1:maxnorbs,1:maxnorbs,1:this%nneighbors))
      this%coulomb_IIJJ = M_ZERO

      call distributed_nullify(dist, this%nneighbors)
#ifdef HAVE_MPI
      if(.not. der%mesh%parallel_in_domains) then
        call distributed_init(dist, this%nneighbors, MPI_COMM_WORLD, 'orbs')
      end if
#endif

      do inn = dist%start, dist%end

        ios = this%map_os(inn)

        !Init a submesh from the union of two submeshes
        call submesh_merge(sm, ions%space, der%mesh, this%sphere, os(ios)%sphere, &
                       shift = this%V_ij(inn, 1:ions%space%dim))

        write(message(1),'(a, i3, a, f6.3, a, i5, a)') 'Neighbor ', inn, ' is located at ', &
                             this%V_ij(inn, ions%space%dim+1), ' Bohr and has ', sm%np, ' grid points.'
        call messages_info(1)

        SAFE_ALLOCATE(orb(1:sm%np, 1:max(this%norbs,os(ios)%norbs),1:2))
        SAFE_ALLOCATE(nn(1:sm%np))
        SAFE_ALLOCATE(vv(1:sm%np))

        do ist = 1, this%norbs
          call datomic_orbital_get_submesh_safe(this%spec, sm, this%ii, this%ll, ist-1-this%ll, &
                          1, orb(1:sm%np, ist,1))
        end do
         
        call submesh_shift_center(sm, ions%space, this%V_ij(inn, 1:ions%space%dim)+os(ios)%sphere%center(1:ions%space%dim))

        do ist = 1, os(ios)%norbs
          call datomic_orbital_get_submesh_safe(os(ios)%spec, sm, os(ios)%ii, os(ios)%ll, ist-1-os(ios)%ll, &
                1, orb(1:sm%np, ist,2))
        end do

        SAFE_ALLOCATE(tmp(1:sm%np))

        !Build information needed for the direct Poisson solver on the submesh
        call submesh_build_global(sm, ions%space)

        select case (sm_poisson)
        case(SM_POISSON_DIRECT)
          call poisson_init_sm(this%poisson, namespace, ions%space, psolver, der, sm, method = POISSON_DIRECT_SUM)
        case(SM_POISSON_ISF)
          call poisson_init_sm(this%poisson, namespace, ions%space, psolver, der, sm, method = POISSON_ISF)
        case(SM_POISSON_PSOLVER)
          call poisson_init_sm(this%poisson, namespace, ions%space, psolver, der, sm, method = POISSON_PSOLVER)
        case(SM_POISSON_FFT)
          call poisson_init_sm(this%poisson, namespace, ions%space, psolver, der, sm, method = POISSON_FFT)
        end select
 
        np_sphere = sm%np

        do ist = 1, this%norbs
            !$omp parallel do
            do ip = 1, np_sphere
              nn(ip) = orb(ip,ist,1)*orb(ip,ist,1)
            end do
            !$omp end parallel do    

            !Here it is important to use a non-periodic poisson solver, e.g. the direct solver
            call dpoisson_solve_sm(this%poisson, sm, vv(1:np_sphere), nn(1:np_sphere))

            do jst = 1, os(ios)%norbs

                !$omp parallel do
                do ip = 1, np_sphere
                  tmp(ip) = vv(ip)*orb(ip, jst, 2)*orb(ip, jst, 2)
                end do
                !$omp end parallel do

                this%coulomb_IIJJ(ist, ist, jst, jst, inn) = dsm_integrate(der%mesh, sm, tmp(1:np_sphere), reduce = .false.)
                if(abs(this%coulomb_IIJJ(ist, ist, jst, jst, inn)) < CNST(1.0e-12)) then
                  this%coulomb_IIJJ(ist, ist, jst, jst, inn) = M_ZERO
                end if

            end do !kst 
        end do !ist

        call poisson_end(this%poisson)

        SAFE_DEALLOCATE_A(nn)
        SAFE_DEALLOCATE_A(vv)
        SAFE_DEALLOCATE_A(tmp)

        call submesh_end_global(sm)
        call submesh_end_cube_map(sm)
        call submesh_end(sm)
        SAFE_DEALLOCATE_A(orb)
      end do !inn

      if(der%mesh%parallel_in_domains) then
        call der%mesh%allreduce(this%coulomb_IIJJ)
      end if
 

      if(dist%parallel) then
        call comm_allreduce(dist%mpi_grp, this%coulomb_IIJJ)
      end if
      
      #ifdef HAVE_MPI
      call distributed_end(dist)
      #endif

    end if

    call messages_print_stress(stdout)

    POP_SUB(orbitalset_init_intersite)
  end subroutine orbitalset_init_intersite



#include "undef.F90"
#include "real.F90"
#include "orbitalset_utils_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "orbitalset_utils_inc.F90"

end module orbitalset_utils_oct_m
