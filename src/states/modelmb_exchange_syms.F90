!! Copyright (C) 2009 N. Helbig and M. Verstraete
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
!! $Id: states.F90 5022 2009-03-03 17:47:58Z nitsche $

#include "global.h"

module modelmb_exchange_syms_m

  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use hypercube_m
  use index_m
  use io_m
  use lalg_adv_m
  use loct_m
  use math_m
  use mesh_function_m
  use messages_m
  use modelmb_particles_m
  use modelmb_density_m
  use modelmb_1part_m
  use mpi_m
  use mpi_lib_m
  use parser_m
  use permutations_m
  use profiling_m
  use states_m
  use young_m

  implicit none

  private

  public :: modelmb_sym_state

contains

  ! project out states with proper symmetry for cases which are of symmetry = unknown
  subroutine modelmb_sym_state(eigenval, iyoungstart, iunit, gr, mm, geo, &
             modelmbparticles, wf, symmetries_satisfied)
    FLOAT,                    intent(in)    :: eigenval
    integer,                  intent(in)    :: iyoungstart
    integer,                  intent(in)    :: iunit
    type(grid_t),             intent(in)    :: gr
    integer,                  intent(in)    :: mm
    type(geometry_t),         intent(in)    :: geo
    type(modelmb_particle_t), intent(in)    :: modelmbparticles
    CMPLX,                    intent(inout) :: wf(1:gr%mesh%np_part_global) ! will be antisymmetrized on output
    logical,                  intent(out)   :: symmetries_satisfied

    integer :: itype, ipart1, ipart2, npptype
    integer :: ip, ipp, iup, idown, iyoung, iyoungstart_eff
    integer :: ikeeppart
    integer :: ndimmb
    integer :: nspindown, nspinup, iperm_up, iperm_down

    type(permutations_t) :: perms_up, perms_down
    type(modelmb_1part_t) :: mb_1part
    type(young_t) :: young

    integer, allocatable :: ix(:), ixp(:), ofst(:)
    integer, allocatable :: symmetries_satisfied_alltypes(:)
 
    FLOAT :: normalizer, norm

    FLOAT, allocatable  :: antisymrho(:)
    CMPLX, allocatable  :: antisymwf(:)
    CMPLX, allocatable  :: antisymwf_swap(:)

    logical :: debug_antisym

    PUSH_SUB(modelmb_sym_states)

    symmetries_satisfied = .false.

    debug_antisym = .false.


    call permutations_nullify(perms_up)
    call permutations_nullify(perms_down)
    call young_nullify (young)

    call modelmb_1part_nullify(mb_1part)

    SAFE_ALLOCATE(ix(1:MAX_DIM))
    SAFE_ALLOCATE(ixp(1:MAX_DIM))

    SAFE_ALLOCATE(symmetries_satisfied_alltypes(1:modelmbparticles%ntype_of_particle))
    symmetries_satisfied_alltypes = 0

    ndimmb=modelmbparticles%ndim

    normalizer = product(gr%mesh%spacing(1:gr%mesh%sb%dim)) !1/units_out%length**gr%mesh%sb%dim

    ! this is in case _none_ of the particles is symetrized,
    !   then the whole itype loop is skipped
    norm = M_ONE

    ! FIXME: could one of these arrays be avoided?
    SAFE_ALLOCATE(antisymwf(1:gr%mesh%np_part_global))
    SAFE_ALLOCATE(antisymwf_swap(1:gr%mesh%np_part_global))

    ! for each particle type
    do itype = 1, modelmbparticles%ntype_of_particle

      ! FIXME: for multiple particle types this needs to be fixed.
      ! Also, for inequivalent spin configurations this should vary, and we get
      ! different 1 body densities, no?
      ikeeppart = modelmbparticles%particles_of_type(1,itype)

      ! if the particle is not fermionic, just cycle to next one
      ! FIXME: boson case is not treated yet
      if (modelmbparticles%bosonfermion(ikeeppart) /= 1) then ! 1 = fermion
        symmetries_satisfied_alltypes(itype) = 1
        cycle
      end if

      call modelmb_1part_init(mb_1part, gr%mesh, ikeeppart, &
             modelmbparticles%ndim, gr%sb%box_offset)

      npptype = modelmbparticles%nparticles_per_type(itype)

      SAFE_ALLOCATE(antisymrho(1:gr%mesh%np_part_global))
      SAFE_ALLOCATE(ofst(1:npptype))
      do ipart1 = 1, npptype
        ofst(ipart1) = (ipart1-1)*ndimmb
      end do

      ! note: use of spin nomenclature is just for visualization, no real spin
      ! here.
      ! FIXME: will not work for several particle types!
      do nspindown = 0, floor(npptype/2.)
        nspinup = npptype - nspindown

        call permutations_init(nspinup,perms_up)
        call permutations_init(nspindown,perms_down)

        ! generate all Young diagrams, decorated, for this distribution of up
        ! and downs
        call young_init (young, nspinup, nspindown)

        ! loop over all Young diagrams for present distribution of spins up and down
        iyoungstart_eff = iyoungstart
        if (iyoungstart_eff > young%nyoung) iyoungstart_eff = 1
        do iyoung = iyoungstart_eff, young%nyoung
          antisymwf = M_z0
          antisymwf(:) = wf(:)

          ! first symmetrize over pairs of particles associated in the present
          ! Young diagram
          do idown = 1, young%ndown
            ipart1 = modelmbparticles%particles_of_type(young%young_down(idown,iyoung), itype)
            ipart2 = modelmbparticles%particles_of_type(  young%young_up(idown,iyoung), itype)
  
            antisymwf_swap=antisymwf
            do ip = 1, gr%mesh%np_part_global
              ! get present position
              call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)
       
              ! invert coordinates of ipart1 and ipart2
              ixp = ix
              ! permutate the particles ipart1 and its spin up partner
              ixp (ofst(ipart1)+1:ofst(ipart1)+ndimmb) = &
                  ix (ofst(ipart2)+1:ofst(ipart2)+ndimmb)
              ixp (ofst(ipart2)+1:ofst(ipart2)+ndimmb) = &
                  ix (ofst(ipart1)+1:ofst(ipart1)+ndimmb)
              
              ! get position of exchanged point
              ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)
  
              antisymwf_swap(ip)=antisymwf_swap(ip) + antisymwf(ipp)
          
            end do ! ip
            antisymwf = antisymwf_swap * M_HALF
          end do
    
          if (debug_antisym) then 
            antisymrho = TOFLOAT(conjg(antisymwf)*antisymwf) * normalizer
            norm = sum(antisymrho)
            write (message(1), '(a,I7,a,I7,a,E20.10)') 'norm of pair-symmetrized-state ',&
                    mm, ' with ', nspindown, ' spins down is ', norm
            call messages_info(1)
          end if
 
          ! for each permutation of particles of this type
          !  antisymmetrize the up and down labeled spins, amongst themselves
          antisymwf_swap = M_z0
          do iperm_up = 1, perms_up%npermutations
            do ip = 1, gr%mesh%np_part_global
              ! get present position
              call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)
              ! initialize coordinates for all particles
              ixp = ix
              ! permute the particles labeled spin up 
              do iup = 1, nspinup
                ! get image of ipart1 under permutation iperm1
                ipart1 = young%young_up(iup,iyoung)
                ipart2 = young%young_up(perms_up%allpermutations(iup,iperm_up),iyoung)
                ixp (ofst(ipart1)+1:ofst(ipart1)+ndimmb) = ix (ofst(ipart2)+1:ofst(ipart2)+ndimmb) ! part1 to 2
              end do
              ! get position of exchanged point
              ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)
              antisymwf_swap(ip)=antisymwf_swap(ip) + perms_up%permsign(iperm_up)*antisymwf(ipp)
            end do ! ip
          end do ! iperm_up

          antisymwf=antisymwf_swap / TOFLOAT(perms_up%npermutations)

          ! the following could be removed for production
          if (debug_antisym) then 
            antisymrho = TOFLOAT(conjg(antisymwf)*antisymwf) * normalizer
            norm = sum(antisymrho)
            write (message(1), '(a,I7,a,I7,a,E20.10)') 'norm of up-antisym+pairsym-state ',&
                    mm, ' with ', nspindown, ' spins down is ', norm
            call messages_info(1)
          end if

          antisymwf_swap = M_z0
          do iperm_down = 1, perms_down%npermutations
            do ip = 1, gr%mesh%np_part_global
              ! get present position
              call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)
              ! initialize coordinates for all particles
              ixp = ix
              ! permute the particles labeled spin down
              do idown = 1, nspindown
                ! get image of ipart1 under permutation iperm1
                ipart1 = young%young_down(idown,iyoung)
                ipart2 = young%young_down(perms_down%allpermutations(idown,iperm_down),iyoung)
                ixp (ofst(ipart1)+1:ofst(ipart1)+ndimmb) = ix (ofst(ipart2)+1:ofst(ipart2)+ndimmb) ! part1 to 2
              end do
              ! get position of exchanged point
              ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)
              antisymwf_swap(ip)=antisymwf_swap(ip) + perms_down%permsign(iperm_down)*antisymwf(ipp)
            end do ! ip
          end do ! iperm_down

          antisymwf=antisymwf_swap / TOFLOAT(perms_down%npermutations)
     
          antisymrho = TOFLOAT(conjg(antisymwf)*antisymwf)
          norm = sum(antisymrho) * normalizer

          if (norm > CNST(1.e-5)) &
            write (iunit, '(I5,3x,E16.6,2x,I4,2x,I7,8x,I3,5x,E14.6)') &
                  mm, eigenval, itype, iyoung, nspindown, norm
   
          ! FIXME: this is probably the problem when running in single precision
          if (norm < CNST(1.e-3)) cycle
  
          ! we have found a valid Young diagram and antisymmetrized the present
          ! state enough. Loop to next type
          symmetries_satisfied_alltypes(itype) = 1
          exit

        end do ! iyoung loop
   
        call permutations_end(perms_up)
        call permutations_end(perms_down)
        call young_end (young)
 
        ! found a good Young diagram for this particle type.
        ! save antisymmetrized copy and go on to next particle type
        if (symmetries_satisfied_alltypes(itype) == 1) then
          wf = antisymwf
          exit
        end if

      end do ! nspindown=0,1...
   
      SAFE_DEALLOCATE_A(antisymrho)
      SAFE_DEALLOCATE_A(ofst)
      call modelmb_1part_end(mb_1part)

    end do ! itype

    ! check if all types of particles have been properly symmetrized
    if (sum(symmetries_satisfied_alltypes) == modelmbparticles%ntype_of_particle) then
      wf = wf / sqrt(norm)
      symmetries_satisfied = .true.
    end if

    SAFE_DEALLOCATE_A(antisymwf_swap)
    SAFE_DEALLOCATE_A(antisymwf)

    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(ixp)

    POP_SUB(modelmb_sym_states)
  end subroutine modelmb_sym_state
  ! ---------------------------------------------------------

end module modelmb_exchange_syms_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
