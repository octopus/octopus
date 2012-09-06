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

  use batch_m
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
  use mesh_batch_m
  use mesh_function_m
  use messages_m
  use modelmb_particles_m
  use modelmb_density_matrix_m
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

  !> project out states with proper symmetry for cases which are of symmetry = unknown
  subroutine modelmb_sym_state(eigenval, iunit, gr, mm, &
             modelmbparticles, ncombo, young_used, wf, symmetries_satisfied)
    FLOAT,                    intent(in)    :: eigenval
    integer,                  intent(in)    :: iunit
    type(grid_t),             intent(in)    :: gr
    integer,                  intent(in)    :: mm
    type(modelmb_particle_t), intent(in)    :: modelmbparticles
    integer,                  intent(in)    :: ncombo
    integer,                  intent(inout) :: young_used(1:ncombo)
    CMPLX,                    intent(inout) :: wf(1:gr%mesh%np) !< will be antisymmetrized on output
    logical,                  intent(out)   :: symmetries_satisfied

    integer :: npptype
    integer :: iyoung
    integer :: itype
    integer :: ikeeppart
    integer :: nspindown, nspinup
    integer :: idiagram_combo

    type(modelmb_1part_t) :: mb_1part
    type(young_t) :: young

    integer, allocatable :: dg_combo_iy(:,:)
    integer, allocatable :: dg_combo_ndown(:,:)
    integer, allocatable :: sym_ok_alltypes(:)
 
    FLOAT :: norm

    CMPLX, allocatable  :: antisymwf(:,:,:)              ! stores progressive steps of antisymmetrized wf


    PUSH_SUB(modelmb_sym_state)

    symmetries_satisfied = .false.

    call young_nullify (young)

    call modelmb_1part_nullify(mb_1part)

    SAFE_ALLOCATE(sym_ok_alltypes(1:modelmbparticles%ntype_of_particle))

! set up combinations of young diagrams (1 for each type)

    SAFE_ALLOCATE(dg_combo_ndown(1:modelmbparticles%ntype_of_particle, 1:ncombo))
    dg_combo_ndown = 0
    SAFE_ALLOCATE(dg_combo_iy(1:modelmbparticles%ntype_of_particle, 1:ncombo))
    dg_combo_iy = 0

    idiagram_combo = 1
    do itype = 1, modelmbparticles%ntype_of_particle
      ikeeppart = modelmbparticles%particles_of_type(1,itype)
      if (modelmbparticles%bosonfermion(ikeeppart) /= 1) then ! 1 is for fermion - might introduce a parameter in modelmb_particles
        dg_combo_ndown(itype, idiagram_combo) = 0
        dg_combo_iy(itype, idiagram_combo) = 1
        idiagram_combo = idiagram_combo + 1
        cycle
      end if
      npptype = modelmbparticles%nparticles_per_type(itype)
      do nspindown = 0, floor(npptype/2.)
        nspinup = npptype - nspindown
        call young_init (young, nspinup, nspindown)
        do iyoung = 1, young%nyoung
          dg_combo_ndown(itype, idiagram_combo) = nspindown
          dg_combo_iy(itype, idiagram_combo) = iyoung
          idiagram_combo = idiagram_combo + 1
        end do
        call young_end (young)
      end do
    end do

    SAFE_ALLOCATE(antisymwf(1:gr%mesh%np,1,1))

    ! index for combination of Young diagrams 
    do idiagram_combo = 1, ncombo
      antisymwf(:,1,1) = wf(1:gr%mesh%np)
      ! skip diagram combinations already used in present degenerate subspace
      if (young_used (idiagram_combo) > 0) cycle

      call modelmb_sym_state_1diag(eigenval, iunit, gr, mm, &
         modelmbparticles, dg_combo_ndown(:, idiagram_combo), &
         dg_combo_iy(:, idiagram_combo), &
         antisymwf, sym_ok_alltypes, norm)
     
      ! test the overall symmetrization (no 0.0 norms for present combination of Young diagrams)
      ! check if all types of particles have been properly symmetrized
      if (sum(sym_ok_alltypes) == modelmbparticles%ntype_of_particle .and. abs(norm) > 1.e-6) then
        wf(:) = antisymwf(:,1,1) / sqrt(norm)
        symmetries_satisfied = .true.
        young_used (idiagram_combo) = 1
        ! eventually exit the combo loop
        exit
      end if

    end do ! idiagram_combo



    SAFE_DEALLOCATE_A(sym_ok_alltypes)
    SAFE_DEALLOCATE_A(antisymwf)
    SAFE_DEALLOCATE_A(dg_combo_ndown)
    SAFE_DEALLOCATE_A(dg_combo_iy)

    POP_SUB(modelmb_sym_state)
  end subroutine modelmb_sym_state


  ! ---------------------------------------------------------
  !> project out states for a single combination of Young diagrams (1 diagram for each particle type)
  subroutine modelmb_sym_state_1diag(eigenval, iunit, gr, mm, &
             modelmbparticles, nspindown_in, iyoung_in, &
             antisymwf, sym_ok_alltypes, norm)
    FLOAT,                    intent(in)    :: eigenval
    integer,                  intent(in)    :: iunit
    type(modelmb_particle_t), intent(in)    :: modelmbparticles
    type(grid_t),             intent(in)    :: gr
    integer,                  intent(in)    :: mm
    integer,                  intent(in)    :: nspindown_in(1:modelmbparticles%ntype_of_particle)
    integer,                  intent(in)    :: iyoung_in(1:modelmbparticles%ntype_of_particle)

    CMPLX,                    intent(inout) :: antisymwf(1:gr%mesh%np,1,1) !< will be antisymmetrized on output
    integer,                  intent(out)   :: sym_ok_alltypes(1:modelmbparticles%ntype_of_particle)
    FLOAT,                    intent(out)   :: norm

! local vars
    integer :: ipart1, npptype
    integer :: ikeeppart, itype
    integer :: ndimmb
    integer :: ofst_so_far
    integer :: nspinup

    type(permutations_t) :: perms_up, perms_down
    type(modelmb_1part_t) :: mb_1part
    type(young_t) :: young

    type(batch_t) :: antisymwfbatch

    integer, allocatable :: ix(:), ofst(:)
    integer, allocatable :: p_of_type_up(:), p_of_type_down(:)
    character(len=500) :: tmpstring, youngstring
 
    FLOAT :: normalizer
    CMPLX :: wfdotp(1,1)

    PUSH_SUB(modelmb_sym_state_1diag)

    call permutations_nullify(perms_up)
    call permutations_nullify(perms_down)
    call young_nullify (young)

    call modelmb_1part_nullify(mb_1part)

    SAFE_ALLOCATE(ix(1:MAX_DIM))

    sym_ok_alltypes = 0

    ndimmb=modelmbparticles%ndim

    normalizer = product(gr%mesh%spacing(1:gr%mesh%sb%dim)) !1/units_out%length**gr%mesh%sb%dim

    ! this is in case _none_ of the particles is symetrized,
    !   then the whole itype loop is skipped
    norm = M_ONE

    youngstring = ""
    tmpstring = ""

    ofst_so_far = 0
    ! for each particle type
    do itype = 1, modelmbparticles%ntype_of_particle

      ! FIXME: for multiple particle types this needs to be fixed.
      ! Also, for inequivalent spin configurations this should vary, and we get
      ! different 1 body densities, no?
      ikeeppart = modelmbparticles%particles_of_type(1,itype)

      ! if the particle is not fermionic, just cycle to next one
      ! FIXME: boson case is not treated yet
      if (modelmbparticles%bosonfermion(ikeeppart) /= 1) then ! 1 = fermion
        sym_ok_alltypes(itype) = 1
        return 
      end if

      call modelmb_1part_init(mb_1part, gr%mesh, ikeeppart, &
             modelmbparticles%ndim, gr%sb%box_offset)

      npptype = modelmbparticles%nparticles_per_type(itype)

      SAFE_ALLOCATE(ofst(1:npptype))
      do ipart1 = 1, npptype
        ofst(ipart1) = ofst_so_far + (ipart1 - 1) * ndimmb
      end do
      ofst_so_far = ofst_so_far + npptype*ndimmb

      ! note: use of spin nomenclature is just for visualization, no real spin
      ! here.
      nspinup = npptype - nspindown_in(itype)

      call permutations_init(nspinup,perms_up)
      call permutations_init(nspindown_in(itype),perms_down)

      ! generate all Young diagrams, decorated, for this distribution of up
      ! and downs
      call young_init (young, nspinup, nspindown_in(itype))

      ! loop over all Young diagrams for present distribution of spins up and down
      SAFE_ALLOCATE(p_of_type_up(nspindown_in(itype)))
      SAFE_ALLOCATE(p_of_type_down(nspindown_in(itype)))
      p_of_type_up   = modelmbparticles%particles_of_type(young%young_up  (1:nspindown_in(itype), iyoung_in(itype)), itype)
      p_of_type_down = modelmbparticles%particles_of_type(young%young_down(1:nspindown_in(itype), iyoung_in(itype)), itype)

      call modelmb_sym_down (ndimmb, npptype, &
             ofst, young, p_of_type_up, p_of_type_down, gr, normalizer, antisymwf)

      SAFE_DEALLOCATE_A(p_of_type_up)
      SAFE_DEALLOCATE_A(p_of_type_down)

      call modelmb_antisym_up (iyoung_in(itype), perms_up, ndimmb, npptype, ofst, &
             young, gr, normalizer, antisymwf)
      call modelmb_antisym_down (iyoung_in(itype), perms_down, ndimmb, npptype, ofst, &
             young, gr, antisymwf)

      call permutations_end(perms_up)
      call permutations_end(perms_down)
      call young_end (young)
 
      SAFE_DEALLOCATE_A(ofst)
      call modelmb_1part_end(mb_1part)

      write (tmpstring, '(3x,I4,1x,I4)') iyoung_in(itype), nspindown_in(itype)
      youngstring = trim(youngstring) // trim(tmpstring)

      if (gr%mesh%parallel_in_domains) then
        call batch_init(antisymwfbatch, 1, 1, 1, antisymwf)
        call zmesh_batch_dotp_self(gr%mesh, antisymwfbatch, wfdotp, reduce=.true.)
        norm = TOFLOAT(wfdotp(1,1))
        call batch_end(antisymwfbatch)
      else
        norm = TOFLOAT(sum(conjg(antisymwf(:,1,1))*antisymwf(:,1,1))) * normalizer
      end if

      if (abs(norm) > 1.e-6) sym_ok_alltypes(itype) = 1
    end do ! itype

    write (iunit, '(a,I5,3x,E16.6,5x,E14.6,2x,a)') &
        "  ", mm, eigenval, norm, trim(youngstring)

    POP_SUB(modelmb_sym_state_1diag)
  end subroutine modelmb_sym_state_1diag
  ! ---------------------------------------------------------

!
!> input 1 wave function, and symetrize wrt spin down labeled particles, according to the given young diagrams
  subroutine modelmb_sym_down(ndimmb, npptype, &
             ofst, young, p_of_type_up, p_of_type_down, gr, normalizer, antisymwf)
    integer, intent(in) :: ndimmb
    integer, intent(in) :: npptype
    integer, intent(in) :: ofst(1:npptype)

    type(young_t), intent(in) :: young
    integer, intent(in) :: p_of_type_up(young%ndown)
    integer, intent(in) :: p_of_type_down(young%ndown)
    type(grid_t), intent(in) :: gr

    FLOAT, intent(in) :: normalizer
    CMPLX, intent(inout) :: antisymwf(1:gr%mesh%np,1,1)

! local vars
    integer :: idown, ipart1, ipart2, ip, ipp
    integer, allocatable :: ix(:), ixp(:)
    integer, allocatable :: forward_map_exchange(:)
    FLOAT :: norm
    CMPLX :: wfdotp(1,1)
    CMPLX, allocatable  :: antisymwf_swap(:,:,:) ! single wf term, with correct permutation of particles
    type(batch_t) :: antisymwfbatch
    logical :: debug_antisym = .false.

    PUSH_SUB(modelmb_sym_down)

    SAFE_ALLOCATE(forward_map_exchange(gr%mesh%np_global))
    SAFE_ALLOCATE(ix(1:MAX_DIM))
    SAFE_ALLOCATE(ixp(1:MAX_DIM))
    SAFE_ALLOCATE(antisymwf_swap(1:gr%mesh%np, 1, 1))

    ! first symmetrize over pairs of particles associated in the present
    ! Young diagram
    do idown = 1, young%ndown
      ipart1 = p_of_type_down(idown)
      ipart2 = p_of_type_up(idown)
  

      ! each processor needs the full map of points for send and recv
      do ip = 1, gr%mesh%np_global
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
        ASSERT (ipp <= gr%mesh%np_global)
        forward_map_exchange(ip) = ipp
      end do ! ip

      if (gr%mesh%parallel_in_domains) then
        antisymwf_swap = antisymwf
        ! set up batch type for global exchange operation: 1 state, the loop over MB states is outside this routine
        call batch_init (antisymwfbatch, 1, 1, 1, antisymwf_swap)
        call zmesh_batch_exchange_points (gr%mesh, antisymwfbatch, forward_map=forward_map_exchange)
        call batch_end(antisymwfbatch)
      else
        antisymwf_swap(:,1,1) = antisymwf(forward_map_exchange(:),1,1)
      end if

      antisymwf = M_HALF * (antisymwf_swap + antisymwf)

    end do ! idown
  
    if (debug_antisym) then 
      if (gr%mesh%parallel_in_domains) then
        call batch_init(antisymwfbatch, 1, 1, 1, antisymwf)
        call zmesh_batch_dotp_self(gr%mesh, antisymwfbatch, wfdotp, reduce=.true.)
        norm = TOFLOAT(wfdotp(1,1))
        call batch_end(antisymwfbatch)
      else
        norm = TOFLOAT(sum(conjg(antisymwf(:,1,1))*antisymwf(:,1,1))) * normalizer
      end if
      write (message(1), '(a,I7,a,E20.10)') 'norm of pair-symmetrized-state with ',&
              young%ndown, ' spins down is ', norm
      call messages_info(1)
    end if
    SAFE_DEALLOCATE_A(forward_map_exchange)
    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(ixp)
    SAFE_DEALLOCATE_A(antisymwf_swap)

    POP_SUB(modelmb_sym_down)
  end subroutine modelmb_sym_down
 
  ! ---------------------------------------------------------

  subroutine modelmb_antisym_up(iyoung_in, perms_up, ndimmb, npptype, ofst, young, gr, normalizer, antisymwf)
    integer, intent(in) :: iyoung_in
    integer, intent(in) :: ndimmb
    integer, intent(in) :: npptype
    integer, intent(in) :: ofst(1:npptype)
    type(young_t), intent(in) :: young
    type(grid_t), intent(in) :: gr
    FLOAT, intent(in) :: normalizer
    type(permutations_t), intent(in) :: perms_up
    CMPLX, intent(inout) :: antisymwf(1:gr%mesh%np,1,1)

! local vars
    integer :: iperm_up, iup
    integer :: ipart1, ipart2, ip
    CMPLX :: wfdotp(1,1)
    FLOAT :: norm
    integer, allocatable :: forward_map_exchange(:)
    integer, allocatable :: ix(:), ixp(:)
    CMPLX, allocatable  :: antisymwf_swap(:,:,:) ! single wf term, with correct permutation of particles
    CMPLX, allocatable  :: antisymwf_acc(:,:,:)
    type(batch_t) :: antisymwfbatch
    logical :: debug_antisym = .false.

    PUSH_SUB(modelmb_antisym_up)

    SAFE_ALLOCATE(forward_map_exchange(gr%mesh%np_global))
    SAFE_ALLOCATE(ix(1:MAX_DIM))
    SAFE_ALLOCATE(ixp(1:MAX_DIM))
    SAFE_ALLOCATE(antisymwf_swap(1:gr%mesh%np, 1, 1))
    SAFE_ALLOCATE(antisymwf_acc(1:gr%mesh%np, 1, 1))
    ! for each permutation of particles of this type
    !  antisymmetrize the up labeled spins, amongst themselves
    antisymwf_acc = M_z0
    do iperm_up = 1, perms_up%npermutations

      do ip = 1, gr%mesh%np_global
        ! get present position
        call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)
        ! initialize coordinates for all particles
        ixp = ix
        ! permute the particles labeled spin up 
        do iup = 1, young%nup
          ! get image of ipart1 under permutation iperm1
          ipart1 = young%young_up(iup,iyoung_in)
          ipart2 = young%young_up(perms_up%allpermutations(iup,iperm_up),iyoung_in)
          ixp (ofst(ipart1)+1:ofst(ipart1)+ndimmb) = ix (ofst(ipart2)+1:ofst(ipart2)+ndimmb) ! part1 to 2
        end do
        ! get position of exchanged point
        forward_map_exchange(ip) = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)
      end do ! ip
    
      if (gr%mesh%parallel_in_domains) then
        antisymwf_swap=antisymwf
        ! set up batch type for global exchange operation: 1 state, the loop over MB states is outside this routine
        call batch_init (antisymwfbatch, 1, 1, 1, antisymwf_swap)
        call zmesh_batch_exchange_points (gr%mesh, antisymwfbatch, forward_map=forward_map_exchange)
        call batch_end(antisymwfbatch)
      else
        antisymwf_swap(:,1,1) = antisymwf(forward_map_exchange(:),1,1)
      end if

      antisymwf_acc = antisymwf_acc + perms_up%permsign(iperm_up)*antisymwf_swap

    end do ! iperm_up

    antisymwf = antisymwf_acc / TOFLOAT(perms_up%npermutations)

    ! the following could be removed for production
    if (debug_antisym) then 
      if (gr%mesh%parallel_in_domains) then
        call batch_init(antisymwfbatch, 1, 1, 1, antisymwf)
        call zmesh_batch_dotp_self(gr%mesh, antisymwfbatch, wfdotp, reduce=.true.)
        norm = TOFLOAT(wfdotp(1,1))
        call batch_end(antisymwfbatch)
      else
        norm = TOFLOAT(sum(conjg(antisymwf(:,1,1))*antisymwf(:,1,1))) * normalizer
      end if
      write (message(1), '(a,I7,a,E20.10)') 'norm of up-antisym+pairsym-state with ',&
              young%ndown, ' spins down is ', norm
      call messages_info(1)
    end if
    SAFE_DEALLOCATE_A(forward_map_exchange)
    SAFE_DEALLOCATE_A(antisymwf_swap)
    SAFE_DEALLOCATE_A(antisymwf_acc)
    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(ixp)

    POP_SUB(modelmb_antisym_up)

  end subroutine modelmb_antisym_up

  ! ---------------------------------------------------------

  subroutine modelmb_antisym_down (iyoung_in, perms_down, ndimmb, npptype, ofst, young, gr, antisymwf)
    integer, intent(in) :: iyoung_in
    integer, intent(in) :: ndimmb
    integer, intent(in) :: npptype
    integer, intent(in) :: ofst(1:npptype)
    type(permutations_t), intent(in) :: perms_down
    type(young_t), intent(in) :: young
    type(grid_t), intent(in) :: gr
    CMPLX, intent(inout) :: antisymwf(1:gr%mesh%np,1,1)

! local vars
    integer :: iperm_down, idown
    integer :: ipart1, ipart2, ip

    integer, allocatable :: forward_map_exchange(:)
    integer, allocatable :: ix(:), ixp(:)
    CMPLX, allocatable  :: antisymwf_swap(:,:,:) !< single wf term, with correct permutation of particles
    CMPLX, allocatable  :: antisymwf_acc(:,:,:)
    type(batch_t) :: antisymwfbatch

    PUSH_SUB(modelmb_antisym_down)

    SAFE_ALLOCATE(forward_map_exchange(gr%mesh%np_global))
    SAFE_ALLOCATE(ix(1:MAX_DIM))
    SAFE_ALLOCATE(ixp(1:MAX_DIM))
    SAFE_ALLOCATE(antisymwf_swap(1:gr%mesh%np, 1, 1))
    SAFE_ALLOCATE(antisymwf_acc(1:gr%mesh%np, 1, 1))


    ! for each permutation of particles of this type
    !  antisymmetrize the down labeled spins, amongst themselves
    antisymwf_acc = M_z0
    do iperm_down = 1, perms_down%npermutations
      do ip = 1, gr%mesh%np_global
        ! get present position
        call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)
        ! initialize coordinates for all particles
        ixp = ix
        ! permute the particles labeled spin down
        do idown = 1, young%ndown
          ! get image of ipart1 under permutation iperm1
          ipart1 = young%young_down(idown,iyoung_in)
          ipart2 = young%young_down(perms_down%allpermutations(idown,iperm_down),iyoung_in)
          ixp (ofst(ipart1)+1:ofst(ipart1)+ndimmb) = ix (ofst(ipart2)+1:ofst(ipart2)+ndimmb) ! part1 to 2
        end do
        ! get position of exchanged point
        forward_map_exchange(ip) = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)
      end do ! ip

      if (gr%mesh%parallel_in_domains) then
        antisymwf_swap = antisymwf
        ! set up batch type for global exchange operation: 1 state, the loop over MB states is outside this routine
        call batch_init (antisymwfbatch, 1, 1, 1, antisymwf_swap)
        call zmesh_batch_exchange_points (gr%mesh, antisymwfbatch, forward_map=forward_map_exchange)
        call batch_end(antisymwfbatch)
      else
        antisymwf_swap(:,1,1) = antisymwf(forward_map_exchange(:),1,1)
      end if

      antisymwf_acc = antisymwf_acc + perms_down%permsign(iperm_down)*antisymwf_swap

    end do ! iperm_down

    antisymwf = antisymwf_acc / TOFLOAT(perms_down%npermutations)

    SAFE_DEALLOCATE_A(forward_map_exchange)
    SAFE_DEALLOCATE_A(antisymwf_swap)
    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(ixp)

    POP_SUB(modelmb_antisym_down)

  end subroutine modelmb_antisym_down

  ! ---------------------------------------------------------

end module modelmb_exchange_syms_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
