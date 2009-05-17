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
  use io_function_m
  use io_m
  use lalg_adv_m
  use loct_m
  use loct_parser_m
  use math_m
  use mesh_function_m
  use messages_m
  use modelmb_particles_m
  use modelmb_1part_m
  use mpi_m
  use mpi_lib_m
  use permutations_m
  use profiling_m
  use states_m
  use units_m
  use young_m

  implicit none

  private

  public :: modelmb_find_exchange_syms, &
            modelmb_sym_state

contains


  ! ---------------------------------------------------------
  !  this routine uses a stochastic method to see if a bunch 
  !  of points give psi(1,2) = +- psi(2,1)
  !  may be inadequate if the symmetry is approximate.
  !  Other option is to force-antisymmetrize (see routine modelmb_sym_states)
  ! ---------------------------------------------------------
  subroutine modelmb_find_exchange_syms(gr, mm, modelMBparticles, wf)

    implicit none

    type(grid_t),             intent(in)   :: gr
    integer,                  intent(in)   :: mm
    CMPLX,                    intent(in)   :: wf(1:gr%mesh%np_part_global)
    type(modelMB_particle_t), intent(inout):: modelMBparticles


    ! local vars
    integer :: itype, ipart1, ipart2, npptype
    integer :: ip, ipp, maxtries, enoughgood
    integer :: ntries, ngoodpoints, nratiocomplaints
    integer :: ofst1, ofst2, ndimMB
    integer :: wf_sign_change, save_sign_change

    integer, allocatable :: ix(:), ixp(:)

    FLOAT :: ratio, wftol, rand0to1, tolonsign
    CMPLX :: wfval1, wfval2

    ! source code

    call push_sub('states.modelmb_find_exchange_syms')

    SAFE_ALLOCATE(ix(1:MAX_DIM))
    SAFE_ALLOCATE(ixp(1:MAX_DIM))

    ndimMB=modelMBparticles%ndim_modelMB

    ! this is the tolerance on the wavefunction in order to use it to determine
    ! sign changes. Should be checked, and might be a function of grid number,
    ! as the norm is diluted over the whole cell...
    ! maybe min(0.001,100*1./ngridpoints) would be better.
    wftol=1.d-5

    ! tolerance on wf(1,2)/wf(2,1) being close to +-1
    tolonsign = 1.d-5

    ! maximum number of times we try to extract the parity
    maxtries = 300
    enoughgood = 20

    ! initialize the container
    modelMBparticles%exchange_symmetry=0

    ! for each particle type
    do itype=1,modelMBparticles%ntype_of_particle_modelMB

      npptype=modelMBparticles%nparticles_per_type(itype)
      nratiocomplaints = 0

      ! for each pair of particles of this type
      do ipart1=1,npptype
        do ipart2=ipart1+1,npptype
 
          ! try a bunch of points in the wf
          ngoodpoints=0
          ntries=0
          save_sign_change=0
          do while (ntries < maxtries .and. ngoodpoints < enoughgood)
            
            ! count real trials of random points (avoids eventual infinite looping)
            ntries = ntries+1

            ! find coordinates of random point
            call random_number(rand0to1)
            ip = floor(gr%mesh%np_part_global*rand0to1)+1
 
            ! wavefunction value at this point ip
            wfval1 = wf(ip)

            ! see if point has enough weight to be used
            if (abs(real(wfval1)) < wftol) cycle

            ! get actual coordinates of MB point
            call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)

            ! offsets in full dimensional space for the particles we care about
            ofst1=(ipart1-1)*ndimMB
            ofst2=(ipart2-1)*ndimMB

            ! invert coordinates of ipart1 and ipart2
            ixp = ix
            ixp (ofst1+1:ofst1+ndimMB) = ix (ofst2+1:ofst2+ndimMB) ! part1 to 2
            ixp (ofst2+1:ofst2+ndimMB) = ix (ofst1+1:ofst1+ndimMB) ! part2 to 1
 
            ! if we are on a diagonal, cycle
            if (all(ix-ixp == 0)) cycle
 
            ! recover MB index for new point
            ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)

            ! wavefunction value at exchanged coordinate point ipp
            wfval2 = wf(ipp)
 
            ! check sign change FIXME: check imaginary part too
            !   normally, as these points are related, wfval2 should also be > wftol
            ratio = real(wfval1) / real(wfval2)
            
            if (abs(ratio+1.0d0) < tolonsign) then
              wf_sign_change = -1
            else if (abs(ratio-1.0d0) < tolonsign) then
              wf_sign_change =  1
            else  ! ratio not clear enough to get parity
              !if (nratiocomplaints < 20) then
              !  write (message(1),'(a,E20.10,a)') 'point skipped for exchange parity: ratio of wf ',ratio, &
              !       ' is not clearly +-1'
              !  call write_info(1)
              !end if
              nratiocomplaints = nratiocomplaints+1
            end if
 
            ! if opposite of previous value, complain
            if (save_sign_change /= 0 .and. save_sign_change /= wf_sign_change) then
              write (message(1),'(a)') ' WARNING: sign change is different between wf points examined'
              call write_info(1)
            end if

            save_sign_change=wf_sign_change
            ! ngoodpoints counts points for which we were able to compare wf values 
            ngoodpoints = ngoodpoints + 1
          end do !  tries and goodpoints
         
          ! if we have found enough points with wf large enough to compare
          !  and the number of points where the ratio was far from +-1 was not
          !  too high
          if (ngoodpoints == enoughgood .and. &
              real(nratiocomplaints)/real(ngoodpoints) < 0.1) then
            modelMBparticles%exchange_symmetry(ipart1,ipart2,itype)=save_sign_change
            modelMBparticles%exchange_symmetry(ipart2,ipart1,itype)=save_sign_change
          else 
            !write (message(1),'(a)') ' WARNING: unable to determine sign change for'
            !write (message(2),'(a,I6,a,I6,a,I6)') ' exchange of ', ipart1, ' and ', ipart2, ' of many-body state ', mm
            !call write_info(2)
          end if 
 
        end do ! ipart2
      end do ! ipart1
    end do ! itype

    ! do some crude output
    do itype=1,modelMBparticles%ntype_of_particle_modelMB
      !write (message(1),'(a)')    ' lower triangle of exchange parity matrix'
      write (message(1),'(a,I6,a,I6)') '   for particles of type ', itype, '   in MB state ', mm
      call write_info(1)

      do ipart1=1,npptype
        do ipart2=1,ipart1
          write (*,'(2x, I6)',ADVANCE='NO') modelMBparticles%exchange_symmetry(ipart2,ipart1,itype)
        end do ! ipart2
        write (*,*)
      end do ! ipart1

      modelMBparticles%bosonfermion(itype) = 'unknown'
      if (sum(modelMBparticles%exchange_symmetry) ==  npptype*(npptype-1)) then
        write (message(1),'(a,I6,a,I6,a)') 'For particles of type ', itype, ',  MB state ', &
               mm, ' is totally symmetric (bosonic)'
        call write_info(1)
        modelMBparticles%bosonfermion(itype) = 'boson'
      else if (sum(modelMBparticles%exchange_symmetry) == -npptype*(npptype-1)) then
        write (message(1),'(a,I6,a,I6,a)') 'For particles of type ', itype, ',  MB state ', &
               mm, ' is totally antisymmetric'
        call write_info(1)
        modelMBparticles%bosonfermion(itype) = 'fermion'
      end if
    end do ! itype

    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(ixp)

    call pop_sub()
  end subroutine modelmb_find_exchange_syms
  ! ---------------------------------------------------------


  ! project out states with proper symmetry for cases which are of symmetry = unknown
  subroutine modelmb_sym_state(dir, gr, mm, geo, modelMBparticles, wf, symmetries_satisfied)

    implicit none

    character(len=*), intent(in)    :: dir
    type(grid_t),     intent(in)    :: gr
    integer,          intent(in)    :: mm
    type(geometry_t), intent(in)    :: geo
    type(modelMB_particle_t), intent(in) :: modelMBparticles
    ! this will be antisymmetrized on output
    CMPLX,                    intent(inout) :: wf(1:gr%mesh%np_part_global)
    logical,          intent(out)   :: symmetries_satisfied

    ! local vars
    integer :: itype, ipart1, ipart2, npptype
    integer :: ip, ipp, iup, idown, iyoung
    integer :: iantisym, ierr, ikeeppart
    integer :: ndimMB
    integer :: nspindown, nspinup, iperm_up, iperm_down
    integer :: iunit, jj, idir

    type(permutations_t) :: perms_up, perms_down
    type(modelmb_1part_t) :: mb_1part
    type(young_t) :: young

    integer, allocatable :: ix(:), ixp(:), ofst(:)
    integer, allocatable :: symmetries_satisfied_alltypes(:)
 
    FLOAT :: normalizer, norm
    CMPLX :: scalprod

    FLOAT, allocatable  :: antisymrho(:)
    FLOAT, allocatable  :: antisymrho_1part(:)
    CMPLX, allocatable  :: antisymwf(:)
    CMPLX, allocatable  :: antisymwf_swap(:)

    character (len=80) :: tmpstring
    character (len=80) :: filename

    logical :: debug_antisym


! source

    call push_sub('states.modelmb_sym_states')

    symmetries_satisfied = .false.

    debug_antisym = .false.

    call permutations_nullify(perms_up)
    call permutations_nullify(perms_down)
    call young_nullify (young)

    call modelmb_1part_nullify(mb_1part)

    SAFE_ALLOCATE(ix(1:MAX_DIM))
    SAFE_ALLOCATE(ixp(1:MAX_DIM))

    SAFE_ALLOCATE(symmetries_satisfied_alltypes(1:modelMBparticles%ntype_of_particle_modelMB))
    symmetries_satisfied_alltypes = 0

    ndimMB=modelMBparticles%ndim_modelMB

    normalizer = product(gr%mesh%h(1:gr%mesh%sb%dim)) !M_ONE/units_out%length%factor**gr%mesh%sb%dim

    if (modelMBparticles%ntype_of_particle_modelMB > 1) then
      write (message(1), '(a)') 'modelmb_sym_state: not coded for several particly types '
      call write_fatal(1)
    end if

    ! FIXME: could one of these arrays be avoided?
    SAFE_ALLOCATE(antisymwf(1:gr%mesh%np_part_global))
    SAFE_ALLOCATE(antisymwf_swap(1:gr%mesh%np_part_global))

    ! for each particle type
    do itype = 1, modelMBparticles%ntype_of_particle_modelMB

      ! FIXME: for multiple particle types this needs to be fixed.
      ! Also, for inequivalent spin configurations this should vary, and we get
      ! different 1 body densities, no?
      ikeeppart = 1

      ! if the particle is not fermionic, just cycle to next one
      ! FIXME: boson case is not treated yet
      if (modelMBparticles%bosonfermion(ikeeppart) /= 'fermion') then
        symmetries_satisfied_alltypes(itype) = 1
        cycle
      end if

      call modelmb_1part_init(mb_1part, gr%mesh, ikeeppart, &
             modelMBparticles%ndim_modelmb, gr%sb%box_offset)

      npptype = modelMBparticles%nparticles_per_type(itype)

      SAFE_ALLOCATE(antisymrho(1:gr%mesh%np_part_global))
      SAFE_ALLOCATE(ofst(1:npptype))
      do ipart1 = 1, npptype
        ofst(ipart1) = (ipart1-1)*ndimMB
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
        do iyoung = 1, young%nyoung
          antisymwf=cmplx(0.0d0,0.0d0)
          antisymwf(:) = wf(:)

          ! first symmetrize over pairs of particles associated in the present
          ! Young diagram
          do idown = 1, young%ndown
            ipart1 = young%young_down(idown,iyoung)
            ipart2 = young%young_up  (idown,iyoung)
  
            antisymwf_swap=antisymwf
            do ip = 1, gr%mesh%np_part_global
              ! get present position
              call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)
       
              ! invert coordinates of ipart1 and ipart2
              ixp = ix
              ! permutate the particles ipart1 and its spin up partner
              ixp (ofst(ipart1)+1:ofst(ipart1)+ndimMB) = &
                  ix (ofst(ipart2)+1:ofst(ipart2)+ndimMB)
              ixp (ofst(ipart2)+1:ofst(ipart2)+ndimMB) = &
                  ix (ofst(ipart1)+1:ofst(ipart1)+ndimMB)
              
              ! get position of exchanged point
              ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)
  
              antisymwf_swap(ip)=antisymwf_swap(ip) + antisymwf(ipp)
          
            end do ! ip
            antisymwf = antisymwf_swap * 0.5d0
          end do
    
          if (debug_antisym) then 
            antisymrho = real(conjg(antisymwf)*antisymwf) * normalizer
            norm = sum(antisymrho)
            write (message(1), '(a,I7,a,I7,a,E20.10)') 'norm of pair-symmetrized-state ',&
                    mm, ' with ', nspindown, ' spins down is ', norm
            call write_info(1)
          end if
 
          ! for each permutation of particles of this type
          !  antisymmetrize the up and down labeled spins, amongst themselves
          antisymwf_swap = cmplx(0.0d0, 0.0d0)
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
                ixp (ofst(ipart1)+1:ofst(ipart1)+ndimMB) = ix (ofst(ipart2)+1:ofst(ipart2)+ndimMB) ! part1 to 2
              end do
              ! get position of exchanged point
              ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)
              antisymwf_swap(ip)=antisymwf_swap(ip) + perms_up%permsign(iperm_up)*antisymwf(ipp)
            end do ! ip
          end do ! iperm_up

          antisymwf=antisymwf_swap / dble(perms_up%npermutations)

          ! the following could be removed for production
          if (debug_antisym) then 
            antisymrho = real(conjg(antisymwf)*antisymwf) * normalizer
            norm = sum(antisymrho)
            write (message(1), '(a,I7,a,I7,a,E20.10)') 'norm of up-antisym+pairsym-state ',&
                    mm, ' with ', nspindown, ' spins down is ', norm
            call write_info(1)
          end if

          antisymwf_swap=cmplx(0.0d0, 0.0d0)
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
                ixp (ofst(ipart1)+1:ofst(ipart1)+ndimMB) = ix (ofst(ipart2)+1:ofst(ipart2)+ndimMB) ! part1 to 2
              end do
              ! get position of exchanged point
              ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)
              antisymwf_swap(ip)=antisymwf_swap(ip) + perms_down%permsign(iperm_down)*antisymwf(ipp)
            end do ! ip
          end do ! iperm_down

          antisymwf=antisymwf_swap / dble(perms_down%npermutations)
     
          antisymrho = real(conjg(antisymwf)*antisymwf)
          norm = sum(antisymrho) * normalizer

          if (norm < 1.d-8) cycle

          call young_write_one (young, iyoung)
          
          antisymwf=antisymwf_swap / sqrt(norm)
  
          scalprod = sum(conjg(antisymwf)*wf)
          write (message(1), '(a,I7,a,I7,a,2E20.10)') 'scalar product to wf of state ', mm, &
                  ' having ', nspindown, ' spins down is ', scalprod
          call write_info(1)

          write (message(1), '(a,I7,a,I7,a,E20.10)') 'norm of state ', mm, ' with ', nspindown, ' spins down is ', norm
          call write_info(1)
   
          SAFE_ALLOCATE(antisymrho_1part(1:mb_1part%npt_1part))
          call zmf_calculate_rho(ikeeppart, mb_1part, npptype, gr%mesh, antisymwf, antisymrho_1part)
          ! FIXME: this will also depend on particle type and idimension
          write(filename,'(a,a,i3.3,a,i3.3,a,i3.3,a,i3.3)') trim(dir),'/asymden_iMB', mm, &
                '_ipar', ikeeppart,'_ndn',nspindown,'_iY',iyoung
          iunit = io_open(filename,action='write')
          do jj = 1, mb_1part%npt_1part
            call hypercube_i_to_x(mb_1part%hypercube_1part, mb_1part%ndim1part, &
                 mb_1part%nr_1part, mb_1part%enlarge_1part(1), jj, ix)
            do idir=1,mb_1part%ndim1part
              write(iunit,'(es11.3)', ADVANCE='no') ix(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
            end do
            write(iunit,'(es18.10)') antisymrho_1part(jj)
          end do
          call io_close(iunit)
          SAFE_DEALLOCATE_A(antisymrho_1part)

          ! we have found a valid Young diagram and antisymmetrized the present
          ! state enough. Loop to next type
          symmetries_satisfied_alltypes(itype) = 1
          exit

        end do ! iyoung loop
   
        call permutations_end(perms_up)
        call permutations_end(perms_down)
        call young_end (young)
 
        if (symmetries_satisfied_alltypes(itype) == 1) exit

      end do ! nspindown=0,1...
   
      SAFE_DEALLOCATE_A(antisymrho)

      call modelmb_1part_end(mb_1part)

    end do ! itype

    ! check if all types of particles have been properly symmetrized
    if (sum(symmetries_satisfied_alltypes) == modelMBparticles%ntype_of_particle_modelMB) then
      wf = antisymwf
      symmetries_satisfied = .true.
    end if

    SAFE_DEALLOCATE_A(antisymwf_swap)
    SAFE_DEALLOCATE_A(antisymwf)

    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(ixp)

    call pop_sub()
  end subroutine modelmb_sym_state
  ! ---------------------------------------------------------


end module modelmb_exchange_syms_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
