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
  use modelmb_particles_m
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
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use profiling_m
  use states_m
  use units_m

  implicit none

  private

  public :: modelmb_find_exchange_syms,&
            modelmb_find_exchange_syms_all

contains

  ! ---------------------------------------------------------
  subroutine modelmb_find_exchange_syms_all(dir, gr, st, geo)
    character(len=*), intent(in) :: dir
    type(grid_t),     intent(in) :: gr
    type(geometry_t),       intent(in)    :: geo
    type(states_t),   intent(inout) :: st

    integer :: mm

    call push_sub('states.modelmb_find_exchange_syms_all')

    do mm=1,st%nst
      call modelmb_find_exchange_syms(gr, mm, st)
! FIXME: separate output from calculation of antisymmetrized states
      call modelmb_sym_states(dir, gr, mm, st, geo)
    end do

    call pop_sub()

  end subroutine modelmb_find_exchange_syms_all



  ! ---------------------------------------------------------
  subroutine modelmb_find_exchange_syms(gr, mm, st)
    type(grid_t),     intent(in) :: gr
    integer,          intent(in) :: mm
    type(states_t),   intent(inout) :: st

    ! local vars
    integer :: itype, ipart1, ipart2, npptype
    integer :: ip, ipp, maxtries
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

    ndimMB=st%modelMBparticles%ndim_modelMB

    ! this is the tolerance on the wavefunction in order to use it to determine
    ! sign changes. Should be checked, and might be a function of grid number,
    ! as the norm is diluted over the whole cell...
    ! maybe min(0.001,100*1./ngridpoints) would be better.
    wftol=1.d-5

    ! tolerance on wf(1,2)/wf(2,1) being close to +-1
    tolonsign = 1.d-4

    ! maximum number of times we try to extract the parity
    maxtries = 100

    ! initialize the container
    st%modelMBparticles%exchange_symmetry=0

    ! for each particle type
    do itype=1,st%modelMBparticles%ntype_of_particle_modelMB

      npptype=st%modelMBparticles%nparticles_per_type(itype)
      nratiocomplaints = 0

      ! for each pair of particles of this type
      do ipart1=1,npptype
        do ipart2=ipart1+1,npptype
 
          ! try a bunch of points in the wf
          ngoodpoints=0
          ntries=0
          save_sign_change=0
          do while (ntries < maxtries .and. ngoodpoints < 10)
            ntries = ntries+1
            
            ! find coordinates of random point
            call random_number(rand0to1)
            ip = nint(gr%mesh%np_global*rand0to1)
 
            ! wavefunction value at this point ip
            if(states_are_real(st)) then
              wfval1 = cmplx(st%dpsi(ip, 1, mm, 1),M_ZERO)
            else
              wfval1 =       st%zpsi(ip, 1, mm, 1)
            end if

            ! see if point has enough weight to be used
            if (abs(real(wfval1)) > wftol) then

              ! get actual coordinates of MB point
              call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)
 
              ! invert coordinates of ipart1 and ipart2
              ixp = ix
              ofst1=(ipart1-1)*ndimMB
              ofst2=(ipart2-1)*ndimMB
              ixp (ofst1+1:ofst1+ndimMB) = ix (ofst2+1:ofst2+ndimMB) ! part1 to 2
              ixp (ofst2+1:ofst2+ndimMB) = ix (ofst1+1:ofst1+ndimMB) ! part2 to 1
 
              ! recover MB index for new point
              ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)

              ! wavefunction value at exchanged coordinate point ipp
              if(states_are_real(st)) then
                wfval2 = cmplx(st%dpsi(ipp, 1, mm, 1),M_ZERO)
              else
                wfval2 =       st%zpsi(ipp, 1, mm, 1)
              end if
 
              ! check sign change FIXME: check imaginary part too
              !   normally, as these points are related, wfval2 should also be > wftol
              ratio = real(wfval1) / real(wfval2)
              
              if (abs(ratio+1.0d0) < tolonsign) then
                wf_sign_change = -1
              else if (abs(ratio-1.0d0) < tolonsign) then
                wf_sign_change =  1
              else  ! ratio not clear enough to get parity
                if (nratiocomplaints < 10) then
                  write (message(1),'(a,E20.10,a)') 'point skipped for exchange parity: ratio of wf ',ratio, &
                       ' is not clearly +-1'
                  call write_info(1)
                  nratiocomplaints = nratiocomplaints+1
                end if
                cycle
              end if
 
              ! if opposite of previous value, complain
              if (save_sign_change /= 0 .and. save_sign_change /= wf_sign_change) then
                write (message(1),'(a)') ' WARNING: sign change is different between wf points examined'
                call write_info(1)
              end if

              save_sign_change=wf_sign_change
              ngoodpoints = ngoodpoints + 1
            end if
          end do !  tries and goodpoints
         
          if (ngoodpoints < 10) then
            !write (message(1),'(a)') ' WARNING: unable to determine sign change for'
            !write (message(2),'(a,I6,a,I6,a,I6)') ' exchange of ', ipart1, ' and ', ipart2, ' of many-body state ', mm
            !call write_info(2)
          else 
            st%modelMBparticles%exchange_symmetry(ipart1,ipart2,itype)=save_sign_change
            st%modelMBparticles%exchange_symmetry(ipart2,ipart1,itype)=save_sign_change
          end if 
 
        end do ! ipart2
      end do ! ipart1
    end do ! itype

    ! do some crude output
    do itype=1,st%modelMBparticles%ntype_of_particle_modelMB
      !write (message(1),'(a)')    ' lower triangle of exchange parity matrix'
      write (message(1),'(a,I6,a,I6)') '   for particles of type ', itype, '   in MB state ', mm
      call write_info(1)

      do ipart1=1,npptype
        do ipart2=1,ipart1
          write (*,'(2x, I6)',ADVANCE='NO') st%modelMBparticles%exchange_symmetry(ipart2,ipart1,itype)
        end do ! ipart2
        write (*,*)
      end do ! ipart1

      st%modelMBparticles%bosonfermion(itype) = 'unknown'
      if (sum(st%modelMBparticles%exchange_symmetry) ==  npptype*(npptype-1)) then
        write (message(1),'(a,I6,a,I6,a)') 'For particles of type ', itype, ',  MB state ', mm, ' is totally symmetric (bosonic)'
        call write_info(1)
        st%modelMBparticles%bosonfermion(itype) = 'boson'
      else if (sum(st%modelMBparticles%exchange_symmetry) == -npptype*(npptype-1)) then
        write (message(1),'(a,I6,a,I6,a)') 'For particles of type ', itype, ',  MB state ', mm, ' is totally antisymmetric'
        call write_info(1)
        st%modelMBparticles%bosonfermion(itype) = 'fermion'
      end if
    end do ! itype

    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(ixp)

    call pop_sub()
  end subroutine modelmb_find_exchange_syms
  ! ---------------------------------------------------------


  ! project out states with proper symmetry for cases which are of symmetry = unknown
  subroutine modelmb_sym_states(dir, gr, mm, st, geo)
    character(len=*), intent(in) :: dir
    type(grid_t),     intent(in)    :: gr
    integer,          intent(in)    :: mm
    type(states_t),   intent(inout) :: st
    type(geometry_t), intent(in)    :: geo

    ! local vars
    integer :: itype, ipart1, ipart2, npptype
    integer :: ip, ipp, maxtries
    integer :: iantisym, ierr
    integer :: ofst1, ofst2, ndimMB

    integer, allocatable :: ix(:), ixp(:)
 
    FLOAT :: normalizer

    CMPLX :: wfval1, wfval2
    CMPLX, allocatable  :: antisymwf(:,:)

    character(len=80) :: fname, tmpstring

! source

    call push_sub('states.modelmb_sym_states')

    SAFE_ALLOCATE(ix(1:MAX_DIM))
    SAFE_ALLOCATE(ixp(1:MAX_DIM))

    ndimMB=st%modelMBparticles%ndim_modelMB

    normalizer = M_ONE/units_out%length%factor**gr%mesh%sb%dim

    ! for each particle type
    do itype = 1,st%modelMBparticles%ntype_of_particle_modelMB
      if (st%modelMBparticles%bosonfermion(itype) /= 'unknown') cycle

      npptype=st%modelMBparticles%nparticles_per_type(itype)
      SAFE_ALLOCATE(antisymwf(1:gr%mesh%np_global,npptype*(npptype-1)/2))

      iantisym = 0
      ! for each pair of particles of this type
      do ipart1 = 1, npptype
        do ipart2 = 1, ipart1-1
          iantisym = iantisym+1
 
          do ip = 1, gr%mesh%np_global
            ! get present position
            call index_to_coords(gr%mesh%idx, gr%sb%dim, ip, ix)

            ! invert coordinates of ipart1 and ipart2
            ixp = ix
            ofst1=(ipart1-1)*ndimMB
            ofst2=(ipart2-1)*ndimMB
            ixp (ofst1+1:ofst1+ndimMB) = ix (ofst2+1:ofst2+ndimMB) ! part1 to 2
            ixp (ofst2+1:ofst2+ndimMB) = ix (ofst1+1:ofst1+ndimMB) ! part2 to 1

            ! get position of exchanged point
            ipp = index_from_coords(gr%mesh%idx, gr%sb%dim, ixp)

            if(states_are_real(st)) then
              wfval1 = cmplx(st%dpsi(ip,  1, mm, 1),M_ZERO)
              wfval2 = cmplx(st%dpsi(ipp, 1, mm, 1),M_ZERO)
            else
              wfval1 =       st%zpsi(ip,  1, mm, 1)
              wfval2 =       st%zpsi(ipp, 1, mm, 1)
            end if

            antisymwf(ip,iantisym)=wfval1-wfval2
      
          end do ! ip

          write(fname,'(a,i4.4,a,i4.4,a,i4.4,a,i4.4)') '/antisymwf.iMB',mm,'_typ',&
             itype,'_asym',ipart1,'_',ipart2

          call zoutput_function(output_xcrysden, dir, fname, gr%mesh, gr%sb, &
                   antisymwf(:,iantisym), sqrt(normalizer), ierr, is_tmp = .false., geo = geo)

        end do ! ipart2
      end do ! ipart1

      ! normalize the bloody thing
      antisymwf=antisymwf*sqrt(0.5d0)

      SAFE_DEALLOCATE_A(antisymwf)

    end do ! itype

    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(ixp)

    call pop_sub()
  end subroutine modelmb_sym_states
  ! ---------------------------------------------------------


end module modelmb_exchange_syms_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
