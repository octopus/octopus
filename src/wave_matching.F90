!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id$

#include "global.h"

module wave_matching
use global
use system
use hamiltonian
use io
use states
use functions
use derivatives
use output
use math
use scf

implicit none


private
public :: wave_matching_run, op

type(der_discr_type), pointer :: der_pointer
type(f_der_type),     pointer :: f_der
FLOAT  :: energy  
  
contains

subroutine op(x, y)
  CMPLX, intent(in)  :: x(:)
  CMPLX, intent(out) :: y(:)
  CMPLX, allocatable :: u(:)
  integer :: n
  integer :: i

  allocate(u(der_pointer%m%np))
  do i = 1, der_pointer%m%np
     u(i) = -energy*x(i)
  enddo
!  call zderivatives_lapl(der_pointer, x, u)
  call zf_laplacian(f_der , x, y)
  do i = 1, der_pointer%m%np
     y(i) = -M_HALF*y(i) + u(i)
  enddo
  deallocate(u)
end subroutine op
  

integer function wave_matching_run(sys, h) result(ierr)
  type(system_type)       :: sys
  type(hamiltonian_type)  :: h

  type(f_der_type),     target :: f_der_current
  type(der_discr_type), target :: der_current
  type(der_discr_type)         :: der_backup  
  type(scf_type)               :: scfv
  CMPLX, allocatable   :: boundary_points(:,:,:,:), b(:,:,:), xout(:,:,:)
  CMPLX                :: zres
  FLOAT, allocatable   :: weights(:,:), occs(:,:), energies(:,:),kpoints(:,:,:)
  FLOAT   :: threshold, match_energy_delta
  integer :: align_axis, boundary_np, nstates
  integer :: ik1, ist1, idim1, ik2, ist2, idim2, j1, j2
  integer :: i, j, s1, s2, iter


  call push_sub('wave_matching_run')

  ierr = 0
  ! consistency checks for the input file we are currently running
  call check_params

  ! write out potentials
  call hamiltonian_output(h, sys%m, trim(current_label)//"static", sys%outp)

  ! determine the number of boundary points we have to consider for the laplacian
  call calc_no_boundary_points
  
  ! since we checked in check_params that all susbystems have the same number of states 
  ! this could also be minval, or simply one of them.
  nstates = maxval(no_of_states)
  allocate(boundary_points(2,2,nstates,boundary_np))
  allocate(occs(2,nstates), energies(2,nstates))
  allocate(kpoints(2,3,sys%st%d%nik))
  allocate(b(nstates,nstates,sys%m%np), xout(nstates,nstates,sys%m%np))
  allocate(weights(sys%f_der%der_discr%op(conf%dim+1)%n,sys%m%np))
  occs = M_ZERO
  energies = M_ZERO

  ! now we read the boundaries of our left and right neighbors respectively
  ! boundaries of left side
  call read_boundary_points(subsys_label(current_subsystem-1), &
       boundary_points(1,:,:,:), occs(1,:), energies(1,:), kpoints(1,:,:))
  ! boundaries of right side
  call read_boundary_points(subsys_label(current_subsystem+1), &
       boundary_points(2,:,:,:), occs(2,:), energies(2,:), kpoints(2,:,:))

  ! keep a backup of the discretization
  allocate(der_backup%op(conf%dim + 1),der_current%op(conf%dim + 1))
  der_backup%grad  => der_backup%op
  der_backup%lapl  => der_backup%op(conf%dim + 1)
  der_current%grad => der_current%op
  der_current%lapl => der_current%op(conf%dim + 1)

  ! some initialization
  weights    = sys%f_der%der_discr%lapl%w
  der_backup = sys%f_der%der_discr
  b          = M_ZERO
  
  call modify_laplacian

  call loct_parse_float(check_inp('SolverThreshold'),  CNST(0.05), threshold)
  call loct_parse_float(check_inp('MatchEnergyDelta'), CNST(1e-3), match_energy_delta)


  iter = 400
  threshold = CNST(1e-9)

  call scf_init(scfv, sys%m, sys%st, sys%geo, h)

  ! SCF cycle
  scf: do iter = 1, 1 !scfv%max_iter

  ! solve linear system for central region
     call linsys_solver_run

    ! occupations
!   this should be replaced by a new function. the occupations are fixed from the
!   leads outside. a "global" state is either occupied or not. there are no states
!   occupied only in the leads and not in the center.
!   what about occupations of localized states?
!    call states_fermi(st, m)

!     call X(states_calc_dens)(st, m%np, st%rho)
     ! compute convergence criteria

    ! are we finished?
!    finish = &
!        (scf%conv_abs_dens > M_ZERO .and. scf%abs_dens <= scf%conv_abs_dens) .or. &
!        (scf%conv_rel_dens > M_ZERO .and. scf%rel_dens <= scf%conv_rel_dens) .or. &
!        (scf%conv_abs_ev   > M_ZERO .and. scf%abs_ev   <= scf%conv_abs_ev) .or. &
!        (scf%conv_rel_ev   > M_ZERO .and. scf%rel_ev   <= scf%conv_rel_ev)

!    call scf_write_iter

    ! mixing
!    select case (scf%what2mix)
!    case (MIXDENS)
!       ....
!    end select

    ! save restart information


  enddo scf

  call scf_end(scfv)

  deallocate(boundary_points, occs, energies)
  deallocate(kpoints, b, xout, weights)

  call pop_sub()
  
  
contains

  !------------------------------------------------------------------------------
  subroutine check_params()
    implicit none

    integer(POINTER_SIZE) :: blk
    character(256)        :: parse_label

    call push_sub('check_params')    

    if( current_subsystem-1 .lt. 1 ) then
       message(1) = 'Error: Missing left neighbor.'
       message(2) = 'Please correct your input file.'
       call write_fatal(2)
    endif
    
    if( current_subsystem+1 .gt. no_syslabels ) then
       message(1) = 'Error: Missing right neighbor.'
       message(2) = 'Please correct your input file.'
       call write_fatal(2)
    endif
    
    if( subsys_runmode(current_subsystem-1) .eq. 10 ) then
       message(1) = 'Error: Left neighbor cannot be in runmode wave_matching.'
       message(2) = 'Please correct your input file.'
       call write_fatal(2)
    endif
    
    if( subsys_runmode(current_subsystem+1) .eq. 10 ) then
       message(1) = 'Error: Right neighbor cannot be in runmode wave_matching.'
       message(2) = 'Please correct your input file.'
       call write_fatal(2)
    endif

    if(no_of_states(current_subsystem).ne.no_of_states(current_subsystem-1).or. &
         no_of_states(current_subsystem).ne.no_of_states(current_subsystem+1) ) then
       message(1) = 'Error: Mismatch of number of states for the different subsystems.'
       message(2) = 'Please correct your input file.'
       call write_fatal(2)     
    endif

    ! check if we have the same number of k-points for both neighbors
    ! left side
    parse_label = trim(subsys_label(current_subsystem-1))// trim('NumberKPoints')
    if(loct_parse_block(trim(parse_label), blk) .eq. 0) then       
       message(1) = 'Error: Only a single block "NumberKPoints" supported for all subsystems.'
       message(2) = 'If the Numbers of KPoints would differ we could not match all states.'
       message(3) = 'Please correct your input file.'
       call write_fatal(3)
    end if

    ! right side
    parse_label = trim(subsys_label(current_subsystem+1))// trim('NumberKPoints')
    if(loct_parse_block(trim(parse_label), blk) .eq. 0) then
       message(1) = 'Error: Only a single block "NumberKPoints" supported for all subsystems.'
       message(2) = 'If the Numbers of KPoints would differ we could not match all states.'
       message(3) = 'Please correct your input file.'
       call write_fatal(3)
    end if
    
    message(1) = 'Info: Starting Wave-Matching'
    message(2) = 'Info: We are                    : '//subsys_label(current_subsystem)
    message(3) = 'Info: Our left neighbor is      : '//subsys_label(current_subsystem-1)
    message(4) = 'Info: And our right neighbor is : '//subsys_label(current_subsystem+1)
    call write_info(4, stress=.true.)
    
    call pop_sub()
  end subroutine check_params

  !------------------------------------------------------------------------------
  subroutine calc_no_boundary_points()
    implicit none
    
    character(16)        :: points

    call push_sub('calc_no_boundary_points')    

    call loct_parse_int('AlignmentAxis', 1, align_axis)
    
    boundary_np = sys%f_der%der_discr%order
    do i = 1, conf%dim
       !if(i.ne.align_axis) 
       boundary_np = boundary_np*(sys%m%nr(2,i)-sys%m%nr(1,i) + 1)
    enddo

    write(points,'(i)') boundary_np
    message(1) = 'Info: Number of boundary points: '//points
    if(conf%verbose > 30) call write_info(1)
    
    call pop_sub()
  end subroutine calc_no_boundary_points

  !------------------------------------------------------------------------------
  subroutine read_boundary_points(sys_label, boundary_p, occs, energies, kpoints)
    implicit none

    character(len=*), intent(in)  :: sys_label
    CMPLX,            intent(out) :: boundary_p(:,:,:)
    FLOAT,            intent(out) :: occs(:), energies(:), kpoints(:,:)

    character(len=10)   :: filenum
    character(len=128)  :: filename, tmp
    CMPLX, allocatable  :: f(:)
    integer             :: iunit1, iunit2, j, ik, ist, idim, k
    integer             :: function_kind, file_kind, file_np

    call push_sub('read_boundary_points')

    call io_assign(iunit1)
    filename = trim(sys_label)//trim('tmp/restart_gs/occs')
    iunit1 = io_open(filename, action='read', status='old')
    ! read two lines
    read(iunit1, *) tmp
    read(iunit1, *) tmp

    j = 1
    do ik = 1, sys%st%d%nik
       do ist = 1, sys%st%nst
          do idim = 1, sys%st%d%dim

             write(filenum,'(i10.10)') j
             filename = trim(sys_label)//trim('tmp/restart_gs/')//trim(filenum)
             call io_assign(iunit2)
             iunit2 = io_open(filename, action='read', status='old', form='unformatted', die=.false.)
             read(unit = iunit2, iostat = ierr) file_kind, file_np
             message(1) = 'Info: Reading boundaries from '//filename
             if(conf%verbose > 30) call write_info(1)             
             
             allocate(f(file_np))
             
             if (ierr == 0) then
                read(unit = iunit2, iostat = ierr) f(1:file_np)
             endif
             
             do k = 1, boundary_np
                boundary_p(1,j,k) = f(k)
                boundary_p(2,j,k) = f(file_np - boundary_np + k)
             enddo
             
             deallocate(f)
             
             call io_close(iunit2)       
             
             read(iunit1,*) occs(j),tmp,energies(j),tmp, &
                  kpoints(1,ik),tmp,kpoints(2,ik),tmp,kpoints(3,ik)
!             write(*,*) occs(j),tmp,energies(j),tmp, &
!                  kpoints(1,ik),tmp,kpoints(2,ik),tmp,kpoints(3,ik)
!             write(*,*) '----------------'

             j = j + 1
          enddo
       enddo
    enddo
    call io_close(iunit1)              
    
    call pop_sub()
  end subroutine read_boundary_points
  
  !------------------------------------------------------------------------------
  subroutine modify_laplacian()
    implicit none

    integer :: ix(3), ixl(3), k, l, qmax, s1, s2

    call push_sub('modify_laplacian')


    ! TODO: multiply boundary_points with plane wave of correspondig kpoint. should the phase shifts of the plane waves
    !       be determined by the box size of the central system?
!    vs1: do s1 = 1, nstates
!       vs2: do s2 = 1, nstates

    j1 = 1
    do ik1 = 1, sys%st%d%nik
       do ist1 = 1, sys%st%nst
          do idim1 = 1, sys%st%d%dim
             j2 = 1
             do ik2 = 1, sys%st%d%nik
                do ist2 = 1, sys%st%nst
                   do idim2 = 1, sys%st%d%dim
          
          sys%f_der%der_discr%lapl%w = weights
          
          vk: do k = 1, sys%m%np ! for all points in the mesh
             
             ix(:) = sys%m%Lxyz(k,:)
             
             vj: do j = 1, sys%f_der%der_discr%op(conf%dim+1)%n   ! for all points in stencil of the laplacian
                l = align_axis
                ixl = ix
                if(ix(l)+sys%f_der%der_discr%op(conf%dim+1)%stencil(l,j) < sys%m%nr(1, l)) then
                   qmax = sys%m%nr(1, l) - ( ix(l)+sys%f_der%der_discr%op(conf%dim+1)%stencil(l,j) )
                   ixl(l) = sys%m%Lxyz(1,l) + sys%f_der%der_discr%order - qmax
                   ! write(*,*) 'L: outside by', k, j, sys%f_der%der_discr%order - qmax, sys%m%Lxyz_inv(ixl(1),ixl(2),ixl(3))
                   ! build up right hand side
                   b(j1,j2,k) = b(j1,j2,k) + M_HALF*boundary_points(1,2,j1,sys%m%Lxyz_inv(ixl(1),ixl(2),ixl(3))) &
                        *sys%f_der%der_discr%lapl%w(j,k)
                   ! write(*,*) 'b',j,k,b(j1,j2,k),sys%f_der%der_discr%lapl%w(j,k)
                   ! set coeff to zero
                   sys%f_der%der_discr%lapl%w(j,k) = M_ZERO
                else if (ix(l)+sys%f_der%der_discr%op(conf%dim+1)%stencil(l,j) > sys%m%nr(2, l)) then
                   
                   qmax = ( ix(l)+sys%f_der%der_discr%op(conf%dim+1)%stencil(l,j) ) - sys%m%nr(2, l)
                   ixl(l) = sys%m%Lxyz(1,l) + qmax - 1
                   !                 write(*,*) 'R: outside by', k, j, qmax, sys%m%Lxyz_inv(ixl(1),ixl(2),ixl(3))
                   ! build up right hand side
                   b(j1,j2,k) = b(j1,j2,k) + M_HALF*boundary_points(2,1,j2,sys%m%Lxyz_inv(ixl(1),ixl(2),ixl(3))) &
                        *sys%f_der%der_discr%lapl%w(j,k)
                   sys%f_der%der_discr%lapl%w(j,k) = M_ZERO    ! set lapl coeff to zero
                endif

             enddo vj
          enddo vk

                      j2 = j2 + 1 
                   enddo
                enddo
             enddo
             j1 = j1 + 1
          enddo
       enddo
    enddo


!       enddo vs2
!    enddo vs1

    call pop_sub()
  end subroutine modify_laplacian


  !------------------------------------------------------------------------------
  subroutine linsys_solver_run()
    implicit none

    character(len=128)   :: filenbase

    call push_sub('linsys_solver_run')
    j1 = 1
    do ik1 = 1, sys%st%d%nik
       do ist1 = 1, sys%st%nst
          do idim1 = 1, sys%st%d%dim
             j2 = 1
             do ik2 = 1, sys%st%d%nik
                do ist2 = 1, sys%st%nst
                   do idim2 = 1, sys%st%d%dim
                      
                      iter = 400
                      threshold = CNST(1e-5)
                      ! will be variable later (right now it doesnt change)
                      der_current = sys%f_der%der_discr
                      f_der_current = sys%f_der
                      der_pointer => der_current
                      f_der => f_der_current
                      
!                      write(*,*) '+_+_+_+_+_+_+_+'
!                      write(*,*) 'o. states',nstates,j1,j2,occs(1,j1),occs(2,j2)
!                      write(*,*) '+_+_+_+_+_+_+_+'                      

!                      write(*,*) ik1,ik2,ist1,ist2,idim1,idim2
                      ! Match only states with the same energy (up to MatchEnergyDelta)
                      if (abs(energies(1,j1)-energies(2,j2)).lt.match_energy_delta) then
!                         write(*,*) 'e. states',j1,j2,energies(1,j1),energies(2,j2)


!                         write(*,'(10F12.5)') kpoints(1,1,ik1),kpoints(1,2,ik1),kpoints(1,3,ik1), &
!                         kpoints(2,1,ik2),kpoints(2,2,ik2),kpoints(2,3,ik2)
                         !        if(s1.eq.s2) then
                         !           write(*,*) energies(1,s1),energies(2,s2)
                         !        endif
                         
                         !                    energy = energies(1,s1)
                         ! TODO: match only states with the same energy 
                         ! probably something like |e1-e2|<crit, alternatively states with the same index.
                         call zconjugate_gradients(sys%m%np, xout(j1,j2,:), b(j1,j2,:), op, &
                              iter, zres, threshold)

                         
                         write(filenbase,'(a,i3.3,a,i3.3)') 'wm-', j1, '-', j2
                         !                       call zoutput_function(output_fill_how("MeshIndex_and_Gnuplot"), & 

                         call zoutput_function(output_fill_how("MeshIndex"), & 
                              trim(current_label)//"static", trim(filenbase), sys%m, b(j1,j2,:), M_ONE, ierr)

                      endif
                      
                      j2 = j2 + 1 
                   enddo
                enddo
             enddo
             j1 = j1 + 1
          enddo
       enddo
    enddo
    
    call pop_sub()
  end subroutine linsys_solver_run

end function wave_matching_run
  

end module wave_matching
