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

implicit none


private
public :: wave_matching_run
  
  
contains
  
integer function wave_matching_run(sys, h) result(ierr)
  type(system_type)       :: sys
  type(hamiltonian_type)  :: h
  
  call push_sub('wave_matching_run')
  
  ierr = 0
  
  call check_params()


  
  call pop_sub()
  
  
contains
  
  subroutine check_params()
    
    if( current_subsystem-1 .lt. 1 ) then
       message(1) = 'Error: Missing left neighbor'
       message(2) = 'Please correct your input file'
       call write_fatal(2)
    endif
    
    if( current_subsystem+1 .gt. no_syslabels ) then
       message(1) = 'Error: Missing right neighbor'
       message(2) = 'Please correct your input file'
       call write_fatal(2)
    endif
    
    if( subsys_runmode(current_subsystem-1) .eq. 10 ) then
       message(1) = 'Error: Left neighbor cannot be in runmode wave_matching'
       message(2) = 'Please correct your input file'
       call write_fatal(2)
    endif
    
    if( subsys_runmode(current_subsystem+1) .eq. 10 ) then
       message(1) = 'Error: Right neighbor cannot be in runmode wave_matching'
       message(2) = 'Please correct your input file'
       call write_fatal(2)
    endif
    
    message(1) = 'Info: Starting Wave-Matching'
    message(2) = 'Info: We are                    : '//subsys_label(current_subsystem)
    message(3) = 'Info: Our left neighbor is      : '//subsys_label(current_subsystem-1)
    message(4) = 'Info: And our right neighbor is : '//subsys_label(current_subsystem+1)
    call write_info(4, stress=.true.)
    
  end subroutine check_params
  
  
end function wave_matching_run
  

end module wave_matching

























































!!$    !  now we check what we have to run in the respective subsystem
!!$    if(loct_parse_block("SystemRunModes", blk) == 0) then
!!$       num_subsys_runmodes = loct_parse_block_cols(blk,0)
!!$    else
!!$       message(1) = "Could not find block SystemRunModes in the input file."
!!$       message(2) = "This block is mandatory for run mode MultiSubsystem."
!!$       call write_fatal(2)
!!$    endif
!!$    
!!$    if(num_subsys_runmodes/=num_syslabels) then
!!$       message(1) = "The blocks SystemLabels and SystemRunModes do not have"
!!$       message(2) = "the same size. Please correct your input file."
!!$       call write_fatal(1)
!!$    endif
!!$  
!!$    do i = 1, num_subsys_runmodes
!!$       call loct_parse_block_int(blk, 0, i-1, subsys_runmodes(i))
!!$    enddo
!!$    call loct_parse_block_end(blk)
!!$    
!!$    
!!$    ! ... and in which order (replace this later by sorted mode priority)
!!$    if(loct_parse_block("SystemRunOrder", blk) == 0) then
!!$       num_subsys_runmodes = loct_parse_block_cols(blk,0)
!!$    else
!!$       message(1) = "Could not find block SystemRunOrder in the input file."
!!$       message(2) = "This block is mandatory for run mode MultiSubsystem."
!!$       call write_fatal(2)
!!$    endif
!!$    
!!$    if(num_subsys_runmodes/=num_syslabels) then
!!$       message(1) = "The blocks SystemLabels and SystemRunOrder do not have"
!!$       message(2) = "the same size. Please correct your input file."
!!$       call write_fatal(1)
!!$    endif
!!$    
!!$    do i = 1, num_subsys_runmodes
!!$       call loct_parse_block_int(blk, 0, i-1, run_order(i))
!!$    enddo
!!$    call loct_parse_block_end(blk)




!!$    commit log: this is a basic implementation for multi system support in octopus. it is the basis 
!!$                   for wave matching
!!$                   there is now a block SysLabels with NumberOfSystem entries. octopus is running
!!$                   all these systems (sequentially) and is for the respective system searching in the 
!!$                   input file for the usual variables with the self defined system labels attached. 
!!$
!!$                   a transport calculation would for example look like
!!$
!!$                   NumberOfSystems = 3
!!$                   %SystemLabels
!!$                   "LeftLead" | "CentralDevice" | "RightLead"
!!$                   %
!!$                   ...
!!$                   PeriodicDimensionsLeftLead = 2
!!$                   DimensionsLeftLead = 2
!!$                   BoxShapeLeftLead = parallelepiped
!!$                   BoxShapeCentralDevice = parallelepiped
!!$
!!$                   %SpacingLeftLead
!!$                   0.4 | 0.4 | 0.4
!!$                   %
!!$                   %LsizeLeftLead
!!$                   10.0 | 8.0 | 1.0
!!$                   %  
!!$
!!$                   The block SysEmbeddingType specifies what kind of boundary conditions ..

!!$    for transport we have then
!!$
!!$    NumberOfSystems = 3
!!$    %SysLabels
!!$      "LeftLead" | "CentralDevice" | "RightLead"
!!$    %
!!$    %SysBoundaryConditions
!!$         0 | 1 | 0                   ??
!!$       "eval" | "linsys" | "eval"    ??
!!$    %





!------------------------------------------------------------------------------------
! old code
!


!!$!  call f_der_init(sys%m, sys%f_der)
!!$!  call f_der_build(sys%f_der)
!!$
!!$!  ASSERT(associated(sys%f_der%der_discr%lapl))

!!$  if (associated(sys%f_der%der_discr%lapl)) then
!!$     write(*,*) 'yes'
!!$  else
!!$     write(*,*) 'no!'
!!$  endif
!!$
!!$
!!$  do k = 1, 4*boundary_np
!!$     write(*,'(I10,10F10.4)') k, sys%f_der%der_discr%lapl%w(:,k)
!!$     write(*,'(20I10)') k, sys%f_der%der_discr%lapl%i(:,k)
!!$     write(*,*)
!!$     write(*,*)
!!$  enddo
!!$
!!$!  write(*,*) sys%f_der%der_discr%lapl%w
!!$!  write(*,*) der_backup%lapl%w







!!$           do q = 1, qmax
!!$              write(*,*) 'L--',q, sys%f_der%der_discr%lapl%w(q+1,k)
!!$              ! store right hand side
!!$              b(k) = b(k) - boundary_points(1,2,sys%m%Lxyz_inv(ixl(1),ixl(2),ixl(3))) &
!!$                   *sys%f_der%der_discr%lapl%w(q+1,k)
!!$              ixl(l) = ixl(l) + 1
!!$              !boundary_points(1,2,k)
!!$           enddo


!!$           do q = 1, qmax
!!$              write(*,*) 'R--',q, sys%m%Lxyz_inv(ixl(1),ixl(2),ixl(3))
!!$! sys%f_der%der_discr%lapl%w(q+1,k)
!!$              ! store right hand side
!!$!              b(k) = b(k) - boundary_points(1,2,sys%m%Lxyz_inv(ixl(1),ixl(2),ixl(3))) &
!!$!                   *sys%f_der%der_discr%lapl%w(q+1,k)
!!$              ixl(l) = ixl(l) + 1
!!$           enddo

!        write(*,*)'1..:',sys%f_der%der_discr%op(conf%dim+1)%stencil(1:conf%dim,j),  &
!             mesh_index(sys%m, ix(1:conf%dim)+sys%f_der%der_discr%op(conf%dim+1)%stencil(1:conf%dim,j), 1), &
!             sys%f_der%der_discr%lapl%w(j,k)
        
        !        do l = 1, conf%dim ! l is actually align_axis




!!$  iter = 400
!!$  threshold = CNST(1.0e-5)
!!$  b5 = M_ONE + M_zI
!!$  num5 = 5
!!$!  call zconjugate_gradients(num5, x5, b5, op5, iter, zres, threshold)
!!$  call zconjugate_gradients(num5, x5, b5, op5, op5t, iter, zres, threshold)
!!$  
!!$  call op5(x5, y5)
!!$  write(*,*) b5
!!$  write(*,*) '--'
!!$  write(*,*) x5
!!$  write(*,*) '--'
!!$  write(*,*) y5

!!$
!!$                 if (.not.rhs) then
!!$                    rhs = .true.
!!$                    k0 = k
!!$                 endif
!!$        rhs = .false.
