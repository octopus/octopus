!! Copyright (C) 2018 N. Tancogne-Dejean 
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
!!

 subroutine X(lowdin_overlap)( basis )
  type(orbitalbasis_t),    intent(in)   :: basis
!  integer,          intent(in)       :: ik
!  logical,          intent(in)       :: has_phase

!  integer :: ios, ios2, im, im2, norbs, np
!  type(orbitalset_t), pointer :: os, os2

  PUSH_SUB(X(loewdin_overlap))

!  np = this%orbsets(1)%sphere%mesh%np

!  do ios = 1, this%norbsets
!    os => this%orbsets(ios)
!    norbs = this%orbsets(ios)%norbs
!    os%X(S)(1:norbs,1:this%maxnorbs,1:this%norbsets) = R_TOTYPE(M_ZERO)
!
!    !TODO: Use symmetry of the overlap matrices
!    norbs = this%orbsets(ios)%norbs
!
!    do im = 1, norbs
!
!      do ios2 = 1, this%norbsets
!        os2 => this%orbsets(ios2)
!
!        if(ios2 == ios) then
!          os%X(S)(im,im,ios2) = M_ONE
!        else
 !         if(this%overlap) then
 !         do im2 = 1, os2%norbs
 !           if(has_phase) then
 !#ifdef R_TCOMPLEX
 !             if(simul_box_is_periodic(os%sphere%mesh%sb) .and. .not. this%submeshforperiodic) then
 !               os%X(S)(im,im2,ios2) = zmf_dotp(os%sphere%mesh, &
 !                                         os%eorb_mesh(1:os%sphere%mesh%np,1,im,ik), &
 !                                         os2%eorb_mesh(1:os%sphere%mesh%np,1,im2,ik))
 !             else
 !               os%X(S)(im,im2,ios2) = zsubmesh_to_submesh_dotp(os2%sphere, &
 !                                        os2%eorb_submesh(1:os2%sphere%np,1,im2,ik), &
 !                                        os%sphere, os%eorb_submesh(1:os%sphere%np,1,im,ik))
 !             end if
 !#endif
 !           else
 !             os%X(S)(im,im2,ios2) = X(submesh_to_submesh_dotp)(os2%sphere, os2%X(orb)(1:os2%sphere%np,1,im2), &
 !                  os%sphere, os%X(orb)(1:os%sphere%np,1,im))
 !           end if
 !                                                                
 !           end do ! im2
 !         end if
 !       end if
 !     end do !im
 !   end do !ios2
 ! end do !ios 

  POP_SUB(X(loewdin_overlap))
 end subroutine X(lowdin_overlap) 
