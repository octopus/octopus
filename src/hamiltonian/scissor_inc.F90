! Copyright (C) 2016 X. Andrade, N. Tancogne-Dejean
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

subroutine X(scissor_apply)(this, mesh, psib, hpsib)
  type(scissor_t),  intent(in)    :: this
  type(mesh_t),     intent(in)    :: mesh
  type(wfs_elec_t), intent(in)    :: psib
  type(wfs_elec_t), intent(inout) :: hpsib

  integer             :: ibatch, ist, idim
  R_TYPE              :: dot
  R_TYPE, allocatable :: psi(:,:)
  R_TYPE, allocatable :: hpsi(:,:)
  R_TYPE, allocatable :: gspsi(:,:)

  PUSH_SUB(scissor_apply)

  !Ading first gap*psi to hpsi
  call batch_axpy(mesh%np, this%gap, psib, hpsib)
 
  SAFE_ALLOCATE(psi(1:mesh%np, 1:this%gs_st%d%dim))
  SAFE_ALLOCATE(hpsi(1:mesh%np,1:this%gs_st%d%dim))
  SAFE_ALLOCATE(gspsi(1:mesh%np, 1:this%gs_st%d%dim))

  do ibatch = 1, psib%nst
    call batch_get_state(psib, ibatch, mesh%np, psi)
    call batch_get_state(hpsib,ibatch, mesh%np, hpsi)
    
    do ist = 1, this%gs_st%nst
      call states_elec_get_state(this%gs_st, mesh, ist, psib%ik, gspsi)

      dot = X(mf_dotp)(mesh, this%gs_st%d%dim, gspsi(:,:), psi) &
         * this%gs_st%occ(ist, psib%ik) / this%gs_st%smear%el_per_state
      do idim = 1, this%gs_st%d%dim
        call lalg_axpy(mesh%np, -this%gap*dot,  gspsi(:, idim), hpsi(:, idim))
      end do !idim
   end do !ist 

   call batch_set_state(hpsib, ibatch, mesh%np, hpsi)
  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(hpsi)
  SAFE_DEALLOCATE_A(gspsi)

  POP_SUB(scissor_apply) 
end subroutine X(scissor_apply)

subroutine X(scissor_commute_r)(this, mesh, ik, psi, gpsi)
   type(scissor_t), intent(in)    :: this
   type(mesh_t),    intent(in)    :: mesh 
   R_TYPE,          intent(in)    :: psi(:,:)
   integer,         intent(in)    :: ik
   R_TYPE,          intent(inout) :: gpsi(:, :, :)

   integer :: ist, idim, idir
   R_TYPE  :: dot
   R_TYPE, allocatable :: gspsi(:,:), tmpstate(:,:), psi_r(:,:)

   PUSH_SUB(scissor_commute_r)

   SAFE_ALLOCATE(gspsi(1:mesh%np, 1:this%gs_st%d%dim))
   SAFE_ALLOCATE(psi_r(1:mesh%np, 1:this%gs_st%d%dim))
   SAFE_ALLOCATE(tmpstate(1:mesh%np, 1:this%gs_st%d%dim))
   
   tmpstate(1:mesh%np, 1:this%gs_st%d%dim) = R_TOTYPE(M_ZERO)
   do ist = 1, this%gs_st%nst
     call states_elec_get_state(this%gs_st, mesh, ist, ik, gspsi )
     !<gpsi|psi>
     dot = X(mf_dotp)(mesh, this%gs_st%d%dim, gspsi, psi) &
           * this%gs_st%occ(ist, ik)/ this%gs_st%smear%el_per_state
     do idim = 1, this%gs_st%d%dim
      call lalg_axpy(mesh%np, dot,  gspsi(1:mesh%np, idim), tmpstate(1:mesh%np, idim)) 
     enddo
   enddo
   ! |gpsi> -= S x|gspsi><gspsi|psi>
   do idim = 1, this%gs_st%d%dim
     do idir = 1, mesh%sb%dim
       gpsi(1:mesh%np, idir, idim) = gpsi(1:mesh%np, idir, idim) &
             - this%gap * mesh%x(1:mesh%np, idir) * tmpstate(1:mesh%np, idim)
     enddo
   enddo

   do idir = 1, mesh%sb%dim
     do idim = 1, this%gs_st%d%dim
       psi_r(1:mesh%np, idim) = mesh%x(1:mesh%np,idir) * psi(1:mesh%np, idim)
     end do
     tmpstate(1:mesh%np,:) = R_TOTYPE(M_ZERO)
     do ist = 1, this%gs_st%nst
       call states_elec_get_state(this%gs_st, mesh, ist, ik, gspsi )
       ! <gspsi|r|psi>
       dot = X(mf_dotp)(mesh, this%gs_st%d%dim, gspsi, psi_r) &
         * this%gs_st%occ(ist, ik) / this%gs_st%smear%el_per_state
       do idim = 1, this%gs_st%d%dim
         call lalg_axpy(mesh%np, dot,  gspsi(1:mesh%np, idim), tmpstate(1:mesh%np, idim))
       end do
     enddo
     do idim = 1, this%gs_st%d%dim
       ! |gpsi> += S |gspsi><gspsi|x|psi>
       gpsi(1:mesh%np, idir, idim) = gpsi(1:mesh%np, idir, idim) &
             + this%gap * tmpstate(1:mesh%np, idim)
     enddo
   enddo

   SAFE_DEALLOCATE_A(gspsi)
   SAFE_DEALLOCATE_A(psi_r)
   SAFE_DEALLOCATE_A(tmpstate) 

   POP_SUB(scissor_commute_r)
end subroutine X(scissor_commute_r)
