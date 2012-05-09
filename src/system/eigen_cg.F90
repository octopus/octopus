!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: eigen_cg.F90 5954 2009-10-17 20:53:52Z xavier $

#include "global.h"

module eigen_cg_m
  use comm_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use lalg_basic_m
  use lalg_adv_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use preconditioners_m
  use profiling_m
  use states_m
  use states_calc_m

  implicit none

  private
  public ::                 &
    deigensolver_cg2,       &
    zeigensolver_cg2,       &
    deigensolver_cg2_new,   &
    zeigensolver_cg2_new,   &
    eigensolver_direct
contains

  subroutine eigensolver_direct(gr, st, hm, pre, tol, niter, converged, ik, diff)
    type(grid_t),           intent(in)    :: gr
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    type(preconditioner_t), intent(in)    :: pre
    FLOAT,                  intent(in)    :: tol
    integer,                intent(inout) :: niter
    integer,                intent(inout) :: converged
    integer,                intent(in)    :: ik
    FLOAT,        optional, intent(out)   :: diff(1:st%nst)
        
    CMPLX, allocatable :: psi(:, :), h_psi(:,:), h_rr(:,:), cL_rr(:,:), cR_rr(:,:), zeigenval(:), manyzeigenval(:)
    FLOAT, allocatable :: eigenval(:), manyeigenval(:), sortkey(:)
    integer, allocatable :: sortindices(:)
    FLOAT :: spacingsquared
    CMPLX :: kinetic_phase, tmp, tmp2
    logical :: fdkinetic

    integer       :: ib, jb, p, errcode
    

    SAFE_ALLOCATE(  psi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(h_psi(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(h_rr(1:gr%mesh%np, 1:gr%mesh%np))
    SAFE_ALLOCATE(cL_rr(1:2, 1:2))
    SAFE_ALLOCATE(cR_rr(1:gr%mesh%np, 1:gr%mesh%np))
    SAFE_ALLOCATE(zeigenval(1:st%nst))
    SAFE_ALLOCATE(manyzeigenval(1:gr%mesh%np))
    SAFE_ALLOCATE(sortkey(1:gr%mesh%np))
    SAFE_ALLOCATE(sortindices(1:gr%mesh%np))
    h_psi=M_z0 !R_TOTYPE(M_ZERO) ! presumably this sets h_psi = 0
    h_rr=M_z0
    psi=M_z0

    kinetic_phase = exp(-M_TWO * M_zI * hm%cmplxscl_th)
    
    print *,"theta",hm%cmplxscl_th
    
    spacingsquared = (gr%mesh%spacing(1) * gr%mesh%spacing(1))

    fdkinetic = .false.
!     fdkinetic = .true.
    
!     if (.not.fdkinetic) then
!        do ib = 1, gr%mesh%np
!           do jb = 1, gr%mesh%np
!              tmp = 0
!              do p = 1, (gr%mesh%np - 1) / 2
!                 tmp2 = p * M_PI / (gr%mesh%np * gr%mesh%spacing(1))
!                 tmp = tmp + cos((p*M_TWO*M_PI*(ib-jb))/gr%mesh%np)*2*tmp2*tmp2
!              end do
!              h_rr(jb, ib) = tmp * M_TWO * kinetic_phase / gr%mesh%np
!              if(jb==ib) then
!                h_rr(ib, ib) = h_rr(ib, ib) + hm%hm_base%potential(ib, 1) + M_zI * hm%hm_base%Impotential(ib, 1)
!              end if
!           end do
!        end do
!     end if

    if (.not.fdkinetic) then
      tmp2 = M_z0
       do ib = 1, gr%mesh%np
          do jb = 1, gr%mesh%np
             tmp = M_z0
             do p = 1, (gr%mesh%np - 1) / 2
               tmp = tmp + cos((p*M_TWO*M_PI*(ib-jb))/gr%mesh%np)*M_TWO*&
                 (((M_PI*p)/(gr%mesh%np * gr%mesh%spacing(1)))**2);
             end do
             h_rr(jb, ib) = tmp * M_TWO * kinetic_phase / gr%mesh%np
             if(jb==ib) then
               h_rr(ib, ib) = h_rr(ib, ib) + hm%hm_base%potential(ib, 1) + M_zI * hm%hm_base%Impotential(ib, 1)
               tmp2= tmp2 + h_rr(ib, ib)
             end if
          end do
       end do
    end if

    if (fdkinetic) then
      do ib = 1, gr%mesh%np
         ! kinetic fd stencil
            h_rr(ib, ib) = M_ONE / spacingsquared * kinetic_phase
            if (ib > 1) then
               h_rr(ib, ib - 1) = -(M_HALF / spacingsquared) * kinetic_phase
               h_rr(ib - 1, ib) = -(M_HALF / spacingsquared) * kinetic_phase
            end if
         h_rr(ib, ib) = h_rr(ib, ib) + hm%hm_base%potential(ib, 1) + M_zI * hm%hm_base%Impotential(ib, 1)
      end do
    end if

   
    cL_rr (1,1) = M_z0 
    cL_rr (2,2) = M_z0
    cL_rr (1,2) = M_z1 * kinetic_phase
    cL_rr (2,1) = M_z1 * kinetic_phase
        
!     call zmout(6, 2, 2, cL_rr, 2, -6, 'pauli scaled')
    
    call lalg_eigensolve_nonh(2, cL_rr, manyzeigenval, errcode, 'R')
    if (errcode.ne.0) then
       print*, 'something went wrong, errcode', errcode
    end if
    
    print *,"eigs", manyzeigenval(1:2)
!     call zmout(6, 2, 2, cL_rr, 2, -6, 'pauli scaled vecs')

!     cL_rr = h_rr
    cR_rr = h_rr
    
!     call zmout(6, st%nst, st%nst, h_rr, gr%mesh%np, -6, 'HRR')
    
    call lalg_eigensolve_nonh(gr%mesh%np, cR_rr, manyzeigenval, errcode, 'R')
    if (errcode.ne.0) then
       print*, 'something went wrong, errcode', errcode
    end if
!     call lalg_eigensolve_nonh(gr%mesh%np, cL_rr, manyzeigenval, errcode, 'L')
!     if (errcode.ne.0) then
!        print*, 'something went wrong, errcode', errcode
!     end if
    
!     call zmout(6, gr%mesh%np, 1, manyzeigenval, gr%mesh%np, -6, 'zeigenval')
     
    tmp = sum(manyzeigenval(:))
    print *, "the trace vals!!", tmp, "the trace H", tmp2

    !sortkey(:) = -imag(manyzeigenval(:))
    !sortkey(:) = real(manyzeigenval(:)) - imag(manyzeigenval(:))
    sortkey(:) = real(manyzeigenval(:))
!     sortkey(:) = abs(manyzeigenval(:))

    call sort(sortkey, sortindices)
    
    do p = 1, st%nst
       zeigenval(p) = manyzeigenval(sortindices(p))
       print *, p, zeigenval(p)
    end do
    
    do p = 1, st%nst
       do ib = 1, gr%mesh%np
!           st%psi%zL(ib, 1, p, 1) = cL_rr(ib, sortindices(p))
           st%psi%zR(ib, 1, p, 1) = cR_rr(ib, sortindices(p))
       end do
       st%zeigenval%Re(p, ik) = real(zeigenval(p))
       st%zeigenval%Im(p, ik) = aimag(zeigenval(p))
    end do

    if (present(diff)) diff = M_ZERO
    
    converged = st%nst
    niter = huge(1)

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(h_psi)
    SAFE_DEALLOCATE_A(h_rr)
    SAFE_DEALLOCATE_A(zeigenval)
    SAFE_DEALLOCATE_A(manyeigenval)
    SAFE_DEALLOCATE_A(manyzeigenval)
    SAFE_DEALLOCATE_A(sortkey)
    SAFE_DEALLOCATE_A(sortindices)
    SAFE_DEALLOCATE_A(cL_rr)
    SAFE_DEALLOCATE_A(cR_rr)
  end subroutine eigensolver_direct


 



#include "real.F90"
#include "eigen_cg_inc.F90"
#include "undef.F90"

#include "complex.F90"
#include "eigen_cg_inc.F90"
#include "undef.F90"

end module eigen_cg_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
