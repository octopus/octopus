!! Copyright (C) 2014 Ask Hjorth Larsen
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

module eigen_feast_m

  use batch_m
  use datasets_m ! parse block and such
  use global_m
  use grid_m
  use hamiltonian_m
  use linear_solver_m
  use math_m ! some linear solvers are in math for some reason
  use messages_m
  use parser_m
  use profiling_m
  use states_m
  
  implicit none
  
  private
  public ::             &
    feast_t,            &
    zeigensolver_feast, &
    feast_init,         &
    feast_end

  type feast_t
    type(linear_solver_t) :: linear_solver

    CMPLX,   allocatable  :: Zedge(:)
    integer, allocatable  :: Tedge(:), Nedge(:)
    integer :: maxiter
    
    !integer :: ijob ! reverse communciation flag
  end type feast_t

contains

  subroutine feast_init(this, gr, nst)
    type(feast_t), intent(inout) :: this
    type(grid_t),  intent(in)    :: gr
    integer,       intent(in)    :: nst

    integer :: ncontourpts, icontour
    FLOAT :: recontour, imcontour
    integer :: nrows
    type(block_t) :: blk

    
    PUSH_SUB(feast_init)

    !%Variable FeastMaxIter
    !%Type integer 
    !%Section SCF::Eigensolver::FEAST
    !%Description 
    !% Maximum number of extra iterations that the FEAST eigensolver
    !% will perform per SCF step.  Must be >= 0.  0 means that only
    !% one iteration will be done.  Default is 20.
    !%End     
    call parse_integer(datasets_check('FeastMaxIter'), 20, this%maxiter)

    call linear_solver_init(this%linear_solver, gr, "FEAST", .false., LS_QMR_SYMMETRIC)

    !%Variable FeastContour
    !%Type block
    !%Section SCF::Eigensolver::FEAST
    !%Description
    !% A list of points tracing out a complex contour within which FEAST will
    !% search for eigenvalues.  Each row is one point.  Real and imaginary
    !% parts are read from first and second column respectively.  For now,
    !% points will be joined by straight line segments.
    !%
    !% Note: Contour must have negative orientation!
    !%End
    if(parse_block(datasets_check('FeastContour'), blk) == 0) then
      ncontourpts = parse_block_n(blk)
      SAFE_ALLOCATE(this%Zedge(1:ncontourpts))
      SAFE_ALLOCATE(this%Nedge(1:ncontourpts))
      SAFE_ALLOCATE(this%Tedge(1:ncontourpts))

      do icontour=1, ncontourpts
        !call parse_block_cmplx(blk, icontour - 1, 0, this%Zedge(icontour))
        call parse_block_float(blk, icontour - 1, 0, recontour)
        call parse_block_float(blk, icontour - 1, 1, imcontour)
        this%Zedge(icontour) = cmplx(recontour, imcontour, REAL_PRECISION)
      end do
    else
      message(1) = 'FeastContour block must be specified when using FEAST eigensolver.'
      call messages_fatal(1)
    end if

    this%Tedge(:) = M_ZERO ! type of each segment
    this%Nedge(:) = M_ONE ! number of interpolation points for each segment

    POP_SUB(feast_init)
  end subroutine feast_init

  subroutine feast_end(this)
    type(feast_t), intent(inout) :: this
    
    PUSH_SUB(feast_end)
    call linear_solver_end(this%linear_solver)

    SAFE_DEALLOCATE_A(this%Tedge)
    SAFE_DEALLOCATE_A(this%Zedge)
    SAFE_DEALLOCATE_A(this%Nedge)

    POP_SUB(feast_end)
  end subroutine feast_end
  
  subroutine zeigensolver_feast(this, gr, st, hm, converged, ik, diff)
    type(feast_t),          intent(inout) :: this
    type(grid_t),           intent(in)    :: gr
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(in)    :: hm
    integer,                intent(inout) :: converged
    integer,                intent(in)    :: ik
    FLOAT,        optional, intent(out)   :: diff(:) !< (1:st%nst)

    integer :: np, M0
    integer, dimension(64) :: feast_fpm
    CMPLX :: feast_Emid
    FLOAT :: feast_r
    CMPLX :: epsout
    integer :: feast_loopvar
    CMPLX, allocatable:: eigenval(:)
    CMPLX, allocatable:: XL(:, :),XR(:, :)
    integer :: nstates_within_contour
    FLOAT, allocatable :: resl(:), resr(:)
    integer :: info
    CMPLX, allocatable :: Zne(:), Wne(:)

    integer :: ii
    integer :: ijob, j,k,s
    CMPLX :: feast_Ze
    CMPLX, allocatable :: workl(:,:),workr(:,:),workc(:,:),zAq(:,:),zSq(:,:)
    FLOAT :: rea, img

    CMPLX, allocatable :: ls_xbuf(:, :), ls_ybuf(:, :)
    !F!LOAT, allocatable :: ls_residue(:)
    FLOAT :: ls_oneresidue
    integer :: ls_niter, ist

    type(batch_t) :: ls_xbatch, ls_ybatch
    integer, allocatable :: ls_iters(:)
    FLOAT, allocatable :: ls_residue(:)
    CMPLX, allocatable :: ls_Zebuf(:)

    CMPLX, allocatable, target :: xbatchbuf(:,:,:), ybatchbuf(:,:,:)

    PUSH_SUB(zeigensolver_feast)

! Let`s just compile an empty subroutine if FEAST is not there.  If it isn`t, we`ll raise an error
#ifdef HAVE_FEAST

    np = gr%mesh%np
    M0 = st%nst ! XXXXX keep name M0.  M0 might be changed by FEAST!
    ASSERT(st%d%dim == 1)
    SAFE_ALLOCATE(ls_xbuf(1:gr%mesh%np_part, 1))
    SAFE_ALLOCATE(ls_ybuf(1:gr%mesh%np, 1))

    ! I am pretty sure it must be possible to do away with some of all the buffers.
    ! Right now we allocate a lot and copy back and forth so that buffers have the correct size.
    ! XXX Try to reduce buffer allocation.
    SAFE_ALLOCATE(xbatchbuf(1:np, 1:1, 1:M0))
    SAFE_ALLOCATE(ybatchbuf(1:np, 1:1, 1:M0))
    xbatchbuf = M_ONE

    SAFE_ALLOCATE(eigenval(1:M0))
    ! Must figure how to distribute XR/XL
    SAFE_ALLOCATE(XR(1:np,1:M0))
    SAFE_ALLOCATE(XL(1:np,1:M0))
    SAFE_ALLOCATE(resr(1:M0))
    SAFE_ALLOCATE(resl(1:M0))

    call feastinit(feast_fpm)
    feast_fpm(1) = 1 ! print out run-time comments on screen
    feast_fpm(2) = sum(this%Nedge) ! total number of contour points.  Default 8 somehow
    feast_fpm(3) = 12 ! stopping convergence criteria for double precision eps=10^(-fpm(3))
    feast_fpm(4) = this%maxiter ! Maximum number of FEAST refinement loop allowed, >= 0 (default 20)
    feast_fpm(5) = 1 ! provide initial guess subspace
    !feast_fpm(6) ! type of convergence criterion (trace vs residual)
    !feast_fpm(7) ! convergence criterion for single precision
    !feast_fpm(9) is the communicator which defaults to mpi world.
    !feast_fpm(14) Return only subspace Q after 1 contour (0: No, 1: Yes)

    SAFE_ALLOCATE(Zne(1:feast_fpm(2))) !! Contains the complex valued contour points 
    SAFE_ALLOCATE(Wne(1:feast_fpm(2))) !! Contains the complex valued integrations weights
    call zfeast_customcontour(feast_fpm,size(this%Nedge, 1),this%Nedge,this%Tedge,this%Zedge,Zne,Wne)

    !open(10,file='contour.txt', status='replace')
    !do i=1,feast_fpm(2)
    !  write(10,*) dble(Zne(i)),aimag(Zne(i))
    !enddo
    !close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! allocate FEAST RCI variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    SAFE_ALLOCATE(zAq(1:M0,1:M0)) !FEAST RCI variables
    SAFE_ALLOCATE(zSq(1:M0,1:M0))
    SAFE_ALLOCATE(workl(1:np,1:M0))
    SAFE_ALLOCATE(workr(1:np,1:M0))
    SAFE_ALLOCATE(workc(1:np,1:M0))

    SAFE_ALLOCATE(ls_iters(1:st%nst))
    SAFE_ALLOCATE(ls_residue(1:st%nst))
    SAFE_ALLOCATE(ls_Zebuf(1:st%nst))

    !zAq = M_ZERO
    !zSq = M_ZERO

    do ist=1, M0
      call states_get_state(st, gr%mesh, ist, ik, XR(:, ist:ist))
    end do
    !XR(:, :) = M_ONE
    !XL(:, :) = M_ONE

    ijob = -1
    do while (ijob/=0) 
      !print *, 'FEEEEEEEAAAAAAST'!, sum(XR), sum(zAq), sum(zSq)
      call zfeast_srcix(&
        ijob,& ! communication flag IO
        np,& ! problem size gr%mesh%np, I
        feast_Ze,& ! IO but apparently O on first call
        workr,& ! N x M0 work array, O
        workl,& ! N x M0 work array, O
        workc,& ! N x M0 work array, O
        zAq,& ! M0 x M0 work array, O
        zSq,& ! M0 x M0 work array, O
        feast_fpm,& ! feast parameters, IO
        epsout,& ! error, O
        feast_loopvar,& ! number of iterations, O
        feast_Emid,& ! just some complex number, O
        feast_r,& ! just a double precision number.... O??
        M0,& ! number of eigenvalues / size of internal space
        eigenval,& ! M0 eigenvalues
        XR,& ! M x M0 right eigenvectors
        XL,& ! M x M0 left eigenvectors
        nstates_within_contour,& ! output
        resR,& ! M0 right residuals
        resL,& ! M0 left residuals
        info,& ! info==0 indicates success in the end.  O
        Zne,& ! contour data
        Wne& ! more contour data
        )
      select case(ijob)

      case(10) !! factorize (ZeB-A)

        ! Nothing for now.  Maybe some preconditioning

      case(11) !!solve the linear system (ZeB-A)x=workc(1:np,1:M0) result in to workc
        
        ybatchbuf(:, 1, :) = workc(1:np, 1:M0)
        
        call batch_init(ls_xbatch, st%d%dim, 1, st%nst, xbatchbuf)
        call batch_init(ls_ybatch, st%d%dim, 1, st%nst, ybatchbuf)

        ls_Zebuf(:) = -feast_Ze

        call zlinear_solver_solve_HXeY_batch(this%linear_solver, hm, gr, st, ik, &
          ls_xbatch, ls_ybatch, ls_Zebuf, CNST(2e-7), &
          ls_residue, ls_iters)
        
        workc(1:np, 1:M0) = xbatchbuf(1:np, 1, 1:M0)

        call batch_end(ls_xbatch)
        call batch_end(ls_ybatch)

        ! Non-batch version below
        !
        !do ist = 1, M0
        !  ls_xbuf(:, :) = M_ONE
        !  ls_ybuf(1:gr%mesh%np, 1) = workc(1:gr%mesh%np, ist)

        !  call zlinear_solver_solve_HXeY(this%linear_solver, hm, gr, st, ist, ik, ls_xbuf, ls_ybuf, &
        !    -feast_Ze, CNST(2e-7), ls_oneresidue, ls_niter)
        !  workc(1:gr%mesh%np, ist) = -ls_xbuf(1:gr%mesh%np, 1)
        !end do
        !workc(1:gr%mesh%np, 1:M0) = M_ZERO
        
      ! case 20 does not occur with symmetric matrices  
      !case(20) !! factorize (zeB-A)^T (needed if transpose of fact cannot be reused) 
      
      ! case 21 doesn`t occur either with symmetric matrices
      !case(21) !!solve the linear system (ZeB-A)x=work2(1:np,1:M0) result in to work2
      !

      case(30) !! perform multiplication A*x(1:np,fpm(24):fpm(24)+fpm(25)-1) result in workr(1:np,fpm(24)+fpm(25)-1)
        do ist=feast_fpm(24), feast_fpm(24) + feast_fpm(25) - 1
          ls_xbuf(:, :) = M_ZERO
          ls_xbuf(1:np, 1) = XR(1:np, ist)
          call zhamiltonian_apply(hm, gr%der, ls_xbuf, ls_ybuf, ist, ik)
          workr(1:np, ist) = ls_ybuf(1:np, 1)
        end do

      case(31) !! perform multiplication A^C*x(1:np,fpm(24):fpm(24)+fpm(25)-1) result in workl(1:np,fpm(24)+fpm(25)-1)
        workl(1:np, feast_fpm(24):feast_fpm(24) + feast_fpm(25) - 1) = M_ZERO

      case(40) !! perform multiplication B*x(1:np,fpm(24):fpm(24)+fpm(25)-1) result in workr(1:np,fpm(24)+fpm(25)-1)
        workr(1:np, feast_fpm(24):feast_fpm(24) + feast_fpm(25) - 1) = XR(1:np, feast_fpm(24):feast_fpm(24) + feast_fpm(25) - 1)

      case(41) !! perform multiplication B^C*x(1:np,fpm(24):fpm(24)+fpm(25)-1) result in workl(1:np,fpm(24)+fpm(25)-1)
        workl(1:np, feast_fpm(24) + feast_fpm(25) - 1) = XL(1:np, feast_fpm(24) + feast_fpm(25) - 1)

      end select
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! FIXME: use messages_info to produce something sensible in parallel.
    print *,'FEAST OUTPUT INFO',info
    if (info==0) then
      print *,'*************************************************'
      print *,'************** REPORT ***************************'
      print *,'*************************************************'
      print *,'# states found/subspace', nstates_within_contour, M0
      print *,'# iterations', feast_loopvar
      print *,'TRACE',sum(eigenval(1:nstates_within_contour))
      print *,'Relative error on the Trace',dble(epsout)
      print *,'Eigenvalues/Max(right/left Residuals)'
      do ist=1, nstates_within_contour
        print *, ist, eigenval(ist),max(resl(ist),resr(ist))
      enddo
    endif

    st%eigenval(:, :) = M_ZERO

    st%eigenval(1:st%nst, ik) = real(eigenval(1:st%nst), REAL_PRECISION)
    ASSERT(associated(st%zeigenval%Im))
    st%zeigenval%Im(1:st%nst, ik) = aimag(eigenval(1:st%nst))

    do ist=1, st%nst
      call states_set_state(st, gr%mesh, ist, ik, XR(:, ist:ist))
    end do

    if(present(diff)) then
      diff(:) = resr(:)
    end if

    SAFE_DEALLOCATE_A(xbatchbuf)
    SAFE_DEALLOCATE_A(ybatchbuf)
    SAFE_DEALLOCATE_A(ls_Zebuf)
    SAFE_DEALLOCATE_A(ls_iters)
    SAFE_DEALLOCATE_A(ls_residue)
    SAFE_DEALLOCATE_A(zAq)
    SAFE_DEALLOCATE_A(zSq)
    SAFE_DEALLOCATE_A(workl)
    SAFE_DEALLOCATE_A(workr)
    SAFE_DEALLOCATE_A(workc)
    SAFE_DEALLOCATE_A(eigenval)
    SAFE_DEALLOCATE_A(XR)
    SAFE_DEALLOCATE_A(XL)
    SAFE_DEALLOCATE_A(resr)
    SAFE_DEALLOCATE_A(resl)
    SAFE_DEALLOCATE_A(ls_xbuf)
    SAFE_DEALLOCATE_A(ls_ybuf)
#else
    message(1) = 'FEAST eigensolver called, but FEAST is not compiled in.  Must compile with FEAST or not call FEAST eigensolver.'
    call messages_fatal(1)
#endif
    POP_SUB(zeigensolver_feast)
  end subroutine zeigensolver_feast
end module

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
