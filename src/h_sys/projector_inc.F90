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
!! $Id$

!------ -----------------------------------------------------------------------
! X(project_psi) calculates the action of a projector on the psi wavefunction.
! The result is summed up to ppsi
subroutine X(project_psi)(mesh, pj, npj, dim, psi, ppsi, reltype, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj(:)
  integer,           intent(in)    :: npj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: ppsi(:, :)  ! ppsi(1:mesh%np, dim)
  integer,           intent(in)    :: reltype
  integer,           intent(in)    :: ik

  integer :: ipj, nreduce, ii, ns, idim, ll, mm
  R_TYPE, allocatable :: reduce_buffer(:), lpsi(:, :)
  integer, allocatable :: ireduce(:, :, :)
#if defined(HAVE_MPI)
  R_TYPE, allocatable   :: reduce_buffer_dest(:)
  type(profile_t), save :: reduce_prof
#endif

  call push_sub('projector_inc.project_psi')

  ALLOCATE(ireduce(1:npj, 0:MAX_L, -MAX_L:MAX_L), npj*(MAX_L + 1)*(2*MAX_L + 1))

  ! generate the reduce buffer and related structures
  nreduce = 0
  do ipj = 1, npj
    if(pj(ipj)%type == M_NONE) cycle
    do ll = 0, pj(ipj)%lmax
      if (ll == pj(ipj)%lloc) cycle
      do mm = -ll, ll
        ireduce(ipj, ll, mm) = 1 + nreduce
        nreduce = nreduce + pj(ipj)%reduce_size
      end do
    end do
  end do

  ALLOCATE(reduce_buffer(1:nreduce), nreduce)
#if defined(HAVE_MPI)
  ALLOCATE(reduce_buffer_dest(1:nreduce), nreduce)
#endif

  ! calculate <p|psi>
  do ipj = 1, npj
    if(pj(ipj)%type == M_NONE) cycle

    ns = pj(ipj)%sphere%ns

    ALLOCATE(lpsi(ns, dim),  ns*dim)

    do idim = 1, dim
      if(simul_box_is_periodic(mesh%sb)) then
        lpsi(1:ns, idim) = psi(pj(ipj)%sphere%jxyz(1:ns), idim)*pj(ipj)%phase(1:ns, ik)
        call profiling_count_operations(ns*R_MUL)
      else
        lpsi(1:ns, idim) = psi(pj(ipj)%sphere%jxyz(1:ns), idim)
      end if
    end do

    do ll = 0, pj(ipj)%lmax
      if (ll == pj(ipj)%lloc) cycle
      do mm = -ll, ll

        ii = ireduce(ipj, ll, mm)

        select case(pj(ipj)%type)
        case(M_KB)
          call X(kb_project_bra)(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(ll, mm), dim, lpsi, reduce_buffer(ii:))
        case(M_RKB)
#ifdef R_TCOMPLEX
          if(ll /= 0) then
            call rkb_project_bra(mesh, pj(ipj)%sphere, pj(ipj)%rkb_p(ll, mm), lpsi, reduce_buffer(ii:))
          else
            call zkb_project_bra(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(1, 1), dim, lpsi, reduce_buffer(ii:))
          end if
#endif
        case(M_HGH)
          call X(hgh_project_bra)(mesh, pj(ipj)%sphere, pj(ipj)%hgh_p(ll, mm), dim, reltype, lpsi, reduce_buffer(ii:))
        end select

      end do ! mm
    end do ! ll

    deallocate(lpsi)

  end do

  ! reduce <p|psi>
#if defined(HAVE_MPI)
  if(mesh%parallel_in_domains) then
    call profiling_in(reduce_prof, "VNLPSI_REDUCE")
    call MPI_Allreduce(reduce_buffer, reduce_buffer_dest, nreduce, R_MPITYPE, MPI_SUM, mesh%vp%comm, mpi_err)
    reduce_buffer = reduce_buffer_dest
    call profiling_out(reduce_prof)
  end if
#endif

  ! calculate |ppsi> += |p><p|psi>
  do ipj = 1, npj
    if(pj(ipj)%type == M_NONE) cycle
    
    ns = pj(ipj)%sphere%ns
    ALLOCATE(lpsi(ns, dim),  ns*dim)

    lpsi = M_ZERO

    do ll = 0, pj(ipj)%lmax
      if (ll == pj(ipj)%lloc) cycle
      do mm = -ll, ll

        ii = ireduce(ipj, ll, mm)

        select case(pj(ipj)%type)
        case(M_KB)
          reduce_buffer(ii:ii + 1) = reduce_buffer(ii:ii + 1)*pj(ipj)%kb_p(ll, mm)%e(1:2)
          call X(kb_project_ket)(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(ll, mm), dim, reduce_buffer(ii:), lpsi)
        case(M_RKB)
#ifdef R_TCOMPLEX
          if(ll /= 0) then
            call rkb_project_ket(mesh, pj(ipj)%sphere, pj(ipj)%rkb_p(ll, mm), reduce_buffer(ii:), lpsi)
          else
            call zkb_project_ket(mesh, pj(ipj)%sphere, pj(ipj)%kb_p(1, 1), dim, reduce_buffer(ii:), lpsi)
          end if
#endif
        case(M_HGH)
          call X(hgh_project_ket)(mesh, pj(ipj)%sphere, pj(ipj)%hgh_p(ll, mm), dim, reltype, reduce_buffer(ii:), lpsi)
        end select

      end do ! mm
    end do ! ll

    do idim = 1, dim
      if(simul_box_is_periodic(mesh%sb)) then
        ppsi(pj(ipj)%sphere%jxyz(1:ns), idim) = &
          ppsi(pj(ipj)%sphere%jxyz(1:ns), idim) + lpsi(1:ns, idim)*conjg(pj(ipj)%phase(1:ns, ik))
        call profiling_count_operations(ns*(R_ADD + R_MUL))
      else
        ppsi(pj(ipj)%sphere%jxyz(1:ns), idim) = ppsi(pj(ipj)%sphere%jxyz(1:ns), idim) + lpsi(1:ns, idim)
        call profiling_count_operations(ns*R_ADD)
      end if
    end do
    
    deallocate(lpsi)

  end do

  deallocate(reduce_buffer, ireduce)

#ifdef HAVE_MPI
  deallocate(reduce_buffer_dest)
#endif

  call pop_sub()
end subroutine X(project_psi)

!------------------------------------------------------------------------------
! X(psia_project_psib) calculates <psia|projector|psib>
R_TYPE function X(psia_project_psib)(mesh, pj, dim, psia, psib, reltype, ik) result(apb)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psia(:, :)  ! psia(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: psib(:, :)  ! psib(1:mesh%np, dim)
  integer,           intent(in)    :: reltype
  integer,           intent(in)    :: ik

  integer ::  ns, idim
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)

  call push_sub('projector_inc.psia_project_psib')

  ns = pj%sphere%ns

  ALLOCATE(lpsi(ns, dim),  ns*dim)
  ALLOCATE(plpsi(ns, dim), ns*dim)

  do idim = 1, dim
    if(simul_box_is_periodic(mesh%sb)) then
      lpsi(1:ns, idim) = psib(pj%sphere%jxyz(1:ns), idim)*pj%phase(1:ns, ik)
    else
      lpsi(1:ns, idim) = psib(pj%sphere%jxyz(1:ns), idim)
    end if
  end do
      
  call X(project_sphere)(mesh, pj, dim, lpsi, plpsi, reltype)

  apb = M_ZERO
  do idim = 1, dim
    if(simul_box_is_periodic(mesh%sb)) then
      plpsi(1:ns, idim) = R_CONJ(psia(pj%sphere%jxyz(1:ns), idim))*plpsi(1:ns, idim)*conjg(pj%phase(1:ns, ik))
    else
      plpsi(1:ns, idim) = R_CONJ(psia(pj%sphere%jxyz(1:ns), idim))*plpsi(1:ns, idim)
    end if
    apb = apb + X(sm_integrate)(mesh, pj%sphere, plpsi(1:ns, idim))
  end do

  deallocate(lpsi, plpsi)

  call pop_sub()
end function X(psia_project_psib)

subroutine X(project_sphere)(mesh, pj, dim, psi, ppsi, reltype)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:ns, dim)
  R_TYPE,            intent(out)   :: ppsi(:, :)  ! ppsi(1:ns, dim)
  integer,           intent(in)    :: reltype

  integer :: ns, ll, mm

  call push_sub('projector_inc.project_sphere')

  ns = pj%sphere%ns

  ppsi = M_ZERO

  do ll = 0, pj%lmax
    if (ll == pj%lloc) cycle
    do mm = -ll, ll
      
      select case (pj%type)
      case (M_HGH)
        call X(hgh_project)(mesh, pj%sphere, pj%hgh_p(ll, mm), dim, psi, ppsi, reltype)
      case (M_KB)
        call X(kb_project)(mesh, pj%sphere, pj%kb_p(ll, mm), dim, psi, ppsi)
      case (M_RKB)
#ifdef R_TCOMPLEX
        if(ll /= 0) then
          call rkb_project(mesh, pj%sphere, pj%rkb_p(ll, mm), psi, ppsi)
        else
          call zkb_project(mesh, pj%sphere, pj%kb_p(1, 1), dim, psi, ppsi)
        end if
#endif
      end select
  
    end do
  end do

  call pop_sub()
end subroutine X(project_sphere)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
