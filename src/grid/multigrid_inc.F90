!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
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

  ! ---------------------------------------------------------
  subroutine X(multigrid_coarse2fine)(tt, coarse_der, fine_mesh, f_coarse, f_fine, order)
    type(transfer_table_t),  intent(in)    :: tt
    type(derivatives_t),     intent(in)    :: coarse_der
    type(mesh_t),            intent(in)    :: fine_mesh
    R_TYPE,                  intent(inout) :: f_coarse(:)
    R_TYPE,                  intent(out)   :: f_fine(:)
    integer, optional,       intent(in)    :: order

    integer :: idir, order_, ii, ifactor
    integer :: ipc, ipf, ipfg, xf(1:3), xc(1:3), dd(1:3)
    FLOAT, allocatable :: factor(:), points(:)

    PUSH_SUB(X(multigrid_coarse2fine))

    call profiling_in(interp_prof, "MG_INTERPOLATION")

    ASSERT(ubound(f_coarse, dim = 1) == coarse_der%mesh%np_part)
    ASSERT(coarse_der%mesh%np == tt%n_coarse)

    order_ = 1
    if(present(order)) order_ = order

    SAFE_ALLOCATE(points(1:2*order_))
    SAFE_ALLOCATE(factor(1:2*order_))

    forall(ii = 1:2*order_) points(ii) = ii

    call interpolation_coefficients(2*order_, points, order_ + M_HALF, factor)

    factor = factor/coarse_der%mesh%sb%dim

    call X(derivatives_set_bc)(coarse_der, f_coarse)

#ifdef HAVE_MPI
    if(coarse_der%mesh%parallel_in_domains) call X(vec_ghost_update)(coarse_der%mesh%vp, f_coarse)
#endif

    do ipf = 1, fine_mesh%np
      
      ipfg = ipf
#ifdef HAVE_MPI
      ! translate to a global index
      if(fine_mesh%parallel_in_domains) ipfg = fine_mesh%vp%local(ipf - 1 + fine_mesh%vp%xlocal(fine_mesh%vp%partno))
#endif 
      xf(1:3) = fine_mesh%idx%lxyz(ipfg, 1:3)

      dd = mod(xf, 2)
      
      f_fine(ipf) = M_ZERO
      
      do idir = 1, coarse_der%mesh%sb%dim
        ifactor = 1
        do ii = -order_, order_
          if(ii == 0) cycle
          xc = xf + (2*ii - sign(1, ii))*dd
          xc = xc/2
          ipc = coarse_der%mesh%idx%lxyz_inv(xc(1), xc(2), xc(3))
#ifdef HAVE_MPI
            ! translate to a local index
          if(coarse_der%mesh%parallel_in_domains) ipc = vec_global2local(coarse_der%mesh%vp, ipc, coarse_der%mesh%vp%partno)
#endif
          f_fine(ipf) = f_fine(ipf) + factor(ifactor)*f_coarse(ipc)
          ifactor = ifactor + 1
        end do
      end do

    end do

    call profiling_out(interp_prof)
    POP_SUB(X(multigrid_coarse2fine))
  end subroutine X(multigrid_coarse2fine)

  ! ---------------------------------------------------------
  subroutine X(multigrid_fine2coarse)(tt, fine_der, coarse_mesh, f_fine, f_coarse, method_p)
    type(transfer_table_t), intent(in)    :: tt
    type(derivatives_t),    intent(in)    :: fine_der
    type(mesh_t),           intent(in)    :: coarse_mesh
    R_TYPE,                 intent(inout) :: f_fine(:)
    R_TYPE,                 intent(out)   :: f_coarse(:)
    integer, optional,      intent(in)    :: method_p

    integer :: method

    PUSH_SUB(X(multigrid_fine2coarse))

    if(present(method_p)) then
      method=method_p
    else
      method=FULLWEIGHT
    end if

    select case(method)
    case(FULLWEIGHT)
      call X(multigrid_restriction)(tt, fine_der, coarse_mesh, f_fine, f_coarse)
    case(INJECTION)
      call X(multigrid_injection)(tt, f_fine, f_coarse)
    case default
      write(message(1), '(a,i2,a)') 'Multigrid: Restriction method  = ', method, ' is not valid.'
      call messages_fatal(1)
    end select

    POP_SUB(X(multigrid_fine2coarse))
  end subroutine X(multigrid_fine2coarse)


  ! ---------------------------------------------------------
  subroutine X(multigrid_injection)(tt, f_fine, f_coarse)
    type(transfer_table_t), intent(in)  :: tt
    R_TYPE,                 intent(in)  :: f_fine(:)
    R_TYPE,                 intent(out) :: f_coarse(:)

    integer :: ii

    PUSH_SUB(X(multigrid_injection))
    call profiling_in(injection_prof, "MG_INJECTION")

    do ii = 1, tt%n_coarse
      f_coarse(ii) = f_fine(tt%to_coarse(ii))
    end do

    call profiling_out(injection_prof)
    POP_SUB(X(multigrid_injection))
  end subroutine X(multigrid_injection)

  ! ---------------------------------------------------------
  subroutine X(multigrid_restriction)(tt, fine_der, coarse_mesh, f_fine, f_coarse)
    type(transfer_table_t), intent(in)    :: tt
    type(derivatives_t),    intent(in)    :: fine_der
    type(mesh_t),           intent(in)    :: coarse_mesh
    R_TYPE,                 intent(inout) :: f_fine(:)
    R_TYPE,                 intent(out)   :: f_coarse(:)

    FLOAT :: weight(-1:1,-1:1,-1:1)

    integer :: nn, fn, di, dj, dk, dd, fi(MAX_DIM)

    PUSH_SUB(X(multigrid_restriction))
    call profiling_in(restrict_prof, "MG_RESTRICTION")

    do di = -1, 1
      do dj = -1, 1
        do dk = -1, 1
          dd = abs(di) + abs(dj) + abs(dk)
          weight(di, dj, dk) = CNST(0.5)**dd
        end do
      end do
    end do

    call X(derivatives_set_bc)(fine_der, f_fine)

#ifdef HAVE_MPI
    if(fine_der%mesh%parallel_in_domains) call X(vec_ghost_update)(fine_der%mesh%vp, f_fine)
#endif

    do nn = 1, tt%n_coarse
      fn = tt%to_coarse(nn)
#ifdef HAVE_MPI
      ! translate to a global index
      if(fine_der%mesh%parallel_in_domains) then
        fn = fine_der%mesh%vp%local(fn - 1 + fine_der%mesh%vp%xlocal(fine_der%mesh%vp%partno))
      end if
#endif
      fi(:) = fine_der%mesh%idx%lxyz(fn, :)

      f_coarse(nn) = M_ZERO

      do di = -1, 1
        do dj = -1, 1
          do dk = -1, 1
            fn = fine_der%mesh%idx%lxyz_inv(fi(1) + di, fi(2) + dj, fi(3) + dk)

#ifdef HAVE_MPI
            ! translate to a local index
            if(fine_der%mesh%parallel_in_domains) fn = vec_global2local(fine_der%mesh%vp, fn, fine_der%mesh%vp%partno)
#endif
            if(fine_der%mesh%use_curvilinear) then
              f_coarse(nn) = f_coarse(nn) + weight(di, dj, dk)*f_fine(fn)*fine_der%mesh%vol_pp(fn)
            else
              f_coarse(nn) = f_coarse(nn) + weight(di, dj, dk)*f_fine(fn)*fine_der%mesh%vol_pp(1)
            end if

          end do
        end do
      end do

      if(fine_der%mesh%use_curvilinear) then
        f_coarse(nn) = f_coarse(nn)/coarse_mesh%vol_pp(nn)
      else
        f_coarse(nn) = f_coarse(nn)/coarse_mesh%vol_pp(1)
      end if
    end do

    call profiling_count_operations(tt%n_coarse*(27*3 + 1))
    call profiling_out(restrict_prof)
    POP_SUB(X(multigrid_restriction))
  end subroutine X(multigrid_restriction)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
