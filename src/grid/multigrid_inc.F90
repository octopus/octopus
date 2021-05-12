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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

  ! ---------------------------------------------------------
  subroutine X(multigrid_coarse2fine)(tt, coarse_der, fine_mesh, f_coarse, f_fine, order, set_bc)
    type(transfer_table_t),  intent(in)    :: tt
    type(derivatives_t),     intent(in)    :: coarse_der
    type(mesh_t),            intent(in)    :: fine_mesh
    R_TYPE,                  intent(inout) :: f_coarse(:)
    R_TYPE,                  intent(out)   :: f_fine(:)
    integer, optional,       intent(in)    :: order
    logical, optional,       intent(in)    :: set_bc

    integer :: idir, order_, ii, ifactor
    integer :: ipc, ipf, xf(coarse_der%dim), xc(coarse_der%dim), dd(coarse_der%dim)
    FLOAT, allocatable :: factor(:), points(:)

    PUSH_SUB(X(multigrid_coarse2fine))

    call profiling_in(interp_prof, TOSTRING(X(MG_INTERPOLATION)))

    ASSERT(ubound(f_coarse, dim = 1) == coarse_der%mesh%np_part)
    ASSERT(coarse_der%mesh%np == tt%n_coarse)

    order_ = 1
    if(present(order)) order_ = order

    SAFE_ALLOCATE(points(1:2*order_))
    SAFE_ALLOCATE(factor(1:2*order_))

    do ii = 1, 2*order_
      points(ii) = ii
    end do

    call interpolation_coefficients(2*order_, points, order_ + M_HALF, factor)

    factor = factor/coarse_der%dim

    if(optional_default(set_bc, .true.)) then
      call boundaries_set(coarse_der%boundaries, f_coarse)

#ifdef HAVE_MPI
      if(coarse_der%mesh%parallel_in_domains) call X(vec_ghost_update)(coarse_der%mesh%vp, f_coarse)
#endif
    end if

    !We perform a trilinear interpolation, see https://en.wikipedia.org/wiki/Trilinear_interpolation
    do ipf = 1, fine_mesh%np
      call mesh_local_index_to_coords(fine_mesh, ipf, xf)

      dd = mod(xf, 2)

      if(all(dd == 0)) then ! This point belongs to the coarse grid
        xc = xf/2
        ipc = mesh_local_index_from_coords(coarse_der%mesh, xc)
        f_fine(ipf) = f_coarse(ipc)
        cycle
      end if

      
      f_fine(ipf) = M_ZERO
      
      do idir = 1, coarse_der%dim
        ifactor = 1
        do ii = -order_, order_
          if(ii == 0) cycle
          xc = xf + (2*ii - sign(1, ii))*dd
          xc = xc/2
          ipc = mesh_local_index_from_coords(coarse_der%mesh, xc)
          f_fine(ipf) = f_fine(ipf) + factor(ifactor)*f_coarse(ipc)
          ifactor = ifactor + 1
        end do
        dd(idir) = -dd(idir)
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
    call profiling_in(injection_prof, TOSTRING(X(MG_INJECTION)))

    do ii = 1, tt%n_coarse
      f_coarse(ii) = f_fine(tt%to_coarse(ii))
    end do

    call profiling_out(injection_prof)
    POP_SUB(X(multigrid_injection))
  end subroutine X(multigrid_injection)

  ! ---------------------------------------------------------
  subroutine X(multigrid_restriction)(tt, fine_der, coarse_mesh, f_fine, f_coarse, set_bc)
    type(transfer_table_t), intent(in)    :: tt
    type(derivatives_t),    intent(in)    :: fine_der
    type(mesh_t),           intent(in)    :: coarse_mesh
    R_TYPE,                 intent(inout) :: f_fine(:)
    R_TYPE,                 intent(out)   :: f_coarse(:)
    logical, optional,      intent(in)    :: set_bc

    FLOAT :: weight(-1:1,-1:1,-1:1)

    integer :: nn, fn, di, dj, dk, dd, fi(MAX_DIM)

    PUSH_SUB(X(multigrid_restriction))
    call profiling_in(restrict_prof, TOSTRING(X(MG_RESTRICTION)))

    !The following code only works in 3D
    ASSERT(fine_der%dim == 3)

    do di = -1, 1
      do dj = -1, 1
        do dk = -1, 1
          dd = abs(di) + abs(dj) + abs(dk)
          weight(di, dj, dk) = CNST(0.5)**dd
        end do
      end do
    end do

    if(optional_default(set_bc, .true.)) then
      call boundaries_set(fine_der%boundaries, f_fine)

#ifdef HAVE_MPI
      if(fine_der%mesh%parallel_in_domains) call X(vec_ghost_update)(fine_der%mesh%vp, f_fine)
#endif
    end if

    do nn = 1, tt%n_coarse
      fn = tt%to_coarse(nn)
      call mesh_local_index_to_coords(fine_der%mesh, fn, fi)

      f_coarse(nn) = M_ZERO

      do di = -1, 1
        do dj = -1, 1
          do dk = -1, 1
            fn = mesh_local_index_from_coords(fine_der%mesh, [fi(1) + di, fi(2) + dj, fi(3) + dk])

            if(fine_der%mesh%use_curvilinear) then
              f_coarse(nn) = f_coarse(nn) + weight(di, dj, dk)*f_fine(fn)*fine_der%mesh%vol_pp(fn)
            else
              f_coarse(nn) = f_coarse(nn) + weight(di, dj, dk)*f_fine(fn)
            end if

          end do
        end do
      end do

      if(fine_der%mesh%use_curvilinear) then
        f_coarse(nn) = f_coarse(nn)/coarse_mesh%vol_pp(nn)
      else
        f_coarse(nn) = f_coarse(nn)*fine_der%mesh%vol_pp(1)/coarse_mesh%vol_pp(1)
      end if
    end do

    call profiling_count_operations(tt%n_coarse*(27*3 + 1))
    call profiling_out(restrict_prof)
    POP_SUB(X(multigrid_restriction))
  end subroutine X(multigrid_restriction)

  ! ---------------------------------------------------------
  subroutine X(multigrid_coarse2fine_batch)(tt, coarse_der, fine_mesh, coarseb, fineb, order)
    type(transfer_table_t),  intent(in)    :: tt
    type(derivatives_t),     intent(in)    :: coarse_der
    type(mesh_t),            intent(in)    :: fine_mesh
    class(batch_t),          intent(inout) :: coarseb
    class(batch_t),          intent(inout) :: fineb
    integer, optional,       intent(in)    :: order

    integer :: idir, order_, ii, ifactor, ist
    integer :: ipc, ipf, xf(coarse_der%dim), xc(coarse_der%dim), dd(coarse_der%dim)
    FLOAT, allocatable :: factor(:), points(:)
    R_TYPE, allocatable :: f_coarse(:), f_fine(:)
    type(pv_handle_batch_t) :: handle

    PUSH_SUB(X(multigrid_coarse2fine_batch))

    call profiling_in(interp_prof, TOSTRING(X(MG_INTERPOLATION_BATCH)))

    ASSERT(coarseb%nst_linear == fineb%nst_linear)

    order_ = optional_default(order, 1)

    SAFE_ALLOCATE(points(1:2*order_))
    SAFE_ALLOCATE(factor(1:2*order_))

    do ii = 1, 2*order_
      points(ii) = ii
    end do

    call interpolation_coefficients(2*order_, points, order_ + M_HALF, factor)

    factor = factor/coarse_der%dim

    call boundaries_set(coarse_der%boundaries, coarseb)
    if(coarse_der%mesh%parallel_in_domains) then
      call X(ghost_update_batch_start)(coarse_der%mesh%vp, coarseb, handle)
      call X(ghost_update_batch_finish)(handle)
    end if

    ASSERT(fineb%status() == coarseb%status())

    select case(fineb%status())
    case(BATCH_PACKED)
      !We perform a trilinear interpolation, see https://en.wikipedia.org/wiki/Trilinear_interpolation
      do ipf = 1, fine_mesh%np
        call mesh_local_index_to_coords(fine_mesh, ipf, xf)

        dd = mod(xf, 2)

        if(all(dd == 0)) then ! This point belongs to the coarse grid
          xc = xf/2
          ipc = mesh_local_index_from_coords(coarse_der%mesh, xc)
          do ist = 1, coarseb%nst_linear
            fineb%X(ff_pack)(ist, ipf) = coarseb%X(ff_pack)(ist, ipc)
          end do
          cycle
        end if

        do ist = 1, coarseb%nst_linear
          fineb%X(ff_pack)(ist, ipf) = M_ZERO
        end do

        do idir = 1, coarse_der%dim
          ifactor = 1
          do ii = -order_, order_
            if(ii == 0) cycle
            xc = xf + (2*ii - sign(1, ii))*dd
            xc = xc/2
            ipc = mesh_local_index_from_coords(coarse_der%mesh, xc)
            do ist = 1, coarseb%nst_linear
              fineb%X(ff_pack)(ist, ipf) = fineb%X(ff_pack)(ist, ipf) + factor(ifactor) * coarseb%X(ff_pack)(ist, ipc)
            end do
            ifactor = ifactor + 1
          end do
          dd(idir) = -dd(idir)
        end do

      end do

    case(BATCH_DEVICE_PACKED, BATCH_NOT_PACKED)
      SAFE_ALLOCATE(f_coarse(1:coarse_der%mesh%np_part))
      SAFE_ALLOCATE(f_fine(1:fine_mesh%np))
      do ist = 1, coarseb%nst_linear
        call batch_get_state(coarseb, ist, coarse_der%mesh%np_part, f_coarse)
        call batch_get_state(fineb, ist, fine_mesh%np, f_fine)

        call X(multigrid_coarse2fine)(tt, coarse_der, fine_mesh, f_coarse, f_fine, order=order, set_bc = .false.) 

        call batch_set_state(fineb, ist, fine_mesh%np, f_fine)
      end do

      SAFE_DEALLOCATE_A(f_coarse)
      SAFE_DEALLOCATE_A(f_fine)

    end select

    call profiling_out(interp_prof)
    POP_SUB(X(multigrid_coarse2fine_batch))
  end subroutine X(multigrid_coarse2fine_batch)

    ! ---------------------------------------------------------
  subroutine X(multigrid_fine2coarse_batch)(tt, fine_der, coarse_mesh, fineb, coarseb, method_p)
    type(transfer_table_t), intent(in)    :: tt
    type(derivatives_t),    intent(in)    :: fine_der
    type(mesh_t),           intent(in)    :: coarse_mesh
    class(batch_t),         intent(inout) :: fineb
    class(batch_t),         intent(inout) :: coarseb
    integer, optional,      intent(in)    :: method_p

    integer :: method, ist
    R_TYPE, allocatable :: f_coarse(:), f_fine(:)

    PUSH_SUB(X(multigrid_fine2coarse_batch))

    ASSERT(coarseb%nst_linear == fineb%nst_linear)

    method = optional_default(method_p, FULLWEIGHT)

    select case(method)
    case(FULLWEIGHT)
     call X(multigrid_restriction_batch)(tt, fine_der, coarse_mesh, fineb, coarseb)
    case(INJECTION)
      SAFE_ALLOCATE(f_coarse(1:coarse_mesh%np))
      SAFE_ALLOCATE(f_fine(1:fine_der%mesh%np_part))
      do ist = 1, coarseb%nst_linear
        call batch_get_state(coarseb, ist, coarse_mesh%np, f_coarse)
        call batch_get_state(fineb, ist, fine_der%mesh%np, f_fine)
        call X(multigrid_injection)(tt, f_fine, f_coarse)
        call batch_set_state(coarseb, ist, coarse_mesh%np, f_coarse)
      end do
      SAFE_DEALLOCATE_A(f_coarse)
      SAFE_DEALLOCATE_A(f_fine)
    case default
      write(message(1), '(a,i2,a)') 'Multigrid: Restriction method  = ', method, ' is not valid.'
      call messages_fatal(1)
    end select
    

    POP_SUB(X(multigrid_fine2coarse_batch))
  end subroutine X(multigrid_fine2coarse_batch)

  ! ---------------------------------------------------------
  subroutine X(multigrid_restriction_batch)(tt, fine_der, coarse_mesh, fineb, coarseb)
    type(transfer_table_t), intent(in)    :: tt
    type(derivatives_t),    intent(in)    :: fine_der
    type(mesh_t),           intent(in)    :: coarse_mesh
    class(batch_t),         intent(inout) :: fineb
    class(batch_t),         intent(inout) :: coarseb

    FLOAT :: weight(-1:1,-1:1,-1:1)
    integer :: nn, fn, di, dj, dk, dd, fi(3)
    R_TYPE, allocatable :: f_coarse(:), f_fine(:)
    integer :: ist
    type(pv_handle_batch_t) :: handle

    PUSH_SUB(X(multigrid_restriction_batch))
    call profiling_in(restrict_prof, TOSTRING(X(MG_RESTRICTION_BATCH)))

    !The following code only works in 3D
    ASSERT(fine_der%dim == 3)

    do di = -1, 1
      do dj = -1, 1
        do dk = -1, 1
          dd = abs(di) + abs(dj) + abs(dk)
          weight(di, dj, dk) = CNST(0.5)**dd
        end do
      end do
    end do

    call boundaries_set(fine_der%boundaries, fineb)
    if(fine_der%mesh%parallel_in_domains) then
      call X(ghost_update_batch_start)(fine_der%mesh%vp, fineb, handle)
      call X(ghost_update_batch_finish)(handle)
    end if

    ASSERT(fineb%status() == coarseb%status())

    select case(fineb%status())
    case(BATCH_PACKED)
      do nn = 1, tt%n_coarse
        fn = tt%to_coarse(nn)
        call mesh_local_index_to_coords(fine_der%mesh, fn, fi)
        do ist = 1, coarseb%nst_linear
          coarseb%X(ff_pack)(ist, nn) = M_ZERO
        end do

        do di = -1, 1
          do dj = -1, 1
            do dk = -1, 1
              fn = mesh_local_index_from_coords(fine_der%mesh, [fi(1) + di, fi(2) + dj, fi(3) + dk])

              if(fine_der%mesh%use_curvilinear) then
                f_coarse(nn) = f_coarse(nn) + weight(di, dj, dk)*f_fine(fn)*fine_der%mesh%vol_pp(fn)
                do ist = 1, coarseb%nst_linear
                  coarseb%X(ff_pack)(ist, nn) = coarseb%X(ff_pack)(ist, nn) &
                           + weight(di, dj, dk)*fine_der%mesh%vol_pp(fn)*fineb%X(ff_pack)(ist, fn)
                end do
              else
                do ist = 1, coarseb%nst_linear
                  coarseb%X(ff_pack)(ist, nn) = coarseb%X(ff_pack)(ist, nn) &
                           + weight(di, dj, dk)*fineb%X(ff_pack)(ist, fn)
                end do
              end if

            end do
          end do
        end do

        if(fine_der%mesh%use_curvilinear) then
          do ist = 1, coarseb%nst_linear
            coarseb%X(ff_pack)(ist, nn) = coarseb%X(ff_pack)(ist, nn)/coarse_mesh%vol_pp(nn)
          end do
        else
          do ist = 1, coarseb%nst_linear
            coarseb%X(ff_pack)(ist, nn) = coarseb%X(ff_pack)(ist, nn)*fine_der%mesh%vol_pp(1)/coarse_mesh%vol_pp(1)
          end do
        end if
      end do

    case(BATCH_DEVICE_PACKED, BATCH_NOT_PACKED)            
      SAFE_ALLOCATE(f_coarse(1:coarse_mesh%np))
      SAFE_ALLOCATE(f_fine(1:fine_der%mesh%np_part))
      do ist = 1, coarseb%nst_linear
        call batch_get_state(coarseb, ist, coarse_mesh%np, f_coarse)
        call batch_get_state(fineb, ist, fine_der%mesh%np_part, f_fine)            

        call X(multigrid_restriction)(tt, fine_der, coarse_mesh, f_fine, f_coarse, set_bc = .false.)

        call batch_set_state(coarseb, ist, coarse_mesh%np, f_coarse)
      end do
      SAFE_DEALLOCATE_A(f_coarse)
      SAFE_DEALLOCATE_A(f_fine)
    end select

    call profiling_count_operations(tt%n_coarse*(27*3 + 1))
    call profiling_out(restrict_prof)
    POP_SUB(X(multigrid_restriction_batch))
  end subroutine X(multigrid_restriction_batch)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
