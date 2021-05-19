!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade
!! Copyright (C) 2021 N. Tancogne-Dejean
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
  subroutine X(multigrid_coarse2fine)(tt, coarse_der, fine_mesh, f_coarse, f_fine)
    type(transfer_table_t),  intent(in)    :: tt
    type(derivatives_t),     intent(in)    :: coarse_der
    type(mesh_t),            intent(in)    :: fine_mesh
    R_TYPE,                  intent(inout) :: f_coarse(:)
    R_TYPE,                  intent(out)   :: f_fine(:)

    FLOAT, allocatable :: weight(:)
    integer, allocatable :: shift(:,:)
    integer :: di, nn, fn, fi(coarse_der%dim), fii(coarse_der%dim) 

    PUSH_SUB(X(multigrid_coarse2fine))

    call profiling_in(interp_prof, TOSTRING(X(MG_PROLONGATION)))

    ASSERT(ubound(f_fine, dim = 1) == fine_mesh%np_part)

    ! We get theiinterpolation weights for the N-linear interpolation
    SAFE_ALLOCATE(weight(3**coarse_der%dim))
    SAFE_ALLOCATE(shift(coarse_der%dim, 3**coarse_der%dim))
    call multigrid_build_stencil(coarse_der%dim, weight, shift)

    f_fine = M_ZERO

    do nn = 1, coarse_der%mesh%np
      call mesh_local_index_to_coords(coarse_der%mesh, nn, fi)

      do di = 1, 3**coarse_der%dim
        fii = 2*fi + shift(:,di)
        fn = mesh_local_index_from_coords(fine_mesh, fii)
        f_fine(fn) = f_fine(fn) + weight(di)*f_coarse(nn)
      end do
    end do

    SAFE_DEALLOCATE_A(weight)
    SAFE_DEALLOCATE_A(shift)

    call profiling_count_operations(coarse_der%mesh%np*3**coarse_der%dim)
    call profiling_out(interp_prof)
    POP_SUB(X(multigrid_coarse2fine))
  end subroutine X(multigrid_coarse2fine)


  ! ---------------------------------------------------------
  subroutine X(multigrid_fine2coarse)(tt, fine_der, coarse_mesh, f_fine, f_coarse, method_p, set_bc)
    type(transfer_table_t), intent(in)    :: tt
    type(derivatives_t),    intent(in)    :: fine_der
    type(mesh_t),           intent(in)    :: coarse_mesh
    R_TYPE,                 intent(inout) :: f_fine(:)
    R_TYPE,                 intent(out)   :: f_coarse(:)
    integer, optional,      intent(in)    :: method_p
    logical, optional,      intent(in)    :: set_bc

    integer :: method

    PUSH_SUB(X(multigrid_fine2coarse))

    if(present(method_p)) then
      method=method_p
    else
      method=FULLWEIGHT
    end if

    select case(method)
    case(FULLWEIGHT)
      call X(multigrid_restriction)(tt, fine_der, coarse_mesh, f_fine, f_coarse, set_bc)
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

    integer :: nn, fn, di, fi(fine_der%dim), fii(fine_der%dim)
    FLOAT, allocatable :: weight(:)
    integer, allocatable :: shift(:,:)

    PUSH_SUB(X(multigrid_restriction))
    call profiling_in(restrict_prof, TOSTRING(X(MG_RESTRICTION)))

    ! We get the interpolation weights for the N-linear interpolation
    SAFE_ALLOCATE(weight(3**fine_der%dim))
    SAFE_ALLOCATE(shift(fine_der%dim, 3**fine_der%dim))
    call multigrid_build_stencil(fine_der%dim, weight, shift)

    ! Normalization
    if(.not. fine_der%mesh%use_curvilinear) then
      weight = weight * fine_der%mesh%vol_pp(1)/coarse_mesh%vol_pp(1)
    end if

    ! We need to set the boundary conditions, as we use values of f_fine outside of the 
    ! np points due to the stencil
    if(optional_default(set_bc, .true.)) then
      call boundaries_set(fine_der%boundaries, f_fine)

#ifdef HAVE_MPI
      if(fine_der%mesh%parallel_in_domains) call X(vec_ghost_update)(fine_der%mesh%vp, f_fine)
#endif
    end if

    do nn = 1, coarse_mesh%np
      call mesh_local_index_to_coords(coarse_mesh, nn, fi)

      f_coarse(nn) = M_ZERO

      do di = 1, 3**fine_der%dim
        fii = 2*fi + shift(:,di)
        fn = mesh_local_index_from_coords(fine_der%mesh, fii)
        if(fine_der%mesh%use_curvilinear) then
          f_coarse(nn) = f_coarse(nn) + weight(di)*f_fine(fn)*fine_der%mesh%vol_pp(fn)
        else
          f_coarse(nn) = f_coarse(nn) + weight(di)*f_fine(fn)
        end if
      end do

      ! Normalization
      if(fine_der%mesh%use_curvilinear) then
        f_coarse(nn) = f_coarse(nn)/coarse_mesh%vol_pp(nn)
      end if
    end do

    SAFE_DEALLOCATE_A(weight)

    call profiling_count_operations(coarse_mesh%np*3**fine_der%dim)
    call profiling_out(restrict_prof)
    POP_SUB(X(multigrid_restriction))
  end subroutine X(multigrid_restriction)

  ! ---------------------------------------------------------
  subroutine X(multigrid_coarse2fine_batch)(tt, coarse_der, fine_mesh, coarseb, fineb)
    type(transfer_table_t),  intent(in)    :: tt
    type(derivatives_t),     intent(in)    :: coarse_der
    type(mesh_t),            intent(in)    :: fine_mesh
    class(batch_t),          intent(inout) :: coarseb
    class(batch_t),          intent(inout) :: fineb

    integer :: ist, nn, ii, ipf
    integer :: xf(coarse_der%dim), xc(coarse_der%dim)
    FLOAT, allocatable :: weight(:)
    integer, allocatable :: shift(:,:)
    R_TYPE, allocatable :: f_coarse(:), f_fine(:)

    PUSH_SUB(X(multigrid_coarse2fine_batch))

    call profiling_in(interp_prof, TOSTRING(X(MG_INTERPOLATION_BATCH)))

    ASSERT(coarseb%nst_linear == fineb%nst_linear)
    ASSERT(fineb%status() == coarseb%status())

    select case(fineb%status())
    case(BATCH_PACKED)

      ! We get the interpolation weights for the N-linear interpolation
      SAFE_ALLOCATE(weight(3**coarse_der%dim))
      SAFE_ALLOCATE(shift(coarse_der%dim, 3**coarse_der%dim))
      call multigrid_build_stencil(coarse_der%dim, weight, shift)

      call batch_set_zero(fineb)

      !We perform a trilinear interpolation, see https://en.wikipedia.org/wiki/Trilinear_interpolation
      do nn = 1, coarse_der%mesh%np
        call mesh_local_index_to_coords(coarse_der%mesh, nn, xf)

        do ii = 1, 3**coarse_der%dim
          xc = 2*xf + shift(:,ii)
          ipf = mesh_local_index_from_coords(fine_mesh, xc)
          do ist = 1, coarseb%nst_linear
            fineb%X(ff_pack)(ist, ipf) = fineb%X(ff_pack)(ist, ipf) + weight(ii) * coarseb%X(ff_pack)(ist, nn)
          end do
        end do
      end do

      SAFE_DEALLOCATE_A(weight)
      SAFE_DEALLOCATE_A(shift)

    case(BATCH_DEVICE_PACKED, BATCH_NOT_PACKED)
      SAFE_ALLOCATE(f_coarse(1:coarse_der%mesh%np))
      SAFE_ALLOCATE(f_fine(1:fine_mesh%np_part))
      do ist = 1, coarseb%nst_linear
        call batch_get_state(coarseb, ist, coarse_der%mesh%np, f_coarse)
        call batch_get_state(fineb, ist, fine_mesh%np, f_fine)

        call X(multigrid_coarse2fine)(tt, coarse_der, fine_mesh, f_coarse, f_fine) 

        call batch_set_state(fineb, ist, fine_mesh%np, f_fine)
      end do

      SAFE_DEALLOCATE_A(f_coarse)
      SAFE_DEALLOCATE_A(f_fine)

    end select

    call profiling_count_operations(coarse_der%mesh%np*3**coarse_der%dim)
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

    integer :: nn, fn, di, fi(fine_der%dim), fii(fine_der%dim)
    R_TYPE, allocatable :: f_coarse(:), f_fine(:)
    integer :: ist
    FLOAT, allocatable   :: weight(:)
    integer, allocatable :: shift(:,:)
    type(pv_handle_batch_t) :: handle

    PUSH_SUB(X(multigrid_restriction_batch))
    call profiling_in(restrict_prof, TOSTRING(X(MG_RESTRICTION_BATCH)))

    call boundaries_set(fine_der%boundaries, fineb)
    if(fine_der%mesh%parallel_in_domains) then
      call X(ghost_update_batch_start)(fine_der%mesh%vp, fineb, handle)
      call X(ghost_update_batch_finish)(handle)
    end if

    ASSERT(fineb%status() == coarseb%status())

    select case(fineb%status())
    case(BATCH_PACKED)
      ! We get the interpolation weights for the N-linear interpolation
      SAFE_ALLOCATE(weight(3**fine_der%dim))
      SAFE_ALLOCATE(shift(fine_der%dim, 3**fine_der%dim))
      call multigrid_build_stencil(fine_der%dim, weight, shift)

      ! Normalization
      if(.not. fine_der%mesh%use_curvilinear) then
        weight = weight * fine_der%mesh%vol_pp(1)/coarse_mesh%vol_pp(1)
      end if

      do nn = 1, coarse_mesh%np
        call mesh_local_index_to_coords(coarse_mesh, nn, fi)
  
        do ist = 1, coarseb%nst_linear
          coarseb%X(ff_pack)(ist, nn) = M_ZERO
        end do

        if(fine_der%mesh%use_curvilinear) then
          do di = 1, 3**fine_der%dim
            fii = 2*fi + shift(:,di)
            fn = mesh_local_index_from_coords(fine_der%mesh, fii)
            do ist = 1, coarseb%nst_linear
              coarseb%X(ff_pack)(ist, nn) = coarseb%X(ff_pack)(ist, nn) &
                     + weight(di)*fine_der%mesh%vol_pp(fn)*fineb%X(ff_pack)(ist, fn)
            end do
          end do

          ! Normalization
          do ist = 1, coarseb%nst_linear
            coarseb%X(ff_pack)(ist, nn) = coarseb%X(ff_pack)(ist, nn)/coarse_mesh%vol_pp(nn)
          end do

        else

          do di = 1, 3**fine_der%dim
            fii = 2*fi + shift(:,di)
            fn = mesh_local_index_from_coords(fine_der%mesh, fii)
            do ist = 1, coarseb%nst_linear
              coarseb%X(ff_pack)(ist, nn) = coarseb%X(ff_pack)(ist, nn) &
                    + weight(di)*fineb%X(ff_pack)(ist, fn)
            end do
          end do

        end if
  
      end do
 
    case(BATCH_DEVICE_PACKED, BATCH_NOT_PACKED)            
      SAFE_ALLOCATE(f_coarse(1:coarse_mesh%np))
      SAFE_ALLOCATE(f_fine(1:fine_der%mesh%np_part))
      do ist = 1, coarseb%nst_linear
        call batch_get_state(coarseb, ist, coarse_mesh%np, f_coarse)
        ! We set the boundary points on fineb, so we copy them
        call batch_get_state(fineb, ist, fine_der%mesh%np_part, f_fine)            

        call X(multigrid_restriction)(tt, fine_der, coarse_mesh, f_fine, f_coarse, set_bc = .false.)

        call batch_set_state(coarseb, ist, coarse_mesh%np, f_coarse)
      end do
      SAFE_DEALLOCATE_A(f_coarse)
      SAFE_DEALLOCATE_A(f_fine)
    end select

    call profiling_count_operations(coarse_mesh%np*3**fine_der%dim)
    call profiling_out(restrict_prof)
    POP_SUB(X(multigrid_restriction_batch))
  end subroutine X(multigrid_restriction_batch)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
