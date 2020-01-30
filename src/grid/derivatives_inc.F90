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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

!> This module calculates the derivatives (gradients, Laplacians, etc.)
!! of a function. Note that the function whose derivative is to be calculated
!! *has* to be defined (1:mesh%np_part), while the (1:mesh%np) values of the derivative
!! are calculated. This was made to simplify the parallel mode, and has to be
!! followed by the rest of the code.

! ---------------------------------------------------------
!> These are the workhorse routines that handle the calculation of derivatives
subroutine X(derivatives_batch_start)(op, der, ff, opff, handle, ghost_update, set_bc, factor)
  type(nl_operator_t),       target, intent(in)    :: op
  type(derivatives_t),       target, intent(in)    :: der
  class(batch_t),            target, intent(inout) :: ff
  class(batch_t),            target, intent(inout) :: opff
  type(derivatives_handle_batch_t),  intent(out)   :: handle
  logical,                 optional, intent(in)    :: ghost_update
  logical,                 optional, intent(in)    :: set_bc
  FLOAT,                   optional, intent(in)    :: factor

  PUSH_SUB(X(derivatives_batch_start))

  handle%ghost_update = optional_default(ghost_update, .true.)

  handle%op   => op
  handle%der  => der
  handle%ff   => ff
  handle%opff => opff

  if(present(factor)) then
    handle%factor_present = .true.
    handle%factor = factor
  else
    handle%factor_present = .false.
  end if

  ASSERT(handle%ff%nst_linear == handle%opff%nst_linear)

  if(optional_default(set_bc, .true.)) then
    call boundaries_set(der%boundaries, ff)
    ASSERT(.not. der%boundaries%spiral)
  end if

#ifdef HAVE_MPI

  if(derivatives_overlap(der) .and. der%mesh%parallel_in_domains .and. handle%ghost_update) then
    call X(ghost_update_batch_start)(der%mesh%vp, ff, handle%pv_h)
  end if
#endif

  POP_SUB(X(derivatives_batch_start))
end subroutine X(derivatives_batch_start)


! ---------------------------------------------------------
subroutine X(derivatives_batch_finish)(handle)
  type(derivatives_handle_batch_t), intent(inout) :: handle

  logical :: done

  PUSH_SUB(X(derivatives_batch_finish))

  done = .false.

#ifdef HAVE_MPI
  if(derivatives_overlap(handle%der) .and. handle%der%mesh%parallel_in_domains .and. handle%ghost_update) then

    done = .true.

    if(handle%factor_present) then
      call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, &
        ghost_update = .false., points = OP_INNER, factor = handle%factor)
    else
      call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update=.false., points=OP_INNER)
    end if

    call X(ghost_update_batch_finish)(handle%pv_h)

    if(handle%factor_present) then
      call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, &
        ghost_update = .false., points = OP_OUTER, factor = handle%factor)
    else
      call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update = .false., points = OP_OUTER)
    end if

  end if
#endif

  if(.not. done) then
    if(handle%factor_present) then
      call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update = .false., factor = handle%factor)
    else
      call X(nl_operator_operate_batch)(handle%op, handle%ff, handle%opff, ghost_update = .false.)
    end if
  end if

  POP_SUB(X(derivatives_batch_finish))
end subroutine X(derivatives_batch_finish)


! ---------------------------------------------------------
subroutine X(derivatives_batch_perform)(op, der, ff, opff, ghost_update, set_bc, factor)
  type(nl_operator_t), intent(in)    :: op
  type(derivatives_t), intent(in)    :: der
  class(batch_t),      intent(inout) :: ff
  class(batch_t),      intent(inout) :: opff
  logical,   optional, intent(in)    :: ghost_update
  logical,   optional, intent(in)    :: set_bc
  FLOAT,     optional, intent(in)    :: factor

  type(derivatives_handle_batch_t) :: handle

  PUSH_SUB(X(derivatives_batch_perform))

  call ff%check_compatibility_with(opff)

  call X(derivatives_batch_start)(op, der, ff, opff, handle, ghost_update, set_bc, factor)
  call X(derivatives_batch_finish)(handle)

  POP_SUB(X(derivatives_batch_perform))

end subroutine X(derivatives_batch_perform)


! ---------------------------------------------------------
!> Now the simplified interfaces
subroutine X(derivatives_perform)(op, der, ff, op_ff, ghost_update, set_bc, factor)
  type(nl_operator_t), target, intent(in)    :: op
  type(derivatives_t),         intent(in)    :: der
  R_TYPE,                      intent(inout) :: ff(:)     !< (der%mesh%np_part)
  R_TYPE,                      intent(inout) :: op_ff(:)  !< (>= der%mesh%np)
  logical, optional,           intent(in)    :: ghost_update
  logical, optional,           intent(in)    :: set_bc
  FLOAT,   optional,           intent(in)    :: factor

  type(batch_t) :: batch_ff, batch_op_ff

  PUSH_SUB(X(derivatives_perform))

  ASSERT(ubound(ff, DIM=1) >= der%mesh%np_part)

  call batch_init(batch_ff, ff)
  call batch_init(batch_op_ff, op_ff)

  ASSERT(batch_ff%is_ok())
  ASSERT(batch_op_ff%is_ok())

  call X(derivatives_batch_perform) (op, der, batch_ff, batch_op_ff, ghost_update, set_bc, factor)

  call batch_ff%end()
  call batch_op_ff%end()

  POP_SUB(X(derivatives_perform))

end subroutine X(derivatives_perform)


! ---------------------------------------------------------
subroutine X(derivatives_lapl)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t),       intent(in)    :: der
  R_TYPE,                    intent(inout) :: ff(:)     !< (der%mesh%np_part)
  R_TYPE,                    intent(out)   :: op_ff(:)  !< (der%mesh%np)
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  PUSH_SUB(X(derivatives_lapl))

  call X(derivatives_perform)(der%lapl, der, ff, op_ff, ghost_update, set_bc)

  POP_SUB(X(derivatives_lapl))
end subroutine X(derivatives_lapl)


! ---------------------------------------------------------
subroutine X(derivatives_grad)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: ff(:)        !< ff(der%mesh%np_part)
  R_TYPE,              intent(out)   :: op_ff(:, :)  !< op_ff(der%mesh%np, der%mesh%sb%dim)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  integer :: idir,ip
  logical :: set_bc_, ghost_update_

  PUSH_SUB(X(derivatives_grad))
  call profiling_in(gradient_prof, "GRADIENT")

  ASSERT(ubound(op_ff, DIM=2) >= der%dim)

  set_bc_ = optional_default(set_bc, .true.)
  ghost_update_ = optional_default(ghost_update, .true.)

  do idir = 1, der%dim
    call X(derivatives_perform) (der%grad(idir), der, ff, op_ff(:, idir), ghost_update_, set_bc_)

    set_bc_       = .false. ! there is no need to update again
    ghost_update_ = .false. ! the boundary or ghost points
  end do

  ! Grad_xyw = Bt Grad_uvw, see Chelikowsky after Eq. 10
  if (simul_box_is_periodic(der%mesh%sb) .and. der%mesh%sb%nonorthogonal ) then
    forall (ip = 1:der%mesh%np)
      op_ff(ip, 1:der%dim) = matmul(der%mesh%sb%klattice_primitive(1:der%dim, 1:der%dim),op_ff(ip, 1:der%dim))
    end forall
  end if

  call profiling_out(gradient_prof)
  POP_SUB(X(derivatives_grad))
end subroutine X(derivatives_grad)


! ---------------------------------------------------------
subroutine X(derivatives_div)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,  target,     intent(inout) :: ff(:,:)   !< ff(der%mesh%np_part, der%mesh%sb%dim)
  R_TYPE,              intent(out)   :: op_ff(:)  !< op_ff(der%mesh%np)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  R_TYPE, allocatable :: tmp(:)
  R_TYPE, pointer     :: ff_uvw(:,:)
  integer             :: idir, ii, ip

  PUSH_SUB(X(derivatives_div))
  call profiling_in(divergence_prof, "DIVERGENCE")

  ASSERT(ubound(ff, DIM=2) >= der%dim)

  ! div_xyw (F)= div_uvw (BF), where B
  if (simul_box_is_periodic(der%mesh%sb) .and. der%mesh%sb%nonorthogonal ) then
    SAFE_ALLOCATE(ff_uvw(1:der%mesh%np_part, 1:der%dim))
    forall (ip = 1:der%mesh%np_part)
      ff_uvw(ip, 1:der%dim) = matmul(transpose(der%mesh%sb%klattice_primitive(1:der%dim, 1:der%dim)),ff(ip, 1:der%dim))
    end forall
  else
    ff_uvw => ff
  end if

  call X(derivatives_perform) (der%grad(1), der, ff_uvw(:, 1), op_ff, ghost_update, set_bc)

  SAFE_ALLOCATE(tmp(1:der%mesh%np))

  do idir = 2, der%dim
    call X(derivatives_perform) (der%grad(idir), der, ff_uvw(:, idir), tmp, ghost_update, set_bc)

    forall(ii = 1:der%mesh%np) op_ff(ii) = op_ff(ii) + tmp(ii)
  end do

  SAFE_DEALLOCATE_A(tmp)
  if(simul_box_is_periodic(der%mesh%sb) .and. der%mesh%sb%nonorthogonal ) then
    SAFE_DEALLOCATE_P(ff_uvw)
  else
    nullify(ff_uvw)
  end if

  call profiling_out(divergence_prof)
  POP_SUB(X(derivatives_div))
end subroutine X(derivatives_div)



! ---------------------------------------------------------
subroutine X(derivatives_curl)(der, ff, op_ff, ghost_update, set_bc)
  type(derivatives_t), intent(in)    :: der
  R_TYPE,              intent(inout) :: ff(:,:)    !< ff(der%mesh%np_part, der%dim)
  R_TYPE,              intent(out)   :: op_ff(:,:)
    !< if dim = 2, op_ff(der%mesh%np, der%dim)
    !! if dim = 1, op_ff(der%mesh%np, 1)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  integer, parameter  :: curl_dim(3) = (/-1, 1, 3/)
  R_TYPE, allocatable :: tmp(:)
  integer             :: ii, np

  PUSH_SUB(X(derivatives_curl))
  call profiling_in(curl_prof, "CURL")

  ASSERT(der%dim==2 .or. der%dim==3)
  ASSERT(ubound(ff,    DIM=2) >= der%dim)
  ASSERT(ubound(op_ff, DIM=2) >= curl_dim(der%dim))

  ASSERT(.not.der%mesh%sb%nonorthogonal)

  SAFE_ALLOCATE(tmp(1:der%mesh%np_part))

  op_ff(:,:) = R_TOTYPE(M_ZERO)
  np = der%mesh%np

  select case(der%dim)
  case(3)
    call X(derivatives_perform) (der%grad(3), der, ff(:,1), tmp, ghost_update, set_bc)
    forall(ii = 1:np) op_ff(ii, 2) = op_ff(ii, 2) + tmp(ii)
    call X(derivatives_perform) (der%grad(2), der, ff(:,1), tmp, .false., .false.)
    forall(ii = 1:np) op_ff(ii, 3) = op_ff(ii, 3) - tmp(ii)

    call X(derivatives_perform) (der%grad(3), der, ff(:,2), tmp, ghost_update, set_bc)
    forall(ii = 1:np) op_ff(ii, 1) = op_ff(ii, 1) - tmp(ii)
    call X(derivatives_perform) (der%grad(1), der, ff(:,2), tmp, .false., .false.)
    forall(ii = 1:np) op_ff(ii, 3) = op_ff(ii, 3) + tmp(ii)

    call X(derivatives_perform) (der%grad(2), der, ff(:,3), tmp, ghost_update, set_bc)
    forall(ii = 1:np) op_ff(ii, 1) = op_ff(ii, 1) + tmp(ii)
    call X(derivatives_perform) (der%grad(1), der, ff(:,3), tmp, .false., .false.)
    forall(ii = 1:np) op_ff(ii, 2) = op_ff(ii, 2) - tmp(ii)

  case(2)
    call X(derivatives_perform) (der%grad(2), der, ff(:,1), tmp, ghost_update, set_bc)
    forall(ii = 1:np) op_ff(ii, 1) = op_ff(ii, 1) - tmp(ii)
    call X(derivatives_perform) (der%grad(1), der, ff(:,2), tmp, .false., .false.)
    forall(ii = 1:np) op_ff(ii, 1) = op_ff(ii, 1) + tmp(ii)
  end select

  SAFE_DEALLOCATE_A(tmp)
  call profiling_out(curl_prof)
  POP_SUB(X(derivatives_curl))
end subroutine X(derivatives_curl)


! ----------------------------------------------------------

subroutine X(derivatives_test)(this, namespace, repetitions, min_blocksize, max_blocksize)
  type(derivatives_t), intent(in) :: this
  type(namespace_t),   intent(in) :: namespace
  integer,             intent(in) :: repetitions
  integer,             intent(in) :: min_blocksize
  integer,             intent(in) :: max_blocksize

  R_TYPE, allocatable :: ff(:), opff(:, :), gradff(:, :)
  R_TYPE :: aa, bb, cc
  integer :: ip, idir, ist
  type(batch_t) :: ffb, opffb
  type(batch_t), allocatable :: gradffb(:)
  integer :: blocksize, itime
  logical :: packstates
  real(8) :: stime, etime
  character(len=20) :: type

  call parse_variable(namespace, 'StatesPack', .true., packstates)

  SAFE_ALLOCATE(ff(1:this%mesh%np_part))
  SAFE_ALLOCATE(opff(1:this%mesh%np, 1:this%mesh%sb%dim))
  SAFE_ALLOCATE(gradff(1:this%mesh%np, 1:this%mesh%sb%dim))

#ifdef R_TREAL
#ifdef SINGLE_PRECISION
  type = 'real single'
#else
  type = 'real'
#endif
#else
#ifdef SINGLE_PRECISION
  type = 'complex single'
#else
  type = 'complex'
#endif
#endif

  ! Note: here we need to use a constant function or anything that
  ! is constant at the borders, since we assume that all boundary
  ! points have equal values to optimize the application of the nl-operator.

  aa = R_TOTYPE(1.0)/this%mesh%sb%lsize(1)
  bb = R_TOTYPE(10.0)
  cc = R_TOTYPE(100.0)

#ifdef R_TCOMPLEX
  ! we make things more "complex"
  aa = aa + M_ZI*R_TOTYPE(0.01)
  bb = bb*exp(M_ZI*R_TOTYPE(0.345))
  cc = cc - M_ZI*R_TOTYPE(50.0)
#endif


  forall(ip = 1:this%mesh%np_part) ff(ip) = bb*exp(-aa*sum(this%mesh%x(ip, 1:this%mesh%sb%dim)**2)) + cc
  do ip = 1, this%mesh%np
    do idir = 1, this%mesh%sb%dim
      gradff(ip, idir) = -M_TWO*aa*bb*this%mesh%x(ip, idir)*exp(-aa*sum(this%mesh%x(ip, :)**2))
    end do
  end do

  message(1) = "Testing Laplacian"
  call messages_info(1)

  ! test Laplacian
  blocksize = min_blocksize
  do

    call batch_init(ffb, 1, blocksize)
    call ffb%X(allocate)(1, blocksize, this%mesh%np_part)

    call batch_init(opffb, 1, blocksize)
    call opffb%X(allocate)(1, blocksize, this%mesh%np)

    forall(ist = 1:blocksize, ip = 1:this%mesh%np_part)
      ffb%X(ff_linear)(ip, ist) = ff(ip)
    end forall

    if(packstates) then
      call ffb%do_pack()
      call opffb%do_pack(copy = .false.)
    end if

    if(repetitions > 1) then
      call X(derivatives_batch_perform)(this%lapl, this, ffb, opffb, set_bc = .false., factor = CNST(0.5))
    end if

    stime = loct_clock()
    do itime = 1, repetitions
      call X(derivatives_batch_perform)(this%lapl, this, ffb, opffb, set_bc = .false., factor = CNST(0.5))
    end do
    etime = (loct_clock() - stime)/dble(repetitions)

    if(packstates) then
      call opffb%do_unpack()
    end if

    forall(ip = 1:this%mesh%np)
      opffb%X(ff_linear)(ip, blocksize) = CNST(2.0)*opffb%X(ff_linear)(ip, blocksize) - &
        (M_FOUR*aa**2*bb*sum(this%mesh%x(ip, :)**2)*exp(-aa*sum(this%mesh%x(ip, :)**2)) &
        - this%mesh%sb%dim*M_TWO*aa*bb*exp(-aa*sum(this%mesh%x(ip, :)**2)))
    end forall

    write(message(1), '(3a,i3,a,es17.10,a,f8.3)') &
      'Laplacian ', trim(type),  &
      ' bsize = ', blocksize,    &
      ' , error = ', X(mf_nrm2)(this%mesh, opffb%X(ff_linear)(:, blocksize)), &
      ' , Gflops = ',  &
#ifdef R_TREAL
      blocksize*this%mesh%np*CNST(2.0)*this%lapl%stencil%size/(etime*CNST(1.0e9))
#else
      blocksize*this%mesh%np*CNST(4.0)*this%lapl%stencil%size/(etime*CNST(1.0e9))
#endif

    call messages_info(1)

    call ffb%end()
    call opffb%end()

    blocksize = 2*blocksize
    if(blocksize > max_blocksize) exit

  end do

  message(1) = "Testing Gradient"
  call messages_info(1)

  ! test Gradient
  SAFE_ALLOCATE(gradffb(1:this%mesh%sb%dim))
  blocksize = min_blocksize
  do
    call batch_init(ffb, 1, blocksize)
    call ffb%X(allocate)(1, blocksize, this%mesh%np_part)

    do ist = 1, blocksize
      do ip = 1,this%mesh%np_part
        ffb%X(ff_linear)(ip, ist) = ff(ip)
      end do
    end do

    if(packstates) then
      call ffb%do_pack()
    end if

    do idir = 1, this%mesh%sb%dim
      call ffb%copy_to(gradffb(idir))
    end do

    if(repetitions > 1) then
      call X(derivatives_batch_grad)(this, ffb, gradffb, set_bc=.false.)
    end if

    stime = loct_clock()
    do itime = 1, repetitions
      call X(derivatives_batch_grad)(this, ffb, gradffb, set_bc=.false.)
    end do
    etime = (loct_clock() - stime)/dble(repetitions)

    if(packstates) then
      do idir = 1, this%mesh%sb%dim
        call gradffb(idir)%do_unpack()
      end do
    end if

    do ip = 1, this%mesh%np
      do idir = 1, this%mesh%sb%dim
        gradffb(idir)%X(ff_linear)(ip, blocksize) = &
          gradffb(idir)%X(ff_linear)(ip, blocksize) - gradff(ip, idir)
      end do
    end do

    write(message(1), '(3a,i3,a,es17.10,a,f8.3)') &
      'Batch gradient ', trim(type),  &
      ' bsize = ', blocksize,    &
      ' , error = ', X(mf_nrm2)(this%mesh, gradffb(1)%X(ff_linear)(:, blocksize)), &
      ' , Gflops = ',  &
#ifdef R_TREAL
      blocksize*this%mesh%np*CNST(2.0)*this%grad(1)%stencil%size*this%mesh%sb%dim/(etime*CNST(1.0e9))
#else
      blocksize*this%mesh%np*CNST(4.0)*this%grad(1)%stencil%size*this%mesh%sb%dim/(etime*CNST(1.0e9))
#endif
    call messages_info(1)

    call ffb%end()
    do idir = 1, this%mesh%sb%dim
      call gradffb(idir)%end()
    end do

    blocksize = 2*blocksize
    if(blocksize > max_blocksize) exit
  end do
  SAFE_DEALLOCATE_A(gradffb)

  call X(derivatives_grad)(this, ff, opff, set_bc = .false.)

  forall(idir = 1:this%mesh%sb%dim, ip = 1:this%mesh%np)
    opff(ip, idir) = opff(ip, idir) - (-M_TWO*aa*bb*this%mesh%x(ip, idir)*exp(-aa*sum(this%mesh%x(ip, :)**2)))
  end forall

  message(1) = ''
  call messages_info(1)


  write(message(1), '(3a, es17.10)') 'Gradient ', trim(type),  &
    ' err = ', X(mf_nrm2)(this%mesh, this%mesh%sb%dim, opff)
  call messages_info(1)

  message(1) = ''
  call messages_info(1)

  SAFE_DEALLOCATE_A(ff)
  SAFE_DEALLOCATE_A(opff)
  SAFE_DEALLOCATE_A(gradff)

end subroutine X(derivatives_test)

! ---------------------------------------------------------
subroutine X(derivatives_partial)(der, ff, op_ff, dir, ghost_update, set_bc)
  type(derivatives_t),       intent(in)   :: der
  R_TYPE,                    intent(inout) :: ff(:)     !< (der%mesh%np_part)
  R_TYPE,                    intent(out)   :: op_ff(:)  !< (der%mesh%np)
  integer,                   intent(in)    :: dir
  logical, optional,         intent(in)    :: ghost_update
  logical, optional,         intent(in)    :: set_bc

  logical :: set_bc_, ghost_update_

  PUSH_SUB(X(derivatives_partial))

  set_bc_ = optional_default(set_bc, .true.)
  ghost_update_ = optional_default(ghost_update, .true.)

  call X(derivatives_perform) (der%grad(dir), der, ff, op_ff, ghost_update_, set_bc_)

  POP_SUB(X(derivatives_partial))
end subroutine X(derivatives_partial)

subroutine X(derivatives_batch_grad)(der, ffb, opffb, ghost_update, set_bc)
  type(derivatives_t), intent(in)    :: der
  class(batch_t),      intent(inout) :: ffb
  class(batch_t),      intent(inout) :: opffb(:)
  logical, optional,   intent(in)    :: ghost_update
  logical, optional,   intent(in)    :: set_bc

  integer :: idir
  logical :: set_bc_, ghost_update_

  PUSH_SUB(X(derivatives_batch_grad))
  call profiling_in(batch_gradient_prof, "BATCH_GRADIENT")

  set_bc_ = optional_default(set_bc, .true.)
  ghost_update_ = optional_default(ghost_update, .true.)

  ASSERT(size(opffb) == der%mesh%sb%dim)

  do idir = 1, der%mesh%sb%dim
    call X(derivatives_batch_perform)(der%grad(idir), der, ffb, opffb(idir), &
      ghost_update=ghost_update_, set_bc=set_bc_)
    set_bc_       = .false. ! there is no need to update again
    ghost_update_ = .false. ! the boundary or ghost points
  end do

  if (simul_box_is_periodic(der%mesh%sb) .and. der%mesh%sb%nonorthogonal ) then
    call X(batch_vector_uvw_to_xyz)(der, opffb)
  end if

  call profiling_out(batch_gradient_prof)
  POP_SUB(X(derivatives_batch_grad))
end subroutine X(derivatives_batch_grad)

subroutine X(batch_vector_uvw_to_xyz)(der, uvw, xyz)
  type(derivatives_t),              intent(in)    :: der
  class(batch_t), target,           intent(inout) :: uvw(:)
  class(batch_t), target, optional, intent(inout) :: xyz(:)

  class(batch_t), pointer:: xyz_(:)
  integer :: ist, ip, idim1, idim2
  R_TYPE, allocatable :: tmp(:)

  PUSH_SUB(X(batch_vector_uvw_to_xyz))
  if(present(xyz)) then
    xyz_ => xyz
  else
    xyz_ => uvw
  end if
  ASSERT(size(uvw) == der%dim)
  ASSERT(size(xyz_) == der%dim)
  do idim1 = 1, der%dim
    ASSERT(uvw(idim1)%status() == xyz_(idim1)%status())
    ASSERT(uvw(idim1)%status() == uvw(1)%status())
  end do

  SAFE_ALLOCATE(tmp(der%dim))

  select case(uvw(1)%status())
  case(BATCH_NOT_PACKED)
    do ist = 1, uvw(1)%nst_linear
      !$omp parallel do private(tmp, ip, idim1, idim2)
      do ip = 1, der%mesh%np
        tmp = R_TOTYPE(M_ZERO)
        ! Grad_xyw = Bt Grad_uvw, see Chelikowsky after Eq. 10
        do idim2 = 1, der%dim
          do idim1 = 1, der%dim
            tmp(idim1) = tmp(idim1) + &
              der%mesh%sb%klattice_primitive(idim1, idim2) * uvw(idim2)%X(ff_linear)(ip, ist)
          end do
        end do
        do idim1 = 1, der%dim
          xyz_(idim1)%X(ff_linear)(ip, ist) = tmp(idim1)
        end do
      end do
    end do
  case(BATCH_PACKED)
    !$omp parallel do private(tmp, ip, ist, idim1, idim2)
    do ip = 1, der%mesh%np
      do ist = 1, uvw(1)%nst_linear
        tmp = R_TOTYPE(M_ZERO)
        ! Grad_xyw = Bt Grad_uvw, see Chelikowsky after Eq. 10
        do idim2 = 1, der%dim
          do idim1 = 1, der%dim
            tmp(idim1) = tmp(idim1) + &
              der%mesh%sb%klattice_primitive(idim1, idim2) * uvw(idim2)%pack%X(psi)(ist, ip)
          end do
        end do
        do idim1 = 1, der%dim
          xyz_(idim1)%pack%X(psi)(ist, ip) = tmp(idim1)
        end do
      end do
    end do
  case(BATCH_DEVICE_PACKED)
    call uvw_to_xyz_opencl()
  end select

  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(X(batch_vector_uvw_to_xyz))

contains
  subroutine uvw_to_xyz_opencl
    integer :: localsize, dim3, dim2
    type(accel_mem_t) :: matrix_buffer

    PUSH_SUB(uvw_to_xyz_opencl)

    call accel_create_buffer(matrix_buffer, ACCEL_MEM_READ_ONLY, TYPE_FLOAT, der%dim**2)
    call accel_write_buffer(matrix_buffer, der%dim**2, der%mesh%sb%klattice_primitive)

    call accel_set_kernel_arg(kernel_uvw_xyz, 0, der%mesh%np)
    call accel_set_kernel_arg(kernel_uvw_xyz, 1, matrix_buffer)
    call accel_set_kernel_arg(kernel_uvw_xyz, 2, uvw(1)%pack%buffer)
    call accel_set_kernel_arg(kernel_uvw_xyz, 3, log2(uvw(1)%pack%size_real(1)))
    call accel_set_kernel_arg(kernel_uvw_xyz, 4, xyz_(1)%pack%buffer)
    call accel_set_kernel_arg(kernel_uvw_xyz, 5, log2(xyz_(1)%pack%size_real(1)))
    if(der%dim > 1) then
      call accel_set_kernel_arg(kernel_uvw_xyz, 6, uvw(2)%pack%buffer)
      call accel_set_kernel_arg(kernel_uvw_xyz, 7, log2(uvw(2)%pack%size_real(1)))
      call accel_set_kernel_arg(kernel_uvw_xyz, 8, xyz_(2)%pack%buffer)
      call accel_set_kernel_arg(kernel_uvw_xyz, 9, log2(xyz_(2)%pack%size_real(1)))
    end if
    if(der%dim > 2) then
      call accel_set_kernel_arg(kernel_uvw_xyz, 10, uvw(3)%pack%buffer)
      call accel_set_kernel_arg(kernel_uvw_xyz, 11, log2(uvw(3)%pack%size_real(1)))
      call accel_set_kernel_arg(kernel_uvw_xyz, 12, xyz_(3)%pack%buffer)
      call accel_set_kernel_arg(kernel_uvw_xyz, 13, log2(xyz_(3)%pack%size_real(1)))
    end if

    localsize = accel_kernel_workgroup_size(kernel_uvw_xyz)/uvw(1)%pack%size_real(1)

    dim3 = der%mesh%np/(accel_max_size_per_dim(2)*localsize) + 1
    dim2 = min(accel_max_size_per_dim(2)*localsize, pad(der%mesh%np, localsize))

    call accel_kernel_run(kernel_uvw_xyz, (/uvw(1)%pack%size_real(1), dim2, dim3/), &
      (/uvw(1)%pack%size_real(1), localsize, 1/))
    call accel_finish()

    call accel_release_buffer(matrix_buffer)

    POP_SUB(uvw_to_xyz_opencl)
  end subroutine uvw_to_xyz_opencl
end subroutine X(batch_vector_uvw_to_xyz)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
