#include "global.h"

module storage_m

#if defined(HAVE_MPI)
  use boundaries_m
#endif
  use global_m
  use grid_m
  use kinds_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use multigrid_m
  use profiling_m
  use simulation_m

  implicit none

  private

  public ::    &
    storage_t

  public ::             &
    storage_init,       &
    storage_start,      &
    storage_update,     &
    storage_stop,       &
    storage_reset,      &
    storage_integrate,  &
    storage_mlt,        &
    storage_add,        &
    storage_sub,        &
    storage_reduce,     &
    storage_get,        &
    storage_transfer,   &
    storage_copy,       &
    storage_end

  type :: storage_t
    private
    type(simulation_t),                pointer :: sim     =>null()
    type(grid_t),                      pointer :: grid    =>null()
    type(mesh_t),                      pointer :: mesh    =>null()
    logical                                    :: alloc   = .true.
    logical                                    :: full    = .true.
    logical                                    :: fine    = .false.
    integer                                    :: ndim    = 0
    integer                                    :: size    = 0
    real(kind=wp)                              :: default = 0.0_wp
    real(kind=wp), dimension(:,:), allocatable :: data
  end type storage_t

  interface storage_init
    module procedure storage_init_simple
    module procedure storage_init_copy
  end interface storage_init

  interface storage_mlt
    module procedure storage_mlt_dim
    module procedure storage_mlt_all
  end interface storage_mlt

  interface storage_integrate
    module procedure storage_integrate_intg_dim
    module procedure storage_integrate_intg_all
    module procedure storage_integrate_dotp
  end interface storage_integrate

  interface storage_get
    module procedure storage_get_info
    module procedure storage_get_sim
    module procedure storage_get_grid
    module procedure storage_get_mesh
    module procedure storage_get_storage_1d
    module procedure storage_get_storage_md
  end interface storage_get

contains

  ! ---------------------------------------------------------
  subroutine storage_init_simple(this, ndim, default, full, allocate)
    type(storage_t),         intent(out) :: this
    integer,       optional, intent(in)  :: ndim
    real(kind=wp), optional, intent(in)  :: default
    logical,       optional, intent(in)  :: full
    logical,       optional, intent(in)  :: allocate

    PUSH_SUB(storage_init_simple)

    nullify(this%sim, this%grid, this%mesh)
    this%full = .true.
    if(present(full)) this%full = full
    this%alloc = .true.
    if(present(allocate)) this%alloc = allocate
    this%ndim = 1
    if(present(ndim))then
      ASSERT(ndim>0)
      this%ndim = ndim
    end if
    this%size = 0
    this%default = 0.0_wp
    if(present(default)) this%default = default

    POP_SUB(storage_init_simple)
  end subroutine storage_init_simple

  ! ---------------------------------------------------------
  subroutine storage_init_copy(this, that)
    type(storage_t), intent(out) :: this
    type(storage_t), intent(in)  :: that

    PUSH_SUB(storage_init_copy)

    ASSERT(that%ndim>0)
    call storage_init_simple(this, that%ndim, that%default, that%full, that%alloc)

    POP_SUB(storage_init_copy)
  end subroutine storage_init_copy

  ! ---------------------------------------------------------
  subroutine storage_start(this, sim, fine)
    type(storage_t),            intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    logical,          optional, intent(in)    :: fine

    PUSH_SUB(storage_start)

    ASSERT(.not.associated(this%sim))
    ASSERT(.not.associated(this%grid))
    ASSERT(.not.associated(this%mesh))
    ASSERT(this%ndim>0)
    ASSERT(.not.allocated(this%data))
    this%sim => sim
    call simulation_get(this%sim, this%grid)
    ASSERT(associated(this%grid))
    this%fine = .false.
    if(this%grid%have_fine_mesh)then
      if(present(fine)) this%fine = fine
    end if
    call simulation_get(this%sim, this%mesh, this%fine)
    ASSERT(associated(this%mesh))
    ASSERT(this%mesh%np>0)
    ASSERT(this%mesh%np_part>0)
    ASSERT(this%mesh%np<this%mesh%np_part)
    if(this%full)then
      this%size = this%mesh%np_part
    else
      this%size = this%mesh%np
    end if
    if(this%alloc)then
      SAFE_ALLOCATE(this%data(this%size,this%ndim))
    end if
    call storage_reset(this)

    POP_SUB(storage_start)
  end subroutine storage_start

  ! ---------------------------------------------------------
  subroutine storage_update(this)
    type(storage_t), intent(inout) :: this

    integer :: indx

    PUSH_SUB(storage_update)

    ASSERT(associated(this%sim))
    ASSERT(associated(this%mesh))
    ASSERT(this%ndim>0)
    if(this%alloc)then
      if(this%full)then
        do indx = this%mesh%np+1, this%size
          this%data(indx,:) = 0.0_wp
        end do
      end if
#if defined(HAVE_MPI)
      do indx = 1, this%ndim
        call dvec_ghost_update(this%mesh%vp, this%data(:,indx))
      end do
#endif
    end if

    POP_SUB(storage_update)
  end subroutine storage_update

  ! ---------------------------------------------------------
  subroutine storage_stop(this)
    type(storage_t), intent(inout) :: this

    PUSH_SUB(storage_stop)

    ASSERT(associated(this%sim))
    ASSERT(associated(this%mesh))
    ASSERT(this%ndim>0)
    ASSERT(this%size>0)
    if(this%alloc)then
      SAFE_DEALLOCATE_A(this%data)
    end if
    this%size = 0

    POP_SUB(storage_stop)
  end subroutine storage_stop

  ! ---------------------------------------------------------
  subroutine storage_reset(this, default)
    type(storage_t),         intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: default

    real(kind=wp) :: dflt

    PUSH_SUB(storage_reset)

    ASSERT(associated(this%sim))
    ASSERT(associated(this%mesh))
    ASSERT(this%ndim>0)
    ASSERT(this%size>0)
    if(this%alloc)then
      dflt = this%default
      if(present(default)) dflt = default
      this%data(1:this%mesh%np,:) = dflt
      call storage_update(this)
    end if

    POP_SUB(storage_reset)
  end subroutine storage_reset

  ! ---------------------------------------------------------
  subroutine storage_integrate_intg_dim(this, dim, value)
    type(storage_t), intent(in)  :: this
    integer,         intent(in)  :: dim
    real(kind=wp),   intent(out) :: value

    PUSH_SUB(storage_integrate_intg_dim)

    ASSERT(associated(this%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(this%mesh))
    ASSERT(dim>0)
    ASSERT(dim<=this%ndim)
    ASSERT(this%size>0)
    value = 0.0_wp
    if(this%alloc)then
      value = dmf_integrate(this%mesh, this%data(:,dim))
    end if

    POP_SUB(storage_integrate_intg_dim)
  end subroutine storage_integrate_intg_dim

  ! ---------------------------------------------------------
  subroutine storage_integrate_intg_all(this, value)
    type(storage_t), intent(in)  :: this
    real(kind=wp),   intent(out) :: value

    integer :: indx

    PUSH_SUB(storage_integrate_intg_all)

    ASSERT(associated(this%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(this%mesh))
    ASSERT(this%ndim>0)
    ASSERT(this%size>0)
    value = 0.0_wp
    if(this%alloc)then
      do indx = 1, this%ndim
        value = value + dmf_integrate(this%mesh, this%data(:,indx))
      end do
    end if

    POP_SUB(storage_integrate_intg_all)
  end subroutine storage_integrate_intg_all

  ! ---------------------------------------------------------
  subroutine storage_integrate_dotp_aux(this, that, value)
    type(storage_t), intent(in)  :: this
    type(storage_t), intent(in)  :: that
    real(kind=wp),   intent(out) :: value

    integer :: indx

    PUSH_SUB(storage_integrate_dotp_aux)

    ASSERT(associated(this%mesh,that%mesh))
    value = 0.0_wp
    do indx = 1, this%ndim
      value = value + dmf_dotp(this%mesh, this%data(:,indx), that%data(:,indx))
    end do

    POP_SUB(storage_integrate_dotp_aux)
  end subroutine storage_integrate_dotp_aux

  ! ---------------------------------------------------------
  subroutine storage_integrate_dotp(this, that, value)
    type(storage_t), target, intent(in)  :: this
    type(storage_t), target, intent(in)  :: that
    real(kind=wp),           intent(out) :: value

    type(storage_t), pointer :: cdat, fdat
    type(storage_t)          :: data

    PUSH_SUB(storage_integrate_dotp)

    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(that%grid))
    ASSERT(associated(this%grid,that%grid))
    ASSERT(associated(this%mesh))
    ASSERT(associated(that%mesh))
    ASSERT(this%ndim>0)
    ASSERT(that%ndim>0)
    ASSERT(this%ndim==that%ndim)
    ASSERT(this%size>0)
    ASSERT(that%size>0)
    value = 0.0_wp
    nullify(cdat, fdat)
    if(this%alloc.and.that%alloc)then
      if(this%fine.eqv.that%fine)then
        call storage_integrate_dotp_aux(this, that, value)
      else
        if(this%fine)then
          cdat => that
          fdat => this
        else
          cdat => this
          fdat => that
        end if
        call storage_init(data, cdat)
        call storage_start(data, cdat%sim)
        call storage_transfer(data, fdat)
        call storage_integrate_dotp_aux(cdat, data, value)
        call storage_end(data)
        nullify(cdat, fdat)
      end if
    end if

    POP_SUB(storage_integrate_dotp)
  end subroutine storage_integrate_dotp

  ! ---------------------------------------------------------
  subroutine storage_mlt_dim(this, dim, mlt)
    type(storage_t), intent(inout) :: this
    integer,         intent(in)    :: dim
    real(kind=wp),   intent(in)    :: mlt

    integer :: indx

    PUSH_SUB(storage_mlt_dim)

    ASSERT(associated(this%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(this%mesh))
    ASSERT(this%ndim>0)
    ASSERT(this%size>0)
    if(this%alloc)then
      do indx = 1, this%mesh%np
        this%data(indx,dim) = mlt * this%data(indx,dim)
      end do
      call storage_update(this)
    end if

    POP_SUB(storage_mlt_dim)
  end subroutine storage_mlt_dim

  ! ---------------------------------------------------------
  subroutine storage_mlt_all(this, mlt)
    type(storage_t), intent(inout) :: this
    real(kind=wp),   intent(in)    :: mlt

    integer :: indx

    PUSH_SUB(storage_mlt_all)

    ASSERT(associated(this%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(this%mesh))
    ASSERT(this%ndim>0)
    ASSERT(this%size>0)
    if(this%alloc)then
      do indx = 1, this%mesh%np
        this%data(indx,:) = mlt * this%data(indx,:)
      end do
      call storage_update(this)
    end if

    POP_SUB(storage_mlt_all)
  end subroutine storage_mlt_all

  ! ---------------------------------------------------------
  subroutine storage_add(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    integer :: indx

    PUSH_SUB(storage_add)

    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(that%grid))
    ASSERT(associated(this%mesh))
    ASSERT(associated(that%mesh))
    ASSERT(associated(this%mesh,that%mesh))
    ASSERT(this%ndim>0)
    ASSERT(that%ndim>0)
    ASSERT(this%ndim==that%ndim)
    ASSERT(this%size>0)
    ASSERT(that%size>0)
    ASSERT(this%size==that%size)
    ASSERT(this%default.equal.that%default)
    if(this%alloc.and.that%alloc)then
      do indx = 1, this%mesh%np
        this%data(indx,:) = this%data(indx,:) + that%data(indx,:)
      end do
      call storage_update(this)
    end if

    POP_SUB(storage_add)
  end subroutine storage_add

  ! ---------------------------------------------------------
  subroutine storage_sub(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    integer :: indx

    PUSH_SUB(storage_sub)

    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(that%grid))
    ASSERT(associated(this%mesh))
    ASSERT(associated(that%mesh))
    ASSERT(associated(this%mesh,that%mesh))
    ASSERT(this%ndim>0)
    ASSERT(that%ndim>0)
    ASSERT(this%ndim==that%ndim)
    ASSERT(this%size>0)
    ASSERT(that%size>0)
    ASSERT(this%size==that%size)
    ASSERT(this%default.equal.that%default)
    if(this%alloc.and.that%alloc)then
      do indx = 1, this%mesh%np
        this%data(indx,:) = this%data(indx,:) - that%data(indx,:)
      end do
      call storage_update(this)
    end if

    POP_SUB(storage_sub)
  end subroutine storage_sub

  ! ---------------------------------------------------------
  subroutine storage_reduce(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    integer :: indx

    PUSH_SUB(storage_reduce)

    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(that%grid))
    ASSERT(associated(this%mesh))
    ASSERT(associated(that%mesh))
    ASSERT(associated(this%mesh,that%mesh))
    ASSERT(this%ndim>0)
    ASSERT(that%ndim>0)
    ASSERT(this%ndim==1)
    ASSERT(this%ndim<=that%ndim)
    ASSERT(this%size>0)
    ASSERT(that%size>0)
    ASSERT(this%size==that%size)
    ASSERT(this%default.equal.that%default)
    if(this%alloc.and.that%alloc)then
      do indx = 1, this%mesh%np
        this%data(indx,1) = sum(that%data(indx,:))
      end do
      call storage_update(this)
    end if

    POP_SUB(storage_reduce)
  end subroutine storage_reduce

  ! ---------------------------------------------------------
  subroutine storage_get_info(this, dim, size, fine, alloc, default)
    type(storage_t), target, intent(in)  :: this
    integer,       optional, intent(out) :: dim
    integer,       optional, intent(out) :: size
    logical,       optional, intent(out) :: fine
    logical,       optional, intent(out) :: alloc
    real(kind=wp), optional, intent(out) :: default

    PUSH_SUB(storage_get_info)

    if(present(dim))then
      dim = 0
      if(this%ndim>0) dim = this%ndim
    end if
    if(present(size))then
      size = 0
      if(associated(this%mesh)) size = this%mesh%np
    end if
    if(present(fine)) fine = this%fine
    if(present(alloc)) alloc = this%alloc
    if(present(default)) default = this%default

    POP_SUB(storage_get_info)
  end subroutine storage_get_info

  ! ---------------------------------------------------------
  subroutine storage_get_sim(this, that)
    type(storage_t),     target, intent(in) :: this
    type(simulation_t), pointer             :: that

    PUSH_SUB(storage_get_sim)

    nullify(that)
    if(associated(this%sim)) that => this%sim

    POP_SUB(storage_get_sim)
  end subroutine storage_get_sim

  ! ---------------------------------------------------------
  subroutine storage_get_grid(this, that)
    type(storage_t), target, intent(in) :: this
    type(grid_t),   pointer             :: that

    PUSH_SUB(storage_get_grid)

    nullify(that)
    if(associated(this%grid)) that => this%grid

    POP_SUB(storage_get_grid)
  end subroutine storage_get_grid

  ! ---------------------------------------------------------
  subroutine storage_get_mesh(this, that)
    type(storage_t), target, intent(in) :: this
    type(mesh_t),   pointer             :: that

    PUSH_SUB(storage_get_mesh)

    nullify(that)
    if(associated(this%mesh)) that => this%mesh

    POP_SUB(storage_get_mesh)
  end subroutine storage_get_mesh

  ! ---------------------------------------------------------
  subroutine storage_get_storage_1d(this, that)
    type(storage_t),              target, intent(in) :: this
    real(kind=wp), dimension(:), pointer             :: that

    PUSH_SUB(storage_get_storage_1d)

    nullify(that)
    if(allocated(this%data))then
      ASSERT(this%alloc)
      ASSERT(this%ndim==1)
      ASSERT(this%size>1)
      that => this%data(:,1)
    end if

    POP_SUB(storage_get_storage_1d)
  end subroutine storage_get_storage_1d

  ! ---------------------------------------------------------
  subroutine storage_get_storage_md(this, that)
    type(storage_t),                target, intent(in) :: this
    real(kind=wp), dimension(:,:), pointer             :: that

    PUSH_SUB(storage_get_storage_md)

    nullify(that)
    if(allocated(this%data))then
      ASSERT(this%alloc)
      ASSERT(this%ndim>0)
      ASSERT(this%size>1)
      that => this%data
    end if

    POP_SUB(storage_get_storage_md)
  end subroutine storage_get_storage_md

  ! ---------------------------------------------------------
  subroutine storage_transfer(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    real(kind=wp), dimension(this%size) :: buff
    integer                             :: indx

    PUSH_SUB(storage_transfer)

    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%grid))
    ASSERT(associated(that%grid))
    ASSERT(associated(this%grid,that%grid))
    ASSERT(associated(this%mesh))
    ASSERT(associated(that%mesh))
    ASSERT(this%ndim>0)
    ASSERT(that%ndim>0)
    ASSERT(this%ndim==that%ndim)
    ASSERT(this%size>0)
    ASSERT(that%size>0)
    ASSERT(this%default.equal.that%default)
    ASSERT(this%alloc.and.that%alloc)
    if(this%fine.eqv.that%fine)then
      ASSERT(associated(this%mesh,that%mesh))
      do indx = 1, this%mesh%np
        this%data(indx,:) = that%data(indx,:)
      end do
    else
      if(this%fine)then
        do indx = 1, this%ndim
          buff = that%data(:,indx)
          call dmultigrid_coarse2fine(this%grid%fine%tt, this%grid%der, this%mesh, &
            buff, this%data(:,indx), order=2)
        end do
      else
        do indx = 1, this%ndim
          buff = that%data(:,indx)
          call dmultigrid_fine2coarse(this%grid%fine%tt, this%grid%fine%der, this%mesh, &
            buff, this%data(:,indx), INJECTION)
        end do
      end if
    end if
    call storage_update(this)

    POP_SUB(storage_transfer)
  end subroutine storage_transfer

  ! ---------------------------------------------------------
  subroutine storage_copy(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    PUSH_SUB(storage_copy)

    call storage_end(this)
    if(that%ndim>0)then
      call storage_init_simple(this, that%ndim, that%default, that%full, that%alloc)
      if(associated(that%sim).and.associated(that%mesh))then
        call storage_start(this, that%sim, that%fine)
        if(that%alloc)then
          this%data(1:this%mesh%np,:) = that%data(1:this%mesh%np,:)
          call storage_update(this)
        end if
      end if
    end if

    POP_SUB(storage_copy)
  end subroutine storage_copy

  ! ---------------------------------------------------------
  subroutine storage_end(this)
    type(storage_t), intent(inout) :: this

    PUSH_SUB(storage_end)

    nullify(this%sim, this%grid, this%mesh)
    this%alloc = .true.
    this%full = .true.
    this%fine = .false.
    this%ndim = 0
    this%size = 0
    this%default = 0.0_wp
    if(allocated(this%data))then
      SAFE_DEALLOCATE_A(this%data)
    end if

    POP_SUB(storage_end)
  end subroutine storage_end

end module storage_m

!! Local Variables:
!! mode: f90
!! End:
