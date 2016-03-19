#include "global.h"

module storage_oct_m

#if defined(HAVE_MPI)
  use boundaries_oct_m
#endif
  use global_oct_m
  use grid_oct_m
  use kinds_oct_m
  use json_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use multigrid_oct_m
  use profiling_oct_m
  use simulation_oct_m

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
    type(json_object_t),               pointer :: config  =>null()
    type(simulation_t),                pointer :: sim     =>null()
    type(grid_t),                      pointer :: grid    =>null()
    type(mesh_t),                      pointer :: mesh    =>null()
    logical                                    :: fine    = .false.
    logical                                    :: alloc   = .false.
    logical                                    :: full    = .true.
    integer                                    :: ndim    = 0
    integer                                    :: size    = 0
    real(kind=wp)                              :: default = 0.0_wp
    real(kind=wp), dimension(:,:), allocatable :: data
  end type storage_t

  interface storage_init
    module procedure storage_init_config
    module procedure storage_init_type
    module procedure storage_init_copy
  end interface storage_init

  interface storage_add
    module procedure storage_add_dim
    module procedure storage_add_all
  end interface storage_add

  interface storage_mlt
    module procedure storage_mlt_dim
    module procedure storage_mlt_all
  end interface storage_mlt

  interface storage_reset
    module procedure storage_reset_dim
    module procedure storage_reset_all
  end interface storage_reset

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
  subroutine storage_init_config(this, ndim, default, fine, full, allocate)
    type(json_object_t),     intent(out) :: this
    integer,       optional, intent(in)  :: ndim
    real(kind=wp), optional, intent(in)  :: default
    logical,       optional, intent(in)  :: fine
    logical,       optional, intent(in)  :: full
    logical,       optional, intent(in)  :: allocate

    PUSH_SUB(storage_init_config)

    call json_init(this)
    if(present(ndim)) call json_set(this, "dimensions", ndim)
    if(present(default)) call json_set(this, "default", default)
    if(present(fine)) call json_set(this, "fine", fine)
    if(present(full)) call json_set(this, "full", full)
    if(present(allocate)) call json_set(this, "allocate", allocate)

    POP_SUB(storage_init_config)
  end subroutine storage_init_config

  ! ---------------------------------------------------------
  subroutine storage_init_type(this, config)
    type(storage_t),             intent(out) :: this
    type(json_object_t), target, intent(in)  :: config

    integer :: ierr

    PUSH_SUB(storage_init_type)

    this%config => config
    call json_get(this%config, "fine", this%fine, ierr)
    if(ierr/=JSON_OK) this%fine = .true.
    call json_get(this%config, "full", this%full, ierr)
    if(ierr/=JSON_OK) this%full = .true.
    call json_get(this%config, "default", this%default, ierr)
    if(ierr/=JSON_OK) this%default = 0.0_wp

    POP_SUB(storage_init_type)
  end subroutine storage_init_type

  ! ---------------------------------------------------------
  subroutine storage_init_copy(this, that, ndim, fine)
    type(storage_t),   intent(out) :: this
    type(storage_t),   intent(in)  :: that
    integer, optional, intent(in)  :: ndim
    logical, optional, intent(in)  :: fine

    integer :: nd
    logical :: fn

    PUSH_SUB(storage_init_copy)

    ASSERT(associated(that%config))
    call storage_init(this, that%config)
    nd = that%ndim
    if(present(ndim)) nd = ndim
    fn = that%fine
    if(present(fine)) fn = fine
    if(associated(that%sim)) call storage_start(this, that%sim, ndim=nd, fine=fn)

    POP_SUB(storage_init_copy)
  end subroutine storage_init_copy

  ! ---------------------------------------------------------
  subroutine storage_start(this, sim, ndim, fine)
    type(storage_t),            intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    integer,          optional, intent(in)    :: ndim
    logical,          optional, intent(in)    :: fine

    integer :: ierr

    PUSH_SUB(storage_start)

    ASSERT(associated(this%config))
    ASSERT(.not.associated(this%sim))
    this%sim => sim
    call simulation_get(this%sim, this%grid)
    ASSERT(associated(this%grid))
    if(present(fine)) this%fine = fine
    if(.not.this%grid%have_fine_mesh) this%fine = .false.
    call simulation_get(this%sim, this%mesh, this%fine)
    ASSERT(associated(this%mesh))
    ASSERT(this%mesh%np>0)
    ASSERT(this%mesh%np_part>0)
    ASSERT(this%mesh%np<this%mesh%np_part)
    call json_get(this%config, "allocate", this%alloc, ierr)
    if(ierr/=JSON_OK) this%alloc = .true.
    this%size = this%mesh%np
    if(this%full) this%size = this%mesh%np_part
    call json_get(this%config, "dimensions", this%ndim, ierr)
    if(ierr/=JSON_OK) this%ndim = 1
    if(present(ndim)) this%ndim = ndim
    ASSERT(this%ndim>0)
    ASSERT(.not.allocated(this%data))
    call storage_alloc(this)
    call storage_reset(this)

    POP_SUB(storage_start)
  end subroutine storage_start

  ! ---------------------------------------------------------
  subroutine storage_update_aux(this, dim)
    type(storage_t), intent(inout) :: this
    integer,         intent(in)    :: dim

    integer :: indx

    PUSH_SUB(storage_update_aux)

    do indx = this%mesh%np+1, this%size
      this%data(indx,dim) = 0.0_wp
    end do
#if defined(HAVE_MPI)
    if(this%mesh%parallel_in_domains) call dvec_ghost_update(this%mesh%vp, this%data(:,dim))
#endif

    POP_SUB(storage_update_aux)
  end subroutine storage_update_aux

  ! ---------------------------------------------------------
  subroutine storage_update(this)
    type(storage_t), intent(inout) :: this

    integer :: indx

    PUSH_SUB(storage_update)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(this%alloc.and.this%full)then
      do indx = 1, this%ndim
        call storage_update_aux(this, indx)
      end do
    end if

    POP_SUB(storage_update)
  end subroutine storage_update

  ! ---------------------------------------------------------
  subroutine storage_stop(this)
    type(storage_t), intent(inout) :: this

    integer :: ierr

    PUSH_SUB(storage_stop)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_dealloc(this)
    nullify(this%sim, this%grid, this%mesh)
    call json_get(this%config, "fine", this%fine, ierr)
    if(ierr/=JSON_OK) this%fine = .true.
    this%alloc = .false.
    this%size = 0
    this%ndim = 0

    POP_SUB(storage_stop)
  end subroutine storage_stop

  ! ---------------------------------------------------------
  subroutine storage_alloc(this)
    type(storage_t), intent(inout) :: this

    PUSH_SUB(storage_alloc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    call storage_dealloc(this)
    if(this%alloc)then
      SAFE_ALLOCATE(this%data(1:this%size,1:this%ndim))
    end if

    POP_SUB(storage_alloc)
  end subroutine storage_alloc

  ! ---------------------------------------------------------
  subroutine storage_dealloc(this)
    type(storage_t), intent(inout) :: this

    PUSH_SUB(storage_dealloc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(allocated(this%data))then
      SAFE_DEALLOCATE_A(this%data)
    end if

    POP_SUB(storage_dealloc)
  end subroutine storage_dealloc

  ! ---------------------------------------------------------
  subroutine storage_realloc(this)
    type(storage_t), intent(inout) :: this

    PUSH_SUB(storage_realloc)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(this%alloc)then
      if(allocated(this%data))then
        if((size(this%data,dim=1)/=this%size).or.&
          (size(this%data,dim=2)/=this%ndim)) call storage_alloc(this)
      else
        call storage_alloc(this)
      end if
    else
      call storage_dealloc(this)
    end if

    POP_SUB(storage_realloc)
  end subroutine storage_realloc

  ! ---------------------------------------------------------
  subroutine storage_reset_aux(this, dim, default)
    type(storage_t),         intent(inout) :: this
    integer,                 intent(in)    :: dim
    real(kind=wp), optional, intent(in)    :: default

    real(kind=wp) :: dflt

    PUSH_SUB(storage_reset_aux)

    dflt = this%default
    if(present(default)) dflt = default
    this%data(1:this%mesh%np,dim) = dflt
    if(this%full) call storage_update_aux(this, dim)

    POP_SUB(storage_reset_aux)
  end subroutine storage_reset_aux

  ! ---------------------------------------------------------
  subroutine storage_reset_dim(this, dim, default)
    type(storage_t),         intent(inout) :: this
    integer,                 intent(in)    :: dim
    real(kind=wp), optional, intent(in)    :: default

    PUSH_SUB(storage_reset_dim)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(this%alloc) call storage_reset_aux(this, dim, default)

    POP_SUB(storage_reset_dim)
  end subroutine storage_reset_dim

  ! ---------------------------------------------------------
  subroutine storage_reset_all(this, default)
    type(storage_t),         intent(inout) :: this
    real(kind=wp), optional, intent(in)    :: default

    integer :: indx

    PUSH_SUB(storage_reset_all)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(this%alloc)then
      do indx = 1, this%ndim
        call storage_reset_aux(this, indx, default)
      end do
    end if

    POP_SUB(storage_reset_all)
  end subroutine storage_reset_all

  ! ---------------------------------------------------------
  subroutine storage_integrate_intg_dim(this, dim, value)
    type(storage_t), intent(in)  :: this
    integer,         intent(in)  :: dim
    real(kind=wp),   intent(out) :: value

    PUSH_SUB(storage_integrate_intg_dim)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    ASSERT(dim<=this%ndim)
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

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
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

    value = 0.0_wp
    do indx = 1, this%ndim
      value = value + dmf_dotp(this%mesh, this%data(:,indx), that%data(:,indx))
    end do

    POP_SUB(storage_integrate_dotp_aux)
  end subroutine storage_integrate_dotp_aux

  ! ---------------------------------------------------------
  recursive subroutine storage_integrate_dotp_dim(this, that, value)
    type(storage_t), intent(in)  :: this
    type(storage_t), intent(in)  :: that
    real(kind=wp),   intent(out) :: value

    type(storage_t) :: data

    PUSH_SUB(storage_integrate_dotp_dim)

    value = 0.0_wp
    if(this%ndim==that%ndim)then
      call storage_integrate_dotp_aux(this, that, value)
    else
      if(this%ndim>1)then
        call storage_init(data, this, ndim=1)
        call storage_reduce_aux(data, this)
        call storage_integrate_dotp_dim(that, data, value)
        call storage_end(data)
      else
        call storage_integrate_dotp_dim(that, this, value)
      end if
    end if

    POP_SUB(storage_integrate_dotp_dim)
  end subroutine storage_integrate_dotp_dim

  ! ---------------------------------------------------------
  subroutine storage_integrate_dotp_f2c(this, that, value)
    type(storage_t), intent(in)  :: this
    type(storage_t), intent(in)  :: that
    real(kind=wp),   intent(out) :: value

    type(storage_t) :: data

    PUSH_SUB(storage_integrate_dotp_f2c)

    value = 0.0_wp
    call storage_init(data, that, fine=.false.)
    call storage_transfer_f2c(data, that)
    call storage_integrate_dotp_dim(this, data, value)
    call storage_end(data)

    POP_SUB(storage_integrate_dotp_f2c)
  end subroutine storage_integrate_dotp_f2c

  ! ---------------------------------------------------------
  subroutine storage_integrate_dotp(this, that, value)
    type(storage_t), intent(in)  :: this
    type(storage_t), intent(in)  :: that
    real(kind=wp),   intent(out) :: value

    PUSH_SUB(storage_integrate_dotp)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%grid,that%grid))
    value = 0.0_wp
    if(this%alloc.and.that%alloc)then
      if(this%fine.eqv.that%fine)then
        call storage_integrate_dotp_dim(this, that, value)
      else
        if(this%fine)then
          call storage_integrate_dotp_f2c(that, this, value)
        else
          call storage_integrate_dotp_f2c(this, that, value)
        end if
      end if
    end if

    POP_SUB(storage_integrate_dotp)
  end subroutine storage_integrate_dotp

  ! ---------------------------------------------------------
  subroutine storage_mlt_aux(this, dim, mlt)
    type(storage_t), intent(inout) :: this
    integer,         intent(in)    :: dim
    real(kind=wp),   intent(in)    :: mlt

    integer :: indx

    PUSH_SUB(storage_mlt_aux)

    do indx = 1, this%mesh%np
      this%data(indx,dim) = mlt * this%data(indx,dim)
    end do
    if(this%full) call storage_update_aux(this, dim)

    POP_SUB(storage_mlt_aux)
  end subroutine storage_mlt_aux

  ! ---------------------------------------------------------
  subroutine storage_mlt_dim(this, dim, mlt)
    type(storage_t), intent(inout) :: this
    integer,         intent(in)    :: dim
    real(kind=wp),   intent(in)    :: mlt

    PUSH_SUB(storage_mlt_dim)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(this%alloc)then
      call storage_mlt_aux(this, dim, mlt)
    end if

    POP_SUB(storage_mlt_dim)
  end subroutine storage_mlt_dim

  ! ---------------------------------------------------------
  subroutine storage_mlt_all(this, mlt)
    type(storage_t), intent(inout) :: this
    real(kind=wp),   intent(in)    :: mlt

    integer :: indx

    PUSH_SUB(storage_mlt_all)

    ASSERT(associated(this%config))
    ASSERT(associated(this%sim))
    if(this%alloc)then
      do indx = 1, this%ndim
        call storage_mlt_aux(this, indx, mlt)
      end do
    end if

    POP_SUB(storage_mlt_all)
  end subroutine storage_mlt_all

  ! ---------------------------------------------------------
  subroutine storage_add_aux(this, that, dim)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that
    integer,         intent(in)    :: dim

    integer :: indx

    PUSH_SUB(storage_add_aux)

    do indx = 1, this%mesh%np
      this%data(indx,dim) = this%data(indx,dim) + that%data(indx,dim)
    end do
    if(this%full) call storage_update_aux(this, dim)

    POP_SUB(storage_add_aux)
  end subroutine storage_add_aux

  ! ---------------------------------------------------------
  subroutine storage_add_dim(this, that, dim)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that
    integer,         intent(in)    :: dim

    PUSH_SUB(storage_add_dim)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%mesh,that%mesh))
    ASSERT(this%ndim==that%ndim)
    ASSERT(this%size==that%size)
    ASSERT(this%default.equal.that%default)
    if(this%alloc.and.that%alloc) call storage_add_aux(this, that, dim)

    POP_SUB(storage_add_dim)
  end subroutine storage_add_dim

  ! ---------------------------------------------------------
  subroutine storage_add_all(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    integer :: indx

    PUSH_SUB(storage_add_all)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%mesh,that%mesh))
    ASSERT(this%ndim==that%ndim)
    ASSERT(this%size==that%size)
    ASSERT(this%default.equal.that%default)
    if(this%alloc.and.that%alloc)then
      do indx = 1, this%ndim
        call storage_add_aux(this, that, indx)
      end do
    end if

    POP_SUB(storage_add_all)
  end subroutine storage_add_all

  ! ---------------------------------------------------------
  subroutine storage_sub(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    integer :: indx, jndx

    PUSH_SUB(storage_sub)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%mesh,that%mesh))
    ASSERT(this%ndim==that%ndim)
    ASSERT(this%size==that%size)
    ASSERT(this%default.equal.that%default)
    if(this%alloc.and.that%alloc)then
      do jndx = 1, this%ndim
        do indx = 1, this%mesh%np
          this%data(indx,jndx) = this%data(indx,jndx) - that%data(indx,jndx)
        end do
        if(this%full) call storage_update_aux(this, jndx)
      end do
    end if

    POP_SUB(storage_sub)
  end subroutine storage_sub

  ! ---------------------------------------------------------
  subroutine storage_reduce_aux(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    integer :: indx, jndx

    PUSH_SUB(storage_reduce_aux)

    do indx = 1, this%mesh%np
      this%data(indx,1) = that%data(indx,1)
    end do
    do jndx = 2, that%ndim
      do indx = 1, this%mesh%np
        this%data(indx,1) = this%data(indx,1) + that%data(indx,jndx)
      end do
    end do
    if(this%full) call storage_update_aux(this, 1)

    POP_SUB(storage_reduce_aux)
  end subroutine storage_reduce_aux

  ! ---------------------------------------------------------
  subroutine storage_reduce(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    PUSH_SUB(storage_reduce)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%mesh,that%mesh))
    ASSERT(this%ndim==1)
    ASSERT(this%ndim<=that%ndim)
    ASSERT(this%size==that%size)
    ASSERT(this%default.equal.that%default)
    if(this%alloc.and.that%alloc)then
      call storage_reduce_aux(this, that)
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
    if(present(alloc)) alloc = this%alloc .and. associated(this%sim)
    if(present(default)) default = this%default

    POP_SUB(storage_get_info)
  end subroutine storage_get_info

  ! ---------------------------------------------------------
  subroutine storage_get_config(this, that)
    type(storage_t),      target, intent(in) :: this
    type(json_object_t), pointer             :: that

    PUSH_SUB(storage_get_config)

    nullify(that)
    if(associated(this%config)) that => this%config

    POP_SUB(storage_get_config)
  end subroutine storage_get_config

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
    if(this%alloc)then
      ASSERT(allocated(this%data))
      ASSERT(this%ndim==1)
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
    if(this%alloc)then
      ASSERT(allocated(this%data))
      that => this%data
    end if

    POP_SUB(storage_get_storage_md)
  end subroutine storage_get_storage_md

  ! ---------------------------------------------------------
  subroutine storage_transfer_c2f(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    real(kind=wp), dimension(that%size) :: buff
    integer                             :: indx, jndx

    PUSH_SUB(storage_transfer_c2f)

    do indx = 1, that%ndim
      do jndx = 1, that%mesh%np
        buff(jndx) = that%data(jndx,indx)
      end do
      do jndx = that%mesh%np+1, that%size
        buff(jndx) = 0.0_wp
      end do
      call dmultigrid_coarse2fine(this%grid%fine%tt, this%grid%der, this%mesh, &
        buff(:), this%data(:,indx), order=2)
      if(this%full) call storage_update_aux(this, indx)
    end do

    POP_SUB(storage_transfer_c2f)
  end subroutine storage_transfer_c2f

  ! ---------------------------------------------------------
  subroutine storage_transfer_f2c(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    real(kind=wp), dimension(that%size) :: buff
    integer                             :: indx, jndx

    PUSH_SUB(storage_transfer_f2c)

    do indx = 1, this%ndim
      do jndx = 1, that%mesh%np
        buff(jndx) = that%data(jndx,indx)
      end do
      do jndx = that%mesh%np+1, that%size
        buff(jndx) = 0.0_wp
      end do
      call dmultigrid_fine2coarse(this%grid%fine%tt, this%grid%fine%der, this%mesh, &
        buff(:), this%data(:,indx), INJECTION)
      if(this%full) call storage_update_aux(this, indx)
    end do

    POP_SUB(storage_transfer_f2c)
  end subroutine storage_transfer_f2c

  ! ---------------------------------------------------------
  subroutine storage_transfer(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    PUSH_SUB(storage_transfer)

    ASSERT(associated(this%config))
    ASSERT(associated(that%config))
    ASSERT(associated(this%sim))
    ASSERT(associated(that%sim))
    ASSERT(associated(this%grid,that%grid))
    ASSERT(this%ndim==that%ndim)
    ASSERT(this%default.equal.that%default)
    ASSERT(this%alloc.and.that%alloc)
    if(this%fine.eqv.that%fine)then
      ASSERT(associated(this%mesh,that%mesh))
      call storage_copy_aux(this, that)
    else
      if(this%fine)then
        call storage_transfer_c2f(this, that)
      else
        call storage_transfer_f2c(this, that)
      end if
    end if

    POP_SUB(storage_transfer)
  end subroutine storage_transfer

  ! ---------------------------------------------------------
  subroutine storage_copy_aux(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    integer :: indx, jndx

    PUSH_SUB(storage_copy_aux)

    do jndx = 1, that%ndim
      do indx = 1, that%mesh%np
        this%data(indx,jndx) = that%data(indx,jndx)
      end do
      if(this%full) call storage_update_aux(this, jndx)
    end do

    POP_SUB(storage_copy_aux)
  end subroutine storage_copy_aux

  ! ---------------------------------------------------------
  subroutine storage_copy(this, that)
    type(storage_t), intent(inout) :: this
    type(storage_t), intent(in)    :: that

    PUSH_SUB(storage_copy)

    this%config => that%config
    this%sim => that%sim
    this%grid => that%grid
    this%mesh => that%mesh
    this%fine = that%fine
    this%alloc = that%alloc
    this%full = that%full
    this%ndim = that%ndim
    this%size = that%size
    this%default = that%default
    if(associated(this%config).and.associated(this%sim))then
      call storage_realloc(this)
      if(this%alloc) call storage_copy_aux(this, that)
    end if

    POP_SUB(storage_copy)
  end subroutine storage_copy

  ! ---------------------------------------------------------
  subroutine storage_end(this)
    type(storage_t), intent(inout) :: this

    PUSH_SUB(storage_end)

    nullify(this%config, this%sim, this%grid, this%mesh)
    this%fine = .false.
    this%alloc = .false.
    this%full = .true.
    this%ndim = 0
    this%size = 0
    this%default = 0.0_wp
    SAFE_DEALLOCATE_A(this%data)

    POP_SUB(storage_end)
  end subroutine storage_end

end module storage_oct_m

!! Local Variables:
!! mode: f90
!! End:
