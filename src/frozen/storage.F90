#include "global.h"

module storage_m

  use global_m
  use messages_m
  use profiling_m

#if defined(HAVE_MPI)
  use boundaries_m, only: dvec_ghost_update
#endif
  use grid_m,       only: grid_t
  use kinds_m,      only: wp
  use mesh_m,       only: mesh_t

  use basis_m, only: &
    basis_t

  use simulation_m, only: &
    simulation_t,         &
    simulation_get

  use intrpl_m, only: &
    intrpl_init

  use intrpl_m, only: &
    intrpl_t,         &
    intrpl_eval,      &
    intrpl_set,       &
    intrpl_copy,      &
    intrpl_end

  use intrpl_m, only:               &
    STORAGE_INTRPL_OK => INTRPL_OK, &
    STORAGE_INTRPL_OD => INTRPL_OD, &
    STORAGE_INTRPL_NI => INTRPL_NI

  implicit none

  private
!!$  public ::      &
!!$    operator(+), &
!!$    operator(-)

  public ::                &
    storage_init,          &
    storage_start,         &
    storage_update,        &
    storage_stop,          &
    storage_eval,          &
    !storage_conform,       &
    storage_set,           &
    storage_get,           &
    storage_get_size,      &
    storage_get_dimension, &
    storage_get_storage,   &
    storage_copy,          &
    storage_end

  public ::            &
    STORAGE_INTRPL_OK, &
    STORAGE_INTRPL_OD, &
    STORAGE_INTRPL_NI

  type, public :: storage_t
    private
    type(simulation_t),                pointer :: sim     =>null()
    type(mesh_t),                      pointer :: mesh    =>null()
    integer                                    :: ndim    = 0
    real(kind=wp)                              :: default = 0.0_wp
    real(kind=wp), dimension(:,:), allocatable :: data
  end type storage_t

  type, public :: storage_intrpl_t
    private
    type(storage_t), pointer :: self =>null()
    type(intrpl_t)           :: intrp
  end type storage_intrpl_t

!!$  interface operator(+)
!!$    module procedure storage_add
!!$  end interface operator(+)
!!$
!!$  interface operator(-)
!!$    module procedure storage_sub
!!$  end interface operator(-)

  interface storage_init
    module procedure storage_init_simple
    module procedure storage_init_copy
    module procedure storage_intrpl_init
  end interface storage_init

  interface storage_eval
    module procedure storage_intrpl_eval_1d
    module procedure storage_intrpl_eval_md
  end interface storage_eval

  interface storage_set
    module procedure storage_intrpl_set
  end interface storage_set

  interface storage_get
    module procedure storage_get_sim
    !module procedure storage_get_storage
  end interface storage_get

  interface storage_get_storage
    module procedure storage_get_storage_1d
    module procedure storage_get_storage_md
  end interface storage_get_storage

  interface storage_copy
    module procedure storage_copy_storage
    module procedure storage_intrpl_copy
  end interface storage_copy

  interface storage_end
    module procedure storage_end_storage
    module procedure storage_intrpl_end
  end interface storage_end

contains

  ! ---------------------------------------------------------
  subroutine storage_init_simple(this, ndim, default)
    type(storage_t),         intent(out) :: this
    integer,       optional, intent(in)  :: ndim
    real(kind=wp), optional, intent(in)  :: default
    !
    integer :: ierr
    !
    PUSH_SUB(storage_init_simple)
    this%ndim=1
    if(present(ndim))then
      ASSERT(ndim>0)
      this%ndim=ndim
    end if
    this%default=0.0_wp
    if(present(default))this%default=default
    POP_SUB(storage_init_simple)
    return
  end subroutine storage_init_simple

  ! ---------------------------------------------------------
  subroutine storage_init_copy(this, that)
    type(storage_t),         intent(out) :: this
    type(storage_t), target, intent(in)  :: that
    !
    PUSH_SUB(storage_init_copy)
    ASSERT(that%ndim>0)
    ASSERT(associated(that%sim))
    call storage_init_simple(this, that%ndim, that%default)
    call storage_start(this, that%sim)
    POP_SUB(storage_init_copy)
    return
  end subroutine storage_init_copy

  ! ---------------------------------------------------------
  subroutine storage_start(this, sim)
    type(storage_t),                      intent(inout) :: this
    type(simulation_t), optional, target, intent(in)    :: sim
    !
    type(grid_t), pointer :: grid
    !
    PUSH_SUB(storage_start)
    nullify(grid)
    if(present(sim))then
      ASSERT(.not.associated(this%sim))
      ASSERT(.not.associated(this%mesh))
      this%sim=>sim
      call simulation_get(sim, grid)
      ASSERT(associated(grid))
      this%mesh=>grid%mesh
      nullify(grid)
      ASSERT(this%mesh%np_part>0)
    end if
    SAFE_ALLOCATE(this%data(this%mesh%np_part,this%ndim))
    this%data=this%default
    if(abs(this%default)>0.0_wp)&
      call storage_update(this)
    POP_SUB(storage_start)
    return
  end subroutine storage_start

  ! ---------------------------------------------------------
  subroutine storage_update(this)
    type(storage_t), intent(inout) :: this
    !
    integer :: i
    !
    PUSH_SUB(storage_update)
    ASSERT(associated(this%sim))
    this%data(this%mesh%np+1:,:)=0.0_wp
#if defined(HAVE_MPI)
    do i = 1, this%ndim
      call dvec_ghost_update(this%mesh%vp, this%data(:,i))
    end do
#endif
    POP_SUB(storage_update)
    return
  end subroutine storage_update

  ! ---------------------------------------------------------
  subroutine storage_stop(this)
    type(storage_t), intent(inout) :: this
    !
    PUSH_SUB(storage_stop)
    ASSERT(associated(this%sim))
    ASSERT(associated(this%mesh))
    SAFE_DEALLOCATE_A(this%data)
    POP_SUB(storage_stop)
    return
  end subroutine storage_stop

!!$  ! ---------------------------------------------------------
!!$  function storage_add(this, that) result(resl)
!!$    type(storage_t), intent(in) :: this
!!$    type(storage_t), intent(in) :: that
!!$    !
!!$    type(storage_t) :: resl
!!$    !
!!$    P!USH_SUB(storage_add)
!!$    A!SSERT(this%ndim==that%ndim)
!!$    A!SSERT(associated(this%grid, that%grid))
!!$    A!SSERT(associated(this%mesh, that%mesh))
!!$    call storage_init_copy(resl, this)
!!$    resl%data(1:this%mesh%np,:)=this%data(1:this%mesh%np,:)+that%data(1:that%mesh%np,:)
!!$    resl%default=this%default+that%default
!!$    call storage_update(resl)
!!$    P!OP_SUB(storage_add)
!!$    return
!!$  end function storage_add
!!$
!!$  ! ---------------------------------------------------------
!!$  function storage_sub(this, that) result(resl)
!!$    type(storage_t), intent(in) :: this
!!$    type(storage_t), intent(in) :: that
!!$    !
!!$    type(storage_t) :: resl
!!$    !
!!$    P!USH_SUB(storage_sub)
!!$    A!SSERT(this%ndim==that%ndim)
!!$    A!SSERT(associated(this%grid, that%grid))
!!$    A!SSERT(associated(this%mesh, that%mesh))
!!$    call storage_init_copy(resl, this)
!!$    resl%data(1:this%mesh%np,:)=this%data(1:this%mesh%np,:)-that%data(1:that%mesh%np,:)
!!$    resl%default=this%default-that%default
!!$    call storage_update(resl)
!!$    P!OP_SUB(storage_sub)
!!$    return
!!$  end function storage_sub

  ! ---------------------------------------------------------
  elemental function storage_get_size(this) result(that)
    type(storage_t), intent(in) :: this
    !
    integer :: that
    !
    that=0
    if(associated(this%mesh))&
      that=this%mesh%np
    return
  end function storage_get_size

  ! ---------------------------------------------------------
  elemental function storage_get_dimension(this) result(that)
    type(storage_t), intent(in) :: this
    !
    integer :: that
    !
    that=0
    if(this%ndim>0)&
      that=this%ndim
    return
  end function storage_get_dimension

  ! ---------------------------------------------------------
  subroutine storage_get_sim(this, that)
    type(storage_t),     target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(storage_get_sim)
    nullify(that)
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(storage_get_sim)
    return
  end subroutine storage_get_sim

!!$  ! ---------------------------------------------------------
!!$  subroutine storage_get_storage(this, that)
!!$    type(storage_t), intent(inout) :: this !Intent should only be in problems in dmultigrid_*2*
!!$    type(storage_t), intent(inout) :: that
!!$    !
!!$    integer :: i
!!$    !
!!$    P!USH_SUB(storage_get_storage)
!!$    A!SSERT(associated(this%config))
!!$    A!SSERT(associated(this%grid))
!!$    A!SSERT(associated(this%mesh))
!!$    A!SSERT(allocated(this%data))
!!$    if(associated(that%config))then
!!$      if(associated(that%grid))then
!!$        A!SSERT(this%ndim==that%ndim)
!!$        A!SSERT(associated(this%grid,that%grid))
!!$        A!SSERT(associated(that%mesh))
!!$        A!SSERT(allocated(that%data))
!!$      else
!!$        A!SSERT(.not.associated(that%mesh))
!!$        A!SSERT(.not.allocated(that%data))
!!$        call storage_start(that, this%grid, this%ndim, this%fine)
!!$      end if
!!$      if(this%fine.eqv.that%fine)then
!!$        that%data(1:this%mesh%np,:)=this%data(1:that%mesh%np,:)
!!$      else
!!$        if(this%fine)then
!!$          do i = 1, this%ndim
!!$            call dmultigrid_fine2coarse(that%grid%fine%tt, that%grid%fine%der, that%mesh, &
!!$              this%data(:,i), that%data(:,i), that%method)
!!$          end do
!!$        else
!!$          do i = 1, this%ndim
!!$            call dmultigrid_coarse2fine(that%grid%fine%tt, that%grid%der, that%mesh, &
!!$              this%data(:,i), that%data(:,i), that%order)
!!$          end do
!!$        end if
!!$      end if
!!$      call storage_update(that)
!!$    else
!!$      call storage_copy(that, this)
!!$    end if
!!$    P!OP_SUB(storage_get_storage)
!!$    return
!!$  end subroutine storage_get_storage

  ! ---------------------------------------------------------
  subroutine storage_get_storage_1d(this, that)
    type(storage_t),              target, intent(in) :: this
    real(kind=wp), dimension(:), pointer             :: that
    !
    PUSH_SUB(storage_get_storage_1d)
    nullify(that)
    if(associated(this%mesh))then
      ASSERT(this%ndim==1)
      that=>this%data(:,1)
    end if
    POP_SUB(storage_get_storage_1d)
    return
  end subroutine storage_get_storage_1d

  ! ---------------------------------------------------------
  subroutine storage_get_storage_md(this, that)
    type(storage_t),                target, intent(in) :: this
    real(kind=wp), dimension(:,:), pointer             :: that
    !
    PUSH_SUB(storage_get_storage_md)
    nullify(that)
    if(associated(this%mesh))then
      that=>this%data
    end if
    POP_SUB(storage_get_storage_md)
    return
  end subroutine storage_get_storage_md

!!$  ! ---------------------------------------------------------
!!$  subroutine storage_set_storage(this, that)
!!$    type(storage_t), intent(inout) :: this
!!$    type(storage_t), intent(inout) :: that !intent should only be in problems in dmultigrid_*2*
!!$    !
!!$    integer :: i
!!$    !
!!$    P!USH_SUB(storage_set_storage)
!!$    A!SSERT(associated(that%config))
!!$    A!SSERT(associated(that%grid))
!!$    A!SSERT(associated(that%mesh))
!!$    A!SSERT(allocated(that%data))
!!$    if(associated(this%config))then
!!$      if(associated(this%grid))then
!!$        A!SSERT(this%ndim==that%ndim)
!!$        A!SSERT(associated(this%grid,that%grid))
!!$        A!SSERT(associated(this%mesh))
!!$        A!SSERT(allocated(this%data))
!!$      else
!!$        A!SSERT(.not.associated(this%mesh))
!!$        A!SSERT(.not.allocated(this%data))
!!$        call storage_start(this, that%grid, that%ndim, that%fine)
!!$      end if
!!$      if(that%fine.eqv.this%fine)then
!!$        this%data(1:this%mesh%np,:)=that%data(1:that%mesh%np,:)
!!$      else
!!$        if(that%fine)then
!!$          do i = 1, this%ndim
!!$            call dmultigrid_fine2coarse(this%grid%fine%tt, this%grid%fine%der, this%mesh, &
!!$              that%data(:,i), this%data(:,i), this%method)
!!$          end do
!!$        else
!!$          do i = 1, this%ndim
!!$            call dmultigrid_coarse2fine(this%grid%fine%tt, this%grid%der, this%mesh, &
!!$              that%data(:,i), this%data(:,i), this%order)
!!$          end do
!!$        end if
!!$      end if
!!$      call storage_update(this)
!!$    else
!!$      call storage_copy(this, that)
!!$    end if
!!$    P!OP_SUB(storage_set_storage)
!!$    return
!!$  end subroutine storage_set_storage

  ! ---------------------------------------------------------
  elemental function storage_conform(this, that) result(rslt)
    type(storage_t), intent(in) :: this
    type(storage_t), intent(in) :: that
    !
    logical :: rslt
    !
    rslt=.false.
    if((this%ndim>0).and.(that%ndim>0))then
      if(this%ndim==that%ndim)then
        if((associated(this%mesh)).and.(associated(that%mesh)))then
          if(associated(this%mesh,that%mesh))then
            rslt=((allocated(this%data)).and.(allocated(that%data)))
          end if
        endif
      end if
    end if
    return
  end function storage_conform

  ! ---------------------------------------------------------
  subroutine storage_copy_storage(this, that)
    type(storage_t),         intent(inout) :: this
    type(storage_t), target, intent(in)    :: that
    !
    PUSH_SUB(storage_copy_storage)
    if(storage_conform(this,that))then
      this%default=that%default
    else
      call storage_end(this)
      if(that%ndim>0)then
        this%ndim=that%ndim
        this%default=that%default
        if(associated(that%mesh))then
          this%mesh=>that%mesh
          SAFE_ALLOCATE(this%data(this%mesh%np_part,this%ndim))
        end if
      end if
    end if
    if(allocated(this%data))then
      this%data(1:this%mesh%np,:)=that%data(1:that%mesh%np,:)
      call storage_update(this)
    end if
    POP_SUB(storage_copy_storage)
    return
  end subroutine storage_copy_storage

  ! ---------------------------------------------------------
  subroutine storage_end_storage(this)
    type(storage_t), intent(inout) :: this
    !
    PUSH_SUB(storage_end_storage)
    SAFE_DEALLOCATE_A(this%data)
    this%default=0.0_wp
    this%ndim=0
    nullify(this%mesh)
    POP_SUB(storage_end_storage)
    return
  end subroutine storage_end_storage

  ! ---------------------------------------------------------
  subroutine storage_intrpl_init(this, that, type)
    type(storage_intrpl_t),  intent(out) :: this
    type(storage_t), target, intent(in)  :: that
    integer,       optional, intent(in)  :: type
    !
    PUSH_SUB(storage_intrpl_init)
    ASSERT(allocated(that%data))
    this%self=>that
    call intrpl_init(this%intrp, that%sim, that%data, type=type, default=that%default)
    POP_SUB(storage_intrpl_init)
    return
  end subroutine storage_intrpl_init

  ! ---------------------------------------------------------
  subroutine storage_intrpl_eval_1d(this, x, val, ierr)
    type(storage_intrpl_t),       intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),                intent(out) :: val
    integer,                      intent(out) :: ierr
    !
    PUSH_SUB(storage_intrpl_eval_1d)
    call intrpl_eval(this%intrp, x, val, ierr)
    POP_SUB(storage_intrpl_eval_1d)
    return
  end subroutine storage_intrpl_eval_1d

  ! ---------------------------------------------------------
  subroutine storage_intrpl_eval_md(this, x, val, ierr)
    type(storage_intrpl_t),       intent(in)  :: this
    real(kind=wp),  dimension(:), intent(in)  :: x
    real(kind=wp),  dimension(:), intent(out) :: val
    integer,                      intent(out) :: ierr
    !
    PUSH_SUB(storage_intrpl_eval_md)
    call intrpl_eval(this%intrp, x, val, ierr)
    POP_SUB(storage_intrpl_eval_md)
    return
  end subroutine storage_intrpl_eval_md

  ! ---------------------------------------------------------
  subroutine storage_intrpl_set(this, that)
    type(storage_intrpl_t), intent(inout) :: this
    type(basis_t),          intent(in)    :: that
    !
    PUSH_SUB(storage_intrpl_set)
    call intrpl_set(this%intrp, that)
    POP_SUB(storage_intrpl_set)
    return
  end subroutine storage_intrpl_set

  ! ---------------------------------------------------------
  subroutine storage_intrpl_copy(this, that)
    type(storage_intrpl_t), intent(inout) :: this
    type(storage_intrpl_t), intent(in)    :: that
    !
    PUSH_SUB(storage_intrpl_copy)
    this%self=>that%self
    call intrpl_copy(this%intrp, that%intrp)
    POP_SUB(storage_intrpl_copy)
    return
  end subroutine storage_intrpl_copy

  ! ---------------------------------------------------------
  subroutine storage_intrpl_end(this)
    type(storage_intrpl_t), intent(inout) :: this
    !
    PUSH_SUB(storage_intrpl_end)
    nullify(this%self)
    call intrpl_end(this%intrp)
    POP_SUB(storage_intrpl_end)
    return
  end subroutine storage_intrpl_end

end module storage_m

!! Local Variables:
!! mode: f90
!! End:
