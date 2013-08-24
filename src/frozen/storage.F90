#include "global.h"

module storage_m

  use global_m
  use messages_m
  use profiling_m

#if defined(HAVE_MPI)
  use boundaries_m, only: dvec_ghost_update
#endif
  use kinds_m,      only: wp

  use grid_m, only: &
    grid_t

  use mesh_m, only: &
    mesh_t

  use simulation_m, only: &
    simulation_t,         &
    simulation_get

  use interpolation_m, only: &
    interpolation_t,         &
    interpolation_init,      &
    interpolation_eval,      &
    interpolation_copy,      &
    interpolation_end

  use interpolation_m, only:    &
    STORAGE_INTRP_OK=>INTRP_OK, &
    STORAGE_INTRP_OD=>INTRP_OD, &
    STORAGE_INTRP_NI=>INTRP_NI

  implicit none

  private
  public ::      &
    operator(+), &
    operator(-)

  public ::                &
    storage_init,          &
    storage_start,         &
    storage_update,        &
    storage_get,           &
    storage_get_size,      &
    storage_get_dimension, &
    storage_get_default,   &
    storage_get_storage,   &
    storage_copy,          &
    storage_end

  public ::           &
    STORAGE_INTRP_OK, &
    STORAGE_INTRP_OD, &
    STORAGE_INTRP_NI

  public ::                     &
    storage_interpolation_init, &
    storage_interpolation_eval, &
    storage_interpolation_copy, &
    storage_interpolation_end

   type, public :: storage_t
    private
    integer                                    :: ndim    = 0
    type(simulation_t),                pointer :: sim     =>null()
    type(mesh_t),                      pointer :: mesh    =>null()
    real(kind=wp)                              :: default = 0.0_wp
    real(kind=wp), dimension(:,:), allocatable :: storage
  end type storage_t

  type, public :: storage_interpolation_t
    private
    type(storage_t), pointer :: self =>null()
    type(interpolation_t)    :: intrp
  end type storage_interpolation_t

  interface operator(+)
    module procedure storage_add
  end interface operator(+)

  interface operator(-)
    module procedure storage_sub
  end interface operator(-)

  interface storage_init
    module procedure storage_init_simple
    module procedure storage_init_copy
  end interface storage_init

  interface storage_get
    module procedure storage_get_simulation
    module procedure storage_get_mesh
  end interface storage_get

  interface storage_get_storage
    module procedure storage_get_storage_1d
    module procedure storage_get_storage_2d
  end interface storage_get_storage

  interface storage_interpolation_eval
    module procedure storage_interpolation_eval_1d
    module procedure storage_interpolation_eval_2d
  end interface storage_interpolation_eval

contains

  ! ---------------------------------------------------------
  subroutine storage_init_simple(this, ndim, default)
    type(storage_t),         intent(out) :: this
    integer,       optional, intent(in)  :: ndim
    real(kind=wp), optional, intent(in)  :: default
    !
    PUSH_SUB(storage_init_simple)
    this%ndim=1
    if(present(ndim))then
      ASSERT(ndim>0)
      this%ndim=ndim
    end if
    this%sim=>null()
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
    call storage_end(this)
    ASSERT(that%ndim>0)
    ASSERT(associated(that%sim))
    ASSERT(associated(that%mesh))
    this%ndim=that%ndim
    this%sim=>that%sim
    this%mesh=>that%mesh
    this%default=that%default
    SAFE_ALLOCATE(this%storage(this%mesh%np_part,this%ndim))
    this%storage=this%default
    call storage_update(this)
    POP_SUB(storage_init_copy)
    return
  end subroutine storage_init_copy

  ! ---------------------------------------------------------
  subroutine storage_start(this, sim, fine)
    type(storage_t),            intent(inout) :: this
    type(simulation_t), target, intent(in)    :: sim
    logical,          optional, intent(in)    :: fine
    !
    type(grid_t), pointer :: grid
    logical               :: fn
    !
    PUSH_SUB(storage_start)
    fn=.false.
    if(present(fine))fn=fine
    ASSERT(this%ndim>0)
    ASSERT(.not.associated(this%sim))
    this%sim=>sim
    nullify(grid)
    call simulation_get(this%sim, grid)
    ASSERT(associated(grid))
    if(fn.and.grid%have_fine_mesh)then
      this%mesh=>grid%fine%mesh
    else
      this%mesh=>grid%mesh
    end if
    ASSERT(this%mesh%np_part>0)
    SAFE_ALLOCATE(this%storage(this%mesh%np_part,this%ndim))
    this%storage=this%default
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
    this%storage(this%mesh%np+1:,:)=0.0_wp
#if defined(HAVE_MPI)
    do i = 1, this%ndim
      call dvec_ghost_update(this%mesh%vp, this%storage(:,i))
    end do
#endif
    POP_SUB(storage_update)
    return
  end subroutine storage_update

  ! ---------------------------------------------------------
  function storage_add(this, that) result(resl)
    type(storage_t), intent(in) :: this
    type(storage_t), intent(in) :: that
    !
    type(storage_t) :: resl
    !
    PUSH_SUB(storage_add)
    ASSERT(this%ndim==that%ndim)
    ASSERT(associated(this%mesh, that%mesh))
    call storage_init_copy(resl, this)
    resl%storage(1:this%mesh%np,:)=this%storage(1:this%mesh%np,:)+that%storage(1:that%mesh%np,:)
    resl%default=this%default+that%default
    call storage_update(resl)
    POP_SUB(storage_add)
    return
  end function storage_add

  ! ---------------------------------------------------------
  function storage_sub(this, that) result(resl)
    type(storage_t), intent(in) :: this
    type(storage_t), intent(in) :: that
    !
    type(storage_t) :: resl
    !
    PUSH_SUB(storage_sub)
    ASSERT(this%ndim==that%ndim)
    ASSERT(associated(this%mesh, that%mesh))
    call storage_init_copy(resl, this)
    resl%storage(1:this%mesh%np,:)=this%storage(1:this%mesh%np,:)-that%storage(1:that%mesh%np,:)
    resl%default=this%default-that%default
    call storage_update(resl)
    POP_SUB(storage_sub)
    return
  end function storage_sub

  ! ---------------------------------------------------------
  elemental function storage_get_size(this) result(that)
    type(storage_t), intent(in) :: this
    !
    integer :: that
    !
    that=0
    if(associated(this%sim))&
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
  elemental function storage_get_default(this) result(that)
    type(storage_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=0.0_wp
    if(this%ndim>0)&
      that=this%default
    return
  end function storage_get_default

  ! ---------------------------------------------------------
  subroutine storage_get_simulation(this, that)
    type(storage_t),     target, intent(in) :: this
    type(simulation_t), pointer             :: that
    !
    PUSH_SUB(storage_get_simulation)
    that=>null()
    if(associated(this%sim))&
      that=>this%sim
    POP_SUB(storage_get_simulation)
    return
  end subroutine storage_get_simulation

  ! ---------------------------------------------------------
  subroutine storage_get_mesh(this, that)
    type(storage_t), target, intent(in) :: this
    type(mesh_t),   pointer             :: that
    !
    PUSH_SUB(storage_get_mesh)
    that=>null()
    if(associated(this%mesh))&
      that=>this%mesh
    POP_SUB(storage_get_mesh)
    return
  end subroutine storage_get_mesh

  ! ---------------------------------------------------------
  subroutine storage_get_storage_1d(this, that)
    type(storage_t),              target, intent(in) :: this
    real(kind=wp), dimension(:), pointer             :: that
    !
    PUSH_SUB(storage_get_storage_1d)
    that=>null()
    if(associated(this%sim))then
      ASSERT(this%ndim==1)
      that=>this%storage(:,1)
    end if
    POP_SUB(storage_get_storage_1d)
    return
  end subroutine storage_get_storage_1d

  ! ---------------------------------------------------------
  subroutine storage_get_storage_2d(this, that)
    type(storage_t),                target, intent(in) :: this
    real(kind=wp), dimension(:,:), pointer             :: that
    !
    PUSH_SUB(storage_get_storage_2d)
    that=>null()
    if(associated(this%sim))then
      that=>this%storage
    end if
    POP_SUB(storage_get_storage_2d)
    return
  end subroutine storage_get_storage_2d

  ! ---------------------------------------------------------
  subroutine storage_copy(this, that)
    type(storage_t),         intent(out) :: this
    type(storage_t), target, intent(in)  :: that
    !
    PUSH_SUB(storage_copy)
    call storage_end(this)
    if(that%ndim>0)then
      this%ndim   = that%ndim
      this%default= that%default
      if(associated(this%sim))then
        this%sim    =>that%sim
        this%mesh   =>that%mesh
        SAFE_ALLOCATE(this%storage(this%mesh%np_part,this%ndim))
        this%storage(1:this%mesh%np,:)=that%storage(1:that%mesh%np,:)
        call storage_update(this)
      end if
    end if
    POP_SUB(storage_copy)
    return
  end subroutine storage_copy

  ! ---------------------------------------------------------
  subroutine storage_end(this)
    type(storage_t), intent(inout) :: this
    !
    PUSH_SUB(storage_end)
    if(associated(this%sim))then
      SAFE_DEALLOCATE_A(this%storage)
    end if
    this%default=0.0_wp
    nullify(this%mesh, this%sim)
    this%ndim=0
    POP_SUB(storage_end)
    return
  end subroutine storage_end

 ! ---------------------------------------------------------
  subroutine storage_interpolation_init(this, that, type)
    type(storage_interpolation_t), intent(out) :: this
    type(storage_t),       target, intent(in)  :: that
    integer,             optional, intent(in)  :: type
    !
    PUSH_SUB(storage_interpolation_init)
    this%self=>that
    if(that%ndim>0)&
      call interpolation_init(this%intrp, that%sim, that%storage, type=type, default=that%default)
    POP_SUB(storage_interpolation_init)
    return
  end subroutine storage_interpolation_init

  ! ---------------------------------------------------------
  subroutine storage_interpolation_eval_1d(this, x, val, ierr)
    type(storage_interpolation_t), intent(in)  :: this
    real(kind=wp),   dimension(:), intent(in)  :: x
    real(kind=wp),                 intent(out) :: val
    integer,                       intent(out) :: ierr
    !
    PUSH_SUB(storage_interpolation_eval_1d)
    if(associated(this%self%sim))then
      call interpolation_eval(this%intrp, x, val, ierr)
    else
      val=this%self%default
    end if
    POP_SUB(storage_interpolation_eval_1d)
    return
  end subroutine storage_interpolation_eval_1d

  ! ---------------------------------------------------------
  subroutine storage_interpolation_eval_2d(this, x, val, ierr)
    type(storage_interpolation_t), intent(in)  :: this
    real(kind=wp),   dimension(:), intent(in)  :: x
    real(kind=wp),   dimension(:), intent(out) :: val
    integer,                       intent(out) :: ierr
    !
    PUSH_SUB(storage_interpolation_eval_2d)
    if(associated(this%self%sim))then
      call interpolation_eval(this%intrp, x, val, ierr)
    else
      val=this%self%default
    end if
    POP_SUB(storage_interpolation_eval_2d)
    return
  end subroutine storage_interpolation_eval_2d

  ! ---------------------------------------------------------
  subroutine storage_interpolation_copy(this, that)
    type(storage_interpolation_t), intent(out) :: this
    type(storage_interpolation_t), intent(in)  :: that
    !
    PUSH_SUB(storage_interpolation_copy)
    this%self=>that%self
    call interpolation_copy(this%intrp, that%intrp)
    POP_SUB(storage_interpolation_copy)
    return
  end subroutine storage_interpolation_copy

  ! ---------------------------------------------------------
  subroutine storage_interpolation_end(this)
    type(storage_interpolation_t), intent(inout) :: this
    !
    PUSH_SUB(storage_interpolation_end)
    call interpolation_end(this%intrp)
    this%self=>null()
    POP_SUB(storage_interpolation_end)
    return
  end subroutine storage_interpolation_end

end module storage_m

!! Local Variables:
!! mode: f90
!! End:
