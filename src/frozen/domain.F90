#include "global.h"

#undef LIST_TEMPLATE_NAME
#undef LIST_TYPE_NAME
#undef LIST_TYPE_MODULE_NAME

#define LIST_TEMPLATE_NAME domain_node

module domain_oct_m

  use basis_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_intrf_oct_m
  use json_oct_m
  use kinds_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m

#define LIST_INCLUDE_PREFIX
#include "tlist_inc.F90"
#undef LIST_INCLUDE_PREFIX

  implicit none

  private

  public ::   &
    domain_t

  public ::        &
    domain_init,   &
    domain_extend, &
    domain_start,  &
    domain_in,     &
    domain_get,    &
    domain_copy,   &
    domain_end

#define LIST_INCLUDE_HEADER
#include "tlist_inc.F90"
#undef LIST_INCLUDE_HEADER

  type :: domain_node_t
    private
    type(json_object_t), pointer :: config =>null()
    type(domain_t),      pointer :: domain =>null()
    type(basis_t)                :: basis
  end type domain_node_t

  type :: domain_t
    private
    type(space_t),      pointer :: space =>null()
    type(grid_intrf_t), pointer :: igrid =>null()
    type(domain_node_list_t)    :: list
  end type domain_t

  interface domain_get
    module procedure domain_get_space
  end interface domain_get

contains
  
#define LIST_INCLUDE_BODY
#include "tlist_inc.F90"
#undef LIST_INCLUDE_BODY

  ! ---------------------------------------------------------
  subroutine domain_node_new(this)
    type(domain_node_t), pointer :: this

    PUSH_SUB(domain_node_new)

    nullify(this)
    SAFE_ALLOCATE(this)

    POP_SUB(domain_node_new)
  end subroutine domain_node_new

  ! ---------------------------------------------------------
  subroutine domain_node_del(this)
    type(domain_node_t), pointer :: this

    PUSH_SUB(domain_node_del)

    if(associated(this))then
      call domain_node_end(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)

    POP_SUB(domain_node_del)
  end subroutine domain_node_del

  ! ---------------------------------------------------------
  subroutine domain_node_init(this, domain, config)
    type(domain_node_t),                   intent(out) :: this
    type(domain_t),                target, intent(in)  :: domain
    type(json_object_t), optional, target, intent(in)  :: config

    type(space_t), pointer :: space

    PUSH_SUB(domain_node_init)

    nullify(this%config, space)
    this%domain => domain
    if(present(config))then
      this%config => config
      call domain_get(domain, space)
      ASSERT(associated(space))
      call basis_init(this%basis, space, config)
      nullify(space)
    end if

    POP_SUB(domain_node_init)
  end subroutine domain_node_init

  ! ---------------------------------------------------------
  function domain_node_in(this, x) result(in)
    type(domain_node_t),         intent(in) :: this
    real(kind=wp), dimension(:), intent(in) :: x

    logical :: in

    real(kind=wp), dimension(size(x)) :: y

    PUSH_SUB(domain_node_in)
    
    in = .false.
    if(associated(this%config))then
      call basis_to_internal(this%basis, x, y)
      in = domain_in(this%domain, y)
    else
      in = domain_in(this%domain, x)
    end if

    POP_SUB(domain_node_in)
  end function domain_node_in

  ! ---------------------------------------------------------
  subroutine domain_node_copy(this, that)
    type(domain_node_t), intent(inout) :: this
    type(domain_node_t), intent(in)    :: that

    PUSH_SUB(domain_node_copy)

    call domain_node_end(this)
    this%config => that%config
    this%domain => that%domain
    call basis_copy(this%basis, that%basis)

    POP_SUB(domain_node_copy)
  end subroutine domain_node_copy

  ! ---------------------------------------------------------
  subroutine domain_node_end(this)
    type(domain_node_t), intent(inout) :: this

    PUSH_SUB(domain_node_end)

    nullify(this%config, this%domain)
    call basis_end(this%basis)

    POP_SUB(domain_node_end)
  end subroutine domain_node_end

  ! ---------------------------------------------------------
  subroutine domain_new_node(this, that)
    type(domain_t),       intent(inout) :: this
    type(domain_node_t), pointer        :: that

    PUSH_SUB(domain_new_node)

    call domain_node_new(that)
    call domain_node_list_push(this%list, that)

    POP_SUB(domain_new_node)
  end subroutine domain_new_node

  ! ---------------------------------------------------------
  subroutine domain_del_node(this, that)
    type(domain_t),       intent(inout) :: this
    type(domain_node_t), pointer        :: that

    PUSH_SUB(domain_del_node)

    if(associated(that))then
      call domain_node_list_del(this%list, that)
      call domain_node_del(that)
    end if

    POP_SUB(domain_del_node)
  end subroutine domain_del_node

  ! ---------------------------------------------------------
  subroutine domain_init(this, space)
    type(domain_t),        intent(out) :: this
    type(space_t), target, intent(in)  :: space

    PUSH_SUB(domain_init)
    
    this%space => space
    call domain_node_list_init(this%list)

    POP_SUB(domain_init)
  end subroutine domain_init

  ! ---------------------------------------------------------
  subroutine domain_start(this, igrid)
    type(domain_t),             intent(inout) :: this
    type(grid_intrf_t), target, intent(in)    :: igrid

    PUSH_SUB(domain_start)
    
    ASSERT(associated(this%space))
    this%igrid => igrid

    POP_SUB(domain_start)
  end subroutine domain_start

  ! ---------------------------------------------------------
  subroutine domain_extend(this, that, config)
    type(domain_t),                intent(inout) :: this
    type(domain_t),                intent(in)    :: that
    type(json_object_t), optional, intent(in)    :: config

    type(domain_node_t), pointer :: node

    PUSH_SUB(domain_extend)

    nullify(node)
    ASSERT(associated(this%space))
    ASSERT(associated(that%space))
    ASSERT(this%space==that%space)
    ASSERT(.not.associated(this%igrid))
    call domain_new_node(this, node)
    call domain_node_init(node, that, config)
    nullify(node)

    POP_SUB(domain_extend)
  end subroutine domain_extend

  ! ---------------------------------------------------------
  function domain_in(this, x) result(in)
    type(domain_t),              intent(in) :: this
    real(kind=wp), dimension(:), intent(in) :: x

    logical :: in

    type(simul_box_t), pointer :: sb
    type(geometry_t),  pointer :: geo

    PUSH_SUB(domain_in)
    
    nullify(sb, geo)
    ASSERT(associated(this%space))
    ASSERT(associated(this%igrid))
    call grid_intrf_get(this%igrid, sb)
    ASSERT(associated(sb))
    call grid_intrf_get(this%igrid, geo)
    ASSERT(associated(geo))
    in = simul_box_in_box(sb, geo, x)
    nullify(sb, geo)

    POP_SUB(domain_in)
  end function domain_in

  ! ---------------------------------------------------------
  subroutine domain_get_space(this, that)
    type(domain_t), target, intent(in) :: this
    type(space_t), pointer             :: that

    PUSH_SUB(domain_get_space)

    nullify(that)
    if(associated(this%space)) that => this%space

    POP_SUB(domain_get_space)
  end subroutine domain_get_space

  ! ---------------------------------------------------------
  subroutine domain__copy__(this, that)
    type(domain_t), intent(inout) :: this
    type(domain_t), intent(in)    :: that

    type(domain_node_list_iterator_t) :: iter
    type(domain_node_t),      pointer :: onod, inod
    integer                           :: ierr

    PUSH_SUB(domain__copy__)

    ASSERT(associated(that%space))
    call domain_init(this, that%space)
    call domain_node_list_init(iter, that%list)
    do
      nullify(onod, inod)
      call domain_node_list_next(iter, inod, ierr)
      if(ierr/=DOMAIN_NODE_LIST_OK)exit
      ASSERT(associated(inod))
      call domain_new_node(this, onod)
      call domain_node_copy(onod, inod)
    end do
    call domain_node_list_end(iter)
    nullify(onod, inod)

    POP_SUB(domain__copy__)
  end subroutine domain__copy__

  ! ---------------------------------------------------------
  subroutine domain_copy(this, that)
    type(domain_t), intent(inout) :: this
    type(domain_t), intent(in)    :: that

    PUSH_SUB(domain_copy)

    call domain_end(this)
    if(associated(that%space))then
      call domain__copy__(this, that)
      if(associated(that%igrid)) call domain_start(this, that%igrid)
    end if

    POP_SUB(domain_copy)
  end subroutine domain_copy

  ! ---------------------------------------------------------
  subroutine domain_end(this)
    type(domain_t), intent(inout) :: this

    type(domain_node_t), pointer :: node

    PUSH_SUB(domain_end)

    do
      nullify(node)
      call domain_node_list_pop(this%list, node)
      if(.not.associated(node))exit
      call domain_node_del(node)
    end do
    nullify(node)
    call domain_node_list_end(this%list)
    nullify(this%space, this%igrid)

    POP_SUB(domain_end)
  end subroutine domain_end

end module domain_oct_m

#undef LIST_TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
