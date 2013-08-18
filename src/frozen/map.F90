#include "global.h"

module map_m

  use global_m
  use kinds_m, only: wp
  use mesh_m, only: mesh_t
  use messages_m
#ifdef HAVE_MPI
  use par_vec_m, only: dvec_gather, dvec_scatter
#endif
  use profiling_m

  implicit none

  private
  public ::       &
    map_init,     &
    operator(==), &
    operator(/=), &
    map_map,      &
    map_copy,     &
    map_end

  interface operator(==)
    module procedure map_equal
  end interface operator(==)

  interface operator(/=)
    module procedure map_not_equal
  end interface operator(/=)

  type, public :: map_t
    private
    type(mesh_t),              pointer :: imsh =>null()
    type(mesh_t),              pointer :: omsh =>null()
    integer, dimension(:), allocatable :: map
  end type map_t

contains
  
  ! ---------------------------------------------------------
  subroutine map_init(this, mesh_out, mesh_in)
    type(map_t),          intent(out) :: this
    type(mesh_t), target, intent(in)  :: mesh_out
    type(mesh_t), target, intent(in)  :: mesh_in
    !
    integer, dimension(MAX_DIM) :: ix
    integer                     :: i, n
    !
    this%imsh=>mesh_in
    this%omsh=>mesh_out
    SAFE_ALLOCATE(this%map(this%imsh%np_part_global))
    this%map=0
    do i = 1, this%imsh%np_part_global
      ix=this%imsh%idx%lxyz(i,:)
      if((all(ix>=this%omsh%idx%nr(1,:))).and.(all(ix<=this%omsh%idx%nr(2,:))))then
        n=this%omsh%idx%lxyz_inv(ix(1),ix(2),ix(3))
        this%map(i)=n
      end if
    end do
    return
  end subroutine map_init

  ! -----------------------------------------------------
  elemental function map_equal(this, that) result(eqv)
    type(map_t), intent(in) :: this
    type(map_t), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=(associated(this%imsh, that%imsh).and.associated(this%omsh, that%omsh))
    return
  end function map_equal
  
  ! -----------------------------------------------------
  elemental function map_not_equal(this, that) result(neqv)
    type(map_t), intent(in) :: this
    type(map_t), intent(in) :: that
    !
    logical :: neqv
    !
    neqv=.not.map_equal(this, that)
    return
  end function map_not_equal

  ! ---------------------------------------------------------
  pure subroutine map_map_s2s(this, array_out, array_in)
    type(map_t),                 intent(in)  :: this
    real(kind=wp), dimension(:), intent(out) :: array_out
    real(kind=wp), dimension(:), intent(in)  :: array_in
    !
    integer :: i, j, np_i, np_o
    !
    np_i=size(array_in)
    np_o=size(array_out)
    do i = 1, np_i
      j=this%map(i)
      if((j>0).and.(j<np_o+1))&
        array_out(j)=array_in(i)
    end do
    return
  end subroutine map_map_s2s

#ifdef HAVE_MPI
  ! ---------------------------------------------------------
  subroutine map_map_s2p(this, array_out, array_in)
    type(map_t),                 intent(in)  :: this
    real(kind=wp), dimension(:), intent(out) :: array_out
    real(kind=wp), dimension(:), intent(in)  :: array_in
    !
    real(kind=wp), dimension(:), allocatable :: arro
    !
    if(this%omsh%vp%rank==this%omsh%vp%root)then
      SAFE_ALLOCATE(arro(this%omsh%np_part_global))
    else
      SAFE_ALLOCATE(arro(1))
    end if
    if(this%imsh%vp%rank==this%imsh%vp%root)call map_map_s2s(this, arro, array_in)
    call dvec_scatter(this%omsh%vp, this%omsh%vp%root, arro, array_out)
    SAFE_DEALLOCATE_A(arro)
    return
  end subroutine map_map_s2p

  ! ---------------------------------------------------------
  subroutine map_map_p2s(this, array_out, array_in)
    type(map_t),                 intent(in)  :: this
    real(kind=wp), dimension(:), intent(out) :: array_out
    real(kind=wp), dimension(:), intent(in)  :: array_in
    !
    real(kind=wp), dimension(:), allocatable :: arri
    !
    if(this%imsh%vp%rank==this%imsh%vp%root)then
      SAFE_ALLOCATE(arri(this%imsh%np_part_global))
    else
      SAFE_ALLOCATE(arri(1))
    end if
    call dvec_gather(this%imsh%vp, this%imsh%vp%root, arri, array_in)
    call map_map_s2s(this, array_out, arri)
    SAFE_DEALLOCATE_A(arri)
    return
  end subroutine map_map_p2s

  ! ---------------------------------------------------------
  subroutine map_map_p2p(this, array_out, array_in)
    type(map_t),                 intent(in)  :: this
    real(kind=wp), dimension(:), intent(out) :: array_out
    real(kind=wp), dimension(:), intent(in)  :: array_in
    !
    real(kind=wp), dimension(:), allocatable :: arri, arro
    !
    if(this%imsh%vp%rank==this%imsh%vp%root)then
      SAFE_ALLOCATE(arri(this%imsh%np_part_global))
    else
      SAFE_ALLOCATE(arri(1))
    end if
    if(this%omsh%vp%rank==this%omsh%vp%root)then
      SAFE_ALLOCATE(arro(this%omsh%np_part_global))
    else
      SAFE_ALLOCATE(arro(1))
    end if
    call dvec_gather(this%imsh%vp, this%imsh%vp%root, arri, array_in)
    if(this%imsh%vp%rank==this%imsh%vp%root)call map_map_s2s(this, arro, arri)
    call dvec_scatter(this%omsh%vp, this%omsh%vp%root, arro, array_out)
    SAFE_DEALLOCATE_A(arri)
    SAFE_DEALLOCATE_A(arro)
    return
  end subroutine map_map_p2p
#endif

  ! ---------------------------------------------------------
  subroutine map_map(this, array_out, array_in)
    type(map_t),                 intent(in)  :: this
    real(kind=wp), dimension(:), intent(out) :: array_out
    real(kind=wp), dimension(:), intent(in)  :: array_in
    !
#ifdef HAVE_MPI
    if(this%imsh%parallel_in_domains)then
      if(this%omsh%parallel_in_domains)then
        call map_map_p2p(this, array_out, array_in)
      else
        call map_map_p2s(this, array_out, array_in)
      end if
    else
      if(this%omsh%parallel_in_domains)then
        call map_map_s2p(this, array_out, array_in)
      else
        call map_map_s2s(this, array_out, array_in)
      end if
    end if
#else
    call map_map_s2s(this, array_out, array_in)
#endif
    return
  end subroutine map_map

  ! ---------------------------------------------------------
  subroutine map_copy(this_out, this_in)
    type(map_t), intent(out) :: this_out
    type(map_t), intent(in)  :: this_in
    !
    call map_end(this_out)
    if(associated(this_in%imsh))then
      this_out%imsh=>this_in%imsh
      this_out%omsh=>this_in%omsh
      SAFE_ALLOCATE(this_out%map(this_out%imsh%np_part_global))
      this_out%map=this_in%map
    end if
    return
  end subroutine map_copy

  ! ---------------------------------------------------------
  subroutine map_end(this)
    type(map_t), intent(inout) :: this
    !
    if(associated(this%imsh))then
      SAFE_DEALLOCATE_A(this%map)
      nullify(this%imsh, this%omsh)
    end if
    return
  end subroutine map_end

end module map_m

module wrap_mesh_m

  use mesh_m, only: mesh_t

  implicit none

  private
  public :: &
    wrap_mesh_init, &
    operator(==),   &
    operator(/=),   &
    wrap_mesh_hash, &
    wrap_mesh_copy, &
    wrap_mesh_end

  interface operator(==)
    module procedure wrap_mesh_equal
  end interface operator(==)

  interface operator(/=)
    module procedure wrap_mesh_not_equal
  end interface operator(/=)

  type, public :: wrap_mesh_t
    private
    type(mesh_t), pointer :: mesh =>null()
  end type wrap_mesh_t

contains

  ! -----------------------------------------------------
  subroutine wrap_mesh_init(this, mesh)
    type(wrap_mesh_t), intent(out) :: this
    type(mesh_t), target, intent(in) :: mesh
    !
    this%mesh=>mesh
    return
  end subroutine wrap_mesh_init

  ! -----------------------------------------------------
  elemental function wrap_mesh_equal(this, that) result(eqv)
    type(wrap_mesh_t), intent(in) :: this
    type(wrap_mesh_t), intent(in) :: that
    !
    logical :: eqv
    !
    eqv=associated(this%mesh, that%mesh)
    return
  end function wrap_mesh_equal
  
  ! ---------------------------------------------------------
  elemental function  wrap_mesh_hash(this, size) result(n)
    type(wrap_mesh_t), intent(in) :: this
    integer,           intent(in) :: size
    !
    integer :: n
    !
    n=1
    return
  end function wrap_mesh_hash

  ! -----------------------------------------------------
  elemental function wrap_mesh_not_equal(this, that) result(neqv)
    type(wrap_mesh_t), intent(in) :: this
    type(wrap_mesh_t), intent(in) :: that
    !
    logical :: neqv
    !
    neqv=.not.associated(this%mesh, that%mesh)
    return
  end function wrap_mesh_not_equal

  ! -----------------------------------------------------
  subroutine wrap_mesh_copy(this_out, this_in)
    type(wrap_mesh_t), intent(out) :: this_out
    type(wrap_mesh_t), intent(in)  :: this_in
    !
    this_out%mesh=>this_in%mesh
    return
  end subroutine wrap_mesh_copy

  ! -----------------------------------------------------
  elemental subroutine wrap_mesh_end(this)
    type(wrap_mesh_t), intent(inout) :: this
    !
    this%mesh=>null()
    return
  end subroutine wrap_mesh_end

end module wrap_mesh_m

#define TEMPLATE_NAME meshmap
#define MODULE_TYPE_KEY wrap_mesh_m
#define TYPE_KEY wrap_mesh_t
#define HASH_FUNCTION wrap_mesh_hash
#define MODULE_TYPE_VAL map_m
#define TYPE_VAL map_t
#include "thash.F90"
#undef MODULE_TYPE_KEY
#undef TYPE_KEY
#undef HASH_FUNCTION
#undef MODULE_TYPE_VAL
#undef TYPE_VAL
#undef TEMPLATE_NAME

!! Local Variables:
!! mode: f90
!! End:
