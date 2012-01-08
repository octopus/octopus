#include "global.h"

module frozen_epot_m
  use frozen_basis_m, only: frozen_basis_t, frozen_basis_to_internal
  use frozen_interpolation_m, only: frozen_interpolation_t, &
    frozen_interpolation_init, frozen_interpolation, frozen_interpolation_copy, &
    frozen_interpolation_end
  use geometry_m, only: geometry_t
  use global_m
  use loct_m, only: loct_pointer_copy
  use mesh_m, only: mesh_t, mesh_init_from_dump, mesh_copy, mesh_end
  use messages_m
  use profiling_m
  use simul_box_m, only: simul_box_in_box
  use species_m, only: species_zval

  implicit none

  private
  public ::           &
    frozen_epot_t,    &
    frozen_epot_init, &
    frozen_epot_vpsl, &
    frozen_epot_copy, &
    frozen_epot_end

  integer, public, parameter :: NONE      = 1
  integer, public, parameter :: CLASSICAL = 2
  integer, public, parameter :: EXTERNAL  = 3
     
  type frozen_epot_t
    private
    integer                                     :: type
    type(frozen_basis_t), pointer               :: basis
    type(geometry_t),     pointer               :: geo
    type(frozen_interpolation_t)                :: intrp
    type(mesh_t)                                :: mesh
    FLOAT,                pointer, dimension(:) :: vpsl
  end type frozen_epot_t

contains

  ! ---------------------------------------------------------
  subroutine frozen_epot_init(this, basis, geo, type, intrp, iunit)
    type(frozen_epot_t),          intent(inout) :: this
    type(frozen_basis_t), target, intent(in)    :: basis
    type(geometry_t),     target, intent(in)    :: geo
    integer,                      intent(in)    :: type
    integer,                      intent(in)    :: intrp
    integer,                      intent(in)    :: iunit
    !
    integer :: np, gb
    !
    this%type=type
    this%basis=>basis
    this%geo=>geo
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    call mesh_init_from_dump(this%mesh, basis%sb, basis%cv, iunit)
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    read(iunit) np
    ASSERT(np==this%mesh%np_global)
    SAFE_ALLOCATE(this%vpsl(np))
    read(iunit) this%vpsl
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    call frozen_interpolation_init(this%intrp, intrp, this%mesh, this%vpsl)
    return
  end subroutine frozen_epot_init

  ! ---------------------------------------------------------
  function frozen_classical(x, y, c) result(v)
    FLOAT, dimension(:), intent(in)  :: x
    FLOAT, dimension(:), intent(in)  :: y
    FLOAT,               intent(in)  :: c
    !
    FLOAT :: v
    !
    FLOAT :: r
    !
    r=sqrt(sum((x-y)**2))
    if(r<r_small) r=r_small
    v=-c/r
    return
  end function frozen_classical

  ! ---------------------------------------------------------
  subroutine frozen_epot_classical(this, x, vpsl)
    type(frozen_epot_t), intent(in)  :: this
    FLOAT, dimension(:), intent(in)  :: x
    FLOAT,               intent(out) :: vpsl
    !
    FLOAT, dimension(size(x)) :: y
    integer                   :: i
    !
    vpsl=M_ZERO
    call frozen_basis_to_internal(this%basis, x, y)
    do i = 1, this%geo%natoms
      vpsl=vpsl+frozen_classical(y, this%geo%atom(i)%x, species_zval(this%geo%atom(i)%spec))
    end do
    do i = 1, this%geo%ncatoms
      vpsl=vpsl+frozen_classical(y, this%geo%catom(i)%x, this%geo%catom(i)%charge)
    end do
    return
  end subroutine frozen_epot_classical

  ! ---------------------------------------------------------
  subroutine frozen_epot_vpsl(this, x, vpsl)
    type(frozen_epot_t), intent(in)  :: this
    FLOAT, dimension(:), intent(in)  :: x
    FLOAT,               intent(out) :: vpsl
    !
    FLOAT, dimension(size(x)) :: y
    !
    select case(this%type)
    case(NONE)
      vpsl=M_ZERO
    case(CLASSICAL)
      call frozen_epot_classical(this, x, vpsl)
    case(EXTERNAL)
      call frozen_basis_to_internal(this%basis, x, y)
      if(simul_box_in_box(this%basis%sb, this%geo, y))then
        call frozen_interpolation(this%intrp, y, vpsl)
      else
        call frozen_epot_classical(this, x, vpsl)
      end if
    case default
      write(message(1), '(a,i2)') "Unknown potential type: ", this%type
      call messages_fatal(1)
    end select
    return
  end subroutine frozen_epot_vpsl

 ! ---------------------------------------------------------
  subroutine frozen_epot_copy(this_out, this_in)
    type(frozen_epot_t), intent(inout) :: this_out
    type(frozen_epot_t), intent(in)    :: this_in
    !
    this_out%type=this_in%type
    this_out%basis=>this_in%basis
    this_out%geo=>this_in%geo
    call frozen_interpolation_copy(this_out%intrp, this_in%intrp)
    call mesh_copy(this_out%mesh, this_in%mesh)
    call loct_pointer_copy(this_out%vpsl, this_in%vpsl)
    return
  end subroutine frozen_epot_copy

  ! ---------------------------------------------------------
  subroutine frozen_epot_end(this)
    type(frozen_epot_t), intent(inout) :: this
    !
    this%type=M_ZERO
    nullify(this%basis)
    nullify(this%geo)
    call frozen_interpolation_end(this%intrp)
    call mesh_end(this%mesh)
    SAFE_DEALLOCATE_P(this%vpsl)
    return
  end subroutine frozen_epot_end

end module frozen_epot_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
