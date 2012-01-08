#include "global.h"

module frozen_states_m
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

  implicit none

  private
  public ::                &
    frozen_states_t,       &
    frozen_states_init,    &
    frozen_states_density, &
    frozen_states_copy,    &
    frozen_states_end

  type frozen_states_t
    private
    integer                        :: nst
    integer                        :: nspin
    type(frozen_basis_t), pointer  :: basis
    type(geometry_t), pointer      :: geo
    type(frozen_interpolation_t)   :: intrp
    type(mesh_t)                   :: mesh
    FLOAT, pointer, dimension(:,:) :: density
  end type frozen_states_t

contains

  ! ---------------------------------------------------------
  subroutine frozen_states_init(this, basis, geo, intrp, iunit)
    type(frozen_states_t),        intent(inout) :: this
    type(frozen_basis_t), target, intent(in)    :: basis
    type(geometry_t), target,     intent(in)    :: geo
    integer,                      intent(in)    :: intrp
    integer,                      intent(in)    :: iunit
    !
    integer :: np, gb
    !
    this%basis=>basis
    this%geo=>geo
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    call mesh_init_from_dump(this%mesh, basis%sb, basis%cv, iunit)
    read(iunit) this%nst
    read(iunit) np
    ASSERT(np==this%mesh%np_part_global)
    read(iunit) this%nspin
    SAFE_ALLOCATE(this%density(np,this%nspin))
    read(iunit) this%density
    read(iunit) gb
    ASSERT(gb==GUARD_BITS)
    call frozen_interpolation_init(this%intrp, intrp, this%mesh, this%density)
    return
  end subroutine frozen_states_init

  ! ---------------------------------------------------------
  subroutine frozen_states_density(this, x, rho)
    type(frozen_states_t), intent(in)  :: this
    FLOAT,   dimension(:), intent(in)  :: x
    FLOAT,   dimension(:), intent(out) :: rho
    !
    FLOAT, dimension(size(x)) :: y
    !
    rho=M_ZERO
    call frozen_basis_to_internal(this%basis, x, y)
    if(simul_box_in_box(this%basis%sb, this%geo, y)) &
      call frozen_interpolation(this%intrp, y, rho)
    return
  end subroutine frozen_states_density

 ! ---------------------------------------------------------
  subroutine frozen_states_copy(this_out, this_in)
    type(frozen_states_t),         intent(inout) :: this_out
    type(frozen_states_t), target, intent(in)    :: this_in
    !
    this_out%nst=this_in%nst
    this_out%nspin=this_in%nspin
    this_out%basis=>this_in%basis
    this_out%geo=>this_in%geo
    call frozen_interpolation_copy(this_out%intrp, this_in%intrp)
    call mesh_copy(this_out%mesh, this_in%mesh)
    call loct_pointer_copy(this_out%density, this_in%density)
    return
  end subroutine frozen_states_copy

  ! ---------------------------------------------------------
  subroutine frozen_states_end(this)
    type(frozen_states_t), intent(inout) :: this
    !
    this%nst=0
    this%nspin=0
    nullify(this%basis)
    nullify(this%geo)
    call frozen_interpolation_end(this%intrp)
    call mesh_end(this%mesh)
    SAFE_DEALLOCATE_P(this%density)
    return
  end subroutine frozen_states_end

end module frozen_states_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
