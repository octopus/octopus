#include "global.h"

module fd_states_m
  use fd_basis_m, only: fd_basis_t, fd_basis_to_internal
  use fd_interpolation_m, only: fd_interpolation_t, &
    fd_interpolation_init, fd_interpolation, fd_interpolation_copy, fd_interpolation_end
  use geometry_m, only: geometry_t
  use global_m
  use loct_m, only: loct_pointer_copy
  use mesh_m, only: mesh_t, mesh_init_from_dump, mesh_copy, mesh_end
  use messages_m
  use profiling_m
  use simul_box_m, only: simul_box_in_box

  implicit none

  private
  public ::            &
    fd_states_t,       &
    fd_states_init,    &
    fd_states_density, &
    fd_states_copy,    &
    fd_states_end

  type fd_states_t
    private
    integer                        :: nst
    integer                        :: nspin
    type(fd_basis_t), pointer      :: basis
    type(geometry_t), pointer      :: geo
    type(fd_interpolation_t)       :: intrp
    type(mesh_t)                   :: mesh
    FLOAT, pointer, dimension(:,:) :: density
  end type fd_states_t

contains

  ! ---------------------------------------------------------
  subroutine fd_states_init(this, basis, geo, intrp, iunit)
    use io_m
    type(fd_states_t),        intent(inout) :: this
    type(fd_basis_t), target, intent(in)    :: basis
    type(geometry_t), target, intent(in)    :: geo
    integer,                  intent(in)    :: intrp
    integer,                  intent(in)    :: iunit
    !
    integer :: np, gb, i, junit
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
    junit=io_open("rho_fd.plt", action="write", is_tmp=.true.)
    do i = 1, this%mesh%np_part_global
      if((this%mesh%x(i,2)==0.0).and.(this%mesh%x(i,3)==0.0))then
        write(junit,*) this%mesh%x(i,1), this%density(i,:)
      end if
    end do
    call io_close(junit)
    call fd_interpolation_init(this%intrp, intrp, this%mesh, this%density)
    return
  end subroutine fd_states_init

  ! ---------------------------------------------------------
  subroutine fd_states_density(this, x, rho)
    type(fd_states_t),   intent(in)  :: this
    FLOAT, dimension(:), intent(in)  :: x
    FLOAT, dimension(:), intent(out) :: rho
    !
    FLOAT, dimension(size(x)) :: y
    !
    rho=M_ZERO
    call fd_basis_to_internal(this%basis, x, y)
    if(simul_box_in_box(this%basis%sb, this%geo, y)) &
      call fd_interpolation(this%intrp, y, rho)
    !if((y(2)==0.0).and.(y(3)==0.0))print *, "###: " y, rho
    !print *, "###: ", x-y
    return
  end subroutine fd_states_density

 ! ---------------------------------------------------------
  subroutine fd_states_copy(this_out, this_in)
    type(fd_states_t),         intent(inout) :: this_out
    type(fd_states_t), target, intent(in)    :: this_in
    !
    this_out%nst=this_in%nst
    this_out%nspin=this_in%nspin
    this_out%basis=>this_in%basis
    this_out%geo=>this_in%geo
    call fd_interpolation_copy(this_out%intrp, this_in%intrp)
    call mesh_copy(this_out%mesh, this_in%mesh)
    call loct_pointer_copy(this_out%density, this_in%density)
    return
  end subroutine fd_states_copy

  ! ---------------------------------------------------------
  subroutine fd_states_end(this)
    type(fd_states_t), intent(inout) :: this
    !
    this%nst=0
    this%nspin=0
    nullify(this%basis)
    nullify(this%geo)
    call fd_interpolation_end(this%intrp)
    call mesh_end(this%mesh)
    SAFE_DEALLOCATE_P(this%density)
    return
  end subroutine fd_states_end

end module fd_states_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
