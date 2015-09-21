#include "global.h"

module fio_grid_m

  use boundaries_m
  use curvilinear_m
  use derivatives_m
  use double_grid_m
  use fio_curvilinear_m
  use fio_mesh_m
  use fio_simul_box_m
  use geometry_m
  use global_m
  use grid_m
  use json_m
  use kinds_m
  use mesh_m
  use messages_m
  use multigrid_m
  use mpi_m
  use profiling_m
  use simul_box_m
  use stencil_m
  use transfer_table_m

  implicit none

  private

  public ::         &
    fio_grid_init,  &
    fio_grid_start, &
    fio_grid_stop,  &
    fio_grid_copy,  &
    fio_grid_end

contains
  
  ! ---------------------------------------------------------
  subroutine transfer_table_nullify(this)
    type(transfer_table_t), intent(out) :: this

    PUSH_SUB(transfer_table_nullify)

    this%n_coarse = 0
    this%n_fine = 0
    this%n_fine1 = 0
    this%n_fine2 = 0
    this%n_fine4 = 0
    this%n_fine8 = 0
    nullify(this%to_coarse, this%to_fine1, this%to_fine2, this%to_fine4, this%to_fine8, this%fine_i)

    POP_SUB(transfer_table_nullify)
  end subroutine transfer_table_nullify

  ! ---------------------------------------------------------
  subroutine multigrid_level_nullify(this)
    type(multigrid_level_t), intent(out) :: this

    PUSH_SUB(multigrid_level_nullify)

    call transfer_table_nullify(this%tt)
    nullify(this%mesh, this%der)

    POP_SUB(multigrid_level_nullify)
  end subroutine multigrid_level_nullify

  ! ---------------------------------------------------------
  subroutine boundaries_nullify(this)
    type(boundaries_t), intent(out) :: this

    PUSH_SUB(boundaries_nullify)

    nullify(this%mesh, this%per_points)
    this%nper = 0
#ifdef HAVE_MPI
    nullify(this%per_send, this%per_recv, this%nsend, this%nrecv)
#endif
    !this%buff_per_points
    !this%buff_per_send
    !this%buff_per_recv
    !this%buff_nsend
    !this%buff_nrecv

    POP_SUB(boundaries_nullify)
  end subroutine boundaries_nullify

  ! ---------------------------------------------------------
  subroutine derivatives_nullify(this)
    type(derivatives_t), intent(out) :: this

    PUSH_SUB(derivatives_nullify)

    call boundaries_nullify(this%boundaries)
    nullify(this%mesh, this%op, this%lapl, this%grad)
    nullify(this%finer, this%coarser, this%to_finer, this%to_coarser)
    this%dim = 0
    this%order = 0
    this%stencil_type = 0
    this%masses = M_ZERO
    this%np_zero_bc = 0
    this%lapl_cutoff = M_ZERO
    this%n_ghost = 0
#if defined(HAVE_MPI)
    this%comm_method = 0
#endif

    POP_SUB(derivatives_nullify)
  end subroutine derivatives_nullify

  !-------------------------------------------------------  
  subroutine stencil_nullify(this)
    type(stencil_t), intent(out) :: this

    PUSH_SUB(stencil_nullify)

    this%center = -1
    this%size = 0
    nullify(this%points)

    POP_SUB(stencil_nullify)
  end subroutine stencil_nullify

  ! ---------------------------------------------------------
  subroutine fio_grid_init(this, geo, config)
    type(grid_t), target, intent(out) :: this
    type(geometry_t),     intent(in)  :: geo
    type(json_object_t),  intent(in)  :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_grid_init)

    nullify(cnfg)
    call json_get(config, "simul_box", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_simul_box_init(this%sb, geo, cnfg)
    nullify(cnfg)
    call json_get(config, "curvilinear", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_curvilinear_init(this%cv, this%sb, geo, cnfg)
    nullify(cnfg)
    call multigrid_level_nullify(this%fine)
    call derivatives_nullify(this%der)
    call double_grid_nullify(this%dgrid)
    call stencil_nullify(this%stencil)
    this%have_fine_mesh = .false.
    this%fine%mesh => this%mesh
    this%fine%der => this%der
    nullify(this%mgrid)

    POP_SUB(fio_grid_init)
  end subroutine fio_grid_init

  ! ---------------------------------------------------------
  subroutine fio_grid_start(this, mpi_grp, config)
    type(grid_t),        intent(inout) :: this
    type(mpi_grp_t),     intent(in)    :: mpi_grp
    type(json_object_t), intent(in)    :: config

    type(json_object_t), pointer :: cnfg
    integer                      :: ierr

    PUSH_SUB(fio_grid_start)

    nullify(cnfg)
    call json_get(config, "mesh", cnfg, ierr)
    ASSERT(ierr==JSON_OK)
    call fio_mesh_init(this%mesh, this%sb, this%cv, mpi_grp, cnfg)
    nullify(cnfg)

    POP_SUB(fio_grid_start)
  end subroutine fio_grid_start

  ! ---------------------------------------------------------
  subroutine fio_grid_stop(this)
    type(grid_t),    intent(inout) :: this

    PUSH_SUB(fio_grid_stop)

    call fio_mesh_end(this%mesh)

    POP_SUB(fio_grid_stop)
  end subroutine fio_grid_stop

  ! ---------------------------------------------------------
  subroutine fio_grid_copy(this, that)
    type(grid_t), target, intent(inout) :: this
    type(grid_t),         intent(in)    :: that

    PUSH_SUB(fio_grid_copy)

    ASSERT(.not.that%have_fine_mesh)
    call fio_simul_box_copy(this%sb, that%sb)
    call fio_mesh_copy(this%mesh, that%mesh)
    call fio_curvilinear_copy(this%cv, that%cv)
    call multigrid_level_nullify(this%fine)
    call derivatives_nullify(this%der)
    call double_grid_nullify(this%dgrid)
    call stencil_nullify(this%stencil)
    this%have_fine_mesh = .false.
    this%fine%mesh => this%mesh
    this%fine%der => this%der
    nullify(this%mgrid)

    POP_SUB(fio_grid_copy)
  end subroutine fio_grid_copy

  !-------------------------------------------------------------------
  subroutine fio_grid_end(this)
    type(grid_t), intent(inout) :: this

    PUSH_SUB(fio_grid_end)

    call fio_mesh_end(this%mesh)
    call fio_curvilinear_end(this%cv)
    call fio_simul_box_end(this%sb)
    call multigrid_level_nullify(this%fine)
    call derivatives_nullify(this%der)
    call double_grid_nullify(this%dgrid)
    call stencil_nullify(this%stencil)
    this%have_fine_mesh = .false.
    nullify(this%fine%mesh, this%fine%der, this%mgrid)

    POP_SUB(fio_grid_end)
  end subroutine fio_grid_end

end module fio_grid_m

!! Local Variables:
!! mode: f90
!! End:
