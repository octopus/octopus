#include "global.h"

module fio_config_m

  implicit none

  private
  public ::              &
    input_mesh_dir,      &
    input_mesh_file,     &
    input_external_file

  character(len=*), parameter :: input_mesh_dir      = "./restart/"//GS_DIR
  character(len=*), parameter :: input_mesh_file     = "mesh"
  character(len=*), parameter :: input_external_file = "v0.obf"

end module fio_config_m

!! Local Variables:
!! mode: f90
!! End:
