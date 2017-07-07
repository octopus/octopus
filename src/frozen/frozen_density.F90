#include "global.h"

module frozen_density_oct_m

  use base_density_oct_m
  use fio_density_oct_m
  use global_oct_m
  use json_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m

  implicit none

  private

  public ::                &
    frozen_density__acc__

contains

  ! ---------------------------------------------------------
  subroutine frozen_density__acc__(this, that, config)
    type(base_density_t), intent(inout) :: this !> frozen
    type(base_density_t), intent(in)    :: that !> fio
    type(json_object_t),  intent(in)    :: config

    type(base_density_t), pointer :: fios
    type(json_object_t),  pointer :: cnfg
    type(json_array_t),   pointer :: list
    type(simulation_t),   pointer :: sim
    type(json_array_iterator_t)   :: iter
    type(fio_density_intrpl_t)    :: intrp
    integer                       :: type, nspin, ierr
    logical                       :: accu

    PUSH_SUB(frozen_density__acc__)

    nullify(fios, cnfg, list, sim)
    call base_density_get(that, use=accu)
    if(accu)then
      call json_get(config, "type", type, ierr)
      if(ierr==JSON_OK)then
        call fio_density_intrpl_init(intrp, that, type)
      else
        call fio_density_intrpl_init(intrp, that)
      end if
      call base_density_get(this, sim)
      ASSERT(associated(sim))
      call base_density_get(this, nspin=nspin)
      ASSERT(nspin>0)
      ASSERT(nspin<3)
      call fio_density_intrpl_start(intrp, sim, nspin)
      nullify(sim)
      call json_get(config, "positions", list, ierr)
      ASSERT(ierr==JSON_OK)
      ASSERT(json_len(list)>0)
      call json_init(iter, list)
      do
        nullify(fios, cnfg)
        call json_next(iter, cnfg, ierr)
        if(ierr/=JSON_OK)exit
        call fio_density_intrpl_eval(intrp, fios, cnfg)
        ASSERT(associated(fios))
        call base_density__acc__(this, fios)
      end do
      call json_end(iter)
      nullify(fios, cnfg, list)
      call fio_density_intrpl_end(intrp)
    end if

    POP_SUB(frozen_density__acc__)
  end subroutine frozen_density__acc__

end module frozen_density_oct_m
 
!! Local Variables:
!! mode: f90
!! End:
