#include "global.h"

module frozen_external_oct_m

  use atom_oct_m
  use base_geometry_oct_m
  use base_potential_oct_m
  use base_system_oct_m
  use global_oct_m
  use intrpl_oct_m
  use json_oct_m
  use kinds_oct_m
  use mesh_oct_m
  use messages_oct_m
  use profiling_oct_m
  use simulation_oct_m
  use species_oct_m
  use storage_oct_m

  implicit none

  private

  public ::                 &
    frozen_external__acc__

contains

  ! ---------------------------------------------------------
  pure function frozen_external_calc(x, y, c) result(v)
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp), dimension(:), intent(in)  :: y
    real(kind=wp),               intent(in)  :: c

    real(kind=wp) :: v

    real(kind=wp) :: r
    integer       :: n

    n = min(size(x), size(y))
    r = sqrt(sum((x(1:n)-y(1:n))**2))
    if(r<r_small) r = r_small
    v = -c / r

  end function frozen_external_calc

  ! ---------------------------------------------------------
  subroutine frozen_external_classical(this, x, v)
    type(base_potential_t),      intent(in)  :: this
    real(kind=wp), dimension(:), intent(in)  :: x
    real(kind=wp),               intent(out) :: v

    type(base_system_t),    pointer :: sys
    type(base_geometry_t),  pointer :: geom
    type(atom_t),           pointer :: atom
    !type(atom_classical_t), pointer :: catom
    type(base_geometry_iterator_t)  :: iter

    PUSH_SUB(frozen_external_classical)

    v = 0.0_wp
    nullify(sys, geom)
    call base_potential_get(this, sys)
    ASSERT(associated(sys))
    call base_system_get(sys, geom)
    ASSERT(associated(geom))
    call base_geometry_init(iter, geom)
    do
      nullify(atom)
      call base_geometry_next(iter, atom)
      if(.not.associated(atom))exit
      ASSERT(associated(atom%species))
      v = v + frozen_external_calc(x, atom%x, species_zval(atom%species))
    end do
    nullify(atom)
    !do
    !  nullify(catom)
    !  call frozen_geom_next(iter, catom)
    !  if(.not.associated(catom))exit
    !  v=v+frozen_external_calc(x, catom%x, catom%charge)
    !end do
    call base_geometry_end(iter)
    !nullify(catom)
    nullify(sys, geom)

    POP_SUB(frozen_external_classical)
  end subroutine frozen_external_classical

  ! ---------------------------------------------------------
  subroutine  frozen_external__acc__intrpl(this, that, intrpl)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    type(intrpl_t),         intent(in)    :: intrpl

    real(kind=wp), dimension(:),     pointer :: potn
    type(simulation_t),              pointer :: sim
    type(mesh_t),                    pointer :: mesh
    real(kind=wp)                            :: pot
    integer                                  :: indx, np, ndim

    PUSH_SUB(frozen_external__acc__intrpl)

    nullify(potn, sim, mesh)
    call base_potential_get(this, sim)
    ASSERT(associated(sim))
    call simulation_get(sim, ndim=ndim)
    ASSERT(ndim>0)
    call simulation_get(sim, mesh)
    ASSERT(associated(mesh))
    nullify(sim)
    call base_potential_get(this, potn)
    ASSERT(associated(potn))
    call base_potential_get(this, size=np)
    do indx = 1, np
      call intrpl_eval(intrpl, mesh%x(indx,1:ndim), pot, func)
      potn(indx) = potn(indx) + pot
    end do
    nullify(potn, mesh)

    POP_SUB(frozen_external__acc__intrpl)

  contains
    
    function func(x) result(val)
      real(kind=wp), dimension(:), intent(in) :: x
      
      real(kind=wp) :: val

      PUSH_SUB(frozen_external_intrpl_eval.func)

      call frozen_external_classical(that, x, val)

      POP_SUB(frozen_external_intrpl_eval.func)
    end function func

  end subroutine frozen_external__acc__intrpl

  ! ---------------------------------------------------------
  subroutine frozen_external__acc__(this, that, config)
    type(base_potential_t), intent(inout) :: this
    type(base_potential_t), intent(in)    :: that
    type(json_object_t),    intent(in)    :: config

    type(json_object_t), pointer :: cnfg
    type(json_array_t),  pointer :: list
    type(storage_t),     pointer :: data
    type(json_array_iterator_t)  :: iter
    type(intrpl_t)               :: intrp
    integer                      :: type, ierr

    PUSH_SUB(frozen_external__acc__)

    nullify(cnfg, list, data)
    call base_potential_get(that, data)
    ASSERT(associated(data))
    call json_get(config, "interpolation", type, ierr)
    if(ierr==JSON_OK)then
      call intrpl_init(intrp, data, type)
    else
      call intrpl_init(intrp, data)
    end if
    nullify(data)
    call json_get(config, "positions", list, ierr)
    ASSERT(ierr==JSON_OK)
    ASSERT(json_len(list)>0)
    call json_init(iter, list)
    do
      nullify(cnfg)
      call json_next(iter, cnfg, ierr)
      if(ierr/=JSON_OK)exit
      call intrpl_attach(intrp, cnfg)
      call frozen_external__acc__intrpl(this, that, intrp)
      call intrpl_detach(intrp)
    end do
    call json_end(iter)
    nullify(cnfg, list)
    call intrpl_end(intrp)

    POP_SUB(frozen_external__acc__)
  end subroutine frozen_external__acc__

end module frozen_external_oct_m

!! Local Variables:
!! mode: f90
!! End:
