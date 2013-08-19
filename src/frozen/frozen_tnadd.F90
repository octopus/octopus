#include "global.h"

module frozen_tnadd_m

  use global_m
  use messages_m
  use profiling_m

  use json_m,          only: JSON_OK, json_object_t, json_get
  use kinds_m,         only: wp
  use mesh_function_m, only: dmf_dotp, dmf_integrate
  use storage_m,       only: storage_t, storage_get_size, storage_init, &
    storage_start, storage_update, storage_get_storage, storage_end

  use functional_m, only: &
    XC_NONE,              &
    XC_KINETIC

  use functional_m, only:                &
    functional_t,                        &
    functional_init,                     &
    functional_start,                    &
    functional_get_kind,                 &
    functional_get_energy,               &
    functional_get_energy_and_potential, &
    functional_get_potential,            &
    functional_copy,                     &
    functional_end

  use frozen_density_m, only:                          &
    density_t           => frozen_density_t,           &
    density_get_nspin   => frozen_density_get_nspin,   &
    density_get_density => frozen_density_get_density

  use frozen_grid_m, only:   &
    grid_t !=> frozen_grid_t

  use frozen_simulation_m, only:             &
    simulation_t   !=> frozen_simulation_t,   &

  use frozen_simulation_m, only:             &
    simulation_get => frozen_simulation_get

  use frozen_system_m, only:     &
    system_t => frozen_system_t

  use frozen_potential_m, only:                                &
    potential_t             => frozen_potential_t,             &
    potential_init          => frozen_potential_init,          &
    potential_start         => frozen_potential_start,         &
    potential_update        => frozen_potential_update,        &
    potential_get           => frozen_potential_get,           &
    potential_get_size      => frozen_potential_get_size,      &
    potential_get_nspin     => frozen_potential_get_nspin,     &
    potential_set_energy    => frozen_potential_set_energy,    &
    potential_get_energy    => frozen_potential_get_energy,    &
    potential_get_potential => frozen_potential_get_potential, &
    potential_copy          => frozen_potential_copy,          &
    potential_end           => frozen_potential_end

  implicit none

  private
  public ::                     &
    frozen_tnadd_init,          &
    frozen_tnadd_start,         &
    frozen_tnadd_update,        &
    frozen_tnadd_calc,          &
    frozen_tnadd_get_energy,    &
    frozen_tnadd_get_potential, &
    frozen_tnadd_copy,          &
    frozen_tnadd_end


  type, public :: frozen_tnadd_t
    private
    type(density_t), pointer :: density =>null()
    real(kind=wp)            :: factor  = 0.0_wp
    real(kind=wp)            :: energy  = 0.0_wp
    type(functional_t)       :: funct
    type(potential_t)        :: potential
  end type frozen_tnadd_t

contains
  
  ! ---------------------------------------------------------
  subroutine frozen_tnadd_init(this, sys, config)
    type(frozen_tnadd_t), intent(out) :: this
    type(system_t),       intent(in)  :: sys
    type(json_object_t),  intent(in)  :: config
    !
    integer :: ierr
    !
    this%density=>null()
    call json_get(config, "factor", this%factor, ierr)
    if(ierr/=JSON_OK)this%factor=1.0_wp
    this%energy=0.0_wp
    call functional_init(this%funct, config)
    call potential_init(this%potential, sys, config)
    return
  end subroutine frozen_tnadd_init

  ! ---------------------------------------------------------
  subroutine frozen_tnadd_start(this, sim)
    type(frozen_tnadd_t), intent(inout) :: this
    type(simulation_t),   intent(in)    :: sim
    !
    real(kind=wp), dimension(:,:), pointer :: density
    !
    call functional_start(this%funct, sim)
    if(functional_get_kind(this%funct)==XC_KINETIC)then
      call potential_start(this%potential, sim)
      call potential_get(this%potential, this%density)
      ASSERT(associated(this%density))
      call density_get_density(this%density, density)
      ASSERT(associated(density))
      this%energy=functional_get_energy(this%funct, density)
    else
      call frozen_tnadd_end(this)
    end if
    return
  end subroutine frozen_tnadd_start

 ! ---------------------------------------------------------
  subroutine frozen_tnadd_update(this)
    type(frozen_tnadd_t), intent(inout) :: this
    !
    type(simulation_t), pointer :: sim
    !
    nullify(sim)
    call potential_get(this%potential, sim)
    ASSERT(associated(sim))
    call potential_update(this%potential)
    return
  end subroutine frozen_tnadd_update

  ! ---------------------------------------------------------
  subroutine frozen_tnadd_calc(this, density, calc_energy)
    type(frozen_tnadd_t),          intent(inout) :: this
    real(kind=wp), dimension(:,:), intent(in)    :: density
    logical,             optional, intent(in)    :: calc_energy
    !
    real(kind=wp), dimension(size(density,dim=1),size(density,dim=2)) :: rho_1
    real(kind=wp), dimension(:,:),allocatable :: pot_1
    real(kind=wp), dimension(:,:), pointer :: pot_t, rho_f
    type(storage_t)                        :: srho, spot
    real(kind=wp)                          :: et, e1
    logical                                :: calc_enrg
    integer                                :: i !, ns, nd
    !
    print *, "***: frozen_tnadd_calc: start"
    calc_enrg=.true.
    if(present(calc_energy))calc_enrg=calc_energy
    et=0.0_wp
    nullify(rho_f, pot_t) !, rho_1, pot_1, sim)
    call potential_get_potential(this%potential, pot_t)
    ASSERT(associated(pot_t))
    !print *, "***: frozen_tnadd_calc: 1"
    if(calc_enrg)then
      call functional_get_energy_and_potential(this%funct, density, et, pot_t)
    else
      call functional_get_potential(this%funct, density, pot_t)
    end if
    call density_get_density(this%density, rho_f)
    ASSERT(associated(rho_f))
    do i = 1, size(density, dim=1)
      rho_1(i,:) = density(i,:) - rho_f(i,:)
    end do
    e1=0.0_wp
    SAFE_ALLOCATE(pot_1(size(pot_t,dim=1),size(pot_t,dim=2)))
    if(calc_enrg)then
      call functional_get_energy_and_potential(this%funct, rho_1, e1, pot_1)
      call potential_set_energy(this%potential, this%factor*(et-e1-this%energy))
    else
      call functional_get_potential(this%funct, rho_1, pot_1)
      call potential_set_energy(this%potential, 0.0_wp)
    end if
    do i = 1, size(pot_t, dim=1)
      pot_t(i,:)=this%factor*(pot_t(i,:)-pot_1(i,:))
    end do
    call frozen_tnadd_update(this)
    SAFE_DEALLOCATE_A(pot_1)
    nullify(rho_f, pot_t)
    !print *, "***: frozen_tnadd_calc: end"
    return
  end subroutine frozen_tnadd_calc

  ! ---------------------------------------------------------
  elemental function frozen_tnadd_get_size(this) result(np)
    type(frozen_tnadd_t), intent(in) :: this
    !
    integer :: np
    !
    np=potential_get_size(this%potential)
    return
  end function frozen_tnadd_get_size

  ! ---------------------------------------------------------
  elemental function frozen_tnadd_get_nspin(this) result(that)
    type(frozen_tnadd_t), intent(in) :: this
    !
    integer :: that
    !
    that=potential_get_nspin(this%potential)
    return
  end function frozen_tnadd_get_nspin

  ! ---------------------------------------------------------
  elemental function frozen_tnadd_get_energy(this) result(that)
    type(frozen_tnadd_t), intent(in) :: this
    !
    real(kind=wp) :: that
    !
    that=potential_get_energy(this%potential)
    return
  end function frozen_tnadd_get_energy
  
  ! ---------------------------------------------------------
  subroutine frozen_tnadd_get_potential(this, that)
    type(frozen_tnadd_t),           intent(in) :: this
    real(kind=wp), dimension(:,:), pointer     :: that
    !
    call potential_get_potential(this%potential, that)
    return
  end subroutine frozen_tnadd_get_potential

  ! ---------------------------------------------------------
  subroutine frozen_tnadd_copy(this, that)
    type(frozen_tnadd_t), intent(out) :: this
    type(frozen_tnadd_t), intent(in)  :: that
    !
    this%density=>that%density
    this%factor=that%factor
    this%energy=that%energy
    call potential_copy(this%potential, that%potential)
    call functional_copy(this%funct, that%funct)
    return
  end subroutine frozen_tnadd_copy

  ! ---------------------------------------------------------
  subroutine frozen_tnadd_end(this)
    type(frozen_tnadd_t), intent(inout) :: this
    !
    call potential_end(this%potential)
    call functional_end(this%funct)
    this%energy=0.0_wp
    this%factor=0.0_wp
    this%density=>null()
    return
  end subroutine frozen_tnadd_end

end module frozen_tnadd_m

!! Local Variables:
!! mode: f90
!! End:

