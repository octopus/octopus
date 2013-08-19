module ispace_m

  use json_m, only: JSON_OK, json_object_t, json_init, json_set, json_get

  use space_m, only: &
    operator(==),    &
    operator(/=)

  use space_m, only: &
    space_t,         &
    space_copy,      &
    space_end

  use space_m, only:                  &
    space_init_octopus => space_init

  implicit none

  private
  public ::       &
    operator(==), &
    operator(/=)

  public ::                   &
    space_t,                  &
    space_init,               &
    space_create_data_object, &
    space_copy,               &
    space_end

  integer, parameter :: default_ndim = 3

  interface space_init
    module procedure space_init_octopus
    module procedure space_init_data_object
  end interface space_init

contains

  ! ---------------------------------------------------------
  subroutine space_init_data_object(this, config)
    type(space_t),       intent(out) :: this
    type(json_object_t), intent(in)  :: config
    !
    integer :: ndim, ierr
    !
    call json_get(config, "dimensions", ndim, ierr=ierr)
    if(ierr/=JSON_OK)ndim=default_ndim
    call space_init_octopus(this, ndim)
    return
  end subroutine space_init_data_object

  ! ---------------------------------------------------------
  subroutine space_create_data_object(this, config)
    type(space_t),       intent(in)  :: this
    type(json_object_t), intent(out) :: config
    !
    call json_init(config)
    call json_set(config, "dimensions", this%dim)
    return
  end subroutine space_create_data_object

end module ispace_m

!! Local Variables:
!! mode: f90
!! End:
