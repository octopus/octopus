#include "global.h"

module fd_geometry_m
  use fd_basis_m, only: fd_basis_t, fd_basis_to_external, fd_basis_rotate_to_external
  use geometry_m, only: geometry_t, atom_classical_t
  use global_m
  use messages_m
  use profiling_m
  use species_m, only: species_zval
  
  implicit none

  private
  public ::                 &
    fd_atoms_get,           &
    fd_atoms_get_classical, &
    fd_atoms_add
    
contains

  ! ---------------------------------------------------------
  subroutine fd_atoms_get(this, basis, atoms)
    type(geometry_t),                              intent(in)    :: this
    type(fd_basis_t),                              intent(in)    :: basis
    type(atom_classical_t), pointer, dimension(:), intent(inout) :: atoms
    !
    integer :: i
    !
    SAFE_ALLOCATE(atoms(this%natoms))
    do i = 1, this%natoms
      atoms(i)%label = this%atom(i)%label
      call fd_basis_to_external(basis, this%atom(i)%x, atoms(i)%x)
      atoms(i)%v = 0.0 !this%atom(i)%v
      !call fd_basis_rotate_to_external(basis, atoms(i)%v)
      atoms(i)%f = 0.0 !this%atom(i)%f
      !call fd_basis_rotate_to_external(basis, atoms(i)%f)
      atoms(i)%charge = species_zval(this%atom(i)%spec)
    end do
    return
  end subroutine fd_atoms_get
  
  ! ---------------------------------------------------------
  subroutine fd_atoms_get_classical(this, basis, atoms)
    type(geometry_t),                              intent(in)    :: this
    type(fd_basis_t),                              intent(in)    :: basis
    type(atom_classical_t), pointer, dimension(:), intent(inout) :: atoms
    !
    integer :: i
    !
    SAFE_ALLOCATE(atoms(this%ncatoms))
    do i = 1, this%ncatoms
      atoms(i)%label = this%catom(i)%label
      call fd_basis_to_external(basis, this%catom(i)%x, atoms(i)%x)
      atoms(i)%v = 0.0 !this%catom(i)%v
      !call fd_basis_rotate_to_external(basis, atoms(i)%v)
      atoms(i)%f = 0.0 !this%catom(i)%f
      !call fd_basis_rotate_to_external(basis, atoms(i)%f)
      atoms(i)%charge = this%catom(i)%charge
    end do
    return
  end subroutine fd_atoms_get_classical

  ! ---------------------------------------------------------
  subroutine fd_atoms_add(n_out, this_out, n_in, this_in)
    integer,                                       intent(inout) :: n_out
    type(atom_classical_t), pointer, dimension(:), intent(inout) :: this_out
    integer,                                       intent(in)    :: n_in
    type(atom_classical_t),          dimension(:), intent(in)    :: this_in
    !
    type(atom_classical_t), pointer, dimension(:) :: this
    !
    if(n_in>0)then
      if(n_out>0)then
        this=>this_out
        nullify(this_out)
        SAFE_ALLOCATE(this_out(n_out))
        this_out(1:n_out)=this
        SAFE_DEALLOCATE_P(this)
        this_out(n_out+1:n_in+n_out)=this_in
        n_out=n_out+n_in
      else
        SAFE_ALLOCATE(this_out(n_in))
        this_out=this_in
        n_out=n_in
      end if
    end if
    return
  end subroutine fd_atoms_add

end module fd_geometry_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
