!! Copyright (C) 2019 N. Tancogne-Dejean, M. Oliveira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with st program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
#include "global.h"

module hamiltonian_abst_oct_m
  use batch_oct_m
  use global_oct_m
  use loct_oct_m
  use mesh_oct_m
  use messages_oct_m
  use namespace_oct_m
  use types_oct_m
  use varinfo_oct_m

  implicit none

  private

  public ::                           &
    hamiltonian_abst_t

  type, abstract :: hamiltonian_abst_t
    !> Spectral range
    FLOAT :: spectral_middle_point
    FLOAT :: spectral_half_span
  contains
    procedure(is_hermitian),              deferred :: is_hermitian
    procedure(hamiltonian_update_span),   deferred :: update_span
    procedure(dhamiltonian_apply),        deferred :: dapply
    procedure(zhamiltonian_apply),        deferred :: zapply
    procedure(dhamiltonian_magnus_apply), deferred :: dmagnus_apply
    procedure(zhamiltonian_magnus_apply), deferred :: zmagnus_apply
  end type hamiltonian_abst_t

  abstract interface
    logical function is_hermitian(hm)
      import
      class(hamiltonian_abst_t), intent(in) :: hm
    end function is_hermitian

    subroutine hamiltonian_update_span(hm, delta, emin)
      import
      class(hamiltonian_abst_t), intent(inout) :: hm
      FLOAT,                     intent(in)    :: delta
      FLOAT,                     intent(in)    :: emin
    end subroutine hamiltonian_update_span

    subroutine dhamiltonian_apply(hm, namespace, mesh, psib, hpsib, terms, set_bc)
      import
      class(hamiltonian_abst_t),   intent(in)    :: hm
      type(namespace_t),           intent(in)    :: namespace
      type(mesh_t),                intent(in)    :: mesh
      class(batch_t),      target, intent(inout) :: psib
      class(batch_t),      target, intent(inout) :: hpsib
      integer,           optional, intent(in)    :: terms
      logical,           optional, intent(in)    :: set_bc
    end subroutine dhamiltonian_apply

    subroutine zhamiltonian_apply(hm, namespace, mesh, psib, hpsib, terms, set_bc)
      import
      class(hamiltonian_abst_t),   intent(in)    :: hm
      type(namespace_t),           intent(in)    :: namespace
      type(mesh_t),                intent(in)    :: mesh
      class(batch_t),      target, intent(inout) :: psib
      class(batch_t),      target, intent(inout) :: hpsib
      integer,           optional, intent(in)    :: terms
      logical,           optional, intent(in)    :: set_bc
    end subroutine zhamiltonian_apply

    subroutine dhamiltonian_magnus_apply(hm, namespace, mesh, psib, hpsib, vmagnus, set_phase)
      import
      class(hamiltonian_abst_t),   intent(in)    :: hm
      type(namespace_t),           intent(in)    :: namespace
      type(mesh_t),                intent(in)    :: mesh
      class(batch_t),              intent(inout) :: psib
      class(batch_t),              intent(inout) :: hpsib
      FLOAT,                       intent(in)    :: vmagnus(:, :, :)
      logical,           optional, intent(in)    :: set_phase
    end subroutine dhamiltonian_magnus_apply

    subroutine zhamiltonian_magnus_apply(hm, namespace, mesh, psib, hpsib, vmagnus, set_phase)
      import
      class(hamiltonian_abst_t),   intent(in)    :: hm
      type(namespace_t),           intent(in)    :: namespace
      type(mesh_t),                intent(in)    :: mesh
      class(batch_t),              intent(inout) :: psib
      class(batch_t),              intent(inout) :: hpsib
      FLOAT,                       intent(in)    :: vmagnus(:, :, :)
      logical,           optional, intent(in)    :: set_phase
    end subroutine zhamiltonian_magnus_apply
  end interface

end module hamiltonian_abst_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
