!! Copyright (C) 2020 M. Oliveira
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
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module interactions_factory_abst_oct_m
  use ghost_interaction_oct_m
  use global_oct_m
  use interaction_abst_oct_m
  use interaction_partner_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use system_abst_oct_m
  implicit none

  private
  public ::                        &
    interactions_factory_abst_t

  type, abstract :: interactions_factory_abst_t
  contains
    procedure :: create_interactions => interactions_factory_abst_create_interactions
    procedure(interactions_factory_abst_create), deferred :: create
  end type interactions_factory_abst_t

  abstract interface
    function interactions_factory_abst_create(this, type, partner) result(interaction)
      import :: interactions_factory_abst_t
      import interaction_partner_t
      import interaction_abst_t
      class(interactions_factory_abst_t),         intent(in)    :: this
      integer,                                    intent(in)    :: type
      class(interaction_partner_t),       target, intent(inout) :: partner
      class(interaction_abst_t),                  pointer       :: interaction
    end function interactions_factory_abst_create
  end interface

contains
  
  ! ---------------------------------------------------------------------------------------
  recursive subroutine interactions_factory_abst_create_interactions(this, system, partners)
    class(interactions_factory_abst_t),    intent(in)    :: this
    class(system_abst_t),                  intent(inout) :: system
    class(partner_list_t),         target, intent(in)    :: partners

    type(integer_iterator_t) :: interaction_iter
    integer :: interaction_type
    class(system_abst_t), pointer :: subsystem
    type(system_iterator_t) :: iter

    PUSH_SUB(interactions_factory_abst_create_interactions)

    if (debug%info) then
      write(message(1), '(a)') "Debug: -- Creating interactions for " + trim(system%namespace%get())
      call messages_info(1)
    end if

    ! Loop over all the interactions supported by the system
    call interaction_iter%start(system%supported_interactions)
    do while (interaction_iter%has_next())
      interaction_type = interaction_iter%get_next()
      call create_interaction_with_partners(this, system%namespace, partners, system%interactions, interaction_type)
    end do

    ! All systems need to be connected to make sure they remain synchronized.
    ! We enforce that be adding a ghost interaction between all systems
    call create_interaction_with_partners(this, system%namespace, partners, system%interactions)

    ! If the system is a multisystem, then we also need to create the interactions for the subsystems
    select type (system)
    class is (multisystem_t)
      call iter%start(system%list)
      do while (iter%has_next())
        subsystem => iter%get_next()
        call this%create_interactions(subsystem, partners)
      end do
    end select

    PUSH_SUB(interactions_factory_abst_create_interactions)
  end subroutine interactions_factory_abst_create_interactions

  ! ---------------------------------------------------------------------------------------
  recursive subroutine create_interaction_with_partners(this, namespace, partners, interactions, interaction_type)
    class(interactions_factory_abst_t), intent(in)    :: this
    type(namespace_t),                  intent(in)    :: namespace
    class(partner_list_t),      target, intent(in)    :: partners
    class(interaction_list_t),          intent(inout) :: interactions
    integer,                  optional, intent(in)    :: interaction_type

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner
    class(interaction_abst_t), pointer :: interaction

    ! Loop over all available partners
    call iter%start(partners)
    do while (iter%has_next())
      partner => iter%get_next()

      select type (partner)
      class is (multisystem_t)
        ! If the partner is a multitsystem, then we need to create the interactions with all the subsystems it contains
        call create_interaction_with_partners(this, namespace, partner%list, interactions, interaction_type)

      class default
        ! No self-interaction
        if (partner%namespace%get() /= namespace%get()) then

          if (debug%info) then
            write(message(1), '(a)') "Debug: ----  with " + trim(partner%namespace%get())
            call messages_info(1)
          end if

          if (present(interaction_type)) then
            ! If the partner also supports this type of interaction, then create the interaction
            if (partner%supported_interactions_as_partner%has(interaction_type)) then
              interaction => this%create(interaction_type, partner)
              call interactions%add(interaction)
            end if
          else
            ! Create a ghost interaction if no interaction type was given
            interaction => ghost_interaction_t(partner)
            call interactions%add(interaction)
          end if
        end if

      end select
    end do

  end subroutine create_interaction_with_partners

end module interactions_factory_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
