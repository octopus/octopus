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
  use interaction_oct_m
  use interaction_partner_oct_m
  use linked_list_oct_m
  use messages_oct_m
  use multisystem_oct_m
  use namespace_oct_m
  use parser_oct_m
  use system_oct_m
  use varinfo_oct_m
  implicit none

  private
  public ::                        &
    interactions_factory_abst_t

  integer, parameter, public ::   &
    NO_PARTNERS   = -1,            &
    ALL_PARTNERS  = -2,            &
    ONLY_PARTNERS = -3,            &
    ALL_EXCEPT    = -4

  type, abstract :: interactions_factory_abst_t
  contains
    procedure :: create_interactions => interactions_factory_abst_create_interactions
    procedure(interactions_factory_abst_create), deferred :: create
    procedure(interactions_factory_abst_default_mode), deferred :: default_mode
    procedure(interactions_factory_abst_block_name), deferred :: block_name
  end type interactions_factory_abst_t

  abstract interface
    function interactions_factory_abst_create(this, type, partner) result(interaction)
      import :: interactions_factory_abst_t
      import interaction_partner_t
      import interaction_t
      class(interactions_factory_abst_t),         intent(in)    :: this
      integer,                                    intent(in)    :: type
      class(interaction_partner_t),       target, intent(inout) :: partner
      class(interaction_t),                       pointer       :: interaction
    end function interactions_factory_abst_create

    integer function interactions_factory_abst_default_mode(this, type)
      import :: interactions_factory_abst_t
      class(interactions_factory_abst_t), intent(in)    :: this
      integer,                            intent(in)    :: type
    end function interactions_factory_abst_default_mode

    character(len=80) function interactions_factory_abst_block_name(this)
      import :: interactions_factory_abst_t
      class(interactions_factory_abst_t), intent(in)    :: this
    end function interactions_factory_abst_block_name
  end interface

contains
  
  ! ---------------------------------------------------------------------------------------
  recursive subroutine interactions_factory_abst_create_interactions(this, system, available_partners)
    class(interactions_factory_abst_t),    intent(in)    :: this
    class(system_t),                       intent(inout) :: system
    class(partner_list_t),         target, intent(in)    :: available_partners

    type(integer_list_t) :: interactions_to_create
    type(integer_iterator_t) :: interaction_iter
    integer :: interaction_type
    type(partner_list_t) :: partners, partners_flat_list
    type(partner_iterator_t) :: partner_iter
    class(interaction_partner_t), pointer :: partner
    type(system_iterator_t) :: iter
    class(system_t), pointer :: subsystem

    integer :: il, ic, mode
    type(block_t) :: blk
    character(len=MAX_NAMESPACE_LEN) :: input_name

    PUSH_SUB(interactions_factory_abst_create_interactions)

    if (debug%info) then
      write(message(1), '(a)') "Debug: -- Creating interactions for " + trim(system%namespace%get())
      call messages_info(1)
    end if

    ! Make a copy of the interactions list so that we can modify it
    interactions_to_create = system%supported_interactions

    ! Get the list of partners as a flat list
    call flatten_partner_list(available_partners, partners_flat_list)

    ! Parse input. The variable name and description should be given by the
    ! factory, as different factories might have different options.
    if (parse_block(system%namespace, this%block_name(), blk) == 0) then

      ! Loop over all interactions specified in the input file
      do il = 0, parse_block_n(blk) - 1
        ! Read the interaction type (first column)
        call parse_block_integer(blk, il, 0, interaction_type)

        ! Sanity check: the interaction type must be known and must not be mistaken for an interaction mode
        if (.not. varinfo_valid_option(this%block_name(), interaction_type) .or. &
          any(interaction_type == (/ALL_PARTNERS, ONLY_PARTNERS, NO_PARTNERS, ALL_EXCEPT/))) then
          call messages_input_error(system%namespace, this%block_name(), details="Unknown interaction type", row=il, column=0)
        end if

        ! Ignore interactions that are not supported by this system
        if (.not. interactions_to_create%has(interaction_type)) cycle

        ! Read how this interaction should be treated (second column)
        call parse_block_integer(blk, il, 1, mode)

        ! Create list of partners for this interaction taking into account the selected mode
        select case (mode)
        case (ALL_PARTNERS)
          ! Use all available partners
          partners = partners_flat_list
        case (NO_PARTNERS)
          ! No partners for this interaction
          call partners%empty()
        case (ONLY_PARTNERS)
          ! Start with an empty list. We will add only the select partners bellow
          call partners%empty()
        case (ALL_EXCEPT)
          ! Start with full list. We will remove the select partners bellow
          partners = partners_flat_list
        case default
          call messages_input_error(system%namespace, this%block_name(), "Unknown interaction mode", row=il, column=1)
        end select

        if (mode == ONLY_PARTNERS .or. mode == ALL_EXCEPT) then
          ! In these two cases we need to read the names of the selected
          ! partners (remaining columns) and handled them appropriatly
          do ic = 2, parse_block_cols(blk, il) - 1
            call parse_block_string(blk, il, ic, input_name)

            ! Loop over available partners and either add them or remove them
            ! from the list depending on the selected mode
            call partner_iter%start(partners_flat_list)
            do while (partner_iter%has_next())
              partner => partner_iter%get_next()
              if (partner%namespace%is_contained_in(input_name)) then
                select case (mode)
                case (ONLY_PARTNERS)
                  call partners%add(partner)
                case (ALL_EXCEPT)
                  call partners%delete(partner)
                end select
              end if
            end do
          end do

        end if

        ! Now actually create the interactions for the selected partners
        call create_interaction_with_partners(this, system%namespace, partners, system%interactions, interaction_type)

        ! Remove this interaction type from the list, as it has just been handled
        call interactions_to_create%delete(interaction_type)
      end do
      call parse_block_end(blk)
    end if

    ! Loop over all the remaining interactions supported by the system
    call interaction_iter%start(interactions_to_create)
    do while (interaction_iter%has_next())
      interaction_type = interaction_iter%get_next()

      ! Check what is the default mode for this interaction type (all or none)
      select case (this%default_mode(interaction_type))
      case (ALL_PARTNERS)
        partners = partners_flat_list
      case (NO_PARTNERS)
        call partners%empty()
      case default
        message(1) = "Default interaction mode can only be all_partners or no_partners."
        call messages_fatal(1, namespace=system%namespace)
      end select

      call create_interaction_with_partners(this, system%namespace, partners, system%interactions, interaction_type)
    end do

    ! All systems need to be connected to make sure they remain synchronized.
    ! We enforce that be adding a ghost interaction between all systems
    call create_interaction_with_partners(this, system%namespace, partners_flat_list, system%interactions)

    ! If the system is a multisystem, then we also need to create the interactions for the subsystems
    select type (system)
    class is (multisystem_t)
      call iter%start(system%list)
      do while (iter%has_next())
        subsystem => iter%get_next()
        call this%create_interactions(subsystem, available_partners)
      end do
    end select

    POP_SUB(interactions_factory_abst_create_interactions)
  contains

    recursive subroutine flatten_partner_list(partners, flat_list)
      class(partner_list_t),  intent(in)    :: partners
      class(partner_list_t),  intent(inout) :: flat_list

      class(interaction_partner_t), pointer :: partner
      type(partner_iterator_t) :: iterator

      PUSH_SUB(interactions_factory_abst_create_interactions.flatten_partner_list)

      call iterator%start(partners)
      do while (iterator%has_next())
        partner => iterator%get_next()

        call flat_list%add(partner)

        select type (partner)
        class is (multisystem_t)
          ! Also incude the subsystems of a multisystem
          call flatten_partner_list(partner%list, flat_list)
        end select

      end do

      POP_SUB(interactions_factory_abst_create_interactions.flatten_partner_list)
    end subroutine flatten_partner_list

  end subroutine interactions_factory_abst_create_interactions

  ! ---------------------------------------------------------------------------------------
  subroutine create_interaction_with_partners(this, namespace, partners, interactions, interaction_type)
    class(interactions_factory_abst_t), intent(in)    :: this
    type(namespace_t),                  intent(in)    :: namespace
    class(partner_list_t),      target, intent(in)    :: partners
    class(interaction_list_t),          intent(inout) :: interactions
    integer,                  optional, intent(in)    :: interaction_type

    type(partner_iterator_t) :: iter
    class(interaction_partner_t), pointer :: partner
    class(interaction_t), pointer :: interaction

    logical :: interaction_used

    PUSH_SUB(create_interaction_with_partners)

    ! Loop over all available partners
    call iter%start(partners)
    do while (iter%has_next())
      partner => iter%get_next()

      interaction_used = .false.

      ! No self-interaction
      if (partner%namespace%get() /= namespace%get()) then

        if (present(interaction_type)) then
          ! If the partner also supports this type of interaction, then create the interaction
          if (partner%supported_interactions_as_partner%has(interaction_type)) then
            interaction => this%create(interaction_type, partner)
            call interactions%add(interaction)
            interaction_used = .true.
          end if
        else
          ! Create a ghost interaction if no interaction type was given
          interaction => ghost_interaction_t(partner)
          call interactions%add(interaction)
          interaction_used = .true.
        end if

        if (debug%info .and. interaction_used) then
          if( associated(interaction) ) then
            write(message(1), '(a)') "Debug: ----  " + trim(interaction%label) + " with " + trim(partner%namespace%get())
          else
            write(message(1), '(a)') "Debug: ----  interaction was not associated."
          endif
          call messages_info(1)
        end if

      end if

    end do

    POP_SUB(create_interaction_with_partners)
  end subroutine create_interaction_with_partners

end module interactions_factory_abst_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
