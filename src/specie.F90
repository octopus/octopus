!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.

#include "global.h"

module specie
use global
use lib_oct_parser
use lib_oct_gsl_spline
use io
use units
use ps
use math

implicit none

type specie_type
  integer :: index                ! just a counter

  character(len=15) :: label      ! Identifier for the species:
                                  ! "jelli" : jellium sphere.
                                  ! "point" : jellium sphere of radius 0.5 a.u.
                                  ! "usdef" : user defined function
                                  !  other  : a pseudopotential.

  FLOAT   :: z          ! charge of the species. 

  FLOAT   :: z_val      ! valence charge of the species -- the total charge
                        ! minus the core charge in the case of the pseudopotentials.

  FLOAT   :: weight     ! mass, in atomic mass units (!= atomic units of mass)

  logical :: local      ! true if the potential is local, which in this case means
                        ! it is *not* a pseudopotential.

  ! for the user defined potential
  character(len=1024) :: user_def

  ! jellium stuff
  FLOAT :: jradius

  ! For the pseudopotential
  character(len=3)       :: ps_flavour
  type(ps_type), pointer :: ps
  logical                :: nlcc       ! true if we have non-local core corrections

  ! the default values for the spacing and atomic radius
  FLOAT :: def_rsize, def_h
end type specie_type

contains

  subroutine specie_init(s, location, line, ispin)
    type(specie_type), intent(inout) :: s
    integer,           intent(in)    :: location, line, ispin

    call push_sub('specie_init')

    ! some defaults
    s%local     = .true.  ! a local potential
    s%nlcc      = .false. ! without non-local core corrections
    s%def_h     = -M_ONE  ! not defined
    s%def_rsize = -M_ONE  ! not defined
  
    ! call dimension specific routines
    if(conf%dim==1.or.conf%dim==2) then
      call specie1D_init(s, location, line)
    else
      call specie3D_init(s, location, line, ispin)
    end if
    
    call pop_sub()
  end subroutine specie_init

subroutine specie_end(ns, s)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)

  integer :: i

  call push_sub('specie_end')

  do i = 1, ns
    if(s(i)%local) cycle

    if(associated(s(i)%ps)) call ps_end(s(i)%ps)
  end do
  
  if(associated(s)) then ! sanity check
    deallocate(s); nullify(s)
  end if

  call pop_sub()
  return
end subroutine specie_end

FLOAT function specie_get_local(s, x) result(l)
  type(specie_type), intent(IN) :: s
  FLOAT, intent(in) :: x(3)

  FLOAT :: a1, a2, Rb2 ! for jellium
  FLOAT :: xx(3), r

  xx = M_ZERO
  xx(1:conf%dim) = x(:)
  r = sqrt(sum(xx(:)**2))

  if(conf%dim.ne.3 .or. s%label(1:5)=='usdef') then
    l = loct_parse_potential(xx(1), xx(2), xx(3), r, s%user_def)

  else
    select case(s%label(1:5))
    case('jelli', 'point')
      a1 = s%Z/(M_TWO*s%jradius**3)
      a2 = s%Z/s%jradius
      Rb2= s%jradius**2
      
      if(r <= s%jradius) then
        l = (a1*(r*r - Rb2) - a2)
      else
        l = - s%Z/r
      end if

    case default
      if(r >= r_small) then
        l = (loct_splint(s%ps%vlocal,  r) - s%Z_val)/r
      else
        l = s%ps%Vlocal_origin
      end if
    end select
  end if

end function specie_get_local

! returns the gradient of the external potential
subroutine specie_get_glocal(s, x, gv)
  type(specie_type), intent(IN) :: s
  FLOAT, intent(in) :: x(conf%dim)
  FLOAT, intent(out) :: gv(conf%dim)

  FLOAT, parameter :: Delta = CNST(1e-4)
  FLOAT :: xx(3), r, l1, l2
  integer :: i

  gv = M_ZERO
  if(conf%dim.ne.3 .or. s%label(1:5)=='usdef') then
    xx = M_ZERO
    do i = 1, conf%dim
      xx(1:conf%dim) = x(1:conf%dim)
      xx(i) = xx(i) - Delta
      l1 = loct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(1:conf%dim)**2)), s%user_def)
      xx(i) = xx(i) + M_TWO*Delta
      l2 = loct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(1:conf%dim)**2)), s%user_def)
      gv(i) = (l2 - l1)/(M_TWO*Delta)
    end do
  else
    xx(:) = x(:)
    r = sqrt(sum(xx(:)**2))

    select case(s%label(1:5))
    case('jelli', 'point')
      l1 = s%Z/(M_TWO*s%jradius**3)
      
      if(r <= s%jradius) then
        gv(:) = l1*x(:)
      else
        gv(:) = s%Z*x(:)/r**3
      end if

    case default
      l1 = loct_splint(s%ps%vlocal, r)
      l2 = loct_splint(s%ps%dvlocal, r)
      
      gv(:) = -(l2 - (l1 - s%Z_val)/r)/r**2 * x(:)
      
    end select
  end if

end subroutine specie_get_glocal

! returns the localized part of the potential in Fourier space
FLOAT function specie_get_local_fourier(s, x) result(l)
  type(specie_type), intent(IN) :: s
  FLOAT, intent(in) :: x(3)

  FLOAT :: gmod

  if(conf%dim /= 3 .or. s%label(1:5) == 'usdef' &
    .or. s%label(1:5) == 'jelli' .or. s%label(1:5) == 'point') then
    message(1) = 'Periodic arrays of usedef, jelli, point,' 
    message(2) = '1D and 2D systems not implemented yet.'
    call write_fatal(2)

  else
    gmod = sqrt(sum(x(:)**2))
    l = loct_splint(s%ps%vlocal_f, gmod)
  end if

end function specie_get_local_fourier


#include "specie1D.F90"
#include "specie3D.F90"

end module specie
