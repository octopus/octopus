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
use units
use ps
use math

implicit none

type specie_type
  integer :: index                ! just a counter

  character(len=10) :: label      ! Identifier for the species:
                                  ! "jelli" : jellium sphere.
                                  ! "point" : jellium sphere of radius 0.5 a.u.
                                  ! "usdef" : user defined function
                                  ! other   : a pseudopotential.

  FLOAT          :: z          ! charge of the species. 

  FLOAT          :: z_val      ! valence charge of the species -- the total charge
                                  ! minus the core charge in the case of the pseudopotentials.

  FLOAT          :: weight     ! mass, in atomic mass units (=! atomic units of mass)

  logical           :: local      ! true if the potential is local, which in this case means
                                  ! it is *not* a pseudopotential.

  ! The rest of the stuff is only used in 3D calculations

  ! for the user defined potential
  character(len=1024) :: user_def

  ! jellium stuff
  FLOAT :: jradius

  ! For the pseudopotential
  character(len=3) :: ps_flavour
  type(ps_type), pointer :: ps
  logical           :: nlcc       ! true if we have non-local core corrections
  
end type specie_type

contains

function specie_init(s)
  integer :: specie_init
  type(specie_type), pointer :: s(:)

  integer :: nspecies, i, j, lmax, lloc
  character(len=80) :: str

  integer :: ispin

  call push_sub('specie_init')

  ! how many do we have?
  str = "Species"
  nspecies = loct_parse_block_n(str)
  if (nspecies < 1) then
    message(1) = "Input: Species block not specified"
    message(2) = '% Species'
    message(3) = '   specie <params>'
    message(4) = '%'
    call write_fatal(4)    
  end if
  allocate(s(nspecies))

  ! Reads the spin components. This is read here, as well as in states_init,
  ! to be able to pass it to the pseudopotential initializations subroutine.
  call loct_parse_int('SpinComponents', 1, ispin)
  if (ispin < 1 .or. ispin > 3) then
    write(message(1),'(a,i4,a)') "Input: '", ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 3)'
    call write_fatal(2)
  end if
  ispin = min(2, ispin)

  ! call dimension specific routines
  if(conf%dim==1.or.conf%dim==2) then
    call specie1D_init(nspecies, str, s)
  else
    call specie3D_init(nspecies, str, s)
  end if

  specie_init = nspecies

  call pop_sub()
end function specie_init

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
  FLOAT :: xx(3), r, vl, dvl, l1, l2
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

#include "specie1D.F90"
#include "specie3D.F90"

end module specie
