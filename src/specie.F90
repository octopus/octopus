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

#include "config_F90.h"

module specie
use global
use units
use ps
use spline
use math

implicit none

type specie_type
  character(len=10) :: label      ! Identifier for the species:
                                  ! "jelli" : jellium sphere.
                                  ! "point" : jellium sphere of radius 0.5 a.u.
                                  ! "usdef" : user defined function
                                  ! other   : a pseudopotential.

  real(r8)          :: z          ! charge of the species. 

  real(r8)          :: z_val      ! valence charge of the species -- the total charge
                                  ! minus the core charge in the case of the pseudopotentials.

  real(r8)          :: weight     ! mass, in atomic mass units (=! atomic units of mass)

  logical           :: local      ! true if the potential is local, which in this case means
                                  ! it is *not* a pseudopotential.

  ! for the user defined potential
  character(len=1024) :: user_def

  ! For the local pseudopotential in Fourier space...
  complex(r8), pointer :: local_fw(:,:,:)

#if defined(THREE_D)
  ! jellium stuff
  real(r8) :: jradius

  ! For the pseudopotential
  character(len=3) :: ps_flavour
  type(ps_type), pointer :: ps
  logical           :: nlcc       ! true if we have non-local core corrections
  
  ! for the core density in Fourier space
  complex(r8), pointer :: rhocore_fw(:,:,:)      
  
  ! For the non-local pp in fourier space
  integer(POINTER_SIZE) :: nl_planb
  integer :: nl_fft_n(3), nl_hfft_n
  complex(r8), pointer :: nl_fw(:,:,:,:,:), nl_dfw(:,:,:,:,:,:)
#endif

end type specie_type

contains

function specie_init(s)
  integer :: specie_init
  type(specie_type), pointer :: s(:)

  integer :: nspecies, i, j, lmax, lloc
  character(len=80) :: str

  integer :: ispin

  sub_name = 'specie_init'; call push_sub()

  ! how many do we have?
  str = C_string("Species")
  nspecies = oct_parse_block_n(str)
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
  call oct_parse_int(C_string('SpinComponents'), 1, ispin)
  if (ispin < 1 .or. ispin > 3) then
    write(message(1),'(a,i4,a)') "Input: '", ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 3)'
    call write_fatal(2)
  end if
  ispin = min(2, ispin)

  ! call dimension specific routines
  call specie_init_dim(nspecies, str, s)

  specie_init = nspecies

  call pop_sub(); return
end function specie_init

#if defined(ONE_D) || defined(TWO_D)
#  include "specie1D.F90"
#elif defined(THREE_D)
#  include "specie3D.F90"
#endif

end module specie
