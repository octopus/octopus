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

integer, parameter :: &
   SPEC_USDEF  = 1, &      ! user defined function
   SPEC_POINT  = 2, &      ! point charge: jellium sphere of radius 0.5 a.u.
   SPEC_JELLI  = 3, &      ! jellium sphere.
   SPEC_PS_TM2 = PS_TM2, & ! Troullier-Martins pseudopotential
   SPEC_PS_HGH = PS_HGH    ! HGH pseudopotential

type specie_type
  integer :: index                ! just a counter

  character(len=15) :: label      ! Identifier for the species
  integer :: type       ! what type of species
  FLOAT   :: z          ! charge of the species
  FLOAT   :: z_val      ! valence charge of the species -- the total charge
                        ! minus the core charge in the case of the pseudopotentials
  FLOAT   :: weight     ! mass, in atomic mass units (!= atomic units of mass)
  logical :: local      ! true if the potential is local, which in this case means
                        ! it is *not* a pseudopotential.

  ! for the user defined potential
  character(len=1024) :: user_def

  ! jellium stuff
  FLOAT :: jradius

  ! For the pseudopotential
  type(ps_type), pointer :: ps
  logical                :: nlcc       ! true if we have non-local core corrections

  ! the default values for the spacing and atomic radius
  FLOAT :: def_rsize, def_h
end type specie_type

contains

  subroutine specie_filter(s, gmax)
    type(specie_type),     intent(inout) :: s
    FLOAT, intent(in) :: gmax
    call push_sub('specie_cutoff')
    call ps_filter(s%ps, gmax)
    call ps_derivatives(s%ps)
    ! This is for debugging. It should not be here, though.
    !if(conf%verbose>=VERBOSE_DEBUG) call ps_debug(s%ps)
    call pop_sub(); return
  end subroutine specie_filter

  ! ---------------------------------------------------------
  subroutine specie_init(s, location, blk, line, ispin)
    type(specie_type),     intent(inout) :: s
    integer,               intent(in)    :: location, line, ispin
    integer(POINTER_SIZE), intent(in)    :: blk

    integer :: lmax, lloc

    call push_sub('specie_init')

    ! some defaults
    s%local     = .true.  ! a local potential
    s%nlcc      = .false. ! without non-local core corrections
    s%def_h     = -M_ONE  ! not defined
    s%def_rsize = -M_ONE  ! not defined
  
    ASSERT(location==1.or.location==2)

    if(location == 1) then
      call from_block()
    else ! from default file
      if(conf%dim==1.or.conf%dim==2) then
        message(1) = "In 1 or 2 dimensions all species must be defined in the %Species block"
        call write_fatal(1)
      end if
      
      call from_default_file()
    end if
      
    ! masses are always in a.u.m, so convert them to a.u.
    s%weight =  units_inp%mass%factor * s%weight

    ! if we are using non-local pseudopotentials, allocate structure
    if(.not.s%local) then
      allocate(s%ps) ! allocate structure
      call ps_init(s%ps, s%label, s%type, s%Z, lmax, lloc, ispin)
      call ps_derivatives(s%ps)
      if(conf%verbose>=VERBOSE_DEBUG) call ps_debug(s%ps)
      
      s%z_val = s%ps%z_val
      s%nlcc = (s%ps%icore /= 'nc  ' )
    end if

    call pop_sub()

  contains
    
    subroutine from_block()
      integer :: j, n

      call loct_parse_block_float (blk, line-1, 1, s%weight)
      call loct_parse_block_int   (blk, line-1, 2, s%type)
      
      select case(s%type)
      case(SPEC_USDEF) ! user defined
        call loct_parse_block_float (blk, line-1, 3, s%Z_val)
        call loct_parse_block_string(blk, line-1, 4, s%user_def)
        
        ! convert to C string
        j = len(trim(s%user_def))
        s%user_def(j+1:j+1) = achar(0) 

      case(SPEC_POINT) ! this is treated as a jellium with radius 0.5
        call loct_parse_block_float(blk, line-1, 3, s%Z)
        s%jradius = M_HALF
        s%Z_val = 0 
        
      case(SPEC_JELLI)
        call loct_parse_block_float(blk, line-1, 3, s%Z)      ! charge of the jellium sphere
        call loct_parse_block_float(blk, line-1, 4, s%jradius)! radius of the jellium sphere
        s%jradius = units_inp%length%factor * s%jradius ! units conversion
        s%Z_val = s%Z
      
      case(SPEC_PS_TM2,SPEC_PS_HGH) ! a pseudopotential file
        s%local = .false.
        n = loct_parse_block_cols(blk, line-1)
        
        call loct_parse_block_float (blk, line-1, 3, s%Z)
        
        lmax = 2 ! default
        if(n>4) call loct_parse_block_int (blk, line-1, 4, lmax)
        
        lloc = 0 ! default
        if(n>5) call loct_parse_block_int (blk, line-1, 5, lloc)
        
        if(n>6) then
          call loct_parse_block_float (blk, line-1, 6, s%def_h)
          s%def_h = s%def_h * units_inp%length%factor
        end if
        
        if(n>7) then
          call loct_parse_block_float (blk, line-1, 7, s%def_rsize)
          s%def_rsize = s%def_rsize * units_inp%length%factor
        end if

      case default
        write(message(1), '(a,i2,a)') "Unknown pseudopotential type: '", s%type, "'"
        call write_fatal(1)

      end select

    end subroutine from_block

    subroutine from_default_file()
      integer :: i, iunit
      character(len=256) :: fname
      character(len=10)  :: label
      
      s%local = .false. ! we have a pseudo-potencial
      
      write(fname, '(2a)') trim(conf%share), "/PP/defaults"
      call io_assign(iunit)
      open(unit=iunit, file=fname, action='read')
      
      ! go to the right line of the file
      do i = 1, line
        read(iunit,*)
      end do
      
      ! read information
      read(iunit,*) label, s%weight, s%type, s%Z, lmax, lloc, s%def_h, s%def_rsize
      s%def_h     = s%def_h     * P_ANG  ! These units are always in Angstrom
      s%def_rsize = s%def_rsize * P_ANG
      call io_close(iunit)

      ! sanity check
      ASSERT(trim(label) == trim(s%label))
      
    end subroutine from_default_file
 
  end subroutine specie_init


  ! ---------------------------------------------------------
  subroutine specie_end(ns, s)
    integer,           intent(in) :: ns
    type(specie_type), pointer    :: s(:)

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
  end subroutine specie_end


  FLOAT function specie_get_local(s, x) result(l)
    type(specie_type), intent(in) :: s
    FLOAT,             intent(in) :: x(3)

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(3), r

    xx = M_ZERO
    xx(1:conf%dim) = x(:)
    r = sqrt(sum(xx(:)**2))

    select case(s%type)
    case(SPEC_USDEF)
      l = loct_parse_potential(xx(1), xx(2), xx(3), r, s%user_def)
      
    case(SPEC_POINT, SPEC_JELLI)
      a1 = s%Z/(M_TWO*s%jradius**3)
      a2 = s%Z/s%jradius
      Rb2= s%jradius**2
      
      if(r <= s%jradius) then
        l = (a1*(r*r - Rb2) - a2)
      else
        l = - s%Z/r
      end if
      
    case(SPEC_PS_TM2, SPEC_PS_HGH)
!!$      if(r >= r_small) then
!!$        l = (loct_splint(s%ps%vlocal,  r) - s%Z_val)/r
!!$        l = loct_splint(s%ps%vl, r)
!!$      else
!!$        l = s%ps%Vlocal_origin
!!$      end if
        l = loct_splint(s%ps%vl, r)
    end select

  end function specie_get_local


  ! ---------------------------------------------------------
  ! returns the gradient of the external potential
  ! ---------------------------------------------------------
  subroutine specie_get_glocal(s, x, gv)
    type(specie_type), intent(IN) :: s
    FLOAT, intent(in) :: x(conf%dim)
    FLOAT, intent(out) :: gv(conf%dim)
    
    FLOAT, parameter :: Delta = CNST(1e-4)
    FLOAT :: xx(3), r, l1, l2
    integer :: i
    
    gv = M_ZERO
    xx = M_ZERO
    
    xx(1:conf%dim) = x(1:conf%dim)
    r = sqrt(sum(xx(:)**2))
    
    select case(s%type)
    case(SPEC_USDEF)
      do i = 1, conf%dim
        xx(1:conf%dim) = x(1:conf%dim)
        xx(i) = xx(i) - Delta
        l1 = loct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(1:conf%dim)**2)), s%user_def)
        xx(i) = xx(i) + M_TWO*Delta
        l2 = loct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(1:conf%dim)**2)), s%user_def)
        gv(i) = (l2 - l1)/(M_TWO*Delta)
      end do
      
    case(SPEC_POINT, SPEC_JELLI)
      l1 = s%Z/(M_TWO*s%jradius**3)
      
      if(r <= s%jradius) then
        gv(:) = l1*x(:)
      else
        gv(:) = s%Z*x(:)/r**3
      end if
      
    case(SPEC_PS_TM2, SPEC_PS_HGH)
      l1 = loct_splint(s%ps%vlocal, r)
      l2 = loct_splint(s%ps%dvlocal, r)
      
      gv(:) = -(l2 - (l1 - s%Z_val)/r)/r**2 * x(:)
      
    end select

  end subroutine specie_get_glocal


  ! ---------------------------------------------------------
  ! returns the localized part of the potential in Fourier space
  ! ---------------------------------------------------------
  FLOAT function specie_get_local_fourier(s, x) result(l)
    type(specie_type), intent(IN) :: s
    FLOAT, intent(in) :: x(3)
    
    FLOAT :: gmod
    
    if(conf%dim /= 3 .or. s%type == SPEC_USDEF &
       .or. s%type == SPEC_POINT .or. s%type == SPEC_JELLI) then
      message(1) = 'Periodic arrays of usedef, jelli, point,' 
      message(2) = '1D and 2D systems not implemented yet.'
      call write_fatal(2)
      
    else
      gmod = sqrt(sum(x(:)**2))
      l = loct_splint(s%ps%vlocal_f, gmod)
    end if
    
  end function specie_get_local_fourier
  

  ! ---------------------------------------------------------
  FLOAT function specie_get_nlcc(s, x) result(l)
    type(specie_type), intent(IN) :: s
    FLOAT, intent(in) :: x(3)
    
    ! only for 3D pseudopotentials, please
    if(s%type==SPEC_PS_TM2.or.s%type==SPEC_PS_HGH) then
      l = loct_splint(s%ps%core, sqrt(sum(x**2)))
    end if
    
  end function specie_get_nlcc


  ! ---------------------------------------------------------
  subroutine specie_get_nl_part(s, x, l, lm, i, uV, duV, so)
    type(specie_type), intent(IN)  :: s
    FLOAT,             intent(in)  :: x(:)        ! (3)
    integer,           intent(in)  :: l, lm, i
    FLOAT,             intent(out) :: uV, duV(:)  ! (3)
    logical, optional, intent(in)  :: so
    
    FLOAT :: r, uVr0, duvr0, ylm, gylm(3)
    FLOAT, parameter :: ylmconst = CNST(0.488602511902920) !  = sqr(3/(4*pi))
    
    r = sqrt(sum(x**2))
    if(present(so)) then
      ASSERT(so)
      
      uVr0  = loct_splint(s%ps%so_kb(l, i), r)
      duvr0 = loct_splint(s%ps%so_dkb(l, i), r)
    else
      uVr0  = loct_splint(s%ps%kb(l, i), r)
      duvr0 = loct_splint(s%ps%dkb(l, i), r)
    endif
    
    call grylmr(x(1), x(2), x(3), l, lm, ylm, gylm)
    uv = uvr0*ylm
    if(r >= r_small) then
      duv(:) = duvr0 * ylm * x(:)/r + uvr0 * gylm(:)
    else
      if(l == 1) then
        duv = M_ZERO
        if(lm == -1) then
          duv(2) = -ylmconst * duvr0
        else if(lm == 0) then
          duv(3) =  ylmconst * duvr0
        else if(lm == 1) then
          duv(1) = -ylmconst * duvr0
        end if
      else
        duv = M_ZERO
      end if
    end if
    
  end subroutine specie_get_nl_part


end module specie
