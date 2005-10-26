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
!!
!! $Id$

#include "global.h"

module specie
use global
use messages
use syslabels
use lib_oct
use lib_oct_parser
use lib_oct_gsl_spline
use io
use units
use ps
use math

implicit none

private
public :: specie_init, &
          specie_end, &
          specie_type, &
          specie_debug, &
          specie_filter, &
          specie_read, &
          specie_get_local, &
          specie_get_glocal, &
          specie_get_local_fourier, &
          specie_get_nl_part, &
          specie_get_nlcc, &
          specie_get_iwf


integer, public, parameter :: &
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

  ! The number of initial wave functions
  integer :: niwfs
  integer, pointer :: iwf_l(:, :), iwf_m(:, :), iwf_i(:, :)

  ! For the TM pseudos, the lmax and lloc.
  integer :: lmax, lloc

  ! The filtering parameters.
  FLOAT :: alpha, beta, rcut, beta2
end type specie_type

contains

  ! ---------------------------------------------------------
  subroutine specie_debug(dir, s)
    character(len=*), intent(in)  :: dir
    type(specie_type), intent(in) :: s

    character(len=256) :: dirname
    integer :: iunit

    if(mpiv%node .ne. 0) then
       call write_debug_newlines(2)
       return
    endif

    call push_sub('specie.specie_debug')

    dirname = trim(dir)//'/'//trim(s%label)

    call io_mkdir(dirname)

    iunit = io_open(trim(dirname)//'/info', action='write')

    write(iunit, '(a,i3)')    'Index  = ', s%index
    write(iunit, '(2a)')      'Label  = ', trim(s%label)
    write(iunit, '(a,i3)')    'Type   = ', s%type
    if (s%type /= SPEC_USDEF ) write(iunit, '(a,f15.2)') 'z      = ', s%z
    write(iunit, '(a,f15.2)') 'z_val  = ', s%z_val
    write(iunit, '(a,f15.2)') 'weight = ', s%weight
    write(iunit, *)           'local  = ', s%local
    write(iunit, '(2a)')      'usdef  = ', trim(s%user_def)
    if (s%type == SPEC_JELLI .or. s%type == SPEC_POINT) then
      write(iunit, '(a,f15.2)') 'jradius= ', s%jradius
    end if
    write(iunit, *)           'nlcc   = ', s%nlcc
    write(iunit, '(a,f15.2)') 'def_rsize = ', s%def_rsize
    write(iunit, '(a,f15.2)') 'def_h = ', s%def_h
    if (s%type /= SPEC_USDEF ) write(iunit, '(a,i3)')    'lmax  = ', s%lmax
    if (s%type /= SPEC_USDEF ) write(iunit, '(a,i3)')    'lloc  = ', s%lloc

    if(.not.s%local) then
      if(in_debug_mode) call ps_debug(s%ps, trim(dirname))
    endif

    call io_close(iunit)
    call pop_sub()
  end subroutine specie_debug


  ! ---------------------------------------------------------
  subroutine specie_filter(s, gmax)
    type(specie_type),     intent(inout) :: s
    FLOAT, intent(in) :: gmax

    call push_sub('specie.specie_filter')

    call ps_filter(s%ps, gmax, s%alpha, s%beta, s%rcut, s%beta2)
    call ps_getradius(s%ps)
    call ps_derivatives(s%ps)

    call pop_sub()

  end subroutine specie_filter


  ! ---------------------------------------------------------
  subroutine specie_read(s, label)
    type(specie_type),     intent(inout) :: s
    character(len=*), intent(in) :: label

    character(len=256) :: fname
    character(len=10) :: lab
    integer :: i, row, n_spec_block, n_spec_def, iunit, read_data
    integer(POINTER_SIZE) :: blk

    call push_sub('specie.specie_read')

    s%local     = .true.  ! a local potential
    s%nlcc      = .false. ! without non-local core corrections
    s%def_h     = -M_ONE  ! not defined
    s%def_rsize = -M_ONE  ! not defined
    s%user_def  = ""
    read_data   = 0

    ! First, find out if there is a Species block.
    n_spec_block = 0
    blk = int(0, POINTER_SIZE)
    if(loct_parse_block(check_inp('Species'), blk) == 0) then
      n_spec_block = loct_parse_block_n(blk)
    end if

    ! Find out if the seeked species is in the block
    row = -1
    block: do i = 1, n_spec_block
       call loct_parse_block_string(blk, i-1, 0, lab)
       if(trim(lab)==trim(label)) then
         row = i - 1
         exit block
       endif
    enddo block

    ! Read whatever may be read from the block
    if(row>=0) then
       call read_from_block(blk, row, s, read_data)
       call loct_parse_block_end(blk)
    endif

    ! Find out if the species is in the defaults file.
    write(fname, '(2a)') trim(conf%share), "/PP/defaults"
    n_spec_def = max(0, loct_number_of_lines(fname))
    if(n_spec_def > 0) n_spec_def = n_spec_def - 1 ! First line is a comment

    iunit = io_open(fname, action='read', status='old', die=.false.)
    if(iunit > 0) then
      read(iunit,*)

      default_file: do i = 1, n_spec_def
        read(iunit,*) lab
        if(trim(lab) == trim(label)) then
          call read_from_default_file(iunit, read_data, s)
          exit default_file
        end if
      end do default_file

      call io_close(iunit)
    end if

    if(read_data == 0) then
      message(1) = 'Species '//trim(label)//' not found.'
      call write_fatal(1)
    endif

    call pop_sub()
  end subroutine specie_read


  ! ---------------------------------------------------------
  subroutine read_from_default_file(iunit, read_data, s)
    integer,           intent(in)    :: iunit
    integer,           intent(inout) :: read_data
    type(specie_type), intent(inout) :: s

    character(len=10) :: label
    FLOAT :: weight, z, def_h, def_rsize, alpha, beta, rcut, beta2
    integer :: lloc, lmax
    integer :: type

    call push_sub('specie.read_from_default_file')

    backspace(iunit)

    ! If it is in the default file, it *has* to be a pseudopotential
    s%local = .false.

    read(iunit,*) label, weight, type, z, lmax, lloc, def_h, def_rsize, alpha, beta, rcut, beta2
    def_h     = def_h     * P_ANG  ! These units are always in Angstrom
    def_rsize = def_rsize * P_ANG
    ASSERT(trim(label) == trim(s%label))

    if(read_data == 0) then ! The Species was not supplied in the block.
      s%weight    = weight
      s%type      = type
    end if
    if(read_data < 4) s%z         = z
    if(read_data < 5) s%lmax      = lmax
    if(read_data < 6) s%lloc      = lloc
    if(read_data < 7) s%def_h     = def_h
    if(read_data < 8) s%def_rsize = def_rsize
    if(read_data < 9) then
      s%alpha     = alpha
      s%beta      = beta
      s%rcut      = rcut
      s%beta2     = beta2
    end if
    read_data = 8

    call pop_sub()
  end subroutine read_from_default_file


  ! ---------------------------------------------------------
  subroutine read_from_block(blk, row, s, read_data)
    integer(POINTER_SIZE), intent(in) :: blk
    integer, intent(in) :: row
    type(specie_type),     intent(inout) :: s
    integer, intent(out) :: read_data
    integer :: j, n, lmax, lloc

    call push_sub('specie.read_from_block')

    read_data = 0

    call loct_parse_block_float (blk, row, 1, s%weight)
    call loct_parse_block_int   (blk, row, 2, s%type)

    select case(s%type)
    case(SPEC_USDEF) ! user defined
     call loct_parse_block_float (blk, row, 3, s%Z_val)
     call loct_parse_block_string(blk, row, 4, s%user_def)
     ! convert to C string
     j = len(trim(s%user_def))
     s%user_def(j+1:j+1) = achar(0)
     read_data = 5

    case(SPEC_POINT) ! this is treated as a jellium with radius 0.5
     call loct_parse_block_float(blk, row, 3, s%Z)
     s%jradius = M_HALF
     s%Z_val = 0
     read_data = 4

    case(SPEC_JELLI)
     call loct_parse_block_float(blk, row, 3, s%Z)      ! charge of the jellium sphere
     call loct_parse_block_float(blk, row, 4, s%jradius)! radius of the jellium sphere
     s%jradius = units_inp%length%factor * s%jradius ! units conversion
     s%Z_val = s%Z
     read_data = 5

    case(SPEC_PS_TM2,SPEC_PS_HGH) ! a pseudopotential file
     s%local = .false.
     n = loct_parse_block_cols(blk, row)

     call loct_parse_block_float (blk, row, 3, s%Z)
     lmax = 2 ! default
     if(n>4) then
        call loct_parse_block_int (blk, row, 4, s%lmax)
        read_data = 5
     endif
     lloc = 0 ! default
     if(n>5) then
        call loct_parse_block_int (blk, row, 5, s%lloc)
        read_data = 6
     endif

     if(n>6) then
       call loct_parse_block_float (blk, row, 6, s%def_h)
       s%def_h = s%def_h * units_inp%length%factor
       read_data = 7
     end if

     if(n>7) then
       call loct_parse_block_float (blk, row, 7, s%def_rsize)
       s%def_rsize = s%def_rsize * units_inp%length%factor
       read_data = 8
     end if

     if(n>8) then
       ! For the moment being, let us assume that this is always in atomic units.
       call loct_parse_block_float (blk, row, 8, s%alpha)
       call loct_parse_block_float (blk, row, 9, s%beta)
       call loct_parse_block_float (blk, row, 10, s%rcut)
       call loct_parse_block_float (blk, row, 11, s%beta2)
       read_data = 10
     else
       ! Reasonable defaults
       s%alpha = CNST(0.7); s%beta = CNST(18.0); s%rcut = CNST(2.5); s%beta2 = CNST(0.4)
     endif

    case default
      write(message(1), '(a,i2,a)') "Unknown pseudopotential type: '", s%type, "'"
      call write_fatal(1)
    end select

    call pop_sub()
  end subroutine read_from_block


  ! ---------------------------------------------------------
  subroutine specie_init(s, ispin)
    type(specie_type), intent(inout) :: s
    integer,           intent(in)    :: ispin

    integer :: i, l

    call push_sub('specie.specie_init')

    ! masses are always in a.u.m, so convert them to a.u.
    s%weight =  units_inp%mass%factor * s%weight

    ! if we are using non-local pseudopotentials, allocate structure
    if(.not.s%local) then
      allocate(s%ps) ! allocate structure
      call ps_init(s%ps, s%label, s%type, s%Z, s%lmax, s%lloc, ispin)
      call ps_getradius(s%ps)
      call ps_derivatives(s%ps)
      s%z_val = s%ps%z_val
      s%nlcc = (s%ps%icore /= 'nc  ' )
    endif
  
    s%niwfs = specie_niwfs(s)
    allocate(s%iwf_l(s%niwfs, ispin), s%iwf_m(s%niwfs, ispin), s%iwf_i(s%niwfs, ispin))
    call specie_iwf_fix_qn(s, ispin)

    call pop_sub()
  end subroutine specie_init


  ! ---------------------------------------------------------
  subroutine specie_end(ns, s)
    integer,           intent(in) :: ns
    type(specie_type), pointer    :: s(:)

    integer :: i

    call push_sub('specie.specie_end')

    do i = 1, ns
      if(s(i)%local) cycle

      if(associated(s(i)%ps)) call ps_end(s(i)%ps)
    end do

    if(associated(s)) then ! sanity check
      deallocate(s); nullify(s)
    end if

    call pop_sub()
  end subroutine specie_end


  ! ---------------------------------------------------------
  FLOAT function specie_get_local(s, x) result(l)
    type(specie_type), intent(in) :: s
    FLOAT,             intent(in) :: x(3)

    FLOAT :: a1, a2, Rb2 ! for jellium
    FLOAT :: xx(3), r

    xx(:) = x(:)
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
        l = loct_splint(s%ps%vl, r)
    end select

  end function specie_get_local


  ! ---------------------------------------------------------
  ! returns the gradient of the external potential
  ! ---------------------------------------------------------
  subroutine specie_get_glocal(s, x, gv)
    type(specie_type), intent(in) :: s
    FLOAT, intent(in) :: x(:)
    FLOAT, intent(out) :: gv(:)

    FLOAT, parameter :: Delta = CNST(1e-4)
    FLOAT :: xx(3), r, l1, l2
    integer :: i

    gv = M_ZERO

    xx(:) = x(:)
    r = sqrt(sum(xx(:)**2))

    select case(s%type)
    case(SPEC_USDEF)
      do i = 1, 3
        xx(:) = x(:)
        xx(i) = xx(i) - Delta
        l1 = loct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(:)**2)), s%user_def)
        xx(i) = xx(i) + M_TWO*Delta
        l2 = loct_parse_potential(xx(1), xx(2), xx(3), sqrt(sum(xx(:)**2)), s%user_def)
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
        gv(:) = M_ZERO
        if(r>CNST(0.00001)) gv(:) = -loct_splint(s%ps%dvl, r)*x(:)/r

    end select

  end subroutine specie_get_glocal


  ! ---------------------------------------------------------
  ! returns the localized part of the potential in Fourier space
  ! ---------------------------------------------------------
  FLOAT function specie_get_local_fourier(dim, s, x) result(l)
    integer,           intent(in) :: dim
    type(specie_type), intent(in) :: s
    FLOAT,             intent(in) :: x(3)

    FLOAT :: gmod

    if(dim /= 3 .or. s%type == SPEC_USDEF &
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
    type(specie_type), intent(in) :: s
    FLOAT, intent(in) :: x(3)

    ! only for 3D pseudopotentials, please
    if(s%type==SPEC_PS_TM2.or.s%type==SPEC_PS_HGH) then
      l = loct_splint(s%ps%core, sqrt(sum(x**2)))
    end if

  end function specie_get_nlcc


  ! ---------------------------------------------------------
  subroutine specie_get_nl_part(s, x, l, lm, i, uV, duV, so)
    type(specie_type), intent(in)  :: s
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

  FLOAT function specie_get_iwf(s, j, dim, is, x) result(phi)
    type(specie_type), intent(in) :: s
    integer,           intent(in) :: j
    integer,           intent(in) :: dim
    integer,           intent(in) :: is
    FLOAT,             intent(in) :: x(3)

    integer :: i, l, m
    FLOAT :: r, b, leng_scale, rprime

    i = s%iwf_i(j, is)
    l = s%iwf_l(j, is)
    m = s%iwf_m(j, is)
    r = sqrt(sum(x*x))

    if(.not.s%local) then
      phi =  loct_splint(s%ps%ur(i, is), r) * loct_ylm(x(1), x(2), x(3), l, m)
    else
      select case(dim)
      case(1)
        b = M_HALF*CNST(0.01) ! b = (1/2)*m*w^2
        leng_scale = (M_HALF/b)**(M_FOURTH)
        rprime = r/leng_scale
        if(rprime > M_ZERO) then
          phi = exp(-rprime**2/M_TWO) * (2**l) * loct_hypergeometric(-M_HALF*l, M_HALF, rprime**2)
        else
          select case(mod(l, 2))
          case(0)
            phi = factorial(l) * (-1)**(l/2) / factorial(l/2)
          case(1)
            phi = M_ZERO
          end select
        endif
      case(2) ! TODO
      case(3) ! TODO
      end select
    endif

  end function specie_get_iwf

  integer function specie_niwfs(s) result(n)
    type(specie_type), intent(in) :: s
    integer :: i, l
    call push_sub('specie.specie_niwfs')
    n = 0
    if(.not.s%local) then
      do i = 1, s%ps%conf%p
         l = s%ps%conf%l(i)
         if(sum(s%ps%conf%occ(i, :)).ne.M_ZERO) n = n + (2*l+1)
      enddo
    else
      n = 4 * s%z_val ! We hard-code the value of 4, rather arbitrarily.
    endif
    call pop_sub()
  end function specie_niwfs

  subroutine specie_iwf_fix_qn(s, ispin)
    type(specie_type), intent(inout) :: s
    integer, intent(in) :: ispin

    integer :: is, n, i, l, m
    call push_sub('specie.specie_iwf_fix_qn')

    if(.not.s%local) then
       do is = 1, ispin
          n = 1
          do i = 1, s%ps%conf%p
             l = s%ps%conf%l(i)
             do m = -l, l
                s%iwf_i(n, is) = i
                s%iwf_l(n, is) = l
                s%iwf_m(n, is) = m
                n = n + 1
             enddo
          enddo
       enddo
    else
     ! TODO
    endif
    call pop_sub()
  end subroutine specie_iwf_fix_qn

end module specie
