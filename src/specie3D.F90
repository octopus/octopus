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

subroutine specie3D_init(s, location, line, ispin)
  type(specie_type), intent(inout) :: s
  integer,           intent(in)    :: location, line, ispin

  integer :: i, j, lmax, lloc

  ASSERT(location==1.or.location==2)

  if(location == 1) then
    call from_block()
  else
    call from_default()
  end if

  ! masses are always in a.u.m, so convert them to a.u.
  s%weight =  units_inp%mass%factor * s%weight

  if(.not.s%local) then
    allocate(s%ps) ! allocate structure
    call ps_init(s%ps, s%label, s%ps_flavour, s%Z, lmax, lloc, ispin)
    if(conf%verbose>999) call ps_debug(s%ps)
    
    s%z_val = s%ps%z_val
    s%nlcc = (s%ps%icore /= 'nc  ' )
  end if

contains
  subroutine from_block()
    integer :: n
    
    call loct_parse_block_float ("Species", line-1, 1, s%weight)
    
    select case(s%label(1:5))
    case('jelli')
      call loct_parse_block_float("Species", line-1, 2, s%Z)      ! charge of the jellium sphere
      call loct_parse_block_float("Species", line-1, 3, s%jradius)! radius of the jellium sphere
      s%jradius = units_inp%length%factor * s%jradius ! units conversion
      s%Z_val = s%Z
      
    case('point') ! this is treated as a jellium with radius 0.5
      call loct_parse_block_float("Species", line-1, 2, s%Z)
      s%jradius = M_HALF
      s%Z_val = 0 
      
    case('usdef') ! user defined
      call loct_parse_block_float ("Species", line-1, 2, s%Z_val)
      call loct_parse_block_string("Species", line-1, 3, s%user_def)

      ! convert to C string
      j = len(trim(s%user_def))
      s%user_def(j+1:j+1) = achar(0) 

    case default ! a pseudopotential file
      s%local = .false.
      n = loct_parse_block_cols("Species", line-1)

      call loct_parse_block_float ("Species", line-1, 2, s%Z)

      s%ps_flavour = "tm2"
      if(n>3) call loct_parse_block_string("Species", line-1, 3, s%ps_flavour)

      lmax = 2 ! default
      if(n>4) call loct_parse_block_int   ("Species", line-1, 4, lmax)

      lloc = 0 ! default
      if(n>5) call loct_parse_block_int   ("Species", line-1, 5, lloc)

      if(n>6) then
        call loct_parse_block_float ("Species", line-1, 6, s%def_h)
        s%def_h = s%def_h * units_inp%length%factor
      end if

      if(n>7) then
        call loct_parse_block_float ("Species", line-1, 7, s%def_rsize)
        s%def_rsize = s%def_rsize * units_inp%length%factor
      end if
    end select
  end subroutine from_block

  subroutine from_default()
    integer :: iunit
    character(len=256) :: fname
    character(len=10)  :: label

    s%local = .false. ! we have a pseudopential

    write(fname, '(2a)') trim(conf%share), "/PP/defaults"
    call io_assign(iunit)
    open(unit=iunit, file=fname, action='read')

    ! go to the right line of the file
    do i = 1, line
      read(iunit,*)
    end do
    
    ! read information
    read(iunit,*) label, s%weight, s%Z, s%ps_flavour, lmax, lloc, s%def_h, s%def_rsize
    s%def_h     = s%def_h     * P_ANG  ! These units are always in Angstrom
    s%def_rsize = s%def_rsize * P_ANG
    call io_close(iunit)

    ! sanity check
    ASSERT(trim(label) == trim(s%label))

  end subroutine from_default

end subroutine specie3D_init

FLOAT function specie_get_nlcc(s, x) result(l)
  type(specie_type), intent(IN) :: s
  FLOAT, intent(in) :: x(3)

  ! only for 3D pseudopotentials, please
  if(conf%dim==3.and.(s%label(1:5).ne.'jelli'.and.s%label(1:5).ne.'point'.and.s%label(1:5).ne.'usdef')) then
    l = loct_splint(s%ps%core, sqrt(sum(x**2)))
  end if

end function specie_get_nlcc

subroutine specie_get_nl_part(s, x, l, lm, i, uV, duV, so)
  type(specie_type), intent(IN)  :: s
  FLOAT,             intent(in)  :: x(3)
  integer,           intent(in)  :: l, lm, i
  FLOAT,             intent(out) :: uV, duV(3)
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

