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

subroutine specie3D_init(nspecies, str, s)
  integer, intent(in) :: nspecies
  character(len=*), intent(in) :: str
  type(specie_type), pointer :: s(:)

  integer :: i, j, lmax, lloc, ispin

  ! Reads the spin components. This is read here, as well as in states_init,
  ! to be able to pass it to the pseudopotential initializations subroutine.
  call loct_parse_int('SpinComponents', 1, ispin)
  if (ispin < 1 .or. ispin > 3) then
    write(message(1),'(a,i4,a)') "Input: '", ispin,"' is not a valid SpinComponents"
    message(2) = '(SpinComponents = 1 | 2 | 3)'
    call write_fatal(2)
  end if
  ispin = min(2, ispin)
  
  do i = 1, nspecies
    s(i)%index = i
    call loct_parse_block_string(str, i-1, 0, s(i)%label)
    call loct_parse_block_float (str, i-1, 1, s(i)%weight)
    s(i)%weight =  units_inp%mass%factor * s(i)%weight ! units conversion
    
    select case(s(i)%label(1:5))
    case('jelli')
      s(i)%local = .true.  ! we only have a local part
      s(i)%nlcc  = .false. ! no non-local core corrections

      call loct_parse_block_float(str, i-1, 2, s(i)%Z)      ! charge of the jellium sphere
      call loct_parse_block_float(str, i-1, 3, s(i)%jradius)! radius of the jellium sphere
      s(i)%jradius = units_inp%length%factor * s(i)%jradius ! units conversion
      s(i)%Z_val = s(i)%Z
      
    case('point') ! this is treated as a jellium with radius 0.5
      s(i)%local = .true.
      s(i)%nlcc  = .false.

      call loct_parse_block_float(str, i-1, 2, s(i)%Z)
      s(i)%jradius = M_HALF
      s(i)%Z_val = 0 
      
    case('usdef') ! user defined
      s(i)%local = .true.
      s(i)%nlcc  = .false.

      call loct_parse_block_float (str, i-1, 2, s(i)%Z_val)
      call loct_parse_block_string(str, i-1, 3, s(i)%user_def)
      ! convert to C string
      j = len(trim(s(i)%user_def))
      s(i)%user_def(j+1:j+1) = achar(0) 

    case default ! a pseudopotential file
      s(i)%local = .false.

      allocate(s(i)%ps) ! allocate structure
      call loct_parse_block_float (str, i-1, 2, s(i)%Z)
      call loct_parse_block_string(str, i-1, 3, s(i)%ps_flavour)
      call loct_parse_block_int(str, i-1, 4, lmax)
      call loct_parse_block_int(str, i-1, 5, lloc)
      call ps_init(s(i)%ps, s(i)%label, s(i)%ps_flavour, s(i)%Z, lmax, lloc, ispin)
      if(conf%verbose>999) call ps_debug(s(i)%ps)

      s(i)%z_val = s(i)%ps%z_val
      s(i)%nlcc = (s(i)%ps%icore /= 'nc  ' )

    end select
  end do
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
  type(specie_type), intent(IN) :: s
  FLOAT, intent(in) :: x(3)
  integer, intent(in) :: l, lm, i
  FLOAT, intent(out) :: uV, duV(3)
  FLOAT :: r, f, uVr0, duvr0, ylm, gylm(3)
  FLOAT, parameter :: ylmconst = CNST(0.488602511902920) !  = sqr(3/(4*pi))
  logical, optional, intent(in) :: so

  r = sqrt(sum(x**2))
  if(present(so)) then
    if(so) then
      uVr0  = loct_splint(s%ps%so_kb(l, i), r)
      duvr0 = loct_splint(s%ps%so_dkb(l, i), r)
    else
      message(1) = 'Internal.'
      call write_fatal(1)
    endif
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

