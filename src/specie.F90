#include "config_F90.h"

module specie
use global
use units
use ps

implicit none

type specie_type
  character(len=10) :: label
  real(r8) :: Z, Z_val
  real(r8) :: weight
  logical :: local   ! true if the potential is local

  ! jellium stuff
  real(r8) :: jradius

  ! for the user defined potential
  character(len=1024) :: user_def

  ! For the pseudopotential
  type(ps_type), pointer :: ps

  ! For the local pseudopotential in Fourier space...
#if defined(THREE_D)
  complex(r8), pointer :: &
       local_fw(:,:,:),    &  ! for the potential
       rhocore_fw(:,:,:)      ! for the core density
#elif defined(ONE_D)
  complex(r8), pointer :: &
       local_fw(:),       &
       rhocore_fw(:, :, :)
#endif

  ! For the non-local pp in fourier space
  integer(POINTER_SIZE) :: nl_planb
  integer :: nl_fft_n(3), nl_hfft_n
  complex(r8), pointer :: nl_fw(:,:,:,:), nl_dfw(:,:,:,:,:)
end type specie_type

contains

function specie_init(s)
  integer :: specie_init
  type(specie_type), pointer :: s(:)

  integer :: nspecies, i, j, lmax, lloc
  character(len=80) :: str

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

  do i = 1, nspecies
    call oct_parse_block_str(str, i-1, 0, s(i)%label)
    call oct_parse_block_double(str, i-1, 1, s(i)%weight)

#if defined(THREE_D)
    select case(s(i)%label(1:5))
    case('jelli')
      s(i)%local = .true.  ! we only have a local part
      ! s(i)%Z       = the charge of the jellium sphere
      ! s(i)%jradius = the radius of the jellium sphere
      call oct_parse_block_double(str, i-1, 2, s(i)%Z)
      call oct_parse_block_double(str, i-1, 3, s(i)%jradius)
      s(i)%jradius = units_inp%length%factor * s(i)%jradius ! units conversion
      s(i)%Z_val = s(i)%Z
      
    case('point') ! this is treated as a jellium with radius 0.5
      s(i)%local = .true.  ! we only have a local part
      call oct_parse_block_double(str, i-1, 2, s(i)%Z)
      s(i)%jradius = 0.5_r8
      s(i)%Z_val = 0 
      
    case('usdef') ! user defined
      s(i)%local = .true.
      call oct_parse_block_double(str, i-1, 2, s(i)%Z_val)
      call oct_parse_block_str   (str, i-1, 3, s(i)%user_def)
      ! convert to C string
      j = len(trim(s(i)%user_def))
      s(i)%user_def(j+1:j+1) = achar(0) 

    case default ! a pseudopotential file
      s(i)%local = .false.
      allocate(s(i)%ps) ! allocate structure
      call oct_parse_block_double(str, i-1, 2, s(i)%Z)
      call oct_parse_block_int(str, i-1, 3, lmax)
      call oct_parse_block_int(str, i-1, 4, lloc)
      call ps_init(s(i)%ps, s(i)%label, s(i)%Z, lmax, lloc, s(i)%Z_val)
      s(i)%nl_planb= int(-1, POINTER_SIZE)
    end select

#elif defined(ONE_D)
    s(i)%local = .true. ! In 1D, potential has to be local.
    select case(s(i)%label(1:5))
    case('usdef') ! user defined
      call oct_parse_block_double(str, i-1, 2, s(i)%Z_val)
      call oct_parse_block_str   (str, i-1, 3, s(i)%user_def)
      ! convert to C string
      j = len(trim(s(i)%user_def))
      s(i)%user_def(j+1:j+1) = achar(0) 
    case default ! built by the program
      allocate(s(i)%ps)
      call oct_parse_block_double(str, i-1, 2, s(i)%z)
      call oct_parse_block_double(str, i-1, 3, s(i)%z_val)
      call oct_parse_block_str   (str, i-1, 4, s(i)%user_def)
      s(i)%user_def = C_string(s(i)%user_def) 
      call ps_init(s(i)%ps, s(i)%label, s(i)%Z, s(i)%Z_val, s(i)%user_def)
    end select
#endif

    s(i)%weight =  units_inp%mass%factor * s(i)%weight ! units conversion
    
  end do

  specie_init = nspecies

  call pop_sub()
  return
end function specie_init

subroutine specie_end(ns, s)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)

  integer :: i

  sub_name = 'specie_end'; call push_sub()

  species: do i = 1, ns

    nlocal: if(.not. s(i)%local) then
      if(associated(s(i)%ps)) then
        if(s(i)%ps%icore /= 'nc  ' .and. associated(s(i)%rhocore_fw)) then
          deallocate(s(i)%rhocore_fw); nullify(s(i)%rhocore_fw)
        end if
        call ps_end(s(i)%ps)
      end if

      if(s(i)%nl_planb.ne. int(-1, POINTER_SIZE)) then
        call fftw_f77_destroy_plan(s(i)%nl_planb)
        deallocate(s(i)%nl_fw, s(i)%nl_dfw)
        nullify(s(i)%nl_fw, s(i)%nl_dfw)
      end if
    end if nlocal

    if(associated(s(i)%local_fw)) then
      deallocate(s(i)%local_fw); nullify(s(i)%local_fw)
    end if

  end do species
  
  if(associated(s)) then ! sanity check
    deallocate(s); nullify(s)
  end if

  call pop_sub()
  return
end subroutine specie_end

end module specie
