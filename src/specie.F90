module specie
use global
use units
use fdf
use ps

implicit none


type specie_type
  character(len=10) :: label
  real(r8) :: Z, Z_val
  real(r8) :: weight
  logical :: local   ! true if the potential is local

  ! jellium stuff
  real(r8) :: jradius

  ! For the pseudopotential
  type(ps_type), pointer :: ps

  ! For the local potential in Fourier space...
  complex(r8), pointer :: local_fw(:,:,:)
  complex(r8), pointer :: rhocore_fw(:,:,:)
end type specie_type

contains

function specie_init(s)
  integer :: specie_init
  type(specie_type), pointer :: s(:)

  integer :: nspecies, iunit, i, lmax, lloc
  character(len=100) :: str

  sub_name = 'specie_init'; call push_sub()

  ! how many do we have?
  nspecies = fdf_integer("NumberSpecies", 0)
  if (nspecies < 1) then
    write(message(1), '(a,i4,a)') "Input: '", nspecies, "' is not a valid NumberSpecies"
    message(2) = '(1 <= NumberSpecies)'
    call write_fatal(2)
  end if
  allocate(s(nspecies))

  if(fdf_block('Species', iunit) ) then
    do i = 1, nspecies
      read(iunit, '(a)') str ! read the complete line
      read(str, *) s(i)%label

      select case(s(i)%label(1:5))
      case('jelli')
        s(i)%local = .true.  ! we only have a local part
        ! s(i)%Z       = the charge of the jellium sphere
        ! s(i)%jradius = the radius of the jellium sphere
        read(str, *) s(i)%label, s(i)%weight, s(i)%Z, s(i)%jradius
        s(i)%jradius = units_inp%length%factor * s(i)%jradius ! units conversion
        s(i)%Z_val = s(i)%Z

      case('point') ! this is treated as a jellium with radius 0.5
        s(i)%local = .true.  ! we only have a local part
        read(str, *) s(i)%label, s(i)%weight, s(i)%Z
        s(i)%jradius = 0.5_r8
        s(i)%Z_val = 0 
        
      case default ! a pseudopotential file
        s(i)%local = .false.
        allocate(s(i)%ps) ! allocate structure
        read(str, *) s(i)%label, s(i)%weight, s(i)%Z, lmax, lloc
        
        call ps_init(s(i)%ps, s(i)%label, s(i)%Z, lmax, lloc, s(i)%Z_val)
      end select

      s(i)%weight =  units_inp%mass%factor * s(i)%weight ! units conversion

    end do
  else
    message(1) = "Input: Species block not specified"
    message(2) = '%block Species'
    message(3) = '   specie <params>'
    message(4) = '%endblock Species'
    call write_fatal(4)    
  end if

  specie_init = nspecies

  call pop_sub()
  return
end function specie_init

subroutine specie_end(ns, s)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)

  integer :: i

  sub_name = 'specie_end'; call push_sub()

  do i = 1, ns
    if(.not. s(i)%local .and. associated(s(i)%ps)) then
      call ps_end(s(i)%ps)
    end if
  end do

  if(associated(s)) then ! sanity check
    deallocate(s); nullify(s)
  end if

  call pop_sub()
  return
end subroutine specie_end

end module specie
