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

module atom
use io
use specie

implicit none

type atom_type
  character(len=10) :: label
  type(specie_type), pointer :: spec ! pointer to specie

  real(r8) :: x(3), v(3), f(3) ! position/velocity/force of atom in real space

  logical :: move              ! should I move this atom in the optimization mode

  ! the mesh around a given atom...
  integer :: Mps
  integer, pointer :: Jxyz(:)
  real(r8), pointer ::    pnts_ps,            &  ! # points in ps sphere
                          duV(:,:,:),         &
                          duVu(:,:,:),        &  ! the Kleinman Bylander projectors
                          dduV(:,:,:,:)
  ! This is for performance reasons.
  complex(r8), pointer :: zpnts_ps,           & 
                          zuV(:,:,:),         &
                          zuVu(:,:,:),        &
                          zduV(:,:,:,:)
  complex(r8), pointer :: so_uv(:, :, :),     &
                          so_uvu(:, :, :),    &
                          so_duv(:, :, :, :), &
                          so_luv(:, :, :, :)
  complex(r8), pointer :: phases(:,:)    ! factors exp(ik*x)
end type atom_type

type atom_classical_type
  real(r8) :: x(3), v(3), f(3)
  real(r8) :: charge

  character(len=4) :: label
end type atom_classical_type

contains

subroutine atom_init(natoms, a, ncatoms, ca, ns, s)
  integer, intent(out) :: natoms, ncatoms
  type(atom_type), pointer :: a(:)
  type(atom_classical_type), pointer :: ca(:)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)

  integer :: iunit, i, j
  integer(POINTER_SIZE) :: random_gen_pointer
  character(len=80) :: str, label
  logical :: l
  real(r8) :: temperature, sigma, x(3), kin1, kin2

  call push_sub('atom_init')

  if(conf%dim == 3.and.oct_parse_isdef("PDBCoordinates").ne.0) then
    call oct_parse_string('PDBCoordinates', 'coords.pdb', label)

    call io_assign(iunit)
    open(iunit, status='unknown', file=trim(label))
    call loadPDB(iunit, natoms, a, ncatoms, ca)
    do i = 1, natoms
       a(i)%spec => s(get_specie(a(i)%label))
    enddo
    call io_close(iunit)

  else
    ! we now load the positions, either from the input, or from a file
    if(conf%dim == 3.and.oct_parse_isdef("XYZCoordinates").ne.0) then ! read a xyz file
      call oct_parse_string('XYZCoordinates', 'coords.xyz', label)

      call io_assign(iunit)
      open(iunit, status='unknown', file=trim(label))
      read(iunit, *) natoms
      read(iunit, *) ! skip comment line
      allocate(a(natoms))
      nullify(ca); ncatoms = 0;
      
      do i = 1, natoms
        read(iunit,*) a(i)%label, a(i)%x(:)
        a(i)%move = .true.
        a(i)%spec => s(get_specie(a(i)%label))
      end do

      call io_close(iunit)
    else
      str = "Coordinates"
      natoms = oct_parse_block_n(str)
      if(natoms <= 0) then
        message(1) = "Input: Coordinates block not specified"
        message(2) = '% Coordinates'
        message(3) = '  specie  x  y  z  move'
        message(4) = '%'
        call write_fatal(4)
      end if
      
      allocate(a(natoms))
      nullify(ca); ncatoms = 0
      do i = 1, natoms
        call oct_parse_block_string(str, i-1, 0, a(i)%label)
        a(i)%spec => s(get_specie(a(i)%label))
        call oct_parse_block_double (str, i-1, 1, a(i)%x(1))
        call oct_parse_block_double (str, i-1, 2, a(i)%x(2))
        call oct_parse_block_double (str, i-1, 3, a(i)%x(3))
        call oct_parse_block_logical(str, i-1, 4, a(i)%move)
      end do
    end if
  end if

  ! units conversion
  do i = 1, natoms
    a(i)%x = a(i)%x * units_inp%length%factor

    ! seems that some compilers do not initilize pointers
    ! so we do it explicitly
    nullify(a(i)%Jxyz, a(i)%phases)
  end do
  do i = 1, ncatoms
    ca(i)%x = ca(i)%x * units_inp%length%factor
  end do

  ! we now load the velocities, either from the temperature, from the input, or from a file
  if(oct_parse_isdef("RandomVelocityTemp").ne.0) then
    
    call oct_ran_init(random_gen_pointer)
    call oct_parse_double("RandomVelocityTemp", 0.0_r8, temperature)
    do i = 1, natoms
       sigma = sqrt( P_Kb*temperature / a(i)%spec%weight )
       do j = 1, 3
          a(i)%v(j) = oct_ran_gaussian(random_gen_pointer, sigma)
       enddo
    enddo
    call oct_ran_end(random_gen_pointer)
    kin1 = kinetic_energy(natoms, a)
    call cm_vel(natoms, a, x)
    do i = 1, natoms
       a(i)%v(:) = a(i)%v(:) - x
    enddo
    kin2 = kinetic_energy(natoms, a)
    do i = 1, natoms
       a(i)%v(:) =  sqrt(kin1/kin2)*a(i)%v(:)
    enddo

    write(message(1),'(a,f10.4,1x,a)') 'Info: Initial velocities ramdomly distributed with T =', & 
                                   temperature, 'K'
    write(message(2),'(a,f8.4,1x,a)') 'Info: <K>       =', &
                                   (kinetic_energy(natoms, a)/natoms)/units_out%energy%factor, &
                                   units_out%energy%abbrev
    write(message(3),'(a,f8.4,1x,a)') 'Info: 3/2 k_B T =', &
                                   (3.0_r8 / 2.0_r8)*P_Kb*temperature/units_out%energy%factor, &
                                   units_out%energy%abbrev
    write(message(4),'(a)')
    call write_info(4)

  elseif(oct_parse_isdef("XYZVelocities").ne.0 .and. conf%dim==3) then ! read a xyz file
    call io_assign(iunit)
    call oct_parse_string('XYZVelocities', 'velocities.xyz', label)
    open(iunit, status='unknown', file=trim(label))
      
    read(iunit, *)
    read(iunit, *) ! skip comment line
      
    do i = 1, natoms
      read(iunit,*) label, a(i)%v(:)
      a(i)%v(:) = units_inp%velocity%factor * a(i)%v(:) !units conversion
    end do

    call io_close(iunit)
  else 
    str = "Velocities"
    if(oct_parse_isdef(str).ne.0) then
      do i = 1, natoms
        call oct_parse_block_double (str, i-1, 1, a(i)%v(1))
        call oct_parse_block_double (str, i-1, 2, a(i)%v(2))
        call oct_parse_block_double (str, i-1, 3, a(i)%v(3))
        a(i)%v = a(i)%v * units_inp%velocity%factor
      end do
    else
      a(:)%v(1) = 0._r8
      a(:)%v(2) = 0._r8
      a(:)%v(3) = 0._r8
    end if
  end if
  
  call pop_sub()

contains

  integer function get_specie(label)
    character(len=*) :: label

    integer :: j
    logical :: l

    l = .false.
    do j = 1, ns
      if( trim(label) == trim(s(j)%label) ) then
        l = .true.
        get_specie = j
        return
      end if
    end do
    
    if(.not. l) then
      message(1) = "Specie '"+trim(label)+"' not found"
      call write_fatal(1)
    end if
    
  end function get_specie

end subroutine atom_init

subroutine loadPDB(iunit, natoms, a, ncatoms, ca)
  integer, intent(in) :: iunit
  integer, intent(inout) :: natoms, ncatoms
  type(atom_type), pointer :: a(:)
  type(atom_classical_type), pointer :: ca(:)
 
  character(len=80) :: record
  character(len=6) :: record_name
  character(len=4) :: atm
  character(len=3) :: res
  integer :: na, nca

    ! First count number of atoms
    rewind(iunit)
    natoms = 0
    ncatoms = 0
    do
      read(iunit, '(a80)', err=990, end=990) record
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(17x,a3)') res
        if(trim(res) == 'QM') then
          natoms = natoms + 1
        else
          ncatoms = ncatoms + 1
        end if
      end if      
    end do
990 continue

  allocate(a(natoms), ca(ncatoms))

  ! read in the data
  rewind(iunit)
  na = 1; nca = 1
  do
      read(iunit, '(a80)', err=991, end=991) record
      read(record, '(a6)') record_name
      if(trim(record_name) == 'ATOM' .or. trim(record_name) == 'HETATOM') then
        read(record, '(12x,a4,1x,a3)') atm, res
        call str_trim(atm)
        if(trim(res) == 'QM') then
          read(record, '(30x,3f8.3)') a(na)%x
          !a(na)%spec => s(get_specie(atm(1:1)))
          a(na)%label = atm(1:1)
          na = na + 1
        else
          ca(nca)%label = atm
          read(record, '(30x,3f8.3,6x,f6.2)') ca(nca)%x, ca(nca)%charge
          nca = nca + 1
        end if
      end if      
    end do
991 continue

end subroutine loadPDB

subroutine atom_end(na, a, nca, ca)
  integer, intent(in) :: na, nca
  type(atom_type), pointer :: a(:)
  type(atom_classical_type), pointer :: ca(:)

  call atom_dealloc(na, a)

  if(associated(a)) then ! sanity check
    deallocate(a); nullify(a)
  end if

  if(nca > 0 .and. associated(ca)) then
    deallocate(ca); nullify(ca)
  end if
end subroutine atom_end

subroutine atom_dealloc(na, a)
  integer, intent(in) :: na
  type(atom_type), pointer :: a(:)

  integer :: ia

  do ia = 1, na
    if(associated(a(ia)%Jxyz)) then
      deallocate(a(ia)%Jxyz, a(ia)%duV,    a(ia)%duVu,   a(ia)%dduV, &
                             a(ia)%zuV,    a(ia)%zuVu,   a(ia)%zduV, &
                             a(ia)%so_uV,  a(ia)%so_uVu, a(ia)%so_duV, &
                             a(ia)%so_luv)
      nullify(a(ia)%Jxyz, a(ia)%duV,    a(ia)%duVu,   a(ia)%duV, &
                          a(ia)%zuV,    a(ia)%zuVu,   a(ia)%zuV, &
                          a(ia)%so_uV,  a(ia)%so_uVu, a(ia)%so_duV, &
                          a(ia)%so_luv)
      if(conf%periodic_dim/=0 .and. associated(a(ia)%phases)) then
        deallocate(a(ia)%phases)
        nullify(a(ia)%phases)
      end if
    end if
  end do

end subroutine atom_dealloc

subroutine atom_adjust(natoms, a, ncatoms, ca)
  integer, intent(in) :: natoms, ncatoms
  type(atom_type), pointer :: a(:)
  type(atom_classical_type), pointer :: ca(:)

  real(r8) :: center(3), from(3), from2(3), to(3)
  character(len=80) :: str

  ! is there something to do
  if(natoms <= 1) return

  ! recenter
  call find_center(center)
  call translate(-center)

  ! get to axis
  str = "MainAxis"
  if(oct_parse_isdef(str) .ne. 0) then
    call oct_parse_block_double(str, 0, 0, to(1))
    call oct_parse_block_double(str, 0, 1, to(2))
    call oct_parse_block_double(str, 0, 2, to(3))
  else
    to(1) = 0._r8; to(2) = 0._r8; to(3) = 1._r8
  end if
  to = to / sqrt(sum(to**2))

  ! rotate to main axis
  call find_axis(from, from2)
  call rotate(from, from2, to)

  ! recenter
  call find_center(center)
  call translate(-center)
  
contains

  subroutine find_center(x)
    real(r8), intent(out) :: x(3)

    real(r8) :: xmin(3), xmax(3)
    integer  :: i, j

    xmin =  1e10_r8
    xmax = -1e10_r8
    do i = 1, natoms
      do j = 1, 3
        if(a(i)%x(j) > xmax(j)) xmax(j) = a(i)%x(j)
        if(a(i)%x(j) < xmin(j)) xmin(j) = a(i)%x(j)
      end do
    end do

    x = (xmax + xmin)/2._r8
  end subroutine find_center

  subroutine translate(x)
    real(r8), intent(in) :: x(3)

    integer  :: i

    do i = 1, natoms
      a(i)%x(:) = a(i)%x(:) + x(:)
    end do
    do i = 1, ncatoms
      ca(i)%x(:) = ca(i)%x(:) + x(:)
    end do
  end subroutine translate

  subroutine find_axis(x, x2)
    real(r8), intent(out) :: x(3), x2(3)

    integer  :: i, j
    real(r8) :: rmax, r, r2

    ! first get the further apart atoms
    rmax = -1e10_r8
    do i = 1, natoms
      do j = 1, natoms/2 + 1
        r = sqrt(sum((a(i)%x-a(j)%x)**2))
        if(r > rmax) then
          rmax = r
          x = a(i)%x - a(j)%x
        end if
      end do
    end do
    x  = x /sqrt(sum(x**2))

    ! now let us find out what is the second most important axis
    rmax = -1e10_r8
    do i = 1, natoms
      r2 = sum(x(:) * a(i)%x(:))
      r = sqrt(sum((a(i)%x(:) - r2*x(:))**2))
      if(r > rmax) then
        rmax = r
        x2 = a(i)%x(:) - r2*x(:)
      end if
    end do
    if(sum(x2**2) == M_ZERO) then ! linear molecule
      if(x(1) == M_ZERO) then
        x2(1) = x(1); x2(2) = -x(3); x2(3) = x(2)
      else if(x(1) == M_ZERO) then
        x2(2) = x(2); x2(1) = -x(3); x2(3) = x(1)
      else
        x2(3) = x(3); x2(1) = -x(2); x2(2) = x(1)
      end if
    end if
    x2 = x2/sqrt(sum(x2**2))

  end subroutine find_axis

  subroutine rotate(from, from2, to)
    real(r8), intent(in) :: from(3), from2(3), to(3) ! assumed to be normalize

    integer :: i
    real(r8) :: m1(3,3), m2(3,3), m3(3,3), f2(3), per(3)
    real(r8) :: alpha, r

    ! initialize matrices
    m1 = 0._r8; m1(1,1) = 1._r8; m1(2,2) = 1._r8; m1(3,3) = 1._r8

    ! rotate the to axis to the z axis
    if(to(2).ne.0._r8) then
      alpha = atan2(to(2), to(1))
      call rotate_z(m1, alpha)
    end if
    alpha = atan2(sqrt(to(1)**2 + to(2)**2), to(3))
    call rotate_y(m1, -alpha)

    ! get perpendicular to z and from
    f2 = matmul(m1, from)
    per(1) = -f2(2)
    per(2) =  f2(1)
    per(3) = 0._r8
    r = sqrt(sum(per**2))
    if(r > 0._r8) then
      per = per/r
    else
      per(2) = 1._r8
    end if

    ! rotate perpendicular axis to the y axis
    m2 = 0._r8; m2(1,1) = 1._r8; m2(2,2) = 1._r8; m2(3,3) = 1._r8
    alpha = atan2(per(1), per(2))
    call rotate_z(m2, -alpha)

    ! rotate from => to (around the y axis)
    m3 = 0._r8; m3(1,1) = 1._r8; m3(2,2) = 1._r8; m3(3,3) = 1._r8
    alpha = acos(sum(from*to))
    call rotate_y(m3, -alpha)
    
    ! join matrices
    m2 = matmul(transpose(m2), matmul(m3, m2))

    ! rotate around the z axis to get the second axis
    per = matmul(m2, matmul(m1, from2))
    alpha = atan2(per(1), per(2))
    call rotate_z(m2, -alpha) ! second axis is now y

    ! get combined transformation
    m1 = matmul(transpose(m1), matmul(m2, m1))

    ! now transform the coordinates
    ! it is written in this way to avoid what I consider a bug in the Intel compiler
    do i = 1, natoms
      f2 = a(i)%x
      a(i)%x = matmul(m1, f2)
    end do

    do i = 1, ncatoms
      f2 = ca(i)%x
      ca(i)%x = matmul(m1, f2)
    end do

  end subroutine rotate

  subroutine rotate_x(m, angle)
    real(r8), intent(inout) :: m(3,3)
    real(r8), intent(in) :: angle

    real(r8) :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = 0._r8
    aux(1, 1) = 1._r8
    aux(2, 2) = ca
    aux(3, 3) = ca
    aux(2, 3) = sa
    aux(3, 2) = -sa

    m = matmul(aux, m)
  end subroutine rotate_x

  subroutine rotate_y(m, angle)
    real(r8), intent(inout) :: m(3,3)
    real(r8), intent(in) :: angle

    real(r8) :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = 0._r8
    aux(2, 2) = 1._r8
    aux(1, 1) = ca
    aux(3, 3) = ca
    aux(1, 3) = sa
    aux(3, 1) = -sa

    m = matmul(aux, m)
  end subroutine rotate_y

  subroutine rotate_z(m, angle)
    real(r8), intent(inout) :: m(3,3)
    real(r8), intent(in) :: angle

    real(r8) :: aux(3,3), ca, sa

    ca = cos(angle)
    sa = sin(angle)

    aux = 0._r8
    aux(3, 3) = 1._r8
    aux(1, 1) = ca
    aux(2, 2) = ca
    aux(1, 2) = sa
    aux(2, 1) = -sa

    m = matmul(aux, m)
  end subroutine rotate_z

end subroutine atom_adjust

real(r8) function ion_ion_energy(natoms, atom)
  integer, intent(in)         :: natoms
  type(atom_type), intent(in) :: atom(natoms)
  
  real(r8) :: r
  integer :: i, j

  ! calculate the ion-ion energy
  ion_ion_energy = 0.0_r8

  do i = 2, natoms
    do j = 1, i - 1
      r = sqrt(sum((atom(i)%x - atom(j)%x)**2))
      ion_ion_energy = ion_ion_energy + &
           atom(i)%spec%Z_val*atom(j)%spec%Z_val/r
    end do
  end do

end function ion_ion_energy

real(r8) function kinetic_energy(natoms, atom)
  integer, intent(in)         :: natoms
  type(atom_type), intent(in) :: atom(natoms)

  integer :: i

  kinetic_energy = 0.0_r8
  do i = 1, natoms
     kinetic_energy = kinetic_energy + 0.5_r8*atom(i)%spec%weight* &
         sum(atom(i)%v(:)**2)
  enddo

end function kinetic_energy

subroutine cm_pos(natoms, atom, pos)
  integer, intent(in)         :: natoms
  type(atom_type), intent(in) :: atom(natoms)
  real(r8), intent(out) :: pos(3)

  real(r8) :: m
  integer :: i

  pos = M_ZERO; m = M_ZERO
  do i = 1, natoms
     pos(:) = pos(:) + atom(i)%spec%weight*atom(i)%x(:)
     m = m + atom(i)%spec%weight
  enddo
  pos = pos/m
end subroutine cm_pos

subroutine cm_vel(natoms, atom, vel)
  integer, intent(in)         :: natoms
  type(atom_type), intent(in) :: atom(natoms)
  real(r8), intent(out) :: vel(3)

  real(r8) :: m
  integer :: i

  vel = M_ZERO; m = M_ZERO
  do i = 1, natoms
     vel(:) = vel(:) + atom(i)%spec%weight*atom(i)%v(:)
     m = m + atom(i)%spec%weight
  enddo
  vel = vel/m
end subroutine cm_vel

subroutine atom_write_xyz(dir, fname, natoms, a, ncatoms, ca)
  character(len=*), intent(in)  :: dir, fname
  integer, intent(in) :: natoms, ncatoms
  type(atom_type), pointer :: a(:)
  type(atom_classical_type), pointer :: ca(:)

  integer i, iunit
  
#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif

    call oct_mkdir(trim(dir))

    call io_assign(iunit)
    open(iunit, file=trim(dir)+"/"+trim(fname)+'.xyz', status='unknown')
    write(iunit, '(i4)') natoms
    write(iunit, '(1x)')
    do i = 1, natoms
      write(iunit, '(6x,a,2x,3f12.6)') a(i)%label, a(i)%x(:)/units_out%length%factor
    end do
    call io_close(iunit)
    
    if(ncatoms > 0) then
      call io_assign(iunit)
      open(iunit, file=trim(dir)+"/"+trim(fname)+'_classical.xyz', status='unknown')
      write(iunit, '(i4)') ncatoms
      write(iunit, '(1x)')
      do i = 1, ncatoms
        write(iunit, '(6x,a1,2x,3f12.6,a,f12.6)') &
             ca(i)%label(1:1), ca(i)%x(:)/units_out%length%factor, &
             " # ", ca(i)%charge
      end do
      call io_close(iunit)
    end if

#ifdef HAVE_MPI
  end if
#endif
  
  return
end subroutine atom_write_xyz

end module atom
