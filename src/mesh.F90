#include "config.h"

module mesh
use global
use geom
  
implicit none

type mesh_type
  integer :: box_shape ! 0->sphere, 1->cilinder, 2->sphere around each atom
  real(r8)    :: h         ! the (constant) spacing between the points
  real(r8)    :: vol_pp    ! the volume per point (normally h**3)
  
  real(r8)    :: rsize     ! the radius of the sphere or of the cilinder
  real(r8)    :: zsize     ! the length of the cilinder in the z direction
  
  integer :: nl        ! number of points in inner mesh
  
  ! return x, y and z for each point
  integer, pointer :: Lx(:), Ly(:), Lz(:)
  ! return points # for each xyz
  integer, pointer :: Lxyz_inv(:,:,:)

  ! mesh elargement... this is currently only used in solving
  ! the poisson equation with conjugated gradients...
  integer :: nk        ! number of points in outer mesh
  integer, pointer :: Kx(:), Ky(:), Kz(:)

  ! some other vars
  integer :: nx ! the # of points in the whole radius
                    ! (rsize/h) + norder

end type mesh_type

private in_sphere
private in_cilinder
private in_atom

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! creates the mesh
! assumes that the box_shape, h, rsize/zsize of m are filled
! if norder is present, then set nk, Kx, Ky and Kz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mesh_init(m, n_extra, nions, ion)
  type(mesh_type), intent(inout) :: m
  integer, intent(in), optional :: n_extra
  integer, intent(in), optional :: nions
  type(geom_type),  pointer, optional :: ion(:)
  
  ! some local stuff
  character(len=40) :: err
  integer :: next = 0
  integer :: il, ik, ix, iy, iz, nr, nx
  integer, pointer :: Lxyz(:,:,:)
  logical :: b

  sub_name = 'mesh_init'; call push_sub()

  ! for the outer mesh
  if(present(n_extra)) next = n_extra

  ! first we do some bound checks...
  if (m%box_shape>3 .or. m%box_shape<1) then
    write(err, *) m%box_shape
    message(1) = "Input: '"//trim(err)//"' is not a valid Box_Shape"
    message(2) = '(1 <= Box_Shape <= 3)'
    call write_fatal(2)
  end if
  if(m%box_shape == 3 .and. (.not.present(nions) .or. .not.present(ion))) then
    call write_fatal(1)
  endif

  if (m%h<0.01_r8 .or. m%h>1.0_r8) then
    write(err, *) m%h
    message(1) = "Input: '"//trim(err)//"' is not a valid spacing"
    message(2) = '(0.01 <= Spacing <= 1)'
    call write_fatal(2)
  end if
  if (m%rsize<1.0_r8 .or. m%rsize>250.0_r8) then
    write(err, *) m%rsize
    message(1) = "Input: '"//trim(err)//"' is not a valid radius"
    message(2) = '(1 <= radius <= 250)'
    call write_fatal(2)
  end if
  if (m%zsize<1.0_r8 .or. m%zsize>250.0_r8) then
    write(err, *) m%zsize
    message(1) = "Input: '"//trim(err)//"' is not a valid zlength"
    message(2) = '(1 <= zlength <= 250)'
    call write_fatal(2)
  end if

  ! set vol_pp and nx
  m%vol_pp = m%h**3
  nr = max(int(m%rsize/m%h), int(m%zsize/m%h))
  nx = nr + n_extra
  m%nx = nx

  ! allocate the xyz arrays
  allocate(m%Lxyz_inv(-Nx:Nx,-Nx:Nx,-Nx:Nx))
  allocate(Lxyz(-Nx:Nx,-Nx:Nx,-Nx:Nx))
  m%Lxyz_inv = 0
  Lxyz = 0

  ! We label 1 the points inside the inner mesh
  ! and 2 for the points in the outer mesh
  do ix = -nr, nr
    do iy = -nr, nr
      do iz = -nr, nr
        select case(m%box_shape)
        case(1)
          b = in_sphere(m, ix, iy, iz)
        case(2)
          b = in_cilinder(m, ix, iy, iz)
        case(3)
          b = in_atom(m, ix, iy, iz, nions, ion)
        end select
        
        if(b) then
          Lxyz(ix-n_extra:ix+n_extra, iy, iz) = 2_i4
          Lxyz(ix, iy-n_extra:iy+n_extra, iz) = 2_i4
          Lxyz(ix, iy, iz-n_extra:iz+n_extra) = 2_i4
        endif
      enddo
    enddo
  enddo

  do ix = -nr, nr
    do iy = -nr, nr
      do iz = -nr, nr
        select case(m%box_shape)
        case(1)
          b = in_sphere(m, ix, iy, iz)
        case(2)
          b = in_cilinder(m, ix, iy, iz)
        case(3)
          b = in_atom(m, ix, iy, iz, nions, ion)
        end select
        
        if(b) then
          Lxyz(ix, iy, iz) = 1_i4
        endif
      enddo
    enddo
  enddo

  ! now, we count the number of internal and external points
  il=0
  ik=0
  do ix = -Nx, Nx
    do iy = -Nx, Nx
      do iz = -Nx, Nx
        if(Lxyz(ix,iy,iz).eq.1) il=il+1
        if(Lxyz(ix,iy,iz).eq.2) ik=ik+1
      enddo
    enddo
  enddo
  m%nl = il
  m%nk = ik

  allocate(m%Lx(il), m%Ly(il), m%Lz(il))
  allocate(m%Kx(ik), m%Ky(ik), m%Kz(ik))

  il=0
  ik=0
  do ix = -Nx, Nx
    do iy = -Nx, Nx
      do iz = -Nx, Nx
        if(Lxyz(ix,iy,iz).eq.1) then
          il = il + 1
          m%Lx(il) = ix
          m%Ly(il) = iy
          m%Lz(il) = iz
          m%Lxyz_inv(ix,iy,iz) = il
        endif
        if(Lxyz(ix,iy,iz).eq.2) then
          ik = ik + 1
          m%Kx(ik) = ix
          m%Ky(ik) = iy
          m%Kz(ik) = iz
        endif
      enddo
    enddo
  enddo

  deallocate(Lxyz); nullify(Lxyz)

  call pop_sub()
  return
end subroutine mesh_init

! The three next functions are self-explanatory!
! They check if a point is within the sphere, cilinder, 
! or atomic spheres
logical function in_sphere(m, ix, iy, iz)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: ix, iy, iz

  in_sphere = (sqrt(real(ix**2+iy**2+iz**2,r8))*m%H <= m%rsize)
  return
end function in_sphere

logical function in_cilinder(m, ix, iy, iz)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: ix, iy, iz

  real(r8) r, z
  r = sqrt(real(ix**2+iy**2, r8))*m%h
  z = iz*m%H

  in_cilinder = (r<=m%rsize .and. abs(z)<=m%zsize)
  return
end function in_cilinder

logical function in_atom(m, ix, iy, iz, nions, ion)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: ix, iy, iz
  integer, intent(in), optional :: nions
  type(geom_type),  pointer, optional :: ion(:)

  integer :: i
  real(r8) :: r

  in_atom = .false.
  ion_loop: do i = 1, nions
    r = sqrt((ix*m%H-ion(i)%x(1))**2+(iy*m%H-ion(i)%x(2))**2+(iz*m%H-ion(i)%x(3))**2)
    if(r <= m%rsize) then
      in_atom = .true.
      exit ion_loop
    end if
  end do ion_loop

  return
end function in_atom

! Deallocates what has to be deallocated ;)
subroutine mesh_end(m)
  type(mesh_type), intent(inout) :: m

  sub_name = 'mesh_end'; call push_sub()

  if(associated(m%Lx)) then
    deallocate(m%Lx, m%Ly, m%Lz, m%Kx, m%Ky, m%Kz)
    nullify(m%Lx, m%Ly, m%Lz, m%Kx, m%Ky, m%Kz)
  end if
  if(associated(m%Lxyz_inv)) then
    deallocate(m%Lxyz_inv)
    nullify(m%Lxyz_inv)
  end if
    
  call pop_sub()
  return
end subroutine mesh_end

subroutine mesh_write_info(m, unit)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: unit

  character(len=15), parameter :: bs(3) = (/ &
      'sphere       ', &
      'cilinder     ', &
      'around nuclei'/)

  sub_name = 'mesh_write_info'; call push_sub()

  write(message(1), '(a,a,1x,a,f7.3)') '  Type = ', bs(m%box_shape), &
      ' Radius = ', m%rsize
  if(m%box_shape == 2) then
    write(message(1), '(a,a, f7.3)') trim(message(1)), ', zlength = ', m%zsize
  end if
  write(message(2),'(a, f6.3, 5x, a, f8.5)') '  Spacing = ', m%h, &
      '   volume/point = ', m%vol_pp
  write(message(3),'(a, i6, a, i6)') '  # inner mesh = ', m%nl, &
      '   # outer mesh = ', m%nk

  call write_info(3, unit)

  call pop_sub()
  return
end subroutine mesh_write_info

! TODO
subroutine mesh_convert(m_from, v_from, m_to, v_to)
  type(mesh_type), intent(IN) :: m_from, m_to
  real(r8), intent(in)  :: v_from
  real(r8), intent(out) :: v_to

end subroutine mesh_convert

end module mesh
