module hamiltonian
use global
use spline
use fft
use system
use hartree
use xc

implicit none

type hamiltonian_type
  integer :: ispin ! How to handle spin (duplicated in states_type)
  integer :: np  ! number of points (duplicated in mesh)

  integer :: soc ! spin-orbit coupling?

  integer :: vpsl_space           ! How should the local potential be calculated
  real(r8), pointer :: Vpsl(:)    ! the external potential
  real(r8), pointer :: Vhxc(:,:)  ! Hartre+xc potential
  real(r8), pointer :: rho_core(:)! core charge for nl core corrections

  ! the energies (total, ion-ion, exchange, correlation)
  real(r8) :: etot, eii, ex, ec

  ! System under the independent particle approximation, or not.
  logical :: ip_app

  ! hartree potential structure
  type(hartree_type) :: hart
  type(xc_type) :: xc

end type hamiltonian_type

contains

subroutine hamiltonian_init(h, sys)
  type(hamiltonian_type), intent(out) :: h
  type(system_type), intent(inout) :: sys

  sub_name = 'hamiltonian_init'; call push_sub()

  ! Duplicate this two variables
  h%ispin = sys%st%ispin
  h%np  = sys%m%np

  ! allocate potentials and density of the cores
  allocate(h%Vpsl(h%np), h%Vhxc(h%np, h%ispin), h%rho_core(h%np))

  ! should we calculate the local pseudopotentials in Fourier space?
  h%vpsl_space = fdf_integer('LocalPotentialSpace', 0)
  if(h%vpsl_space < 0 .or. h%vpsl_space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", h%vpsl_space, &
         "' is not a valid LocalPotentialSpace"
    message(2) = '(0 <= LocalPotentialSpace <=1)'
    call write_fatal(2)
  end if

  if(h%vpsl_space == 1) then
    call specie_fourier_init(sys%nspecies, sys%specie, sys%m)
  end if

  call hartree_init(h%hart, sys%m)
  call xc_init(h%xc, sys%m, h%ispin)

  call pop_sub()
end subroutine hamiltonian_init

subroutine hamiltonian_end(h)
  type(hamiltonian_type) :: h

  sub_name = 'hamiltonian_end'; call push_sub()

  if(associated(h%Vpsl)) then
    deallocate(h%Vpsl, h%Vhxc)
    nullify(h%Vpsl, h%Vhxc)
  end if

  call hartree_end(h%hart)
  call xc_end(h%xc)

  call pop_sub()
end subroutine hamiltonian_end

! This subroutine should be in specie.F90, but due to the limitations
! of f90 to handle circular dependences it had to come here!
subroutine specie_fourier_init(ns, s, m)
  integer, intent(in) :: ns
  type(specie_type), pointer :: s(:)
  type(mesh_type), intent(IN) :: m

  integer :: i, j, ix, iy, iz
  real(r8) :: r, vl
  real(r8), allocatable :: fr(:,:,:)

  sub_name = 'specie_init_fourier'; call push_sub()

  allocate(fr(m%fft_n,  m%fft_n, m%fft_n))
  do i = 1, ns

    allocate(s(i)%local_fw(m%hfft_n, m%fft_n, m%fft_n))
    fr = 0.0_r8
    do j = 1, m%np
      call mesh_r(m, j, r)
      ix = m%lx(j) + m%fft_n/2 + 1;
      iy = m%ly(j) + m%fft_n/2 + 1;
      iz = m%lz(j) + m%fft_n/2 + 1;
      vl  = splint(s(i)%ps%vlocal, r)
      if(r >= r_small) then
        fr(ix, iy, iz) = (vl - s(i)%Z_val)/r
      else
        fr(ix, iy, iz) = s(i)%ps%vlocal_origin
      endif
    enddo
    call rfftwnd_f77_one_real_to_complex(m%dplanf, fr, s(i)%local_fw)
    call zscal(m%fft_n**2*m%hfft_n, cmplx(1.0_r8/m%fft_n**3, 0.0_r8, r8), s(i)%local_fw, 1)

    if(s(i)%ps%icore /= 'nc  ') then
      fr = 0._r8
      do j = 1, m%np
        call mesh_r(m, j, r)
        ix = m%lx(j) + m%fft_n/2 + 1;
        iy = m%ly(j) + m%fft_n/2 + 1;
        iz = m%lz(j) + m%fft_n/2 + 1;
        vl  = splint(s(i)%ps%core, r)
        fr(ix, iy, iz) = vl
      enddo
      call rfftwnd_f77_one_real_to_complex(m%dplanf, fr, s(i)%rhocore_fw)
      call zscal(m%fft_n**2*m%hfft_n, cmplx(1.0_r8/m%fft_n**3, 0.0_r8, r8), s(i)%rhocore_fw, 1)
    end if
  end do
  
  call pop_sub()
end subroutine specie_fourier_init

subroutine generate_external_pot(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  integer :: ia
  type(specie_type), pointer :: s
  type(atom_type),   pointer :: a
  complex(r8), allocatable :: fw(:,:,:), fwc(:,:,:)
  real(r8), allocatable :: fr(:,:,:)

  ! WARNING DEBUG
!!$  integer :: i, j
!!$  real(r8), allocatable :: f(:,:,:)

  sub_name = 'generate_external_pot'; call push_sub()

  if(h%vpsl_space == 1) then
    allocate(fw(sys%m%hfft_n, sys%m%fft_n, sys%m%fft_n), &
         fwc(sys%m%hfft_n, sys%m%fft_n, sys%m%fft_n))
    fw = M_z0; fwc = M_z0
  end if

  h%Vpsl = 0._r8
  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec

    select case(s%label(1:5))
    case('jelli', 'point')
      call from_jellium(sys%m)
    case default
      call from_pseudopotential(sys%m)
    end select
  end do
  
  if(h%vpsl_space == 1) then
    allocate(fr(sys%m%fft_n, sys%m%fft_n, sys%m%fft_n))
    call rfftwnd_f77_one_complex_to_real(sys%m%dplanb, fw, fr)
    call dcube_to_mesh(sys%m, fr, h%vpsl)
    call rfftwnd_f77_one_complex_to_real(sys%m%dplanb, fwc, fr)
    call dcube_to_mesh(sys%m, fr, h%rho_core)

    deallocate(fw, fwc, fr)
  end if

  ! WARNING DEBUG
!!$  allocate(f(sys%m%fft_n, sys%m%fft_n, sys%m%fft_n))
!!$  call dmesh_to_cube(sys%m, h%Vpsl, f)
!!$  do i = 1, sys%m%fft_n
!!$    do j = 1, sys%m%fft_n
!!$      print *, i, j, f(i, j, sys%m%fft_n/2 + 1)
!!$    end do
!!$    print *
!!$  end do
  
  call pop_sub()
contains
  !***************************************************
  !  jellium stuff
  !***************************************************
  subroutine from_jellium(m)
    type(mesh_type), intent(in) :: m

    integer :: i
    real(r8) :: a1, a2, Rb2, r

    a1 = s%Z/(2._r8*s%jradius**3)
    a2 = s%Z/s%jradius
    Rb2= s%jradius**2
    do i = 1, h%np
      call mesh_r(m, i, r, a=a%x)
      if(r <= s%jradius) then
        h%Vpsl(i) = h%Vpsl(i) + (a1*(r*r - Rb2) - a2)
      else
        h%Vpsl(i) = h%Vpsl(i) - s%Z/r
      end if
    end do
    
  end subroutine from_jellium

  !***************************************************
  !  pseudopotential stuff
  !***************************************************
  subroutine from_pseudopotential(m)
    type(mesh_type), intent(in) :: m

    call build_local_part(m)
    call build_kb_sphere(m)

  end subroutine from_pseudopotential

  ! this actually adds to outp
  subroutine phase_factor(m, vec, inp, outp)
    implicit none
    type(mesh_type), intent(IN) :: m
    real(r8), intent(IN)        :: vec(3)
    complex(r8), intent(IN)     :: inp(m%hfft_n, m%fft_n, m%fft_n)
    complex(r8), intent(inout)  :: outp(m%hfft_n, m%fft_n, m%fft_n)
    
    complex(r8) :: k
    integer     :: ix, iy, iz, ixx, iyy, izz
    
    k = M_zI * ((2.0_r8*M_Pi)/(m%fft_n*m%h))
    do iz = 1, m%fft_n
      izz = pad_feq(iz, m%fft_n, .true.)
      do iy = 1, m%fft_n
        iyy = pad_feq(iy, m%fft_n, .true.)
        do ix = 1, m%hfft_n
          ixx = pad_feq(ix, m%fft_n, .true.)
          outp(ix, iy, iz) = outp(ix, iy, iz) + &
               exp( k * (vec(1)*ixx + vec(2)*iyy + vec(3)*izz) ) * inp(ix, iy, iz)
        end do
      end do
    end do
  end subroutine phase_factor

  subroutine build_local_part(m)
    type(mesh_type), intent(in) :: m

    integer :: i
    real(r8) :: r, vl

    if(h%vpsl_space == 0) then
      do i = 1, m%np
        call mesh_r(m, i, r, a=a%x)
        vl  = splint(s%ps%vlocal, r)
        if(r >= r_small) then
          h%Vpsl(i) = h%Vpsl(i) + (vl - s%Z_val)/r
        else
          h%Vpsl(i) = h%Vpsl(i) + s%ps%Vlocal_origin
        end if
        if(s%ps%icore /= 'nc  ' ) then
          h%rho_core(i) = h%rho_core(i) + splint(s%ps%core, r)
        end if
      end do
    else
      call phase_factor(m, a%x, s%local_fw, fw)
      if(s%ps%icore /= 'nc  ') then
        call phase_factor(m, a%x, s%rhocore_fw, fwc)
      end if
    end if

  end subroutine build_local_part
  
  subroutine build_kb_sphere(m)
    type(mesh_type), intent(in) :: m
    
    integer :: j, k
    real(r8) :: r
    
    j = 0
    do k = 1, h%np
      call mesh_r(m, k, r, a=a%x)
      if(r < s%ps%rc_max) j = j + 1
    end do
    a%Mps = j
    
    allocate(a%Jxyz(j), a%uV(j, (s%ps%L_max+1)**2), &
         a%duV(j, (s%ps%L_max+1)**2, 3), a%uVu((s%ps%L_max+1)**2))
    
    a%uV  = 0.0_r8; a%duV = 0.0_r8; a%uVu = 0.0_r8
    
    j = 0
    do k = 1, h%np
      call mesh_r(m, k, r, a=a%x)
      if(r < s%ps%rc_max) then
        j = j + 1
        a%Jxyz(j) = k
      end if
    end do

  end subroutine build_kb_sphere

end subroutine generate_external_pot

end module hamiltonian
