module hamiltonian
use global
use spline
use fft
use system

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

end type hamiltonian_type

contains

subroutine hamiltonian_init(h, ispin, np)
  type(hamiltonian_type), intent(out) :: h
  integer, intent(in) :: ispin, np

  h%ispin = ispin
  h%np  = np

  allocate(h%Vpsl(np), h%rho_core(np))
  if(ispin == 1) then
    allocate(h%Vhxc(np, 1))
  else
    allocate(h%Vhxc(np, 2))
  end if

  h%vpsl_space = fdf_integer('LocalPotentialSpace', 0)
  if(h%vpsl_space < 0 .or. h%vpsl_space > 1) then
    write(message(1), '(a,i5,a)') "Input: '", h%vpsl_space, &
         "' is not a valid LocalPotentialSpace"
    message(2) = '(0 <= LocalPotentialSpace <=1)'
    call write_fatal(2)
  end if

end subroutine hamiltonian_init

subroutine hamiltonian_end(h)
  type(hamiltonian_type) :: h

  if(associated(h%Vpsl)) then
    deallocate(h%Vpsl, h%Vhxc)
    nullify(h%Vpsl, h%Vhxc)
  end if
end subroutine hamiltonian_end

subroutine generate_external_pot(h, sys)
  type(hamiltonian_type), intent(inout) :: h
  type(system_type), intent(inout) :: sys

  integer :: ia
  type(specie_type), pointer :: s
  type(atom_type),   pointer :: a
  complex(r8), allocatable :: fw(:,:,:), fwc(:,:,:)
  real(r8), allocatable :: fr(:,:,:)

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
    call rfftwnd_f77_one_complex_to_real(sys%m%zplanb, fw, fr)
    call dcube_to_mesh(sys%m, fr, h%vpsl)
    call rfftwnd_f77_one_complex_to_real(sys%m%zplanb, fwc, fr)
    call dcube_to_mesh(sys%m, fr, h%rho_core)

    deallocate(fw, fwc, fr)
  end if
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
      call mesh_r(m, i, a%x, r)
      if(r <= s%jradius) then
        h%Vpsl(i) = h%Vpsl(i) + (a1*(r*r - Rb2) - a2) * P_E2
      else
        h%Vpsl(i) = h%Vpsl(i) - s%Z/r * P_E2
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
        call mesh_r(m, i, a%x, r)
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
      call mesh_r(m, k, a%x, r)
      if(r < s%ps%rc_max) j = j + 1
    end do
    a%Mps = j
    
    allocate(a%Jxyz(j), a%uV(j, (s%ps%L_max+1)**2), &
         a%duV(j, (s%ps%L_max+1)**2, 3), a%uVu((s%ps%L_max+1)**2))
    
    a%uV  = 0.0_r8; a%duV = 0.0_r8; a%uVu = 0.0_r8
    
    j = 0
    do k = 1, h%np
      call mesh_r(m, k, a%x, r)
      if(r < s%ps%rc_max) then
        j = j + 1
        a%Jxyz(j) = k
      end if
    end do

  end subroutine build_kb_sphere

end subroutine generate_external_pot

end module hamiltonian
