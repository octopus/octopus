#include "config.h"

module timedep
use system
use states
use hamiltonian
use mix
use lasers

implicit none

type td_type
  complex(r8), pointer :: zpsi(:,:,:) ! the complex wavefunctions

  integer :: max_iter  ! maximum number of iterations to perform
  integer :: save_iter ! save every save_iter iterations
  integer :: iter      ! the actual iteration
  integer :: evolution_method ! which evolution method to use

  real(r8) :: dt            ! time step
  integer  :: move_ions     ! how do we move the ions?

  complex(r8), pointer :: kin_2(:,:,:) ! for split operator

  real(r8) :: delta_strength ! strength of the delta excitation
  real(r8) :: pol(3)         ! direction of the polarization of the efield

  integer :: lmax        ! maximum multipole moment to write

  integer :: gauge ! in which gauge shall we work in
                   ! 1 = length gauge
                   ! 2 = velocity gauge

  ! lasers stuff
  integer :: no_lasers ! number of laser pulses used
  logical :: output_laser ! write laser field
  type(laser_type), pointer :: lasers(:)

  ! absorbing boundaries
  integer  :: ab         ! do we have absorbing boundaries?
  real(r8) :: ab_width   ! width of the absorbing boundary
  real(r8) :: ab_height  ! height of the absorbing boundary
  real(r8), pointer :: ab_pot(:) ! where we store the ab potential

  ! occupational analysis
  logical :: occ_analysis ! do we perform occupational analysis?

!#ifndef NO_PES
!  logical :: calc_PES_rc    ! PES using rc method
!  type(PES_rc_type) :: PES_rc
!  logical :: calc_PES_mask  ! PES using mask method
!  type(PES_mask_type) :: PES_mask
!#endif

  real(r8), pointer :: v_old1(:,:), v_old2(:,:)
end type td_type

contains

subroutine td_run(td, sys, h)
  type(td_type), intent(inout) :: td
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h

  integer :: i, ii
  real(r8), allocatable :: dipole(:,:), multipole(:,:,:)

  sub_name = 'systm_td'; call push_sub()

  if(td%iter == 0) call td_run_zero_iter(sys%m)
  td%iter = td%iter + 1

  write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Dipole(is) '
  call write_info(1)
  
  allocate(dipole(sys%st%nspin, td%save_iter))
  allocate(multipole((td%lmax + 1)**2, sys%st%nspin, td%save_iter))

  ii = 1
  do i = td%iter, td%max_iter
    if(clean_stop()) exit

    ! time iterate wavefunctions
    call td_rti(sys, h, td, i*td%dt)

    ! update density
    call zcalcdens(sys%st, sys%m%np, sys%st%rho, reduce=.true.)

    ! update hamiltonian and eigenvalues (fermi is *not* called)
    call hamiltonian_setup(h, sys)
    call zhamiltonian_eigenval (h, sys, 1, sys%st%nst) ! eigenvalues
    call hamiltonian_energy(h, sys, -1, reduce=.true.)

    ! measuring
    call states_calculate_multipoles(sys%m, sys%st, td%pol, td%lmax, &
         dipole(:,ii), multipole(:,:,ii))

    write(message(1), '(i7,1x,f14.5,f14.6,4e17.6)') i, &
         i*td%dt       / units_out%time%factor, &
         h%etot        / units_out%energy%factor, &
         dipole(:, ii) / units_out%length%factor
    call write_info(1)
  end do

  deallocate(dipole, multipole)
contains

  subroutine td_run_zero_iter(m)
    type(mesh_type), intent(IN) :: m

    integer :: i, iunit
    real(r8) :: x(3)    
    complex(r8) :: c
    real(r8), allocatable :: dipole(:), multipole(:,:)
    
    do i = 1, min(2, sys%st%ispin)
      td%v_old1(:, i) = h%Vhartree(:) + h%Vxc(:, i)
      td%v_old2(:, i) = td%v_old1(:, i)
    end do
      
    ! we now apply the delta(0) impulse to the wf
    if(td%delta_strength .ne. 0._r8) then
      do i = 1, m%np
        call mesh_xyz(m, i, x)
        c = exp(M_zI * td%delta_strength * sum(x(:)*td%pol(:)))

        sys%st%zpsi(i,:,:,:) = c * sys%st%zpsi(i,:,:,:)
      end do
    end if

    ! output static dipole (iter 0)
    
    allocate(dipole(sys%st%nspin), multipole((td%lmax + 1)**2, sys%st%nspin))
    call states_calculate_multipoles(sys%m, sys%st, td%pol, td%lmax, &
         dipole, multipole)
    call io_assign(iunit)
    open(iunit, status='unknown', file=trim(sys%sysname)//'.mult')
    call td_write_multipole(iunit, 0_i4, 0._r8, dipole, multipole, .true.)
    call io_close(iunit)

    if(td%output_laser) then
      call io_assign(iunit)
      open(unit=iunit, file='laser.out', status='unknown')
      write(iunit, '(a,e17.10)') '# dt = ', td%dt
      call io_close(iunit)
    end if

    ! output initial positions and velocities
!!$    if(td%move_ions > 0) then
!!$      call io_assign(iunit)
!!$      open(iunit, file=trim(sys%sysname)//'.nbo')
!!$      write(iunit, trim(nbo_file_format)) 0.0_r8, ekin, sys%etot, sys%etot + ekin, &
!!$           (sys%ion(j)%x, j=1, sys%nions),                    &
!!$           (sys%ion(j)%v, j=1, sys%nions),                    &
!!$           (forces(1:3, j), j=1, sys%nions)
!!$      call io_close(iunit)
!!$    end if
    
  end subroutine td_run_zero_iter

  subroutine td_write_multipole(iunit, iter, t, dipole, multipole, header)
    integer, intent(in) :: iunit, iter
    real(r8), intent(IN) :: t, dipole(sys%st%nspin), multipole((td%lmax+1)**2, sys%st%nspin)
    logical, intent(in) :: header

    integer :: is, l, m, add_lm
    
    if(header) then
      write(iunit, '(a10,i2,a8,i2)') '# nspin = ', sys%st%nspin, ' lmax = ', td%lmax
      write(iunit, '(a8, a20)', advance='no') 'Iter', 't [eV^-1]    '
      do is = 1, sys%st%nspin
        write(iunit, '(a13,i1,a6)', advance='no') 'dipole(', is, ')     '
      end do
      do is = 1, sys%st%nspin
        do l = 0, td%lmax
          do m = -l, l
            write(iunit, '(a6,i2,a4,i2,a2,i1,a3)', advance='no') 'l=', l, ', m=', m, ' (', is,')   '
            add_lm = add_lm + 1
          end do
        end do
      end do
    
    write(iunit, '(1x)', advance='yes')
  end if

  write(iunit, '(i8, e20.12)', advance='no') iter, t
  do is = 1, sys%st%nspin
    write(iunit, '(e20.12)', advance='no') dipole(is)
  end do

  do is = 1, sys%st%nspin
    add_lm = 1
    do l = 0, td%lmax
      do m = -l, l
        write(iunit, '(e20.12)', advance='no') multipole(add_lm, is)
        add_lm = add_lm + 1
      end do
    end do
  end do
  write(iunit, '(1x)', advance='yes')

  return
end subroutine td_write_multipole
  
end subroutine td_run

#include "td_init.F90"
#include "td_rti.F90"

end module timedep
