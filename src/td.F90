#include "config.h"

module timedep
use global
use io
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
  character(len=100) :: filename ! name of the continuation file
end type td_type

contains

subroutine td_run(td, u_st, sys, h)
  type(td_type), intent(inout) :: td
  type(states_type), intent(IN) :: u_st
  type(system_type), intent(inout) :: sys
  type(hamiltonian_type), intent(inout) :: h

  integer :: i, ii, idim, ist, ik
  real(r8), allocatable :: dipole(:,:), multipole(:,:,:)
  complex(r4), allocatable :: projections(:,:,:,:)
  character(len=100) :: proj_filename

  sub_name = 'systm_td'; call push_sub()

  write(message(1), '(a7,1x,a14,a14,a17)') 'Iter ', 'Time ', 'Energy ', 'Dipole(is) '
  call write_info(1)
  
  allocate(dipole(sys%st%nspin, td%save_iter))
  allocate(multipole((td%lmax + 1)**2, sys%st%nspin, td%save_iter))

  ! occupational analysis stuff
  if(td%occ_analysis) then
    write(proj_filename, '(a,a,i3.3,a)') trim(sys%sysname), &
         '.', mpiv%node, '.proj'
    allocate(projections(u_st%nst, sys%st%st_start:sys%st%st_end, sys%st%nik, td%save_iter))
  end if

  if(td%iter == 0) call td_run_zero_iter(sys%m)
  td%iter = td%iter + 1

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
    if(td%occ_analysis) then
      call td_calc_projection(projections(:,:,:,ii))
    end if

    ! mask function?
    if(td%ab == 2) then
      do ik = 1, sys%st%nik
        do ist = sys%st%st_start, sys%st%st_end
          do idim = 1, sys%st%dim
            sys%st%zpsi(1:sys%m%np, idim, ist, ik) = sys%st%zpsi(1:sys%m%np, idim, ist, ik) * &
                 (1._r8 - td%ab_pot(1:sys%m%np))
          end do
        end do
      end do
    end if

    ! write info
    write(message(1), '(i7,1x,f14.5,f14.6,4es17.6)') i, &
         i*td%dt       / units_out%time%factor, &
         h%etot        / units_out%energy%factor, &
         dipole(:, ii) / units_out%length%factor
    call write_info(1)

    ! write down data
    ii = ii + 1
    save: if(ii==td%save_iter+1 .or. i == td%max_iter) then ! output
      if(i == td%max_iter) td%save_iter = ii - 1
      ii = 1

      ! first resume file
      call zstates_write_restart(trim(td%filename), sys%m, sys%st, &
           iter=i, v1=td%v_old1, v2=td%v_old2)

      if(mpiv%node == 0) call td_write_data()
    end if save

  end do

  deallocate(dipole, multipole)
  if(td%occ_analysis) then
    deallocate(projections)
  end if
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
      call td_write_laser(iunit, 0, 0._r8, .true.)
      call io_close(iunit)
    end if
    
    if(td%occ_analysis) then
      call io_assign(iunit)
      open(iunit, form='unformatted', status='unknown', file=proj_filename)
      write(iunit) sys%st%nik, sys%st%st_start, sys%st%st_end, u_st%nst
      call io_close(iunit)
    end if
    ! output initial positions and velocities
!TODO
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

  subroutine td_write_data()
    integer :: iunit, j, jj, ist, ik, uist

    ! output multipoles
    call io_assign(iunit)
    open(iunit, position='append', file=trim(sys%sysname)//".mult")
    do j = 1, td%save_iter
      jj = i - td%save_iter + j
      call td_write_multipole(iunit, jj, jj*td%dt, &
           dipole(:, j), multipole(:,:, j), .false.)
    end do
    call io_close(iunit)

    ! and now we should output the projections
    if(td%occ_analysis) then
      call io_assign(iunit)
      open(iunit, form='unformatted', position='append', file=proj_filename)
      do j = 1, td%save_iter
        do ik = 1, sys%st%nik
          do ist = 1, sys%st%st_start, sys%st%st_end
            write(iunit) (projections(uist, ist, ik, j), uist = 1, u_st%nst)
          end do
        end do
      end do
      call io_close(iunit)
    end if

    ! output the laser field
    if(td%output_laser) then
      call io_assign(iunit)
      open(iunit, position='append', file='laser.out')
      do j = 1, td%save_iter
        jj = i - td%save_iter + j
        call td_write_laser(iunit, jj, jj*td%dt, .false.)
      end do
      call io_close(iunit)
    end if

  end subroutine td_write_data

  subroutine td_write_laser(iunit, iter, t, header)
    integer, intent(in) :: iunit, iter
    real(r8), intent(IN) :: t
    logical, intent(in) :: header

    integer :: i
    real(r8) :: l_field(1:3), l_vector_field(1:3)
    character(len=20) :: aux

    ! TODO -> confirm these stupid units, especially for the vector field
    if(header) then
      ! first line
      write(iunit, '(a7,e20.12,3a)') '# dt = ', td%dt/units_out%time%factor, &
           "[", trim(units_out%time%abbrev), "]"

      ! second line
      write(iunit, '(a8,a20)', advance='no') '# Iter  ', str_center('t', 20)
      do i = 1, 3
        write(aux, '(a,i1,a)') 'E(', i, ')'
        write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
      end do
      do i = 1, 3
        write(aux, '(a,i1,a)') 'A(', i, ')'
        write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
      end do
      write(iunit, '(1x)', advance='yes')

      ! third line
      write(iunit, '(a8,a20)', advance='no') '#       ', &
           str_center('['//trim(units_out%time%abbrev)//']', 20)
      aux = str_center('['//trim(units_out%energy%abbrev) // ' / ' // &
           trim(units_inp%length%abbrev) // ']', 20)
      write(iunit, '(3a20)', advance='no') aux, aux, aux
      aux = str_center('[1/'// trim(units_inp%length%abbrev) // ']', 20)
      write(iunit, '(3a20)', advance='no') aux, aux, aux
      write(iunit, '(1x)', advance='yes')
    end if

    call laser_field(td%no_lasers, td%lasers, t, l_field)
    call laser_vector_field(td%no_lasers, td%lasers, t, l_vector_field)
    write(iunit,'(i8,7es20.12)') iter, t/units_out%time%factor, &
         l_field(1:3) * units_inp%length%factor / units_inp%energy%factor, &
         l_vector_field(1:3) * units_inp%length%factor
    
  end subroutine td_write_laser

  subroutine td_write_multipole(iunit, iter, t, dipole, multipole, header)
    integer, intent(in) :: iunit, iter
    real(r8), intent(IN) :: t, dipole(sys%st%nspin), multipole((td%lmax+1)**2, sys%st%nspin)
    logical, intent(in) :: header

    integer :: is, l, m, add_lm
    character(len=50) :: aux
    
    if(header) then
      ! first line
      write(iunit, '(a10,i2,a8,i2)') '# nspin = ', sys%st%nspin, ' lmax = ', td%lmax

      ! second line -> columns name
      write(iunit, '(a8,a20)', advance='no') '# Iter  ', str_center('t', 20)
      do is = 1, sys%st%nspin
        write(aux, '(a,i1,a)') 'dipole(', is, ')'
        write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
      end do
      do is = 1, sys%st%nspin
        do l = 0, td%lmax
          do m = -l, l
            write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', l, ', m=', m, ' (', is,')'
            write(iunit, '(a20)', advance='no') str_center(trim(aux), 20)
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')

      ! third line -> units
      write(iunit, '(a8,a20)', advance='no') '#       ', &
           str_center('['//trim(units_out%time%abbrev)//']', 20)
      do is = 1, sys%st%nspin
        write(iunit, '(a20)', advance='no') &
             str_center('['//trim(units_out%length%abbrev)//']', 20)
      end do
      do is = 1, sys%st%nspin
        do l = 0, td%lmax
          do m = -l, l
            select case(l)
            case(0)
              write(iunit, '(a20)', advance='no') ' '
            case(1)
              write(iunit, '(a20)', advance='no') &
                   str_center('['//trim(units_out%length%abbrev)//']', 20)
            case default
              write(aux, '(a,a2,i1)') trim(units_out%length%abbrev), "**", l
              write(iunit, '(a20)', advance='no') str_center('['//trim(aux)//']', 20)
            end select
            add_lm = add_lm + 1
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')
      
    end if
    
    write(iunit, '(i8, es20.12)', advance='no') iter, t/units_out%time%factor
    do is = 1, sys%st%nspin
      write(iunit, '(es20.12)', advance='no') dipole(is)/units_out%length%factor
    end do
    
    do is = 1, sys%st%nspin
      add_lm = 1
      do l = 0, td%lmax
        do m = -l, l
          write(iunit, '(es20.12)', advance='no') multipole(add_lm, is)/units_out%length%factor**l
          add_lm = add_lm + 1
        end do
      end do
    end do
    write(iunit, '(1x)', advance='yes')
    
  end subroutine td_write_multipole

  subroutine td_calc_projection(p)
    complex(r4), intent(out) :: p(u_st%nst, sys%st%st_start:sys%st%st_end, sys%st%nik)

    integer :: uist, uik, ist, ik

    do ik = 1, sys%st%nik
      do ist = sys%st%st_start, sys%st%st_end
        do uist = 1, u_st%nst
          p(uist, ist, ik) = cmplx(sum(sys%st%zpsi(1:sys%m%np,:, ist, ik)* &
               u_st%R_FUNC(psi) (1:sys%m%np,:, uist, ik)), kind=r4)*sys%m%vol_pp
        end do
      end do
    end do
  end subroutine td_calc_projection

end subroutine td_run

#include "td_init.F90"
#include "td_rti.F90"

end module timedep
