#include "config_F90.h"

module spectrum
use global
use io
use units

implicit none

type spec_type
  real(r8) :: start_time  ! start time for the transform
  real(r8) :: end_time    ! when to stop the transform
  real(r8) :: energy_step ! step in energy mesh
  real(r8) :: min_energy  ! maximum of energy mesh
  real(r8) :: max_energy  ! maximum of energy mesh
end type spec_type

! For the strength function
type spec_sf
  ! input
  integer  :: transform        ! The Fourier transform to perform (sin, cos)
  integer  :: damp             ! Damp type (none, exp or pol)
  real(r8) :: damp_factor      ! factor used in damping
  real(r8) :: delta_strength   ! strength of the delta perturbation

  ! output
  integer :: no_e, nspin ! dimensions of sp
  real(r8), pointer :: sp(:,:) ! do not forget to deallocate this
  real(r8) :: ewsum  ! electronic sum rule
  real(r8) :: alpha  ! Polariz. (sum rule)
  real(r8) :: alpha2 ! Polariz. (F.T.)
end type spec_sf

type spec_sh
  ! input
  character :: pol

  ! output
  integer :: no_e ! dimensions of sp
  real(r8), pointer :: sp(:) ! do not forget to deallocate this
end type spec_sh

contains

subroutine spectrum_strength_function(sysname, out_file, s, sf, print_info)
  character(len=*), intent(in) :: sysname, out_file
  type(spec_type), intent(inout) :: s
  type(spec_sf), intent(inout) :: sf
  logical, intent(in) :: print_info

  integer :: iunit, i, is, ie, &
      ntiter, j, jj, k, isp, time_steps
  real(r8) :: dump, dt, x
  real(r8), allocatable :: dipole(:,:)

  call spectrum_mult_info(sysname, iunit, sf%nspin, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)

  ! load dipole from file
  allocate(dipole(0:time_steps, sf%nspin))
  do i = 0, time_steps
    read(iunit, *) j, dump, dipole(i,:)
    dipole(i,:) = dipole(i,:) * units_out%length%factor
  end do
  call io_close(iunit)

  ! subtract static dipole
  do i = 1, sf%nspin
     dipole(:, i) = dipole(:, i) - dipole(0, i)
  end do

  sf%no_e = s%max_energy / s%energy_step
  allocate(sf%sp(0:sf%no_e, sf%nspin))
  sf%sp = 0._r8
  sf%alpha = 0._r8; sf%alpha2 = 0._r8; sf%ewsum = 0._r8

  do k = 0, sf%no_e
    do j = is, ie
      
      jj = j - is
      select case(sf%damp)
      case(0)
        dump = 1._r8
      case(1)
        dump = exp(-jj*dt*sf%damp_factor)
      case(2)
        dump = 1.0 - 3.0*(real(jj)/ntiter)**2                          &
            + 2.0*(real(jj)/ntiter)**3
      end select
      
      select case(sf%transform)
      case(1)
        x = sin(k*s%energy_step*jj*dt)
      case(2)
        x = cos(k*s%energy_step*jj*dt)
      end select

      do isp = 1, sf%nspin
        sf%sp(k, isp) = sf%sp(k, isp) + x*dump*dipole(j, isp)

        ! polarizability sum rule
        if(k == 0) then
          sf%alpha2 = sf%alpha2 + dump*dipole(j, isp)
        end if
      end do

    end do
    sf%sp(k, :) = sf%sp(k, :)*dt

    ! calculate strength function
    sf%sp(k, :) = (sf%sp(k, :)*s%energy_step*k*2._r8)/(M_pi*sf%delta_strength)

    do isp = 1, sf%nspin
      if(k.ne.0) then
        sf%alpha = sf%alpha + (sf%sp(k, isp)/(k*s%energy_step)**2)*s%energy_step
      endif
      sf%ewsum = sf%ewsum + sf%sp(k, isp)*s%energy_step
    end do
  end do

  sf%alpha2 = sf%alpha2 * dt / sf%delta_strength
  deallocate(dipole)

  ! output
  if(trim(out_file) .ne. '-') then
    call io_assign(iunit)
    open(iunit, file=out_file, status='unknown')
    ! should output units, etc...
    do i = 0, sf%no_e
      write(iunit,*) i*s%energy_step / units_out%energy%factor, &
           sf%sp(i, :) * units_out%energy%factor
    end do
    call io_close(iunit)
  end if

  ! print some info
  if(print_info) then
    write(message(1),'(a,i8)')    'Number of spin       = ', sf%nspin
    write(message(2),'(a,i8)')    'Number of time steps = ', ntiter
    write(message(3),'(a,i4)')    'SpecTransformMode    = ', sf%transform
    write(message(4),'(a,i4)')    'SpecDampMode         = ', sf%damp
    write(message(5),'(a,f10.4)') 'SpecDampFactor       = ', sf%damp_factor * units_out%time%factor
    write(message(6),'(a,f10.4)') 'SpecStartTime        = ', s%start_time   / units_out%time%factor
    write(message(7),'(a,f10.4)') 'SpecEndTime          = ', s%end_time     / units_out%time%factor
    write(message(8),'(a,f10.4)') 'SpecMaxEnergy        = ', s%max_energy   / units_inp%energy%factor
    write(message(9),'(a,f10.4)') 'SpecEnergyStep       = ', s%energy_step  / units_inp%energy%factor
    call write_info(9)

    message(1) = ""
    write(message(2),'(a,f16.6)') 'Electronic sum rule  = ', sf%ewsum
    write(message(3),'(a,f16.6)') 'Polariz. (sum rule)  = ', sf%alpha  / units_inp%length%factor**3
    write(message(4),'(a,f16.6)') 'Polariz. (F.T.)      = ', sf%alpha2 / units_inp%length%factor**3
    call write_info(4)
  end if

  return
end subroutine spectrum_strength_function

subroutine spectrum_hs_from_mult(sysname, out_file, s, sh, print_info)
  character(len=*), intent(in) :: sysname, out_file
  type(spec_type), intent(inout) :: s
  type(spec_sh), intent(inout) :: sh
  logical, intent(in) :: print_info

  integer :: i, j, iunit, nspin, time_steps, is, ie, ntiter
  real(r8) :: dt, dump
  real(r8), allocatable :: d(:,:)
  complex(r8) :: c
  complex(r8), allocatable :: dipole(:), ddipole(:)

  call spectrum_mult_info(sysname, iunit, nspin, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)
  
  ! load dipole from file
  allocate(dipole(0:time_steps))
  allocate(d(3, nspin))
  do i = 1, time_steps
    read(iunit, *) j, dump, dump, dump, d
    select case(sh%pol)
    case('x')
      dipole(i) = -sum(d(3, :))
    case('y')
      dipole(i) = -sum(d(1, :))
    case('z')
      dipole(i) = sum(d(2, :))
    case('+')
      dipole(i) = -sum(d(3, :) + M_zI*d(1, :)) / sqrt(2._r8)
    case('-')
      dipole(i) = -sum(d(3, :) - M_zI*d(1, :)) / sqrt(2._r8)
    end select
    dipole(i) = dipole(i) * units_out%length%factor * sqrt(4*M_PI/3)
  end do
  deallocate(d)

  ! we now calculate the first time derivative
  allocate(ddipole(0:time_steps))
  ddipole(0) = (dipole(1) - dipole(0))/dt
  do i = 1, time_steps - 1
    ddipole(i) = (dipole(i + 1) - dipole(i - 1))/(2._r8*dt)
  end do
  ddipole(time_steps) = (dipole(time_steps) - dipole(time_steps - 1))/dt

  ! and the second time derivative
  dipole(0) = (ddipole(1) - ddipole(0))/dt
  do i = 1, time_steps - 1
    dipole(i) = (ddipole(i + 1) - ddipole(i - 1))/(2._r8*dt)
  end do
  dipole(time_steps) = (ddipole(time_steps) - ddipole(time_steps - 1))/dt
  deallocate(ddipole)

  ! now we Fourier transform
  sh%no_e = s%max_energy / s%energy_step
  allocate(sh%sp(0:sh%no_e))
  sh%sp = 0._r8
  
  do i = 0, sh%no_e
    c = M_z0
    do j = is, ie
      c = c + exp(M_zI*i*s%energy_step*j*dt)*dipole(j)
    end do
    sh%sp(i) = abs(c)**2
  end do
  deallocate(dipole)

  ! output
  if(trim(out_file) .ne. '-') then
    call io_assign(iunit)
    open(iunit, file=trim(out_file)//"."//trim(sh%pol), status='unknown')
    ! should output units, etc...
    do i = 0, sh%no_e
      write(iunit,*) i*s%energy_step / units_out%energy%factor, &
           sh%sp(i) * units_out%energy%factor
    end do
    call io_close(iunit)
  end if

end subroutine spectrum_hs_from_mult

subroutine spectrum_hs_from_acc(sysname, out_file, s, sh, print_info)
  character(len=*), intent(in) :: sysname, out_file
  type(spec_type), intent(inout) :: s
  type(spec_sh), intent(inout) :: sh
  logical, intent(in) :: print_info

  integer :: i, j, iunit, time_steps, is, ie, ntiter
  real(r8) :: dt, dummy, a(3)
  complex(r8), allocatable :: acc(:)
  complex(r8) :: c
  
  call spectrum_acc_info(sysname, iunit, time_steps, dt)
  call spectrum_fix_time_limits(time_steps, dt, s%start_time, s%end_time, is, ie, ntiter)
  
  ! load dipole from file
  allocate(acc(0:time_steps))
  acc = 0._r8
  do i = 1, time_steps
    read(iunit, *) j, dummy, a
    select case(sh%pol)
    case('x')
      acc(i) = a(1)
    case('y')
      acc(i) = a(2)
    case('z')
      acc(i) = a(3)
    case('+')
      acc(i) = (a(1) + M_zI*a(2)) / sqrt(2._r8)
    case('-')
      acc(i) = (a(1) - M_zI*a(2)) / sqrt(2._r8)
    end select
    acc(i) = acc(i) * units_out%acceleration%factor
  end do

  ! now we Fourier transform
  sh%no_e = s%max_energy / s%energy_step
  allocate(sh%sp(0:sh%no_e))
  sh%sp = 0._r8
  
  do i = 0, sh%no_e
    c = M_z0
    do j = is, ie
      c = c + exp(M_zI*i*s%energy_step*j*dt)*acc(j)
    end do
    sh%sp(i) = abs(c)**2
  end do
  deallocate(acc)

  ! output
  if(trim(out_file) .ne. '-') then
    call io_assign(iunit)
    open(iunit, file=trim(out_file)//"."//trim(sh%pol), status='unknown')
    ! should output units, etc...
    do i = 0, sh%no_e
      write(iunit,*) i*s%energy_step / units_out%energy%factor, &
           sh%sp(i) * units_out%energy%factor
    end do
    call io_close(iunit)
  end if

end subroutine spectrum_hs_from_acc

subroutine spectrum_mult_info(sysname, iunit, nspin, time_steps, dt)
  character(len=*), intent(IN) :: sysname
  integer, intent(out) :: iunit, nspin, time_steps
  real(r8), intent(out) :: dt

  integer :: i, j
  real(r8) :: t1, t2, dummy

  ! open files
  call io_assign(iunit)
  open(iunit, file=trim(sysname)//'.mult', status='old', iostat=i)
  if(i.ne.0) then
    write(message(1),'(3a)') "Could not open file '", trim(sysname), ".mult'"
    call write_fatal(1)
  endif
  
  ! read in dipole
  read(iunit, '(10x,i2)') nspin
  read(iunit, *); read(iunit, *) ! skip header

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
100 continue
  dt = (t2 - t1) * units_out%time%factor ! units_out is OK
  time_steps = time_steps - 1
  
  if(time_steps < 3) then
    message(1) = "Empty multipole file?"
    call write_fatal(1)
  end if

  rewind(iunit)
  read(iunit, *); read(iunit, *); read(iunit, *) ! skip header

end subroutine spectrum_mult_info

subroutine spectrum_acc_info(sysname, iunit, time_steps, dt)
  character(len=*), intent(IN) :: sysname
  integer, intent(out) :: iunit, time_steps
  real(r8), intent(out) :: dt

  integer :: i, j
  real(r8) :: t1, t2, dummy

  ! open files
  call io_assign(iunit)
  open(iunit, file=trim(sysname)//'.acc', status='old', iostat=i)
  if(i.ne.0) then
    write(message(1),'(3a)') "Could not open file '", trim(sysname), ".acc'"
    call write_fatal(1)
  endif
  
  ! read in dipole
  read(iunit, *); read(iunit, *) ! skip header

  ! count number of time_steps
  time_steps = 0
  do
    read(iunit, *, end=100) j, dummy
    time_steps = time_steps + 1
    if(time_steps == 1) t1 = dummy
    if(time_steps == 2) t2 = dummy
  end do
100 continue
  dt = (t2 - t1) * units_out%time%factor ! units_out is OK
  time_steps = time_steps - 1
  
  if(time_steps < 3) then
    message(1) = "Empty multipole file?"
    call write_fatal(1)
  end if

  rewind(iunit)
  read(iunit, *); read(iunit, *) ! skip header

end subroutine spectrum_acc_info

subroutine spectrum_fix_time_limits(time_steps, dt, start_time, end_time, is, ie, ntiter)
  integer, intent(in) :: time_steps
  real(r8), intent(in) :: dt
  real(r8), intent(inout) :: start_time, end_time
  integer, intent(out) :: is, ie, ntiter

  real(r8) :: ts, te, dummy

  ts = 0._r8; te = time_steps*dt
  if(start_time < ts) start_time = ts
  if(start_time > te) start_time = te
  if(end_time   > te .or. end_time <= 0._r8) end_time   = te
  if(end_time   < ts) end_time   = ts

  if(end_time < start_time) then
    dummy = end_time ! swap
    end_time = start_time
    start_time = dummy
  end if
  is = int(start_time/dt)
  ie = int(end_time/dt)
  ntiter = ie - is + 1

end subroutine spectrum_fix_time_limits

end module spectrum
