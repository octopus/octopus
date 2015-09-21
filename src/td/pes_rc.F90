!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$

#include "global.h"

module pes_rc_m
  use global_m
  use io_m
  use mesh_m
  use messages_m
  use mpi_m
  use comm_m
  use parser_m
  use profiling_m
  use restart_m
  use simul_box_m
  use states_m
  use unit_m
  use unit_system_m
  use hamiltonian_m
  use lasers_m
  use varinfo_m

  private

  public ::                             &
    pes_rc_t,                           &
    pes_rc_init,                        &
    pes_rc_init_write,                  &
    pes_rc_output,                      &
    pes_rc_calc,                        &
    pes_rc_end,                         &
    pes_rc_dump,                        &
    pes_rc_load

  type PES_rc_t
    integer          :: npoints                   !< how many points we store the wf
    integer, pointer :: points(:)                 !< which points to use (local index)
    integer, pointer :: points_global(:)          !< global index of the points
    FLOAT, pointer   :: coords(:,:)               !< coordinates of the sample points
    CMPLX, pointer   :: wf(:,:,:,:,:)   => NULL() !< wavefunctions at sample points
    integer, pointer :: rankmin(:)                !< partition of the mesh containing the points
    FLOAT, pointer   :: dq(:,:)         => NULL() !< part 1 of Volkov phase (recipe phase) 
    FLOAT, pointer   :: domega(:)       => NULL() !< part 2 of Volkov phase (recipe phase)
    integer          :: recipe                    !< type of calculation (RAW/PHASE)
    CMPLX, pointer   :: wfft(:,:,:,:,:) => NULL() !< Fourier transform of wavefunction
    FLOAT            :: omegamax                  !< maximum frequency of the spectrum
    FLOAT            :: delomega                  !< frequency spacing of the spectrum
    integer          :: nomega                    !< number of frequencies of the spectrum
    logical          :: onfly                     !< spectrum is calculated on-the-fly when true
    integer          :: save_iter                 !< output interval and size of arrays 
  end type PES_rc_t

  integer, parameter :: &
    M_RAW   = 1,        &
    M_PHASE = 2

contains

  ! ---------------------------------------------------------
  subroutine PES_rc_init(pesrc, mesh, st, save_iter)
    type(PES_rc_t), intent(out) :: pesrc
    type(mesh_t),   intent(in)  :: mesh
    type(states_t), intent(in)  :: st
    integer,        intent(in)  :: save_iter

    type(block_t) :: blk
    integer       :: ip
    FLOAT         :: xx(MAX_DIM)
    FLOAT         :: dmin
    integer       :: rankmin

    PUSH_SUB(PES_rc_init)

    message(1) = 'Info: Calculating PES using rc technique.'
    call messages_info(1)

    !%Variable PhotoElectronSpectrumPoints
    !%Type block
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% List of points at which to calculate the photoelectron spectrum by Suraud method.
    !% Required when <tt>PhotoElectronSpectrum = pes_rc</tt>.
    !% The exact syntax is:
    !%
    !% <tt>%PhotoElectronSpectrumPoints
    !% <br>&nbsp;&nbsp;x1 | y1 | z1
    !% <br>%
    !% </tt>
    !%End
    call messages_obsolete_variable('PES_rc_points', 'PhotoElectronSpectrumPoints')
    if (parse_block('PhotoElectronSpectrumPoints', blk) < 0) then
      message(1) = 'The PhotoElectronSpectrumPoints block is required when PhotoElectronSpectrum = pes_rc'
      call messages_fatal(1)
    end if

    !%Variable PES_rc_recipe
    !%Type integer
    !%Default raw
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The type for calculating the photoelectron spectrum in the sample point method.
    !%Option raw 1
    !% Calculate the photoelectron spectrum according to A. Pohl, P.-G. Reinhard, and 
    !% E. Suraud, <i>Phys. Rev. Lett.</i> <b>84</b>, 5090 (2000).
    !%Option phase 2
    !% Calculate the photoelectron spectrum by including the Volkov phase (approximately), see
    !% P. M. Dinh, P. Romaniello, P.-G. Reinhard, and E. Suraud, <i>Phys. Rev. A.</i> <b>87</b>, 032514 (2013).
    !%End
    call parse_variable('PES_rc_recipe', M_RAW, pesrc%recipe)
    if(.not.varinfo_valid_option('PES_rc_recipe', pesrc%recipe, is_flag = .true.)) &
      call messages_input_error('PES_rc_recipe')
    call messages_print_var_option(stdout, "PES_rc_recipe", pesrc%recipe)

    !%Variable PES_rc_OmegaMax
    !%Type float
    !%Default 0.0
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% If PES_rc_OmegaMax > 0, the photoelectron spectrum is directly calculated during 
    !% time-propagation, evaluated by the PES_rc method. PES_rc_OmegaMax is then the maximum frequency 
    !% (approximate kinetic energy) and PES_rc_DelOmega the spacing in frequency domain of the spectrum.
    !%End
    call parse_variable('PES_rc_OmegaMax', units_to_atomic(units_inp%energy, M_ZERO), pesrc%omegamax)
    pesrc%onfly = .false.
    if(pesrc%omegamax > M_ZERO) then
      pesrc%onfly = .true.
      message(1) = 'Info: Calculating PES during time propagation.'
      call messages_info(1)
      call messages_print_var_value(stdout, "PES_rc_OmegaMax", pesrc%omegamax)
    end if
 
    !%Variable PES_rc_DelOmega
    !%Type float
    !%Default 0.005
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !% The spacing in frequency domain for the photoelectron spectrum (if PES_rc_OmegaMax > 0).
    !%End
    call parse_variable('PES_rc_DelOmega', units_to_atomic(units_inp%energy, CNST(0.005)), pesrc%delomega)
    if(pesrc%onfly) then
      if(pesrc%delomega <= M_ZERO) call messages_input_error('PES_rc_DelOmega')
      call messages_print_var_value(stdout, "PES_rc_DelOmega", pesrc%delomega)
    end if

    pesrc%npoints = parse_block_n(blk)

    call messages_print_var_value(stdout, "Number of PhotoElectronSpectrumPoints", pesrc%npoints)

    SAFE_ALLOCATE(pesrc%coords(1:mesh%sb%dim, 1:pesrc%npoints))

    SAFE_ALLOCATE(pesrc%points(1:pesrc%npoints))
    SAFE_ALLOCATE(pesrc%points_global(1:pesrc%npoints))
    SAFE_ALLOCATE(pesrc%rankmin(1:pesrc%npoints))
    pesrc%points_global = 0

   ! read points from input file
    do ip = 1, pesrc%npoints
      call parse_block_float(blk, ip - 1, 0, xx(1), units_inp%length)
      call parse_block_float(blk, ip - 1, 1, xx(2), units_inp%length)
      call parse_block_float(blk, ip - 1, 2, xx(3), units_inp%length)
      pesrc%coords(1:mesh%sb%dim, ip) = xx(1:mesh%sb%dim)
    end do
    call parse_block_end(blk)

    message(1) = 'Info: Calculating nearest points.'
    call messages_info(1)

    do ip = 1, pesrc%npoints
      ! nearest point
      pesrc%points(ip) = mesh_nearest_point(mesh, pesrc%coords(1:mesh%sb%dim, ip), dmin, rankmin)
      pesrc%rankmin(ip)= rankmin

      if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
        if(mesh%mpi_grp%rank == rankmin) &
          pesrc%points_global(ip) = mesh%vp%local(mesh%vp%xlocal + pesrc%points(ip) - 1)
        call comm_allreduce(mesh%mpi_grp%comm, pesrc%points_global)
#endif
      else
        pesrc%points_global(ip) = pesrc%points(ip)
      end if
    end do

    SAFE_ALLOCATE(pesrc%wf(1:st%nst, 1:st%d%dim, 1:st%d%nik, 1:pesrc%npoints, 0:save_iter-1))

    if(pesrc%recipe == M_PHASE) then
      SAFE_ALLOCATE(pesrc%dq(1:pesrc%npoints, 0:save_iter-1))
      SAFE_ALLOCATE(pesrc%domega(0:save_iter-1))
      pesrc%dq = M_ZERO
      pesrc%domega = M_ZERO
    end if

    if(pesrc%onfly) then
      pesrc%nomega = nint(pesrc%omegamax/pesrc%delomega)
      SAFE_ALLOCATE(pesrc%wfft(1:st%nst, 1:st%d%dim, 1:st%d%nik, 1:pesrc%npoints, 1:pesrc%nomega))
      pesrc%wfft = M_z0
    end if

    pesrc%save_iter = save_iter

    POP_SUB(PES_rc_init)
  end subroutine PES_rc_init


  ! ---------------------------------------------------------
  subroutine PES_rc_end(pesrc)
    type(PES_rc_t), intent(inout) :: pesrc

    PUSH_SUB(PES_rc_end)

    SAFE_DEALLOCATE_P(pesrc%points)
    SAFE_DEALLOCATE_P(pesrc%points_global)
    SAFE_DEALLOCATE_P(pesrc%wf)
    SAFE_DEALLOCATE_P(pesrc%rankmin)
    SAFE_DEALLOCATE_P(pesrc%coords)

    SAFE_DEALLOCATE_P(pesrc%wfft)

    SAFE_DEALLOCATE_P(pesrc%dq)
    SAFE_DEALLOCATE_P(pesrc%domega)

    POP_SUB(PES_rc_end)
  end subroutine PES_rc_end


  ! ---------------------------------------------------------
  subroutine PES_rc_calc(pesrc, st, mesh, dt, iter, hm)
    type(PES_rc_t),      intent(inout) :: pesrc
    type(states_t),      intent(in)    :: st
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter
    type(hamiltonian_t), intent(in)    :: hm

    CMPLX, allocatable :: psi(:,:,:,:), wfftact(:,:,:,:,:)
    integer            :: ip, ii
    integer            :: dim, stst, stend, kptst, kptend
    logical            :: contains_ip
    CMPLX              :: rawfac
    CMPLX, allocatable :: phasefac(:)
#if defined(HAVE_MPI)
    integer            :: isdim
#endif
    FLOAT              :: omega
    integer            :: iom

    PUSH_SUB(PES_rc_calc)

    ii = mod(iter-1, pesrc%save_iter)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    dim    = st%d%dim

    SAFE_ALLOCATE(psi(1:st%nst, dim, 1:1, 1:st%d%nik))

    if(pesrc%onfly) then
      SAFE_ALLOCATE(wfftact(1:st%nst, dim, 1:st%d%nik, 1:pesrc%npoints, 1:pesrc%nomega))
      wfftact = M_z0
    endif

    if(pesrc%recipe == M_PHASE) then
      SAFE_ALLOCATE(phasefac(1:pesrc%npoints))
    end if

    ! needed for allreduce, otherwise it will take values from previous cycle
    pesrc%wf(:,:,:,:,ii) = M_z0

    contains_ip = .true.

    do ip = 1, pesrc%npoints

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        if(mesh%mpi_grp%rank  ==  pesrc%rankmin(ip)) then ! needed if mesh%parallel_in_domains is true
          contains_ip = .true.
        else
          contains_ip = .false.
        end if
      end if
#endif

      if(contains_ip) then
        call states_get_points(st, pesrc%points(ip), pesrc%points(ip), psi)
        pesrc%wf(stst:stend, dim, kptst:kptend, ip, ii) = psi(stst:stend, dim, 1, kptst:kptend)
      end if
    end do
#if defined(HAVE_MPI)
    call comm_allreduce(mpi_world%comm, pesrc%wf(:,:,:,:,ii))
#endif

    if(pesrc%recipe == M_PHASE) then
      call pes_rc_calc_rcphase(pesrc, mesh, iter, dt, hm, ii)
    end if

    if(pesrc%onfly) then
      do iom = 1, pesrc%nomega
        omega = iom*pesrc%delomega
        rawfac = exp(M_zI * omega * iter * dt) * dt

        select case(pesrc%recipe)
        case(M_RAW)
          wfftact(stst:stend, 1:dim, kptst:kptend, :, iom) = pesrc%wf(stst:stend, 1:dim, kptst:kptend, :, ii) * rawfac
        case(M_PHASE)
          phasefac(:) = exp(M_zI * (pesrc%domega(ii) - sqrt(M_TWO * omega) * pesrc%dq(:, ii)))

          do ik = kptst, kptend
            do ist = stst, stend
              do idim = 1, dim
                wfftact(ist, idim, ik, :, iom) = pesrc%wf(ist, idim, ik, :, ii) * phasefac(:) * rawfac
              end do
            end do
          end do
        end select
#if defined(HAVE_MPI)
        call comm_allreduce(mpi_world%comm, wfftact(:,:,:,:,iom))
#endif
      end do

      pesrc%wfft = pesrc%wfft + wfftact
    end if

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(phasefac)

    SAFE_DEALLOCATE_A(wfftact)

    POP_SUB(PES_rc_calc)
  end subroutine PES_rc_calc


  ! ---------------------------------------------------------
  subroutine PES_rc_output(pesrc, st, iter, dt)
    type(PES_rc_t), intent(in) :: pesrc
    type(states_t), intent(in) :: st
    integer,        intent(in) :: iter
    FLOAT,          intent(in) :: dt

    integer            :: ip, iunit, ii, jj, ik, ist, idim, iom
    CMPLX              :: vfu
    character(len=4)   :: filenr
    FLOAT              :: wfu, omega
    FLOAT, allocatable :: wffttot(:)
    integer            :: save_iter

    PUSH_SUB(PES_rc_output)

    save_iter = pesrc%save_iter

    if(pesrc%onfly) then
      SAFE_ALLOCATE(wffttot(1:pesrc%nomega))
    end if

    do ip = 1, pesrc%npoints
      write(filenr, '(i4.4)') ip

      iunit = io_open('td.general/'//'PES_rc.'//filenr//'.wavefunctions.out', action='write', position='append')
      do ii = 1, save_iter - mod(iter, save_iter)
        jj = iter - save_iter + ii + mod(save_iter - mod(iter, save_iter), save_iter)
        write(iunit, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            do idim = 1, st%d%dim
              vfu = units_from_atomic(sqrt(units_out%length**(-3)), pesrc%wf(ist, idim, ik, ip, ii-1))
              write(iunit, '(1x,e18.10E3,1x,e18.10E3)', advance='no') &
                real(vfu),  aimag(vfu) 
            end do
          end do
        end do
        write(iunit, '(1x)', advance='yes')
      end do
      call io_close(iunit)

      if(pesrc%recipe == M_PHASE) then
        iunit = io_open('td.general/'//'PES_rc.'//filenr//'.phase.out', action='write', position='append')
        do ii = 1, save_iter - mod(iter, save_iter)
          jj = iter - save_iter + ii + mod(save_iter - mod(iter, save_iter), save_iter)
          write(iunit, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
          write(iunit, '(1x,e18.10E3,1x,e18.10E3)', advance='no') &
            pesrc%domega(ii-1), pesrc%dq(ip, ii-1)
          write(iunit, '(1x)', advance='yes')
        end do
        call io_close(iunit)
      end if

      if(pesrc%onfly) then
        wffttot = M_ZERO
        do iom = 1, pesrc%nomega
          do ik = 1, st%d%nik
            do ist = 1, st%nst
              do idim = 1, st%d%dim
                wffttot(iom) = real(pesrc%wfft(ist, idim, ik, ip, iom))**2 + aimag(pesrc%wfft(ist, idim, ik, ip, iom))**2 &
                  + wffttot(iom)
              end do
            end do
          end do
        end do

        iunit = io_open('td.general/'//'PES_rc.'//filenr//'.spectrum.out', action='write', position='rewind')
        write(iunit, '(a44)') '# frequency, total spectrum, orbital spectra'
        do iom = 1, pesrc%nomega 
          omega = iom*pesrc%delomega
          write(iunit, '(e17.10, 1x, e17.10)', advance='no') omega, wffttot(iom)
          do ik = 1, st%d%nik
            do ist = 1, st%nst 
              do idim = 1, st%d%dim
                wfu = real(pesrc%wfft(ist, idim, ik, ip, iom))**2 + aimag(pesrc%wfft(ist, idim, ik, ip, iom))**2
                write(iunit,'(1x,e18.10e3)', advance='no') wfu
              end do
            end do
          end do
        write(iunit,'(1x)', advance='yes')
        end do
        call io_close(iunit)
      end if

    end do

    SAFE_DEALLOCATE_A(wffttot)

    POP_SUB(PES_rc_output)
  end subroutine PES_rc_output

  ! ---------------------------------------------------------
  subroutine PES_rc_init_write(pesrc, mesh, st)
    type(PES_rc_t), intent(in) :: pesrc
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(in) :: st

    integer          :: ip,ik,ist,idim,iunit
    FLOAT            :: xx(MAX_DIM)
    character(len=4) :: filenr

    PUSH_SUB(PES_rc_init_write)

    xx = M_ZERO
    if(mpi_grp_is_root(mpi_world)) then
      do ip = 1, pesrc%npoints
        write(filenr, '(i4.4)') ip

        iunit = io_open('td.general/'//'PES_rc.'//filenr//'.wavefunctions.out', action='write')
        xx(1:mesh%sb%dim) = pesrc%coords(1:mesh%sb%dim, ip) 
        write(iunit,'(a1)') '#'
        write(iunit, '(a7,f17.6,a1,f17.6,a1,f17.6,5a)') &
          '# R = (',xx(1),' ,',xx(2),' ,',xx(3), &
          ' )  [', trim(units_abbrev(units_inp%length)), ']'

        write(iunit,'(a1)') '#'  
        write(iunit, '(a3,14x)', advance='no') '# t' 
        do ik = 1, st%d%nik
          do ist = 1, st%nst
            do idim = 1, st%d%dim
              write(iunit, '(3x,a8,i3,a7,i3,a8,i3,3x)', advance='no') &
                "ik = ", ik, " ist = ", ist, " idim = ", idim
            end do
          end do
        end do
        write(iunit, '(1x)', advance='yes')

        call io_close(iunit)

        if(pesrc%recipe == M_PHASE) then
          iunit = io_open('td.general/'//'PES_rc.'//filenr//'.phase.out', action='write')
          write(iunit,'(a24)') '# time, dq(t), dOmega(t)'
          call io_close(iunit)
        end if
      end do
    end if

    POP_SUB(PES_rc_init_write)
  end subroutine PES_rc_init_write

  ! ---------------------------------------------------------
  subroutine pes_rc_dump(restart, pesrc, st, ierr)
    type(restart_t), intent(in)  :: restart    
    type(pes_rc_t),  intent(in)  :: pesrc
    type(states_t),  intent(in)  :: st
    integer,         intent(out) :: ierr
    
    integer :: err, iunit
    
    PUSH_SUB(pes_rc_dump)

    ierr = 0
    
    if (restart_skip(restart)) then
      POP_SUB(pes_rc_dump)
      return
    end if
    
    if (in_debug_mode) then
      message(1) = "Debug: Writing pes_rc restart."
      call messages_info(1)
    end if

    if(pesrc%onfly) then 
      call zrestart_write_binary(restart, 'pesrc', pesrc%npoints*st%d%dim*st%nst*st%d%nik*pesrc%nomega, &
        pesrc%wfft, err) 
    end if

    if(pesrc%recipe == M_PHASE) then
      iunit = restart_open(restart, 'rcphase')
      write(iunit, *)  pesrc%domega(:), pesrc%dq(:,:)
      call restart_close(restart, iunit)
    end if

    if (err /= 0) ierr = ierr + 1
    
    if (in_debug_mode) then
      message(1) = "Debug: Writing pes_rc restart done."
      call messages_info(1)
    end if
    
    POP_SUB(pes_rc_dump)
  end subroutine pes_rc_dump

  ! ---------------------------------------------------------
  subroutine pes_rc_load(restart, pesrc, st, ierr)
    type(restart_t), intent(in)    :: restart    
    type(pes_rc_t),  intent(inout) :: pesrc
    type(states_t),  intent(inout) :: st
    integer,         intent(out)   :: ierr
    
    integer :: err, iunit
    
    PUSH_SUB(pes_rc_load)
    
    ierr = 0
    
    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(pes_rc_load)
      return
    end if
    
    if (in_debug_mode) then
      message(1) = "Debug: Reading pes_rc restart."
      call messages_info(1)
    end if

    if(pesrc%onfly) then
      call zrestart_read_binary(restart, 'pesrc', pesrc%npoints*st%d%dim*st%nst*st%d%nik*pesrc%nomega, &
        pesrc%wfft, err)
    end if

    if(pesrc%recipe == M_PHASE) then
      iunit = restart_open(restart, 'rcphase')
      read(iunit, *)  pesrc%domega(:), pesrc%dq(:,:)
      call restart_close(restart, iunit)
    end if

    if (err /= 0) ierr = ierr + 1
    
    if(in_debug_mode) then
      message(1) = "Debug: Reading pes_rc restart done."
      call messages_info(1)
    end if
    
    POP_SUB(pes_rc_load)
  end subroutine pes_rc_load 

  ! ---------------------------------------------------------
  subroutine PES_rc_calc_rcphase(pesrc, mesh, iter, dt, hm, ii)
    type(PES_rc_t),      intent(inout) :: pesrc
    type(mesh_t),        intent(in)    :: mesh
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: ii
    
    integer :: dim, il, ip
    FLOAT   :: vp(1:MAX_DIM)
    FLOAT   :: xx(MAX_DIM), er(MAX_DIM)

    PUSH_SUB(PES_rc_calc_rcphase)

    dim = mesh%sb%dim

    vp = M_ZERO
    do il = 1, hm%ep%no_lasers
      call laser_field(hm%ep%lasers(il), vp(1:dim), iter*dt)
    end do
    vp(1:dim) = -vp(1:dim)

    do ip = 1, pesrc%npoints
      xx(1:dim) = pesrc%coords(1:dim, ip)
      er(1:dim) = xx(1:dim)/sqrt(dot_product(xx(1:dim), xx(1:dim)))
      pesrc%dq(ip, ii) = pesrc%dq(ip, ii) + dot_product(er(1:dim), vp(1:dim))/P_C * dt
    end do

    pesrc%domega(ii) = pesrc%domega(ii) + dot_product(vp(1:dim), vp(1:dim)) / (M_TWO * P_C**M_TWO) * dt

    POP_SUB(PES_rc_calc_rcphase)
  end subroutine PES_rc_calc_rcphase

end module pes_rc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
