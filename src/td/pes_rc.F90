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
    integer          :: npoints                 !< how many points we store the wf
    integer, pointer :: points(:)               !< which points to use (local index)
    integer, pointer :: points_global(:)        !< global index of the points
    CMPLX, pointer   :: wf(:,:,:,:,:)           !< wavefunctions at sample points
    integer, pointer :: rankmin(:)              !< partition of the mesh containing the points
    FLOAT, pointer   :: dq(:,:)                 !< part 1 of Volkov phase (type PHASES) 
    FLOAT, pointer   :: domega(:)               !< part 2 of Volkov phase (type PHASES)
    integer          :: recipe                  !< type of calculation (RAW/PHASES)
  end type PES_rc_t

  integer, parameter :: &
    M_RAW    = 1,       &
    M_PHASES = 2

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
    integer       :: buf
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
    !%Option phases 2
    !% Calculate the photoelectron spectrum by including the Volkov phase (approximately), see
    !% P. M. Dinh, P. Romaniello, P.-G. Reinhard, and E. Suraud, <i>Phys. Rev. A.</i> <b>87</b>, 032514 (2013).
    !%End
    call parse_variable('PES_rc_recipe', M_RAW, pesrc%recipe)
    if(.not.varinfo_valid_option('PES_rc_recipe', pesrc%recipe, is_flag = .true.)) then 
      call messages_input_error('PES_rc_recipe')
    end if
    call messages_print_var_option(stdout, "PES_rc_recipe", pesrc%recipe)

    pesrc%npoints = parse_block_n(blk)

    ! read points
    SAFE_ALLOCATE(pesrc%points   (1:pesrc%npoints))
    SAFE_ALLOCATE(pesrc%points_global(1:pesrc%npoints))
    SAFE_ALLOCATE(pesrc%rankmin   (1:pesrc%npoints))

    do ip = 1, pesrc%npoints
      call parse_block_float(blk, ip - 1, 0, xx(1), units_inp%length)
      call parse_block_float(blk, ip - 1, 1, xx(2), units_inp%length)
      call parse_block_float(blk, ip - 1, 2, xx(3), units_inp%length)

      pesrc%points(ip) = mesh_nearest_point(mesh, xx, dmin, rankmin)
      pesrc%rankmin(ip)= rankmin

      if(mesh%parallel_in_domains) then
#if defined(HAVE_MPI)
        if(mesh%mpi_grp%rank == rankmin) then
          buf = mesh%vp%local(mesh%vp%xlocal + pesrc%points(ip) - 1)
        else
          buf = 0
        end if
        call MPI_Allreduce(buf, pesrc%points_global(ip), 1, MPI_INTEGER, mpi_sum, mpi_comm_world, mpi_err)
#endif
      else
        pesrc%points_global(ip) = pesrc%points(ip)
      end if
    end do

    call parse_block_end(blk)

    SAFE_ALLOCATE(pesrc%wf(1:st%nst, 1:st%d%dim, 1:st%d%nik, 1:pesrc%npoints, 0:save_iter-1))
    pesrc%wf = M_z0

    if(pesrc%recipe == M_PHASES) then
      SAFE_ALLOCATE(pesrc%dq(1:pesrc%npoints, 0:save_iter-1))
      SAFE_ALLOCATE(pesrc%domega(0:save_iter-1))
      pesrc%dq = M_ZERO
      pesrc%domega = M_ZERO
    end if

    POP_SUB(PES_rc_init)
  end subroutine PES_rc_init


  ! ---------------------------------------------------------
  subroutine PES_rc_end(pesrc)
    type(PES_rc_t), intent(inout) :: pesrc

    PUSH_SUB(PES_rc_end)

    if(associated(pesrc%points)) then
      SAFE_DEALLOCATE_P(pesrc%points)
      SAFE_DEALLOCATE_P(pesrc%points_global)
      SAFE_DEALLOCATE_P(pesrc%wf)
      SAFE_DEALLOCATE_P(pesrc%rankmin)

      if(pesrc%recipe == M_PHASES) then
        SAFE_DEALLOCATE_P(pesrc%dq)
        SAFE_DEALLOCATE_P(pesrc%domega)
      end if
    end if

    POP_SUB(PES_rc_end)
  end subroutine PES_rc_end


  ! ---------------------------------------------------------
  subroutine PES_rc_calc(pesrc, st, mesh, ii, dt, iter, hm)
    type(PES_rc_t),      intent(inout) :: pesrc
    type(states_t),      intent(in)    :: st
    integer,             intent(in)    :: ii
    type(mesh_t),        intent(in)    :: mesh
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: iter
    type(hamiltonian_t), intent(in)    :: hm

    CMPLX, allocatable :: psi(:,:,:,:), wfact(:,:,:,:)
    integer            :: ip, isdim
    integer            :: dim, stst, stend, kptst, kptend
    logical            :: contains_ip
    CMPLX              :: cfac
    FLOAT, allocatable :: buf(:)

    PUSH_SUB(PES_rc_calc)

    stst   = st%st_start
    stend  = st%st_end
    kptst  = st%d%kpt%start
    kptend = st%d%kpt%end
    dim    = st%d%dim

    SAFE_ALLOCATE(psi(1:st%nst, dim, 1:1, 1:st%d%nik))

    SAFE_ALLOCATE(wfact(1:st%nst, dim, 1:st%d%nik, 1:pesrc%npoints))
    wfact = M_z0

    if(pesrc%recipe == M_PHASES) then
      SAFE_ALLOCATE(buf(1:pesrc%npoints))
      buf = M_ZERO
    end if

    contains_ip = .true.

    do ip = 1, pesrc%npoints

#if defined(HAVE_MPI)
      if(mesh%parallel_in_domains) then
        if(mesh%mpi_grp%rank  ==  pesrc%rankmin(ip)) then !needed if mesh%parallel_in_domains is true
          contains_ip = .true.
        else
          contains_ip = .false.
        end if
      end if
#endif

      if(contains_ip) then
        call states_get_points(st, pesrc%points(ip), pesrc%points(ip), psi)
        wfact(stst:stend, dim, kptst:kptend, ip) = psi(stst:stend, dim, 1, kptst:kptend)

        if(pesrc%recipe == M_PHASES) then
          call PES_rc_calc_rcphase(pesrc, st, mesh, iter, dt, hm, ip, ii)
        end if
      end if
    end do

#if defined(HAVE_MPI)
    isdim = st%nst * dim * st%d%nik * pesrc%npoints
    call MPI_Allreduce(wfact(:,:,:,:), pesrc%wf(:,:,:,:,ii), isdim, MPI_CMPLX, mpi_sum, mpi_comm_world, mpi_err)

    if(pesrc%recipe == M_PHASES) then
      call MPI_Allreduce(pesrc%dq(:,ii), buf(:), pesrc%npoints, MPI_FLOAT, mpi_sum, mpi_comm_world, mpi_err)
      pesrc%dq(:,ii) = buf(:)
    end if
#else
    pesrc%wf(:,:,:,:,ii) = wfact(:,:,:,:)
#endif

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(wfact)

    POP_SUB(PES_rc_calc)
  end subroutine PES_rc_calc


  ! ---------------------------------------------------------
  subroutine PES_rc_output(pesrc, st, iter, save_iter, dt)
    type(PES_rc_t), intent(in) :: pesrc
    type(states_t), intent(in) :: st
    integer,        intent(in) :: iter, save_iter
    FLOAT,          intent(in) :: dt

    integer          :: ip, iunit, ii, jj, ik, ist, idim
    CMPLX            :: vfu
    character(len=2) :: filenr

    PUSH_SUB(PES_rc_output)

    do ip = 1, pesrc%npoints
      write(filenr, '(i2.2)') ip

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

      if(pesrc%recipe == M_PHASES) then
        iunit = io_open('td.general/'//'PES_rc.'//filenr//'.phases.out', action='write', position='append')
        do ii = 1, save_iter - mod(iter, save_iter)
          jj = iter - save_iter + ii + mod(save_iter - mod(iter, save_iter), save_iter)
          write(iunit, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
          write(iunit, '(1x,e18.10E3,1x,e18.10E3)', advance='no') &
            pesrc%domega(ii-1), pesrc%dq(ip, ii-1)
          write(iunit, '(1x)', advance='yes')
        end do
        call io_close(iunit)
      end if

    end do

    POP_SUB(PES_rc_output)
  end subroutine PES_rc_output

  ! ---------------------------------------------------------
  subroutine PES_rc_init_write(pesrc, mesh, st)
    type(PES_rc_t), intent(in) :: pesrc
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(in) :: st

    integer          :: ip,ik,ist,idim,iunit
    FLOAT            :: xx(MAX_DIM)
    character(len=2) :: filenr

    PUSH_SUB(PES_rc_init_write)

    xx = M_ZERO
    if(mpi_grp_is_root(mpi_world)) then
      do ip = 1, pesrc%npoints
        write(filenr, '(i2.2)') ip

        iunit = io_open('td.general/'//'PES_rc.'//filenr//'.wavefunctions.out', action='write')
        xx = mesh_x_global(mesh, pesrc%points_global(ip))
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

        if(pesrc%recipe == M_PHASES) then
          iunit = io_open('td.general/'//'PES_rc.'//filenr//'.phases.out', action='write')
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
   
    if(pesrc%recipe == M_PHASES) then
      iunit = restart_open(restart, 'rcphases')
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

    if(pesrc%recipe == M_PHASES) then
      iunit = restart_open(restart, 'rcphases')
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
  subroutine PES_rc_calc_rcphase(pesrc, st, mesh, iter, dt, hm, ip, ii)
    type(PES_rc_t),      intent(inout) :: pesrc
    type(states_t),      intent(in)    :: st
    type(mesh_t),        intent(in)    :: mesh
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: iter
    integer,             intent(in)    :: ip
    FLOAT,               intent(in)    :: dt
    integer,             intent(in)    :: ii
    
    integer :: dim, il
    FLOAT   :: vp(1:MAX_DIM)
    FLOAT   :: xx(MAX_DIM), er(MAX_DIM)

    PUSH_SUB(PES_rc_calc_rcphase)

    dim = mesh%sb%dim

    xx(1:dim) = mesh%x(pesrc%points(ip), 1:dim)
    er(1:dim) = xx(1:dim)/sqrt(dot_product(xx(1:dim), xx(1:dim)))

    vp = M_ZERO
    do il = 1, hm%ep%no_lasers
      call laser_field(hm%ep%lasers(il), vp(1:dim), iter*dt)
    end do
    vp(1:dim) = -vp(1:dim)

    pesrc%domega(ii) = pesrc%domega(ii) + dot_product(vp(1:dim), vp(1:dim)) / (M_TWO * P_C**M_TWO) * dt
    pesrc%dq(ip, ii) = pesrc%dq(ip, ii) + dot_product(er(1:dim), vp(1:dim))/P_C * dt

    POP_SUB(PES_rc_calc_rcphase)
  end subroutine PES_rc_calc_rcphase

end module pes_rc_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
