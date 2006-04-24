! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!!
!! $Id$

#include "global.h"

module scf_m
  use global_m
  use mpi_m
  use messages_m
  use datasets_m
  use lib_oct_parser_m
  use units_m
  use geometry_m
  use simul_box_m
  use mesh_m
  use mesh_function_m
  use functions_m
  use states_m
  use output_m
  use restart_m
  use v_ks_m
  use hamiltonian_m
  use external_pot_m
  use xc_m
  use eigen_solver_m
  use mix_m
  use lcao_m
  use io_m
  use grid_m
  use profiling_m
  use varinfo_m

  implicit none

  private
  public ::             &
    scf_t,              &
    scf_init,           &
    scf_run,            &
    scf_end

  integer, parameter :: &
    MIXPOT  = 1,        &
    MIXDENS = 2

  type scf_t      ! some variables used for the scf cycle
    integer :: max_iter   ! maximum number of scf iterations

    FLOAT :: lmm_r

    ! several convergence criteria
    FLOAT :: conv_abs_dens, conv_rel_dens, conv_abs_ev, conv_rel_ev
    FLOAT :: abs_dens, rel_dens, abs_ev, rel_ev

    integer :: what2mix
    logical :: lcao_restricted

    type(mix_t) :: smix
    type(eigen_solver_t) :: eigens
  end type scf_t


contains

  ! ---------------------------------------------------------
  subroutine scf_init(gr, scf, st, h)
    type(grid_t),        intent(in)    :: gr
    type(scf_t),         intent(inout) :: scf
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(in)    :: h

    integer :: dim
    FLOAT :: rmin

    call push_sub('scf.scf_init')

    !%Variable MaximumIter
    !%Type integer
    !%Default 200
    !%Section SCF::Convergence
    !%Description
    !% Maximum number of SCF iterations. The code will stop even if convergence
    !% has not been achieved. <tt>0</tt> means unlimited.
    !%End
    call loct_parse_int  (check_inp('MaximumIter'), 200, scf%max_iter)

    !%Variable ConvAbsDens
    !%Type float
    !%Default 1e-5
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the density: 
    !% <math>\epsilon = \int {\rm d}^3r (\rho^{out}(\bf r) -\rho^{inp}(\bf r))^2</math>.
    !% A zero value means do not use this criterion.
    !%End
    call loct_parse_float(check_inp('ConvAbsDens'), CNST(1e-5), scf%conv_abs_dens)

    !%Variable ConvRelDens
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the density:
    !% <math>\epsilon = {1\over N} \int {\rm d}^3r (\rho^{out}(\bf r) -\rho^{inp}(\bf r))^2</math>.
    !% <i>N</i> is the total number of electrons in the problem.
    !% A zero value means do not use this criterion.
    !%End
    call loct_parse_float(check_inp('ConvRelDens'), M_ZERO, scf%conv_rel_dens)

    !%Variable ConvAbsEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the eigenvalues:
    !% <math>\epsilon = \sum_{j=1}^{N_{occ}} \vert \epsilon_j^{out}-\epsilon_j^{inp}\vert</math>.
    !% A zero value means do not use this criterion.
    !%End
    call loct_parse_float(check_inp('ConvAbsEv'), M_ZERO, scf%conv_abs_ev)

    !%Variable ConvRelEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the eigenvalues:
    !% <math>\epsilon = {1 \over N} \sum_{j=1}^{N_{occ}} \vert \epsilon_j^{out}-\epsilon_j^{inp}\vert</math>.
    !% <i>N</i> is the total number of electrons. A zero value means do not use this criterion._m
    !%End
    call loct_parse_float(check_inp('ConvRelEv'), M_ZERO, scf%conv_rel_ev)

    if(scf%max_iter <= 0 .and. &
      scf%conv_abs_dens <= M_ZERO .and. scf%conv_rel_dens <= M_ZERO .and. &
      scf%conv_abs_ev <= M_ZERO .and. scf%conv_rel_ev <= M_ZERO) then
      message(1) = "Input: Not all convergence criteria can be <= 0"
      message(2) = "Please set one of the following:"
      message(3) = "MaximumIter | ConvAbsDens | ConvRelDens | ConvAbsEv | ConvRelEv"
      call write_fatal(3)
    end if

    if(scf%max_iter <= 0) scf%max_iter = huge(scf%max_iter)

    !%Variable What2Mix
    !%Type integer
    !%Default density
    !%Section SCF::Mixing
    !%Description
    !% Selects what should be mixed during the SCF cycle.
    !%Option potential 1
    !% The Kohn-Sham potential
    !%Option density 2
    !% The density
    !%End
    call loct_parse_int(check_inp('What2Mix'), MIXDENS, scf%what2mix)
    if(.not.varinfo_valid_option('What2Mix', scf%what2mix)) call input_error('What2Mix')
    call messages_print_var_option(stdout, "What2Mix", scf%what2mix, "what to mix during SCF cycles")

    if (scf%what2mix == MIXPOT.and.h%ip_app) then
      message(1) = "Input: Cannot mix the potential with non-interacting electrons."
      call write_fatal(1)
    end if

    ! Handle mixing now...
    dim = 1
    if (h%d%cdft) dim = 1 + NDIM
    call mix_init(scf%smix, gr%m, dim, st%d%nspin)

    ! now the eigen solver stuff
    call eigen_solver_init(gr, scf%eigens, st, 25)

    !%Variable SCFinLCAO
    !%Type logical
    !%Default no
    !%Section SCF
    !%Description
    !% Performs all the SCF cycle restricting the calculation to the LCAO subspace.
    !% This may be useful for systems with convergence problems (first do a 
    !% calculation within the LCAO subspace, then restart from that point for
    !% an unrestricted calculation).
    !%End
    call loct_parse_logical(check_inp('SCFinLCAO'), .false., scf%lcao_restricted)
    if(scf%lcao_restricted) then
      message(1) = 'Info: SCF restricted to LCAO subspace'
      call write_info(1)
    end if

    call geometry_min_distance(gr%geo, rmin)
    if(gr%geo%natoms == 1) rmin = CNST(100.0)
    call loct_parse_float(check_inp('LocalMagneticMomentsSphereRadius'), rmin*M_HALF/units_inp%length%factor, scf%lmm_r)
    scf%lmm_r = scf%lmm_r * units_inp%length%factor

    call pop_sub()
  end subroutine scf_init


  ! ---------------------------------------------------------
  subroutine scf_end(scf)
    type(scf_t), intent(inout) :: scf

    call push_sub('scf.scf_end')

    call eigen_solver_end(scf%eigens)
    call mix_end(scf%smix)

    call pop_sub()
  end subroutine scf_end


  ! ---------------------------------------------------------
  subroutine scf_run(scf, gr, st, ks, h, outp)
    type(scf_t),         intent(inout) :: scf
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: st
    type(v_ks_t),        intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: h
    type(output_t),      intent(in)    :: outp

    type(lcao_t) :: lcao_data

    integer :: iter, iunit, is, idim, nspin, dim, err
    FLOAT :: evsum_out, evsum_in
    FLOAT, allocatable :: rhoout(:,:,:), rhoin(:,:,:), rhonew(:,:,:)
    FLOAT, allocatable :: vout(:,:,:), vin(:,:,:), vnew(:,:,:)
    FLOAT, allocatable :: tmp(:)
    character(len=8) :: dirname
    logical :: finish, forced_finish

    call push_sub('scf.scf_run')

    if(scf%lcao_restricted) then
      call lcao_init(gr, lcao_data, st, h)
      if(.not.lcao_data%state == 1) then
        message(1) = 'Nothing to do'
        call write_fatal(1)
      end if
    end if

    nspin = st%d%nspin

    dim = 1
    if (h%d%cdft) dim = 1 + NDIM

    ALLOCATE(rhoout(NP, dim, nspin), NP*dim*nspin)
    ALLOCATE(rhoin (NP, dim, nspin), NP*dim*nspin)

    rhoin(1:NP, 1, 1:nspin) = st%rho(1:NP, 1:nspin)
    rhoout = M_ZERO
    if (st%d%cdft) then
      rhoin(1:NP, 2:dim, 1:nspin) = st%j(1:NP, 1:NDIM, 1:nspin)
    end if

    if (scf%what2mix == MIXPOT) then
      ALLOCATE(vout(NP, dim, nspin), NP*dim*nspin)
      ALLOCATE( vin(NP, dim, nspin), NP*dim*nspin)
      ALLOCATE(vnew(NP, dim, nspin), NP*dim*nspin)

      vin(:, 1, :) = h%vhxc
      vout = M_ZERO
      if (st%d%cdft) vin(:, 2:dim, :) = h%axc(:,:,:)
    else
      ALLOCATE(rhonew(NP, dim, nspin), NP*dim*nspin)
    end if
    evsum_in = states_eigenvalues_sum(st)

    ! SCF cycle
    do iter = 1, scf%max_iter
      call profiling_in(C_PROFILING_SCF_CYCLE)

      if(scf%lcao_restricted) then
        call lcao_wf(lcao_data, st, gr%m, gr%sb, h)
      else
        scf%eigens%converged = 0
        call eigen_solver_run(scf%eigens, gr, st, h, iter)
      end if

      ! occupations
      call states_fermi(st, gr%m)

      ! compute output density, potential (if needed) and eigenvalues sum
      call X(states_calc_dens)(st, NP, st%rho)
      rhoout(1:NP, 1, :) = st%rho(1:NP, :)
      if (h%d%cdft) then
        call states_calc_physical_current(gr, st, st%j)
        rhoout(1:NP, 2:dim, :) = st%j(1:NP, 1:NDIM, :)
      end if
      if (scf%what2mix == MIXPOT) then
        call X(v_ks_calc) (gr, ks, h, st)
        vout(:, 1, :) = h%vhxc
        if (h%d%cdft) vout(:, 2:dim, :) = h%axc(:,:,:)
      end if
      evsum_out = states_eigenvalues_sum(st)

      ! recalculate total energy
      call hamiltonian_energy(h, st, gr%geo%eii, 0)

      ! compute convergence criteria
      scf%abs_dens = M_ZERO
      ALLOCATE(tmp(NP), NP)
      do is = 1, nspin
        do idim = 1, dim
          tmp = (rhoin(:, idim, is) - rhoout(:, idim, is))**2
          scf%abs_dens = scf%abs_dens + dmf_integrate(gr%m, tmp)
        end do
      end do
      deallocate(tmp)

      scf%abs_dens = sqrt(scf%abs_dens)
      scf%rel_dens = scf%abs_dens / st%qtot
      scf%abs_ev = abs(evsum_out - evsum_in)
      scf%rel_ev = scf%abs_ev / abs(evsum_out)

      ! are we finished?
      finish = &
        (scf%conv_abs_dens > M_ZERO .and. scf%abs_dens <= scf%conv_abs_dens) .or. &
        (scf%conv_rel_dens > M_ZERO .and. scf%rel_dens <= scf%conv_rel_dens) .or. &
        (scf%conv_abs_ev   > M_ZERO .and. scf%abs_ev   <= scf%conv_abs_ev)   .or. &
        (scf%conv_rel_ev   > M_ZERO .and. scf%rel_ev   <= scf%conv_rel_ev)

      call scf_write_iter

      ! mixing
      select case (scf%what2mix)
      case (MIXDENS)
        ! mix input and output densities and compute new potential
        call mixing(scf%smix, gr%m, iter, dim, nspin, rhoin, rhoout, rhonew)
        st%rho(1:NP,:) = rhonew(1:NP, 1, :)
        if (h%d%cdft) st%j(1:NP,:,:) = rhonew(1:NP, 2:dim, :)
        call X(v_ks_calc) (gr, ks, h, st)
      case (MIXPOT)
        ! mix input and output potentials
        call mixing(scf%smix, gr%m, iter, dim, nspin, vin, vout, vnew)
        h%vhxc = vnew(:, 1, :)
        if (h%d%cdft) h%axc(:,:,:) = vnew(:,2:dim,:)
      end select

      ! Are we asked to stop? (Whenever Fortran is ready for signals, this should go away)
      forced_finish = clean_stop()
 
      ! save restart information
      if(finish.or.(modulo(iter, 3) == 0).or.iter==scf%max_iter.or.forced_finish) then
        call X(restart_write) (trim(tmpdir)//'restart_gs', st, gr, err, iter=iter)
        if(err.ne.0) then
          message(1) = 'Unsuccesfull write of "'//trim(tmpdir)//'restart_gs"'
          call write_fatal(1)
        end if
      end if

      if(finish) then
        write(message(1), '(a, i4, a)')'Info: SCF converged in ', iter, ' iterations'
        call write_info(1)
        if(scf%lcao_restricted) call lcao_end(lcao_data, st%nst)
        call profiling_out(C_PROFILING_SCF_CYCLE)
        exit
      end if

      if(outp%duringscf) then
        write(dirname,'(a,i4.4)') "scf.",iter
        call X(states_output) (st, gr, dirname, outp)
        call hamiltonian_output(h, gr%m, gr%sb, dirname, outp)
      end if

      ! save information for the next iteration
      rhoin(:, 1, :) = st%rho(1:NP, :)
      if (h%d%cdft) rhoin(1:NP, 2:dim, :) = st%j(1:NP, 1:NDIM, :)
      if (scf%what2mix == MIXPOT) then
        vin(:, 1, :) = h%vhxc
        if (h%d%cdft) vin(:,2:dim,:) = h%axc(:,:,:)
      end if
      evsum_in = evsum_out

      if(forced_finish) then
        call profiling_out(C_PROFILING_SCF_CYCLE)
        exit
      end if

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()


      call profiling_out(C_PROFILING_SCF_CYCLE)
    end do

    if (scf%what2mix == MIXPOT) then
      call X(v_ks_calc) (gr, ks, h, st)
      deallocate(vout, vin, vnew)
    else
      deallocate(rhonew)
    end if
    deallocate(rhoout, rhoin)

    if(.not.finish) then
      message(1) = 'SCF *not* converged!'
      call write_warning(1)
    end if

    ! calculate forces
    call X(epot_forces)(gr, h%ep, st)

    ! output final information
    call scf_write_static("static", "info")
    call X(states_output) (st, gr, "static", outp)
    if(iand(outp%what, output_geometry).ne.0) &
      call atom_write_xyz("static", "geometry", gr%geo)
    call hamiltonian_output(h, gr%m, gr%sb, "static", outp)

    if(simul_box_is_periodic(gr%sb).and.st%d%nik > st%d%nspin) then
      iunit = io_open('static/bands.dat', action='write')
      call states_write_bands(iunit, st%nst, st, gr%sb)
      call io_close(iunit)
    end if

    call pop_sub()

  contains


    ! ---------------------------------------------------------
    subroutine scf_write_iter
      character(len=50) :: str

      call push_sub('scf.scf_write_iter')

      write(str, '(a,i5)') 'SCF CYCLE ITER #' ,iter
      call messages_print_stress(stdout, trim(str))

      write(message(1),'(a,es15.8,2(a,es9.2))') ' etot = ', h%etot/units_out%energy%factor, &
        ' abs_ev   = ', scf%abs_ev/units_out%energy%factor, ' rel_ev   = ', scf%rel_ev
      write(message(2),'(23x,2(a,es9.2))') &
        ' abs_dens = ', scf%abs_dens, ' rel_dens = ', scf%rel_dens
      call write_info(2)

      if(.not.scf%lcao_restricted) then
        write(message(1),'(a,i6)') 'Matrix vector products: ', scf%eigens%matvec
        write(message(2),'(a,i6)') 'Converged eigenvectors: ', scf%eigens%converged
        call write_info(2)
        call states_write_eigenvalues(stdout, st%nst, st, gr%sb, scf%eigens%diff)
      else
        call states_write_eigenvalues(stdout, st%nst, st, gr%sb)
      end if

      if(st%d%ispin > UNPOLARIZED) then
        call write_magnetic_moments(stdout, gr%m, st)
      end if

      call messages_print_stress(stdout)

      call pop_sub()
    end subroutine scf_write_iter


    ! ---------------------------------------------------------
    subroutine scf_write_static(dir, fname)
      character(len=*), intent(in) :: dir, fname

      FLOAT :: e_dip(4, st%d%nspin), n_dip(MAX_DIM)
      FLOAT, parameter :: ATOMIC_TO_DEBYE = CNST(2.5417462)
      integer :: iunit, i, j

      call push_sub('scf.scf_write_static')

      if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
        call io_mkdir(dir)
        iunit = io_open(trim(dir) // "/" // trim(fname), action='write')

        write(iunit, '(a,a)') 'System name: ', gr%geo%sysname

        call grid_write_info(gr, iunit)

        if(simul_box_is_periodic(gr%sb)) then
          call kpoints_write_info(st%d, iunit)
          write(iunit,'(1x)')
        end if

        if(.not. h%ip_app) then
          call v_ks_write_info(ks, iunit)
        else
          write(iunit, '(a)') 'Independent Particles'
        end if

        ! scf information
        if(finish) then
          write(iunit, '(a, i4, a)')'SCF converged in ', iter, ' iterations'
        else
          write(iunit, '(a)') 'SCF *not* converged!'
        end if
        write(iunit, '(1x)')

        call states_write_eigenvalues(iunit, st%nst, st, gr%sb)
        write(iunit, '(1x)')

        write(iunit, '(a)') 'Energy:'
      end if

      call hamiltonian_energy(h, st, gr%geo%eii, iunit)

      if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
      if(st%d%ispin > UNPOLARIZED) then
        call write_magnetic_moments(iunit, gr%m, st)
        if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
      end if


      ! Next lines of code calculate the dipole of the molecule, summing the electronic and
      ! ionic contributions.
      call states_calculate_multipoles(gr, st, 1, e_dip(:, :))
      do j = 1, 3
        e_dip(j+1, 1) = sum(e_dip(j+1, :))
      end do
      call geometry_dipole(gr%geo, n_dip)
      n_dip(1:NDIM) = n_dip(1:NDIM) - e_dip(2:NDIM+1, 1)

      if(mpi_grp_is_root(mpi_world)) then
        write(iunit, '(3a)') 'Dipole [', trim(units_out%length%abbrev), ']:                    [Debye]'
        do j = 1, NDIM
          write(iunit, '(6x,a,i1,a,es14.5,3x,2es14.5)') '<x', j, '> = ', n_dip(j) / units_out%length%factor, &
            n_dip(j)*ATOMIC_TO_DEBYE
        end do
        write(iunit,'(a)')
      end if

      ! Next is the angular momentum. Only applies to 2D and 3D.
      if(NDIM.ne.1) call write_angular_momentum(iunit)
      if(mpi_grp_is_root(mpi_world)) write(iunit, '(a)')

      if(mpi_grp_is_root(mpi_world)) then
        write(iunit, '(a)') 'Convergence:'
        write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'abs_dens = ', scf%abs_dens, &
          ' (', scf%conv_abs_dens, ')'
        write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'rel_dens = ', scf%rel_dens, &
          ' (', scf%conv_rel_dens, ')'
        write(iunit, '(6x, a, es14.8,a,es14.8,4a)') 'abs_ev = ', scf%abs_ev, &
          ' (', scf%conv_abs_ev / units_out%energy%factor, ')', &
          ' [',  trim(units_out%energy%abbrev), ']'
        write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'rel_ev = ', scf%rel_ev, &
          ' (', scf%conv_rel_ev, ')'
        write(iunit,'(1x)')

        write(iunit,'(3a)') 'Forces on the ions [', trim(units_out%force%abbrev), "]"
        write(iunit,'(a,10x,14x,a,14x,a,14x,a)') ' Ion','x','y','z'
        do i = 1, gr%geo%natoms
          write(iunit,'(i4,a10,3f15.6)') i, trim(gr%geo%atom(i)%spec%label), &
            gr%geo%atom(i)%f(:) / units_out%force%factor
        end do

        call io_close(iunit)
      end if

      call pop_sub()
    end subroutine scf_write_static


    ! ---------------------------------------------------------
    subroutine write_angular_momentum(iunit)
      integer,        intent(in) :: iunit

      integer :: ik, ist, ns, j
      character(len=80) tmp_str(MAX_DIM), cspin
      FLOAT :: angular(3), lsquare, o, oplus, ominus
      FLOAT, allocatable :: ang(:, :, :), ang2(:, :)

      ns = 1
      if(st%d%nspin == 2) ns = 2

      if(mpi_grp_is_root(mpi_world)) then
        write(iunit,'(a)') 'Angular Momentum of the KS states [adimensional]:'
        if (st%d%nik > ns) then
          message(1) = 'Kpoints [' // trim(units_out%length%abbrev) // '^-1]'
          call write_info(1, iunit)
        end if
      end if

      ALLOCATE(ang (st%st_start:st%st_end, st%d%nik, 3), (st%st_end - st%st_start + 1)*st%d%nik*3)
      ALLOCATE(ang2(st%st_start:st%st_end, st%d%nik), (st%st_end - st%st_start + 1)*st%d%nik)
      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          call X(states_angular_momentum)(gr, st%X(psi)(:, :, ist, ik), ang(ist, ik, :), ang2(ist, ik))
        end do
      end do
      angular(1) =  states_eigenvalues_sum(st, ang(:, :, 1))
      angular(2) =  states_eigenvalues_sum(st, ang(:, :, 2))
      angular(3) =  states_eigenvalues_sum(st, ang(:, :, 3))
      lsquare    =  states_eigenvalues_sum(st, ang2)

      do ik = 1, st%d%nik, ns
        if(st%d%nik > ns) then
          write(message(1), '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
            st%d%kpoints(1, ik)*units_out%length%factor, ',',            &
            st%d%kpoints(2, ik)*units_out%length%factor, ',',            &
            st%d%kpoints(3, ik)*units_out%length%factor, ')'
          call write_info(1, iunit)
        end if

        write(message(1), '(a4,1x,a5,1x,4a12,4x,a12,1x)')       &
          '#st',' Spin','        <Lx>', '        <Ly>', '        <Lz>', '        <L2>', 'Occupation '
        call write_info(1, iunit)

        do j = 1, st%nst
          do is = 0, ns-1

            if(j > st%nst) then
              o = M_ZERO
              if(st%d%ispin == SPINORS) oplus = M_ZERO; ominus = M_ZERO
            else
              o = st%occ(j, ik+is)
              if(st%d%ispin == SPINORS) then
                oplus  = st%mag(j, ik+is, 1)
                ominus = st%mag(j, ik+is, 2)
              end if
            end if

            if(is.eq.0) cspin = 'up'
            if(is.eq.1) cspin = 'dn'
            if(st%d%ispin.eq.UNPOLARIZED.or.st%d%ispin.eq.SPINORS) cspin = '--'

            write(tmp_str(1), '(i4,3x,a2)') j, trim(cspin)
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,4f12.6,3x,f5.2,a1,f5.2)') &
                ang(j, 1:st%d%nik, ik), ang2(j, ik), oplus, '/', ominus
            else
              write(tmp_str(2), '(1x,4f12.6,3x,f12.6)') &
                ang(j, 1:st%d%nik, ik+is), ang2(j, ik+is), o
            end if
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))
            call write_info(1, iunit)
          end do
        end do

      end do

      write(message(1),'(a)') 'Total Angular Momentum L [adimensional]'
      write(message(2),'(10x,4f12.6)') angular(1:3), lsquare
      call write_info(2, iunit)

      deallocate(ang, ang2)
    end subroutine write_angular_momentum


    ! ---------------------------------------------------------
    subroutine write_magnetic_moments(iunit, m, st)
      integer,        intent(in) :: iunit
      type(mesh_t),   intent(in) :: m
      type(states_t), intent(in) :: st

      integer :: i
      FLOAT :: mm(MAX_DIM)
      FLOAT, allocatable :: lmm(:,:)

      call push_sub('scf.write_magnetic_moments')

      call states_magnetic_moment(m, st, st%rho, mm)
      ALLOCATE(lmm(3, gr%geo%natoms), 3*gr%geo%natoms)
      call states_local_magnetic_moments(m, st, gr%geo, st%rho, scf%lmm_r, lmm)

      if(mpi_grp_is_root(mpi_world)) then

        write(iunit, '(a)') 'Total Magnetic Moment:'
        if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
          write(iunit, '(a,f10.6)') ' mz = ', mm(3)
        else if(st%d%ispin == SPINORS) then ! non-collinear
          write(iunit, '(1x,3(a,f10.6,3x))') 'mx = ',mm(1),'my = ',mm(2),'mz = ',mm(3)
        end if

        write(iunit, '(a,a,a,f7.3,a)') 'Local Magnetic Moments (sphere radius [', &
             trim(units_out%length%abbrev),'] = ', scf%lmm_r/units_out%length%factor, '):'
        if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
          write(iunit,'(a,6x,14x,a)') ' Ion','mz'
          do i = 1, gr%geo%natoms
            write(iunit,'(i4,a10,f15.6)') i, trim(gr%geo%atom(i)%spec%label), lmm(3, i)
          end do
        else if(st%d%ispin == SPINORS) then ! non-collinear
          write(iunit,'(a,8x,13x,a,13x,a,13x,a)') ' Ion','mx','my','mz'
          do i = 1, gr%geo%natoms
            write(iunit,'(i4,a10,3f15.6)') i, trim(gr%geo%atom(i)%spec%label), lmm(1, i), lmm(2, i), lmm(3, i)
          end do
        end if

      end if
      
      deallocate(lmm)

      call pop_sub()
    end subroutine write_magnetic_moments

  end subroutine scf_run

end module scf_m
