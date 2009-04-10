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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"

module scf_m
  use datasets_m
  use energy_m
  use eigensolver_m
  use external_pot_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use lcao_m
  use loct_m
  use loct_parser_m
  use magnetic_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mix_m
  use mpi_m
  use mpi_lib_m
  use multigrid_m
  use h_sys_output_m
  use ob_lippmann_schwinger_m
  use preconditioners_m
  use profiling_m
  use restart_m
  use simul_box_m
  use solids_m
  use states_m
  use states_dim_m
  use states_calc_m
  use units_m
  use v_ks_m
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

  integer, public, parameter :: &
    VERB_NO      = 0,   &
    VERB_COMPACT = 1,   &
    VERB_FULL    = 3
  
  type scf_t      ! some variables used for the scf cycle
    integer :: max_iter   ! maximum number of scf iterations

    FLOAT :: lmm_r

    ! several convergence criteria
    FLOAT :: conv_abs_dens, conv_rel_dens, conv_abs_ev, conv_rel_ev, conv_abs_force
    FLOAT :: abs_dens, rel_dens, abs_ev, rel_ev, abs_force

    integer :: what2mix
    logical :: lcao_restricted

    type(mix_t) :: smix
    type(eigensolver_t) :: eigens
  end type scf_t

contains

  ! ---------------------------------------------------------
  subroutine scf_init(scf, gr, geo, st, hm)
    type(scf_t),         intent(inout) :: scf
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(in)    :: geo
    type(states_t),      intent(in)    :: st
    type(hamiltonian_t), intent(inout) :: hm

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
    call loct_parse_int  (datasets_check('MaximumIter'), 200, scf%max_iter)

    !%Variable ConvAbsDens
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the density: 
    !%
    !% <math>\epsilon = \int {\rm d}^3r \vert \rho^{out}(\bf r) -\rho^{inp}(\bf r) \vert</math>.
    !%
    !% A zero value (the default) means do not use this criterion.
    !%End
    call loct_parse_float(datasets_check('ConvAbsDens'), M_ZERO, scf%conv_abs_dens)

    !%Variable ConvRelDens
    !%Type float
    !%Default 1e-3
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the density: 
    !%
    !% <math>\epsilon = {1\over N} ConvAbsDens</math>.
    !% 
    !% <i>N</i> is the total number of electrons in the problem.  A
    !% zero value means do not use this criterion. By default this
    !% value is set to 1e-3.
    !%End
    call loct_parse_float(datasets_check('ConvRelDens'), CNST(1e-3), scf%conv_rel_dens)

    !%Variable ConvAbsEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the eigenvalues:
    !%
    !% <math>\epsilon = \sum_{j=1}^{N_{occ}} \vert \epsilon_j^{out}-\epsilon_j^{inp}\vert</math>.
    !%
    !% A zero value (the default) means do not use this criterion.
    !%End
    call loct_parse_float(datasets_check('ConvAbsEv'), M_ZERO, scf%conv_abs_ev)
    scf%conv_abs_ev = scf%conv_abs_ev * units_inp%energy%factor

    !%Variable ConvRelEv
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Relative convergence of the eigenvalues:
    !%
    !% <math>\epsilon = {1 \over E} \sum_{j=1}^{N_{occ}} \vert \epsilon_j^{out}-\epsilon_j^{inp}\vert</math>.
    !%
    !% <i>E</i> is the sum of the eigenvalues. A zero value (the default) means do not use this criterion.
    !%End
    call loct_parse_float(datasets_check('ConvRelEv'), M_ZERO, scf%conv_rel_ev)

    call obsolete_variable("ConvAbsForce", "ConvForce")
    call obsolete_variable("ConvRelForce", "ConvForce")

    !%Variable ConvForce
    !%Type float
    !%Default 0.0
    !%Section SCF::Convergence
    !%Description
    !% Absolute convergence of the forces: 
    !% maximum variation of any component of the ionic forces in consecutive iterations.
    !% A zero value (the default) means do not use this criterion.
    !%End
    call loct_parse_float(datasets_check('ConvForce'), M_ZERO, scf%conv_abs_force)
    scf%conv_abs_force = scf%conv_abs_force * units_inp%force%factor

    if(scf%max_iter <= 0 .and. &
      scf%conv_abs_dens <= M_ZERO .and. scf%conv_rel_dens <= M_ZERO .and. &
      scf%conv_abs_ev <= M_ZERO .and. scf%conv_rel_ev <= M_ZERO .and. &
      scf%conv_abs_force <= M_ZERO) then
      message(1) = "Input: Not all convergence criteria can be <= 0"
      message(2) = "Please set one of the following:"
      message(3) = "MaximumIter | ConvAbsDens | ConvRelDens | ConvAbsEv | ConvRelEv | ConvForce "
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
    call loct_parse_int(datasets_check('What2Mix'), MIXDENS, scf%what2mix)
    if(.not.varinfo_valid_option('What2Mix', scf%what2mix)) call input_error('What2Mix')
    call messages_print_var_option(stdout, "What2Mix", scf%what2mix, "what to mix during SCF cycles")

    if (scf%what2mix == MIXPOT.and.hm%theory_level==INDEPENDENT_PARTICLES) then
      message(1) = "Input: Cannot mix the potential with non-interacting particles."
      call write_fatal(1)
    end if

    ! Handle mixing now...
    dim = 1
    if (hm%d%cdft) dim = 1 + gr%mesh%sb%dim
    call mix_init(scf%smix, gr%mesh%np, dim, st%d%nspin)
    call mesh_init_mesh_aux(gr%mesh)

    ! now the eigensolver stuff
    call eigensolver_init(gr, scf%eigens, st)

    if(scf%eigens%es_type == RS_MG .or. preconditioner_is_multigrid(scf%eigens%pre)) then
      if(.not. associated(gr%mgrid)) then
        ALLOCATE(gr%mgrid, 1)
        call multigrid_init(gr%mgrid, geo, gr%cv,gr%mesh, gr%der, gr%stencil)
      end if
      call hamiltonian_mg_init(hm, gr)
    end if

    !%Variable SCFinLCAO
    !%Type logical
    !%Default no
    !%Section SCF
    !%Description
    !% Performs the SCF cycle with the calculation restricted to the LCAO subspace.
    !% This may be useful for systems with convergence problems (first do a 
    !% calculation within the LCAO subspace, then restart from that point for
    !% an unrestricted calculation).
    !%End
    call loct_parse_logical(datasets_check('SCFinLCAO'), .false., scf%lcao_restricted)
    if(scf%lcao_restricted) then
      message(1) = 'Info: SCF restricted to LCAO subspace'
      call write_info(1)
    end if

    call geometry_min_distance(geo, rmin)
    if(geo%natoms == 1) rmin = CNST(100.0)
    call loct_parse_float(datasets_check('LocalMagneticMomentsSphereRadius'), rmin*M_HALF/units_inp%length%factor, scf%lmm_r)
    scf%lmm_r = scf%lmm_r * units_inp%length%factor

    call pop_sub()
  end subroutine scf_init


  ! ---------------------------------------------------------
  subroutine scf_end(scf)
    type(scf_t), intent(inout) :: scf

    call push_sub('scf.scf_end')

    call eigensolver_end(scf%eigens)
    call mix_end(scf%smix)

    call pop_sub()
  end subroutine scf_end


  ! ---------------------------------------------------------
  subroutine scf_run(scf, gr, geo, st, ks, hm, outp, gs_run, verbosity)
    type(scf_t),         intent(inout) :: scf
    type(grid_t),        intent(inout) :: gr
    type(geometry_t),    intent(inout) :: geo
    type(states_t),      intent(inout) :: st
    type(v_ks_t),        intent(inout) :: ks
    type(hamiltonian_t), intent(inout) :: hm
    type(h_sys_output_t),intent(in)    :: outp
    logical, optional,   intent(in)    :: gs_run
    integer, optional,   intent(in)    :: verbosity 

    type(lcao_t) :: lcao
    type(profile_t), save :: prof

    integer :: iter, is, idim, iatom, nspin, dim, err
    FLOAT :: evsum_out, evsum_in, forcetmp
    real(8) :: etime, itime
    FLOAT, allocatable :: rhoout(:,:,:), rhoin(:,:,:), rhonew(:,:,:)
    FLOAT, allocatable :: vout(:,:,:), vin(:,:,:), vnew(:,:,:)
    FLOAT, allocatable :: forceout(:,:), forcein(:,:), forcediff(:), tmp(:)
    character(len=8) :: dirname
    logical :: finish, forced_finish, gs_run_
    integer :: verbosity_

    call push_sub('scf.scf_run')

    if(present(gs_run)) then 
      gs_run_ = gs_run
    else 
      gs_run_ = .true.
    end if
    
    verbosity_ = VERB_FULL
    if(present(verbosity)) verbosity_ = verbosity

    if(scf%lcao_restricted) then
      call lcao_init(lcao, gr, geo, st)
      if(.not. lcao_is_available(lcao)) then
        message(1) = 'Nothing to do'
        call write_fatal(1)
      end if
    end if

    nspin = st%d%nspin

    dim = 1
    if (hm%d%cdft) dim = 1 + gr%mesh%sb%dim

    ALLOCATE(rhoout(gr%mesh%np, dim, nspin), gr%mesh%np*dim*nspin)
    ALLOCATE(rhoin (gr%mesh%np, dim, nspin), gr%mesh%np*dim*nspin)

    rhoin(1:gr%mesh%np, 1, 1:nspin) = st%rho(1:gr%mesh%np, 1:nspin)
    rhoout = M_ZERO
    if (st%d%cdft) then
      rhoin(1:gr%mesh%np, 2:dim, 1:nspin) = st%j(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
    end if

    if (scf%what2mix == MIXPOT) then
      ALLOCATE(vout(gr%mesh%np, dim, nspin), gr%mesh%np*dim*nspin)
      ALLOCATE( vin(gr%mesh%np, dim, nspin), gr%mesh%np*dim*nspin)
      ALLOCATE(vnew(gr%mesh%np, dim, nspin), gr%mesh%np*dim*nspin)

      vin(1:gr%mesh%np, 1, 1:nspin) = hm%vhxc(1:gr%mesh%np, 1:nspin)
      vout = M_ZERO
      if (st%d%cdft) vin(1:gr%mesh%np, 2:dim, 1:nspin) = hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
    else
      ALLOCATE(rhonew(gr%mesh%np, dim, nspin), gr%mesh%np*dim*nspin)
    end if
    evsum_in = states_eigenvalues_sum(st)

    ! allocate and compute forces only if they are used as convergence criteria
    if (scf%conv_abs_force > M_ZERO) then
      ALLOCATE(forcein(geo%natoms, gr%mesh%sb%dim), geo%natoms*gr%mesh%sb%dim)
      ALLOCATE(forceout(geo%natoms, gr%mesh%sb%dim), geo%natoms*gr%mesh%sb%dim)
      ALLOCATE(forcediff(gr%mesh%sb%dim), gr%mesh%sb%dim)
      call epot_forces(gr, geo, hm%ep, st)
      do iatom = 1, geo%natoms
        forcein(iatom,1:gr%mesh%sb%dim) = geo%atom(iatom)%f(1:gr%mesh%sb%dim)
      end do
    endif

    if ( verbosity_ /= VERB_NO ) then
      write(message(1),'(a)') 'Info: Starting SCF iteration'
      call write_info(1)
    end if

    ! SCF cycle
    itime = loct_clock()
    do iter = 1, scf%max_iter
      call profiling_in(prof, "SCF_CYCLE")

      if(scf%lcao_restricted) then
        call lcao_wf(lcao, st, gr, geo, hm)
      else
        ! FIXME: Currently, only the eigensolver or the
        ! Lippmann-Schwinger approach can be used (exclusively),
        ! i. e. no bound states for open boundaries.
        if(gr%sb%open_boundaries) then
          call lippmann_schwinger(scf%eigens, hm, gr, st)
        else
          scf%eigens%converged = 0
          call eigensolver_run(scf%eigens, gr, st, hm, iter)
        end if
      end if

      ! occupations
      call states_fermi(st, gr%mesh)

      ! compute output density, potential (if needed) and eigenvalues sum
      call states_calc_dens(st, gr%mesh%np)
      rhoout(1:gr%mesh%np, 1, 1:nspin) = st%rho(1:gr%mesh%np, 1:nspin)
      if (hm%d%cdft) then
        call calc_physical_current(gr, st, st%j)
        rhoout(1:gr%mesh%np, 2:dim, 1:nspin) = st%j(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      end if
      if (scf%what2mix == MIXPOT) then
        call v_ks_calc(gr, ks, hm, st)
        vout(1:gr%mesh%np, 1, 1:nspin) = hm%vhxc(1:gr%mesh%np, 1:nspin)
        if (hm%d%cdft) vout(1:gr%mesh%np, 2:dim, 1:nspin) = hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      end if
      evsum_out = states_eigenvalues_sum(st)

      ! recalculate total energy
      call total_energy(hm, gr, st, iunit = 0)

      ! compute convergence criteria
      scf%abs_dens = M_ZERO
      ALLOCATE(tmp(gr%mesh%np), gr%mesh%np)
      do is = 1, nspin
        do idim = 1, dim
          tmp = abs(rhoin(1:gr%mesh%np, idim, is) - rhoout(1:gr%mesh%np, idim, is))
          scf%abs_dens = scf%abs_dens + dmf_integrate(gr%mesh, tmp)
        end do
      end do
      deallocate(tmp)

      ! compute forces only if they are used as convergence criteria
      if (scf%conv_abs_force > M_ZERO) then
        call epot_forces(gr, geo, hm%ep, st)
        scf%abs_force = M_ZERO
        do iatom = 1, geo%natoms
          forceout(iatom,1:gr%mesh%sb%dim) = geo%atom(iatom)%f(1:gr%mesh%sb%dim)
          forcediff(1:gr%mesh%sb%dim) = abs( forceout(iatom,1:gr%mesh%sb%dim) - forcein(iatom,1:gr%mesh%sb%dim) )
          forcetmp = maxval( forcediff )
          if ( forcetmp > scf%abs_force ) then
            scf%abs_force = forcetmp
          end if
        end do
      end if

      scf%abs_dens = sqrt(scf%abs_dens)
      scf%rel_dens = scf%abs_dens / st%qtot
      scf%abs_ev = abs(evsum_out - evsum_in)
      scf%rel_ev = scf%abs_ev / abs(evsum_out)

      ! are we finished?
      finish = &
        (scf%conv_abs_dens  <= M_ZERO .or. scf%abs_dens  <= scf%conv_abs_dens)  .and. &
        (scf%conv_rel_dens  <= M_ZERO .or. scf%rel_dens  <= scf%conv_rel_dens)  .and. &
        (scf%conv_abs_force <= M_ZERO .or. scf%abs_force <= scf%conv_abs_force) .and. &
        (scf%conv_abs_ev    <= M_ZERO .or. scf%abs_ev    <= scf%conv_abs_ev)    .and. &
        (scf%conv_rel_ev    <= M_ZERO .or. scf%rel_ev    <= scf%conv_rel_ev)

      etime = loct_clock() - itime
      itime = etime + itime
      call scf_write_iter

      ! mixing
      select case (scf%what2mix)
      case (MIXDENS)
        ! mix input and output densities and compute new potential
        call dmixing(scf%smix, iter, rhoin, rhoout, rhonew, dmf_dotp_aux)
        st%rho(1:gr%mesh%np,1:nspin) = rhonew(1:gr%mesh%np, 1, 1:nspin)
        if (hm%d%cdft) st%j(1:gr%mesh%np,1:gr%mesh%sb%dim,1:nspin) = rhonew(1:gr%mesh%np, 2:dim, 1:nspin)
        call v_ks_calc(gr, ks, hm, st)
      case (MIXPOT)
        ! mix input and output potentials
        call dmixing(scf%smix, iter, vin, vout, vnew, dmf_dotp_aux)
        hm%vhxc(1:gr%mesh%np, 1:nspin) = vnew(1:gr%mesh%np, 1, 1:nspin)
        call hamiltonian_update_potential(hm, gr%mesh)
        if (hm%d%cdft) hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin) = vnew(1:gr%mesh%np, 2:dim, 1:nspin)
      end select

      ! Are we asked to stop? (Whenever Fortran is ready for signals, this should go away)
      forced_finish = clean_stop()

      if(gs_run_) then 
        ! save restart information
        if(finish.or.(modulo(iter, 3) == 0).or.iter==scf%max_iter.or.forced_finish) then
          call restart_write(trim(tmpdir)//'gs', st, gr, err, iter=iter)
          if(err.ne.0) then
            message(1) = 'Unsuccessful write of "'//trim(tmpdir)//'gs"'
            call write_fatal(1)
          end if
        end if
      end if

      if(finish) then
        if(verbosity_ >= VERB_COMPACT) then
          write(message(1), '(a, i4, a)') 'Info: SCF converged in ', iter, ' iterations'
          write(message(2), '(a)')        '' 
          call write_info(2)
        end if
        if(scf%lcao_restricted) call lcao_end(lcao)
        call profiling_out(prof)
        exit
      end if

      if(outp%duringscf .and. gs_run_) then
        write(dirname,'(a,i4.4)') "scf.",iter
        call h_sys_output_all(outp, gr, geo, st, hm, dirname)
      end if

      ! save information for the next iteration
      rhoin(1:gr%mesh%np, 1, 1:nspin) = st%rho(1:gr%mesh%np, 1:nspin)
      if (hm%d%cdft) rhoin(1:gr%mesh%np, 2:dim, 1:nspin) = st%j(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      if (scf%what2mix == MIXPOT) then
        vin(1:gr%mesh%np, 1, 1:nspin) = hm%vhxc(1:gr%mesh%np, 1:nspin)
        if (hm%d%cdft) vin(1:gr%mesh%np, 2:dim, 1:nspin) = hm%axc(1:gr%mesh%np, 1:gr%mesh%sb%dim, 1:nspin)
      end if
      evsum_in = evsum_out
      if (scf%conv_abs_force > M_ZERO) then
        forcein(1:geo%natoms, 1:gr%mesh%sb%dim) = forceout(1:geo%natoms, 1:gr%mesh%sb%dim)
      end if

      if(forced_finish) then
        call profiling_out(prof)
        exit
      end if

      ! check if debug mode should be enabled or disabled on the fly
      call io_debug_on_the_fly()


      call profiling_out(prof)
    end do

    if (scf%what2mix == MIXPOT) then
      call v_ks_calc(gr, ks, hm, st)
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
    call epot_forces(gr, geo, hm%ep, st)

    if (st%wfs_type == M_REAL) then
      call dstates_calc_momentum(gr, st)
    else
      call zstates_calc_momentum(gr, st)
    end if

    if(gs_run_) then 
      ! output final information
      call scf_write_static(STATIC_DIR, "info")
      call h_sys_output_all(outp, gr, geo, st, hm, STATIC_DIR)
    end if

    if(simul_box_is_periodic(gr%sb).and.st%d%nik > st%d%nspin) then
      call states_write_bands(STATIC_DIR, st%nst, st, gr%sb)
      call states_write_fermi_energy(STATIC_DIR, st, gr%mesh, gr%sb)
      call states_degeneracy_matrix(st)
    end if

    call pop_sub()

  contains


    ! ---------------------------------------------------------
    subroutine scf_write_iter
      character(len=50) :: str
      FLOAT :: mem
#ifdef HAVE_MPI
      FLOAT :: mem_tmp
#endif

      call push_sub('scf.scf_write_iter')

      if ( verbosity_ == VERB_FULL ) then

        write(str, '(a,i5)') 'SCF CYCLE ITER #' ,iter
        call messages_print_stress(stdout, trim(str))

        write(message(1),'(a,es15.8,2(a,es9.2))') ' etot = ', hm%etot/units_out%energy%factor, &
             ' abs_ev   = ', scf%abs_ev/units_out%energy%factor, ' rel_ev   = ', scf%rel_ev
        write(message(2),'(23x,2(a,es9.2))') &
             ' abs_dens = ', scf%abs_dens, ' rel_dens = ', scf%rel_dens
        ! write info about forces only if they are used as convergence criteria
        if (scf%conv_abs_force > M_ZERO) then
          write(message(3),'(23x,a,es9.2)') &
             ' force    = ', scf%abs_force/units_out%force%factor
          call write_info(3)
        else
          call write_info(2)
        end if

        if(.not.scf%lcao_restricted) then
          write(message(1),'(a,i6)') 'Matrix vector products: ', scf%eigens%matvec
          write(message(2),'(a,i6)') 'Converged eigenvectors: ', sum(scf%eigens%converged(1:st%d%nik))
          call write_info(2)
          call states_write_eigenvalues(stdout, st%nst, st, gr%sb, scf%eigens%diff)
        else
          call states_write_eigenvalues(stdout, st%nst, st, gr%sb)
        end if

        if(st%d%ispin > UNPOLARIZED) then
          call write_magnetic_moments(stdout, gr%mesh, st)
        end if

        write(message(1),'(a)') ''
        write(message(2),'(a,f14.2)') 'Elapsed time for SCF step :', etime
        call write_info(2)

        if(conf%report_memory) then
          mem = get_memory_usage()/(CNST(1024.0)**2)
#ifdef HAVE_MPI
          call MPI_Allreduce(mem, mem_tmp, 1, MPI_FLOAT, MPI_SUM, mpi_world%comm, mpi_err)
          mem = mem_tmp
#endif
          write(message(1),'(a,f14.2)') 'Memory usage [Mbytes]     :', mem
          call write_info(1)
        end if

        call messages_print_stress(stdout)
        
      end if

      if ( verbosity_ == VERB_COMPACT ) then
        ! write info about forces only if they are used as convergence criteria
        if (scf%conv_abs_force > M_ZERO) then
        write(message(1),'(a,i4,a,es15.8, 2(a,es9.2), a, f7.1, a)') &
             'iter ', iter, &
             ' : etot ', hm%etot/units_out%energy%factor, &
             ' : abs_dens', scf%abs_dens, &
             ' : force ', scf%abs_force/units_out%force%factor, &
             ' : etime ', etime, 's'
        else
        write(message(1),'(a,i4,a,es15.8, a,es9.2, a, f7.1, a)') &
             'iter ', iter, &
             ' : etot ', hm%etot/units_out%energy%factor, &
             ' : abs_dens', scf%abs_dens, &
             ' : etime ', etime, 's'
        end if
        call write_info(1)
      end if

      call pop_sub()
    end subroutine scf_write_iter


    ! ---------------------------------------------------------
    subroutine scf_write_static(dir, fname)
      character(len=*), intent(in) :: dir, fname

      FLOAT :: e_dip(4, st%d%nspin), n_dip(MAX_DIM)
      integer :: iunit, ispin, idir, iatom, nquantumpol

      call push_sub('scf.scf_write_static')

      if(mpi_grp_is_root(mpi_world)) then ! this the absolute master writes
        call io_mkdir(dir)
        iunit = io_open(trim(dir) // "/" // trim(fname), action='write')

        call grid_write_info(gr, iunit)

        if(simul_box_is_periodic(gr%sb)) then
          call kpoints_write_info(st%d, gr%mesh%sb%dim, iunit)
          write(iunit,'(1x)')
        end if

        call v_ks_write_info(ks, iunit)

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
      else
        iunit = 0
      end if

      call total_energy(hm, gr, st, iunit, full = .true.)

      if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
      if(st%d%ispin > UNPOLARIZED) then
        call write_magnetic_moments(iunit, gr%mesh, st)
        if(mpi_grp_is_root(mpi_world)) write(iunit, '(1x)')
      end if

      ! Next lines of code calculate the dipole of the molecule, summing the electronic and
      ! ionic contributions.

      do ispin = 1, st%d%nspin
        call dmf_multipoles(gr%mesh, st%rho(:, ispin), 1, e_dip(:, ispin))
      end do

      call geometry_dipole(geo, n_dip)

      do idir = 1, 3
        ! in periodic directions use single-point Berry`s phase calculation
        if(idir .le. gr%sb%periodic_dim) then
          n_dip(idir) = n_dip(idir) - epot_dipole_periodic(st, gr, idir)
          
          ! use quantum of polarization to reduce to smallest possible magnitude
          nquantumpol = FLOOR(n_dip(idir)/(2 * gr%sb%lsize(idir)))
          if(n_dip(idir) .lt. M_ZERO) then
             nquantumpol = nquantumpol + 1
             ! this makes that if n_dip = -1.1 R, it becomes -0.1 R, not 0.9 R
          endif

          n_dip(idir) = n_dip(idir) - nquantumpol * (2 * gr%sb%lsize(idir))

        ! in aperiodic directions use normal dipole formula
        else
          e_dip(idir + 1, 1) = sum(e_dip(idir + 1, :))
          n_dip(idir) = n_dip(idir) - e_dip(idir + 1, 1)         
        endif
      end do

      if(mpi_grp_is_root(mpi_world)) then
        call io_output_dipole(iunit, n_dip, gr%mesh%sb%dim)
        
        if (simul_box_is_periodic(gr%sb)) then
           write(iunit, '(a)') "Defined only up to quantum of polarization (e * lattice vector)."
           write(iunit, '(a)') "Single-point Berry's phase method only accurate for large supercells."

           if (st%d%nik * st%smear%el_per_state .ne. 2) then
              write(iunit, '(a)') &
                   "WARNING: Single-point Berry's phase method for dipole should not be used when there is more than one k-point."
           endif
        endif
        
        write(iunit, *)
      endif

      ! Output expectation values of the momentum operator
      call write_momentum(iunit)

      ! Next is the angular momentum. Only applies to 2D and 3D.
      if(gr%mesh%sb%dim == 2 .or. gr%mesh%sb%dim == 3) call write_angular_momentum(iunit)

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
        do iatom = 1, geo%natoms
          write(iunit,'(i4,a10,10f15.6)') iatom, trim(geo%atom(iatom)%spec%label), &
            geo%atom(iatom)%f(1:gr%mesh%sb%dim) / units_out%force%factor
        end do

        call io_close(iunit)
      end if

      call pop_sub()
    end subroutine scf_write_static


    ! ---------------------------------------------------------
    subroutine write_angular_momentum(iunit)
      integer, intent(in) :: iunit

      integer            :: ik, ist, ns, j
      character(len=80)  :: tmp_str(MAX_DIM), cspin
      FLOAT              :: angular(3), lsquare, o
      FLOAT, allocatable :: ang(:, :, :), ang2(:, :)
#if defined(HAVE_MPI)
      integer            :: tmp
      FLOAT, allocatable :: lang(:, :)
      integer            :: kstart, kend, kn
#endif

      call push_sub('scf.write_angular_momentum')

      ns = 1
      if(st%d%nspin == 2) ns = 2

      if(mpi_grp_is_root(mpi_world)) then
        write(iunit,'(a)') 'Angular Momentum of the KS states [dimensionless]:'
        if (st%d%nik > ns) then
          message(1) = 'Kpoints [' // trim(units_out%length%abbrev) // '^-1]'
          call write_info(1, iunit)
        end if
      end if

      ALLOCATE(ang (1:st%nst, st%d%nik, MAX_DIM), st%nst*st%d%nik*MAX_DIM)
      ALLOCATE(ang2(1:st%nst, st%d%nik), st%nst*st%d%nik)
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          if (st%wfs_type == M_REAL) then
            call dstates_angular_momentum(gr, st%dpsi(:, :, ist, ik), ang(ist, ik, :), ang2(ist, ik))
          else
            call zstates_angular_momentum(gr, st%zpsi(:, :, ist, ik), ang(ist, ik, :), ang2(ist, ik))
          end if
        end do
      end do

      angular(1) =  states_eigenvalues_sum(st, ang (st%st_start:st%st_end, :, 1))
      angular(2) =  states_eigenvalues_sum(st, ang (st%st_start:st%st_end, :, 2))
      angular(3) =  states_eigenvalues_sum(st, ang (st%st_start:st%st_end, :, 3))
      lsquare    =  states_eigenvalues_sum(st, ang2(st%st_start:st%st_end, :))

      do ik = 1, st%d%nik, ns
        if(st%d%nik > ns) then
          write(message(1), '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
            st%d%kpoints(1, ik)*units_out%length%factor, ',',            &
            st%d%kpoints(2, ik)*units_out%length%factor, ',',            &
            st%d%kpoints(3, ik)*units_out%length%factor, ')'
          call write_info(1, iunit)
        end if

        ! Exchange ang and ang2.
#if defined(HAVE_MPI)
        if(st%d%kpt%parallel) then
          kstart = st%d%kpt%start
          kend = st%d%kpt%end
          kn = st%d%kpt%nlocal

          ASSERT(.not. st%parallel_in_states)

          ALLOCATE(lang(1:st%lnst, 1:kn), st%lnst*kn)
          do j = 1, 3
            lang(1:st%lnst, 1:kn) = ang(st%st_start:st%st_end, kstart:kend, j)
            call MPI_Allgatherv(lang, st%nst*kn, MPI_FLOAT, &
                 ang(:, :, j), st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
                 st%d%kpt%mpi_grp%comm, mpi_err)
          end do
          lang(1:st%lnst, 1:kn) = ang2(st%st_start:st%st_end, kstart:kend)
            call MPI_Allgatherv(lang, st%nst*kn, MPI_FLOAT, &
                 ang2(:, :), st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
                 st%d%kpt%mpi_grp%comm, mpi_err)
          deallocate(lang)
        end if

        if(st%parallel_in_states) then
          ALLOCATE(lang(1:st%lnst, 1), st%lnst)
          do j = 1, 3
            lang(1:st%lnst, 1) = ang(st%st_start:st%st_end, ik, j)
            call lmpi_gen_allgatherv(st%lnst, lang(:, 1), tmp, ang(:, ik, j), st%mpi_grp)
          end do
          lang(1:st%lnst, 1) = ang2(st%st_start:st%st_end, ik)
          call lmpi_gen_allgatherv(st%lnst, lang(:, 1), tmp, ang2(:, ik), st%mpi_grp)
          deallocate(lang)
        end if
#endif
        write(message(1), '(a4,1x,a5,4a12,4x,a12,1x)')       &
          '#st',' Spin','        <Lx>', '        <Ly>', '        <Lz>', '        <L2>', 'Occupation '
        call write_info(1, iunit)

        if(mpi_grp_is_root(mpi_world)) then
          do j = 1, st%nst
            do is = 0, ns-1

              if(j > st%nst) then
                o = M_ZERO
              else
                o = st%occ(j, ik+is)
              end if

              if(is.eq.0) cspin = 'up'
              if(is.eq.1) cspin = 'dn'
              if(st%d%ispin.eq.UNPOLARIZED.or.st%d%ispin.eq.SPINORS) cspin = '--'

              write(tmp_str(1), '(i4,3x,a2)') j, trim(cspin)
              write(tmp_str(2), '(1x,4f12.6,3x,f12.6)') &
                ang(j, ik+is, 1:3), ang2(j, ik+is), o
              message(1) = trim(tmp_str(1))//trim(tmp_str(2))
              call write_info(1, iunit)
            end do
          end do
        end if
        write(message(1),'(a)') ''
        call write_info(1, iunit)

      end do

      write(message(1),'(a)') 'Total Angular Momentum L [dimensionless]'
      write(message(2),'(10x,4f12.6)') angular(1:3), lsquare
      call write_info(2, iunit)

      deallocate(ang, ang2)

      call pop_sub()
    end subroutine write_angular_momentum


    ! ---------------------------------------------------------
    subroutine write_momentum(iunit)
      integer,        intent(in) :: iunit

      integer           :: ik, j, is, ns, iunit2
      character(len=80) :: cspin
      FLOAT             :: o

      call push_sub('scf.write_momentum')   

      ns = 1
      if(st%d%nspin == 2) ns = 2

      write(message(1),'(a)') 'Momentum of the KS states [a.u.]:'
      call write_info(1, iunit)      
      if (st%d%nik > ns) then
        message(1) = 'Kpoints [' // trim(units_out%length%abbrev) // '^-1]'
        call write_info(1, iunit)
      end if

      do ik = 1, st%d%nik, ns
        if(st%d%nik > ns) then
          write(message(1), '(a,i4,3(a,f12.6),a)') '#k =',ik,', k = (',  &
            st%d%kpoints(1, ik)*units_out%length%factor, ',',            &
            st%d%kpoints(2, ik)*units_out%length%factor, ',',            &
            st%d%kpoints(3, ik)*units_out%length%factor, ')'
          call write_info(1, iunit)
        end if

        write(message(1), '(a4,1x,a5,3a12,4x,a12,1x)')       &
          '#st',' Spin','       <px>', '        <py>', '        <pz>', 'Occupation '
        call write_info(1, iunit)

        do j = 1, st%nst
          do is = 0, ns-1

            if(j > st%nst) then
              o = M_ZERO
            else
              o = st%occ(j, ik+is)
            end if
            
            if(is.eq.0) cspin = 'up'
            if(is.eq.1) cspin = 'dn'
            if(st%d%ispin.eq.UNPOLARIZED.or.st%d%ispin.eq.SPINORS) cspin = '--'

            write(message(1), '(i4,3x,a2,1x,3f12.6,3x,f12.6)')        &
                j, trim(cspin), st%momentum(1:min(gr%sb%dim, 3), j, ik), o
            call write_info(1, iunit)

          end do
        end do
        
        write(message(1),'(a)') ''
        call write_info(1, iunit)      

      end do

      ! also write to disk for further processing
      iunit2 = io_open(STATIC_DIR//'/momentum', action='write')

      do ik = 1, st%d%nik
        do j = 1, st%nst
          write(message(1), '(2i4,3x, 3f12.6,3x,f12.6)') &
            j, ik, st%momentum(1:min(gr%sb%dim, 3), j, ik), st%occ(j, ik)
          call write_info(1, iunit2)
        end do
      end do

      call io_close(iunit2)

      call pop_sub()
    end subroutine write_momentum


    ! ---------------------------------------------------------
    subroutine write_magnetic_moments(iunit, m, st)
      integer,        intent(in) :: iunit
      type(mesh_t),   intent(in) :: m
      type(states_t), intent(in) :: st

      integer :: i
      FLOAT :: mm(MAX_DIM)
      FLOAT, allocatable :: lmm(:,:)

      call push_sub('scf.write_magnetic_moments')

      call magnetic_moment(m, st, st%rho, mm)
      ALLOCATE(lmm(MAX_DIM, geo%natoms), MAX_DIM*geo%natoms)
      call magnetic_local_moments(m, st, geo, st%rho, scf%lmm_r, lmm)

      if(mpi_grp_is_root(mpi_world)) then

        write(iunit, '(a)') 'Total Magnetic Moment:'
        if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
          write(iunit, '(a,f10.6)') ' mz = ', mm(3)
        else if(st%d%ispin == SPINORS) then ! non-collinear
          write(iunit, '(1x,3(a,f10.6,3x))') 'mx = ', mm(1),'my = ', mm(2),'mz = ', mm(3)
        end if

        write(iunit, '(a,a,a,f7.3,a)') 'Local Magnetic Moments (sphere radius [', &
             trim(units_out%length%abbrev),'] = ', scf%lmm_r/units_out%length%factor, '):'
        if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
          write(iunit,'(a,6x,14x,a)') ' Ion','mz'
          do i = 1, geo%natoms
            write(iunit,'(i4,a10,f15.6)') i, trim(geo%atom(i)%spec%label), lmm(3, i)
          end do
        else if(st%d%ispin == SPINORS) then ! non-collinear
          write(iunit,'(a,8x,13x,a,13x,a,13x,a)') ' Ion','mx','my','mz'
          do i = 1, geo%natoms
            write(iunit,'(i4,a10,9f15.6)') i, trim(geo%atom(i)%spec%label), lmm(1:m%sb%dim, i)
          end do
        end if

      end if
      
      deallocate(lmm)

      call pop_sub()
    end subroutine write_magnetic_moments

  end subroutine scf_run

end module scf_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
