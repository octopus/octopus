!!  Copyright (C) 2012-2019 I. Theophilou, N. Helbig
!!  Copyright (C) 2019 F. Buchholz, M. Oliveira
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

#include "global.h"

module rdmft_oct_m
  use density_oct_m
  use derivatives_oct_m
  use eigen_cg_oct_m
  use eigensolver_oct_m
  use energy_oct_m
  use exchange_operator_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use hamiltonian_elec_base_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use lalg_adv_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use minimizer_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use output_oct_m
  use parser_oct_m
  use photon_mode_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use space_oct_m
  use species_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_restart_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use v_ks_oct_m
  use xc_oep_oct_m
 
  implicit none

  private
  public ::                   &
    rdm_t,                    &
    rdmft_init,               &
    rdmft_end,                &
    scf_rdmft

  type rdm_t
    private
    type(eigensolver_t) :: eigens
    integer  :: max_iter !< maximum number of scf iterations
    integer  :: iter
    integer  :: nst !< number of states
    integer  :: n_twoint !number of unique two electron integrals
    logical  :: do_basis
    logical  :: hf
    FLOAT    :: mu
    FLOAT    :: occsum
    FLOAT    :: qtot
    FLOAT    :: scale_f
    FLOAT    :: toler
    FLOAT    :: conv_ener
    FLOAT    :: maxFO
    FLOAT    :: tolerFO

    FLOAT, allocatable   :: eone(:)
    FLOAT, allocatable   :: eone_int(:,:)
    FLOAT, allocatable   :: twoint(:)
    FLOAT, allocatable   :: hartree(:,:)
    FLOAT, allocatable   :: exchange(:,:)
    FLOAT, allocatable   :: evalues(:)
    FLOAT, allocatable   :: vecnat(:,:)
    FLOAT, allocatable   :: Coul(:,:,:)
    FLOAT, allocatable   :: Exch(:,:,:)

    integer, allocatable :: i_index(:,:)
    integer, allocatable :: j_index(:,:)
    integer, allocatable :: k_index(:,:)
    integer, allocatable :: l_index(:,:)
  end type rdm_t

  type(rdm_t), pointer :: rdm_ptr
  
contains

  ! ---------------------------------------------------------
  subroutine rdmft_init(rdm, namespace, gr, st, mc, space, fromScratch)
    type(rdm_t),         intent(out)   :: rdm
    type(namespace_t),   intent(in)    :: namespace
    type(grid_t),        intent(inout) :: gr  !< grid
    type(states_elec_t), intent(in)    :: st  !< States
    type(multicomm_t),   intent(in)    :: mc
    type(space_t),       intent(in)    :: space
    logical,             intent(in)    :: fromScratch

    integer :: ist

    PUSH_SUB(rdmft_init)

    if(st%nst < st%qtot + 1) then   
      message(1) = "Too few states to run RDMFT calculation"
      message(2) = "Number of states should be at least the number of electrons plus one"
      call messages_fatal(2)
    end if
   
    if (states_are_complex(st)) then
      call messages_not_implemented("Complex states for RDMFT")
    end if

    ! The documentation for the variable is found in scf_init.
    call parse_variable(namespace, 'MaximumIter', 200, rdm%max_iter)

    !%Variable RDMTolerance
    !%Type float
    !%Default 1e-7 Ha
    !%Section SCF::RDMFT
    !%Description
    !% Convergence criterion for stopping the occupation numbers minimization. Minimization is
    !% stopped when all derivatives of the energy wrt. each occupation number 
    !% are smaller than this criterion. The bisection for finding the correct mu that is needed
    !% for the occupation number minimization also stops according to this criterion.
    !%End
    call parse_variable(namespace, 'RDMTolerance', CNST(1.0e-7), rdm%toler)

    !%Variable RDMToleranceFO
    !%Type float
    !%Default 1e-4 Ha
    !%Section SCF::RDMFT
    !%Description
    !% Convergence criterion for stopping the diagonalization of the Fock matrix in the Piris method.
    !% Orbital minimization is stopped when all off-diagonal ellements of the Fock matrix 
    !% are smaller than this criterion.
    !%End
    call parse_variable(namespace, 'RDMToleranceFO', CNST(1.0e-4), rdm%tolerFO)

    !%Variable RDMConvEner
    !%Type float
    !%Default 1e-6 Ha
    !%Section SCF::RDMFT
    !%Description
    !% Convergence criterion for stopping the overall minimization of the energy with
    !% respect to occupation numbers and the orbitals. The minimization of the 
    !% energy stops when the total energy difference between two subsequent 
    !% minimizations of the energy with respect to the occupation numbers and the
    !% orbitals is smaller than this criterion. It is also used to exit the orbital minimization.
    !%End
    call parse_variable(namespace, 'RDMConvEner', CNST(1.0e-7), rdm%conv_ener)

    !%Variable RDMBasis
    !%Type logical
    !%Default yes 
    !%Section SCF::RDMFT
    !%Description
    !% If true, all the energy terms and corresponding derivatives involved in RDMFT will
    !% not be calculated on the grid but on the basis of the initial orbitals
    !%End
    call parse_variable(namespace, 'RDMBasis',.true., rdm%do_basis)
    
    if (rdm%do_basis .and. fromScratch) then
      call messages_write("RDMFT calculations with RDMBasis = yes cannot be started FromScratch", new_line=.true.)
      call messages_write("Run a calculation for independent particles first")
      call messages_fatal()
    end if

    !%Variable RDMHartreeFock
    !%Type logical
    !%Default no
    !%Section SCF::RDMFT
    !%Description
    !% If true, the code simulates a HF calculation, by omitting the occ.num. optimization
    !% can be used for test reasons
    !%End
    call parse_variable(namespace, 'RDMHartreeFock',.false., rdm%hf)

    rdm%nst = st%nst
    if (rdm%do_basis) then
      rdm%n_twoint = rdm%nst*(rdm%nst + 1)*(rdm%nst**2 + rdm%nst + 2)/8
      SAFE_ALLOCATE(rdm%eone_int(1:rdm%nst, 1:rdm%nst))
      SAFE_ALLOCATE(rdm%twoint(1:rdm%n_twoint))
      SAFE_ALLOCATE(rdm%i_index(1:2,1:rdm%n_twoint))
      SAFE_ALLOCATE(rdm%j_index(1:2,1:rdm%n_twoint))
      SAFE_ALLOCATE(rdm%k_index(1:2,1:rdm%n_twoint))
      SAFE_ALLOCATE(rdm%l_index(1:2,1:rdm%n_twoint))
      SAFE_ALLOCATE(rdm%vecnat(1:rdm%nst, 1:rdm%nst))
      SAFE_ALLOCATE(rdm%Coul(1:rdm%nst, 1:rdm%nst, 1:rdm%nst))
      SAFE_ALLOCATE(rdm%Exch(1:rdm%nst, 1:rdm%nst, 1:rdm%nst))
      rdm%eone_int = M_ZERO
      rdm%twoint = M_ZERO
      rdm%vecnat = M_ZERO
      rdm%i_index = M_ZERO
      rdm%j_index = M_ZERO
      rdm%k_index = M_ZERO
      rdm%l_index = M_ZERO
      rdm%Coul = M_ZERO
      rdm%Exch = M_ZERO
      do ist = 1, rdm%nst
        rdm%vecnat(ist, ist) = M_ONE
      end do
    else
      ! initialize eigensolver. 
      call eigensolver_init(rdm%eigens, namespace, gr, st, mc, space)
      if (rdm%eigens%additional_terms) call messages_not_implemented("CG Additional Terms with RDMFT.")
    end if

    SAFE_ALLOCATE(rdm%eone(1:rdm%nst))
    SAFE_ALLOCATE(rdm%hartree(1:rdm%nst, 1:rdm%nst))
    SAFE_ALLOCATE(rdm%exchange(1:rdm%nst, 1:rdm%nst))
    SAFE_ALLOCATE(rdm%evalues(1:rdm%nst))

    rdm%eone = M_ZERO
    rdm%hartree = M_ZERO
    rdm%exchange = M_ZERO
    rdm%evalues = M_ZERO
    rdm%mu = M_TWO*st%eigenval(int(st%qtot*M_HALF), 1)
    rdm%qtot = st%qtot
    rdm%occsum = M_ZERO
    rdm%scale_f = CNST(1e-2)
    rdm%maxFO = M_ZERO
    rdm%iter = 0

    POP_SUB(rdmft_init)
  end subroutine rdmft_init

  ! ----------------------------------------

  subroutine rdmft_end(rdm)
    type(rdm_t),  intent(inout) :: rdm

    PUSH_SUB(rdmft_end)

    SAFE_DEALLOCATE_A(rdm%evalues)
    SAFE_DEALLOCATE_A(rdm%eone)
    SAFE_DEALLOCATE_A(rdm%hartree)
    SAFE_DEALLOCATE_A(rdm%exchange)

    if (rdm%do_basis) then
      SAFE_DEALLOCATE_A(rdm%eone_int)
      SAFE_DEALLOCATE_A(rdm%twoint)
      SAFE_DEALLOCATE_A(rdm%i_index)
      SAFE_DEALLOCATE_A(rdm%j_index)
      SAFE_DEALLOCATE_A(rdm%k_index)
      SAFE_DEALLOCATE_A(rdm%l_index)
      SAFE_DEALLOCATE_A(rdm%vecnat)
      SAFE_DEALLOCATE_A(rdm%Coul)
      SAFE_DEALLOCATE_A(rdm%Exch)
    else
      call eigensolver_end(rdm%eigens)
    end if

    POP_SUB(rdmft_end)
  end subroutine rdmft_end

  ! ----------------------------------------

  ! scf for the occupation numbers and the natural orbitals
  subroutine scf_rdmft(rdm, namespace, space, gr, ions, st, ks, hm, outp, restart_dump)
    type(rdm_t),              intent(inout) :: rdm
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr  !< grid
    type(ions_t),             intent(in)    :: ions !< geometry
    type(states_elec_t),      intent(inout) :: st  !< States
    type(v_ks_t),             intent(inout) :: ks  !< Kohn-Sham
    type(hamiltonian_elec_t), intent(inout) :: hm  !< Hamiltonian
    type(output_t),           intent(in)    :: outp !< output
    type(restart_t),          intent(in)    :: restart_dump

    type(states_elec_t) :: states_save
    integer :: iter, icount, ip, ist, ierr, maxcount, iorb
    integer(8) :: what_i
    FLOAT :: energy, energy_dif, energy_old, energy_occ, xpos, xneg, rel_ener
    FLOAT, allocatable :: dpsi(:,:), dpsi2(:,:)
    logical :: conv
    character(len=MAX_PATH_LEN) :: dirname    

    PUSH_SUB(scf_rdmft)

    if (hm%d%ispin /= 1) then
      call messages_not_implemented("RDMFT exchange function not yet implemented for spin_polarized or spinors")
    end if

    ! problem is about k-points for exchange
    if (space%is_periodic()) then
      call messages_not_implemented("Periodic system calculations for RDMFT", namespace=namespace)
    end if

    ! exchange routine needs all states on each processor currently
    if(st%parallel_in_states) then
      call messages_not_implemented("RDMFT parallel in states", namespace=namespace)
    end if

    call messages_print_stress(stdout, 'RDMFT Calculation', namespace=namespace)
    call messages_print_var_value(stdout, 'RDMBasis', rdm%do_basis)

    !set initial values
    energy_old = CNST(1.0e20)
    xpos = M_ZERO 
    xneg = M_ZERO
    energy = M_ZERO 
    if (.not. rdm%do_basis) then
      maxcount = 1 !still needs to be checked
    else
      maxcount = 50
      !precalculate matrix elements in basis
      write(message(1),'(a)') 'Calculating Coulomb and exchange matrix elements in basis'
      write(message(2),'(a)') '--this may take a while--'
      call messages_info(2)

      call dstates_elec_me_two_body(st, namespace, space, gr, hm%kpoints, hm%exxop%psolver, 1, &
                                      st%nst, rdm%i_index, rdm%j_index, rdm%k_index, &
        rdm%l_index, rdm%twoint)
      call rdm_integrals(rdm, namespace, hm, st, gr%mesh)
      call sum_integrals(rdm)
    endif

    ! Start the actual minimization, first step is minimization of occupation numbers
    ! Orbital minimization is according to Piris and Ugalde, Vol. 30, No. 13, J. Comput. Chem. (scf_orb) or
    ! using conjugated gradient (scf_orb_cg)
    do iter = 1, rdm%max_iter
      rdm%iter = rdm%iter + 1
      write(message(1), '(a)') '**********************************************************************'
      write(message(2),'(a, i4)') 'Iteration:', iter
      call messages_info(2)
      ! occupation number optimization unless we are doing Hartree-Fock
      if (rdm%hf) then
        call scf_occ_NO(rdm, namespace, gr, hm, space, st, energy_occ)
      else
        call scf_occ(rdm, namespace, gr, hm, space, st, energy_occ)
      end if
      ! orbital optimization
      write(message(1), '(a)') 'Optimization of natural orbitals'
      call messages_info(1)
      do icount = 1, maxcount 
        if (rdm%do_basis) then
          call scf_orb(rdm, namespace, gr, st, hm, space, energy)
        else
          call scf_orb_cg(rdm, namespace, space, gr, ions, st, ks, hm, energy)
        end if
        energy_dif = energy - energy_old
        energy_old = energy
        if (rdm%do_basis) then
          if (abs(energy_dif)/abs(energy) < rdm%conv_ener .and. rdm%maxFO < rdm%tolerFO)  exit
          if (energy_dif < M_ZERO) then 
            xneg = xneg + 1
          else
            xpos = xpos + 1
          end if
          if (xneg > CNST(1.5e0)*xpos) then
            rdm%scale_f = CNST(1.01)*rdm%scale_f
          elseif (xneg < CNST(1.1e0)*xpos) then 
            rdm%scale_f = CNST(0.95)* rdm%scale_f 
          end if
        endif !rdm%do_basis
      end do !icount
      xneg = M_ZERO
      xpos = M_ZERO
      
      rel_ener = abs(energy_occ-energy)/abs(energy)

      write(message(1),'(a,11x,es20.10)') 'Total energy:', units_from_atomic(units_out%energy,energy + hm%ep%eii) 
      write(message(2),'(a,1x,es20.10)') 'Rel. energy difference:', rel_ener
      call messages_info(2)

      if (.not. rdm%hf .and. rdm%do_basis) then
        write(message(1),'(a,18x,es20.10)') 'Max F0:', rdm%maxFO
        call messages_info(1)
      end if


      if (rdm%do_basis) then
        conv = (rel_ener < rdm%conv_ener) .and. rdm%maxFO < rdm%tolerFO
      else
        conv = rel_ener < rdm%conv_ener
      endif
 
      if (rdm%toler > CNST(1e-4)) rdm%toler = rdm%toler*CNST(1e-1) !Is this still okay or does it restrict the possible convergence? FB: Does this makes sense at all?

      ! save restart information
      if ((conv .or. (modulo(iter, outp%restart_write_interval) == 0) .or. iter == rdm%max_iter)) then
        if (rdm%do_basis) then
          call states_elec_copy(states_save, st)
          SAFE_ALLOCATE(dpsi(1:gr%mesh%np, 1:st%d%dim))
          SAFE_ALLOCATE(dpsi2(1:gr%mesh%np, 1:st%d%dim))
          do iorb = 1, st%nst
            dpsi = M_ZERO
            do ist = 1, st%nst
              call states_elec_get_state(st, gr%mesh, ist, 1, dpsi2)
              do ip = 1, gr%mesh%np
                dpsi(ip,1) = dpsi(ip,1) + rdm%vecnat(ist, iorb)*dpsi2(ip,1)
              end do
            end do
            call states_elec_set_state(states_save, gr%mesh, iorb, 1, dpsi)
          end do
          call density_calc(states_save, gr, states_save%rho)
          ! if other quantities besides the densities and the states are needed they also have to be recalculated here!
          call states_elec_dump(restart_dump, space, states_save, gr%mesh, hm%kpoints, ierr, iter=iter) 

          if (conv .or. iter == rdm%max_iter) then
            call states_elec_end(st)
            call states_elec_copy(st, states_save)
          end if
        
          call states_elec_end(states_save)
      
          SAFE_DEALLOCATE_A(dpsi)
          SAFE_DEALLOCATE_A(dpsi2)
        else
          call states_elec_dump(restart_dump, space, st, gr%mesh, hm%kpoints, ierr, iter=iter) 
          
          ! calculate maxFO for cg-solver 
          if (.not. rdm%hf) then
            call calc_maxFO (namespace, hm, st, gr, rdm)
            write(message(1),'(a,18x,es20.10)') 'Max F0:', rdm%maxFO
            call messages_info(1)
          end if 
        endif

        if (ierr /= 0) then
          message(1) = 'Unable to write states wavefunctions.'
          call messages_warning(1)
        end if
        
      endif

      ! write output for iterations if requested
      if (any(outp%what) .and. outp%duringscf) then
        do what_i = lbound(outp%what, 1), ubound(outp%what, 1)
          if (outp%what_now(what_i, iter)) then
            write(dirname,'(a,a,i4.4)') trim(outp%iter_dir), "scf.", iter
            call output_all(outp, namespace, space, dirname, gr, ions, iter, st, hm, ks)
            call scf_write_static(dirname, "info")
            exit
          end if
        end do
      end if

      if (conv) exit
    end do 

    if(conv) then 
      write(message(1),'(a,i3,a)')  'The calculation converged after ',rdm%iter,' iterations'
      write(message(2),'(a,9x,es20.10)')  'The total energy is ', units_from_atomic(units_out%energy,energy + hm%ep%eii)
      call messages_info(2)
    else
      write(message(1),'(a,i3,a)')  'The calculation did not converge after ', iter-1, ' iterations '
      write(message(2),'(a,es15.5)') 'Relative energy difference between the last two iterations ', rel_ener
      write(message(3),'(a,es15.5)') 'The maximal non-diagonal element of the Hermitian matrix F is ', rdm%maxFO
      call messages_info(3)
    end if

    call scf_write_static(STATIC_DIR, "info")
    call output_all(outp, namespace, space, STATIC_DIR, gr, ions, -1, st, hm, ks)

    POP_SUB(scf_rdmft) 

  contains
    ! ---------------------------------------------------------
    subroutine scf_write_static(dir, fname)
      character(len=*), intent(in) :: dir, fname

      integer :: iunit, ist
      FLOAT, allocatable :: photon_number_state (:), ekin_state (:), epot_state (:)

      PUSH_SUB(scf_rdmft.scf_write_static)

      SAFE_ALLOCATE(photon_number_state(1:st%nst))
      SAFE_ALLOCATE(ekin_state(1:st%nst))
      SAFE_ALLOCATE(epot_state(1:st%nst))

      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir(dir, namespace)
        iunit = io_open(trim(dir) // "/" // trim(fname), namespace, action='write')

        call grid_write_info(gr, iunit)

        call v_ks_write_info(ks, iunit, namespace)
        
        if (rdm%do_basis) then
          write(iunit, '(a)')'Orbital optimization with [basis set]'
        else
          write(iunit, '(a)')'Orbital optimization with [conjugated gradients]'
        end if
        write(iunit, '(1x)')
        
        if (rdm%hf) then
          write(iunit, '(a)')'Hartree Fock calculation'
          write(iunit, '(1x)')
        end if
        
        if (hm%psolver%is_dressed) then
          write(iunit, '(a)')'Dressed state calculation'
          call photon_mode_write_info(hm%psolver%photons, iunit)
          write(iunit, '(1x)')
        end if

        ! scf information
        if(conv) then
          write(iunit, '(a, i4, a)')'SCF converged in ', iter, ' iterations'
        else
          write(iunit, '(a)') 'SCF *not* converged!'
        end if
        write(iunit, '(1x)')

        write(iunit, '(3a,es20.10)') 'Total Energy [', trim(units_abbrev(units_out%energy)), ']:', &
          units_from_atomic(units_out%energy, energy + hm%ep%eii) 
        write(iunit,'(a,1x,f16.12)') 'Sum of occupation numbers:', rdm%occsum
      else
        iunit = 0
      end if

      if (hm%psolver%is_dressed) then
        call calc_photon_number(gr, st, hm%psolver%photons, photon_number_state, ekin_state, epot_state)
        if(mpi_grp_is_root(mpi_world)) then
          write(iunit,'(a,1x,f14.12)') 'Total mode occupation:', hm%psolver%photons%number(1)
        end if
      end if

      if(mpi_grp_is_root(mpi_world)) then
        if (rdm%max_iter > 0) then
          write(iunit, '(a)') 'Convergence:'
          write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'maxFO = ', rdm%maxFO
          write(iunit, '(6x, a, es15.8,a,es15.8,a)') 'rel_ener = ', rel_ener
          write(iunit,'(1x)')
        end if
        ! otherwise, these values are uninitialized, and unknown.
      end if

      if (mpi_grp_is_root(mpi_world)) then
        ! Write header
        write(iunit,'(a)') 'Natural occupation numbers:'
        write(iunit,'(a4,5x,a12)', advance='no') '#st', 'Occupation'
        if (.not. rdm%do_basis) write(iunit,'(5x,a12)', advance='no') 'conv'
        if (hm%psolver%is_dressed) write(iunit,'(3(5x,a12))', advance='no') 'Mode Occ.', '-1/2d^2/dq^2', '1/2w^2q^2'
        write(iunit,*)

        ! Write values
        do ist = 1, st%nst
          write(iunit,'(i4,3x,f14.12)', advance='no') ist, st%occ(ist, 1)
          if (.not. rdm%do_basis) write(iunit,'(3x,f14.12)', advance='no') rdm%eigens%diff(ist, 1)
          if (hm%psolver%is_dressed) then
            write(iunit,'(3(3x,f14.12))', advance='no') photon_number_state(ist), ekin_state(ist), epot_state(ist)
          end if
          write(iunit,*)
        end do
      end if
      
      if (mpi_grp_is_root(mpi_world)) then
        call io_close(iunit)
      end if
      
      SAFE_DEALLOCATE_A(photon_number_state)
      SAFE_DEALLOCATE_A(ekin_state)
      SAFE_DEALLOCATE_A(epot_state)

      POP_SUB(scf_rdmft.scf_write_static)
    end subroutine scf_write_static 
  end subroutine scf_rdmft
    
  ! ---------------------------------------------------------
  subroutine calc_maxFO (namespace, hm, st, gr, rdm)
    type(namespace_t),         intent(in)    :: namespace
    type(rdm_t),               intent(inout) :: rdm
    type(grid_t),              intent(in)    :: gr
    type(hamiltonian_elec_t),  intent(inout) :: hm
    type(states_elec_t),       intent(inout) :: st

    FLOAT, allocatable ::  lambda(:,:), FO(:,:)
    integer :: ist, jst

    PUSH_SUB(calc_maxFO)

    SAFE_ALLOCATE(lambda(1:st%nst,1:st%nst))
    SAFE_ALLOCATE(FO(1:st%nst, 1:st%nst))

    ! calculate FO operator to check Hermiticity of lagrange multiplier matrix (lambda)
    lambda = M_ZERO
    FO = M_ZERO
    call construct_lambda(namespace, hm, st, gr, lambda, rdm)

    !Set up FO matrix to check maxFO
    do ist = 1, st%nst
      do jst = 1, ist - 1
        FO(jst, ist) = - (lambda(jst, ist) - lambda(ist ,jst))
      end do
    end do
    rdm%maxFO = maxval(abs(FO))

    SAFE_DEALLOCATE_A(lambda)
    SAFE_DEALLOCATE_A(FO)

    POP_SUB(calc_maxFO)
  end subroutine calc_maxFO

  ! ---------------------------------------------------------
  subroutine calc_photon_number(gr, st, photons, photon_number_state, ekin_state, epot_state)
    type(grid_t),                intent(in)    :: gr
    type(states_elec_t),         intent(in)    :: st
    type(photon_mode_t),         intent(inout) :: photons
    FLOAT,                       intent(out)   :: photon_number_state(:)
    FLOAT,                       intent(out)   :: ekin_state(:)
    FLOAT,                       intent(out)   :: epot_state(:)

    integer :: ist, dim_photon
    FLOAT   :: q2_exp, laplace_exp
    FLOAT, allocatable :: psi(:, :), psi_q2(:), dpsidq(:), d2psidq2(:)

    PUSH_SUB(calc_photon_number)

    ! The photon dimension is always the last
    dim_photon = gr%sb%dim

    SAFE_ALLOCATE(psi(1:gr%mesh%np_part, 1))
    SAFE_ALLOCATE(psi_q2(1:gr%mesh%np))
    SAFE_ALLOCATE(dpsidq(1:gr%mesh%np_part))
    SAFE_ALLOCATE(d2psidq2(1:gr%mesh%np))

    photons%number(1) = M_ZERO

    do ist = 1, st%nst
      call states_elec_get_state(st, gr%mesh, ist, 1, psi)

      ! <phi(ist)|d^2/dq^2|phi(ist)> ~= <phi(ist)| d/dq (d/dq|phi(ist)>)
      call dderivatives_partial(gr%der, psi(:, 1), dpsidq(:), dim_photon, ghost_update = .true., set_bc = .true.)
      call dderivatives_partial(gr%der, dpsidq(1:gr%mesh%np_part), d2psidq2(:), dim_photon, ghost_update = .true., set_bc = .true.)
      laplace_exp = dmf_dotp(gr%mesh, psi(:, 1), d2psidq2(:))
      ekin_state(ist) = -M_HALF*laplace_exp

      ! <phi(ist)|q^2|psi(ist)>= |q|psi(ist)>|^2
      psi_q2(1:gr%mesh%np) = psi(1:gr%mesh%np, 1) * gr%mesh%x(1:gr%mesh%np, dim_photon)**2
      q2_exp = dmf_dotp(gr%mesh, psi(:, 1), psi_q2(:))
      epot_state(ist) = M_HALF * photons%omega(1)**2 * q2_exp

      !! N_phot(ist)=( <phi_i|H_ph|phi_i>/omega - 0.5 ) / N_elec
      !! with <phi_i|H_ph|phi_i>=-0.5* <phi(ist)|d^2/dq^2|phi(ist)> + 0.5*omega <phi(ist)|q^2|psi(ist)>
      photon_number_state(ist) = -M_HALF*laplace_exp / photons%omega(1) + M_HALF * photons%omega(1) * q2_exp
      photon_number_state(ist) = photon_number_state(ist) - M_HALF

      !! N_phot_total= sum_ist occ_ist*N_phot(ist)
      photons%number(1) = photons%number(1) + (photon_number_state(ist) + M_HALF)*st%occ(ist, 1) ! 0.5 must be added again to do the normalization due to the total charge correctly
    end do

    photons%number(1) =  photons%number(1) - st%qtot/M_TWO

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(psi_q2)
    SAFE_DEALLOCATE_A(dpsidq)
    SAFE_DEALLOCATE_A(d2psidq2)

    POP_SUB(calc_photon_number)
  end subroutine calc_photon_number

  ! ---------------------------------------------------------
  
  ! reset occ.num. to 2/0
  subroutine set_occ_pinning(st)
    type(states_elec_t), intent(inout) :: st

    FLOAT, allocatable ::  occin(:,:)

    PUSH_SUB(set_occ_pinning)    

    SAFE_ALLOCATE(occin(1:st%nst, 1:st%d%nik))
    
    occin(1:st%nst, 1:st%d%nik) = st%occ(1:st%nst, 1:st%d%nik)
    where(occin(:,:) < M_ONE) occin(:,:) = M_ZERO
    where(occin(:,:) > M_ONE) occin(:,:) = st%smear%el_per_state
    
    st%occ = occin

    SAFE_DEALLOCATE_A(occin)
    
    POP_SUB(set_occ_pinning)
  end subroutine set_occ_pinning


  ! ---------------------------------------------------------
  ! dummy routine for occupation numbers which only calculates the necessary variables for further use
  ! used in Hartree-Fock mode
  subroutine scf_occ_NO(rdm, namespace, gr, hm, space, st, energy)
    type(rdm_t),              intent(inout) :: rdm
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(space_t),            intent(in)    :: space
    type(states_elec_t),      intent(inout) :: st
    FLOAT,                    intent(out)   :: energy
     
    integer :: ist

    PUSH_SUB(scf_occ_NO)

    write(message(1),'(a)') 'SKIP Optimization of occupation numbers'
    call messages_info(1)
    
    call set_occ_pinning(st)
    
    energy = M_ZERO

    call rdm_derivatives(rdm, namespace, hm, st, gr, space)

    call total_energy_rdm(rdm, st%occ(:,1), energy)

    rdm%occsum = sum(st%occ(1:st%nst, 1:st%d%nik))

    write(message(1),'(a4,5x,a12)')'#st','Occupation'
    call messages_info(1)   

    do ist = 1, st%nst
      write(message(1),'(i4,3x,f11.9)') ist, st%occ(ist, 1)
      call messages_info(1)
    end do

    write(message(1),'(a,1x,f13.9)') 'Sum of occupation numbers', rdm%occsum
    write(message(2),'(a,es20.10)') 'Total energy occ', units_from_atomic(units_out%energy,energy + hm%ep%eii) 
    call messages_info(2)   
    
    POP_SUB(scf_occ_NO)
  end subroutine scf_occ_NO

  ! scf for the occupation numbers 
  subroutine scf_occ(rdm, namespace, gr, hm, space, st, energy)
    type(rdm_t), target,      intent(inout) :: rdm
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(in)    :: gr
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(space_t),            intent(in)    :: space
    type(states_elec_t),      intent(inout) :: st
    FLOAT,                    intent(out)   :: energy

    integer :: ist, icycle, ierr
    FLOAT ::  sumgi1, sumgi2, sumgim, mu1, mu2, mum, dinterv, thresh_occ
    FLOAT, allocatable ::  occin(:,:)
    FLOAT, parameter :: smallocc = CNST(0.00001) 
    FLOAT, allocatable ::   theta(:)
    FLOAT :: objective
    type(profile_t), save :: prof_occ

    PUSH_SUB(scf_occ)
    call profiling_in(prof_occ, "SCF_OCC")

    write(message(1),'(a)') 'Optimization of occupation numbers'
    call messages_info(1)
    
    SAFE_ALLOCATE(occin(1:st%nst, 1:st%d%nik))
    SAFE_ALLOCATE(theta(1:st%nst))

    occin = M_ZERO
    theta  = M_ZERO
    energy = M_ZERO
    
    ! Defines a threshold on occ nums to avoid numerical instabilities.
    ! Needs to be changed consistently with the same variable in objective_rdmft
    thresh_occ = CNST(1e-14)
    
    !Initialize the occin. Smallocc is used for numerical stability
    occin(1:st%nst, 1:st%d%nik) = st%occ(1:st%nst, 1:st%d%nik)
    where(occin(:,:) < smallocc) occin(:,:) = smallocc
    where(occin(:,:) > st%smear%el_per_state - smallocc) occin(:,:) = st%smear%el_per_state - smallocc

    !Renormalize the occupation numbers 
    rdm%occsum = st%qtot

    st%occ = occin
    
    call rdm_derivatives(rdm, namespace, hm, st, gr, space)

    !finding the chemical potential mu such that the occupation numbers sum up to the number of electrons
    !bisection to find the root of rdm%occsum-st%qtot=M_ZERO
    mu1 = rdm%mu   !initial guess for mu 
    mu2 = -CNST(1.0e-6)
    dinterv = M_HALF

    ! Set pointer to rdm, so that it is available in the functions called by the minimizer
    rdm_ptr => rdm

    !use n_j=sin^2(2pi*theta_j) to treat pinned states, minimize for both intial mu
    theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
    call minimize_multidim(MINMETHOD_BFGS, st%nst, theta, CNST(0.05), CNST(0.01), &
        CNST(1e-12), CNST(1e-12), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
    sumgi1 = rdm%occsum - st%qtot
    rdm%mu = mu2
    theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
    call minimize_multidim(MINMETHOD_BFGS, st%nst, theta, CNST(0.05), CNST(0.01), &
        CNST(1e-12), CNST(1e-12), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
    sumgi2 = rdm%occsum - st%qtot

    ! Adjust the interval between the initial mu to include the root of rdm%occsum-st%qtot=M_ZERO
    do while (sumgi1*sumgi2 > M_ZERO) 
      if (sumgi2 > M_ZERO) then
        mu2 = mu1
        sumgi2 = sumgi1
        mu1 = mu1 - dinterv
        rdm%mu = mu1
        theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
        call minimize_multidim(MINMETHOD_BFGS, st%nst, theta, CNST(0.05), CNST(0.01), &
            CNST(1e-12), CNST(1e-12), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
        sumgi1 = rdm%occsum - st%qtot 
      else
        mu1 = mu2
        sumgi1 = sumgi2
        mu2 = mu2 + dinterv
        rdm%mu = mu2
        theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
        call minimize_multidim(MINMETHOD_BFGS, st%nst, theta, CNST(0.05), CNST(0.01), &
            CNST(1e-12), CNST(1e-12), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
        sumgi2 = rdm%occsum - st%qtot 
      end if
    end do

    do icycle = 1, 50
      mum = (mu1 + mu2)*M_HALF
      rdm%mu = mum
      theta(:) = asin(sqrt(occin(:, 1)/st%smear%el_per_state))*(M_HALF/M_PI)
      call minimize_multidim(MINMETHOD_BFGS, st%nst, theta, CNST(0.05), CNST(0.0001), &
          CNST(1e-12), CNST(1e-12), 200, objective_rdmft, write_iter_info_rdmft, objective, ierr)
      sumgim = rdm%occsum - st%qtot

      if (sumgi1*sumgim < M_ZERO) then
        mu2 = mum
      else
        mu1 = mum
        sumgi1 = sumgim
      end if

      ! check occ.num. threshold again after minimization
      do ist = 1, st%nst
        st%occ(ist,1) = M_TWO*sin(theta(ist)*M_PI*M_TWO)**2
        if (st%occ(ist,1) <= thresh_occ ) st%occ(ist,1) = thresh_occ
      end do

      if (abs(sumgim) < rdm%toler .or. abs((mu1-mu2)*M_HALF) < rdm%toler)  exit
    end do

    nullify(rdm_ptr)

    if (icycle >= 50) then
      write(message(1),'(a,1x,f11.4)') 'Bisection ended without finding mu, sum of occupation numbers:', rdm%occsum
      call messages_fatal(1)
    end if

    do ist = 1, st%nst
      st%occ(ist, 1) = st%smear%el_per_state*sin(theta(ist)*M_PI*M_TWO)**2
    end do

    objective = objective + rdm%mu*(rdm%occsum - rdm%qtot)
    energy = objective

    write(message(1),'(a4,5x,a12)')'#st','Occupation'
    call messages_info(1)   

    do ist = 1, st%nst
      write(message(1),'(i4,3x,f14.12)') ist, st%occ(ist, 1)
      call messages_info(1)
    end do

    write(message(1),'(a,3x,f11.9)') 'Sum of occupation numbers: ', rdm%occsum
    write(message(2),'(a,11x,es20.10)') 'Total energy: ', units_from_atomic(units_out%energy, energy + hm%ep%eii)
    call messages_info(2)   

    SAFE_DEALLOCATE_A(occin)
    SAFE_DEALLOCATE_A(theta)

    call profiling_out(prof_occ)
    POP_SUB(scf_occ)
  end subroutine scf_occ

  ! ---------------------------------------------------------
  subroutine objective_rdmft(size, theta, objective, getgrad, df)
    integer,     intent(in)    :: size
    REAL_DOUBLE, intent(in)    :: theta(size)
    REAL_DOUBLE, intent(inout) :: objective
    integer,     intent(in)    :: getgrad
    REAL_DOUBLE, intent(inout) :: df(size)

    integer :: ist
    FLOAT   :: thresh_occ, thresh_theta
    FLOAT, allocatable :: dE_dn(:),occ(:)
 
    PUSH_SUB(objective_rdmft)

    ASSERT(size == rdm_ptr%nst)

    SAFE_ALLOCATE(dE_dn(1:size))
    SAFE_ALLOCATE(occ(1:size))

    occ = M_ZERO
      
    ! Defines a threshold on occ nums to avoid numerical instabilities.
    ! Needs to be changed consistently with the same variable in scf_occ
    thresh_occ = CNST(1e-14)
    thresh_theta = asin(sqrt(thresh_occ/M_TWO))*(M_HALF/M_PI)
    
    do ist = 1, size
      occ(ist) = M_TWO*sin(theta(ist)*M_PI*M_TWO)**2
      if (occ(ist) <= thresh_occ ) occ(ist) = thresh_occ
    end do
      
    rdm_ptr%occsum = sum(occ(1:size))
      
    !calculate the total energy without nuclei interaction and the energy
    !derivatives with respect to the occupation numbers

    call total_energy_rdm(rdm_ptr, occ, objective, dE_dn)
    do ist = 1, size
      if (occ(ist) <= thresh_occ ) then
        df(ist) = M_FOUR*M_PI*sin(M_FOUR*thresh_theta*M_PI)*(dE_dn(ist) - rdm_ptr%mu)
      else
        df(ist) = M_FOUR*M_PI*sin(M_FOUR*theta(ist)*M_PI)*(dE_dn(ist) - rdm_ptr%mu)
      end if
    end do
    objective = objective - rdm_ptr%mu*(rdm_ptr%occsum - rdm_ptr%qtot)

    SAFE_DEALLOCATE_A(dE_dn)
    SAFE_DEALLOCATE_A(occ)

    POP_SUB(objective_rdmft)
  end subroutine objective_rdmft

  ! ---------------------------------------------------------
  subroutine write_iter_info_rdmft(iter, size, energy, maxdr, maxdf, theta)
    integer,     intent(in) :: iter
    integer,     intent(in) :: size
    FLOAT,       intent(in) :: energy, maxdr, maxdf
    FLOAT,       intent(in) :: theta(size)

    PUSH_SUB(write_iter_info_rdmft)

    ! Nothing to do.

    POP_SUB(write_iter_info_rdmft)
  end subroutine write_iter_info_rdmft

  ! scf for the natural orbitals
  subroutine scf_orb(rdm, namespace, gr, st, hm, space, energy)
    type(rdm_t),              intent(inout) :: rdm
    type(namespace_t),        intent(in)    :: namespace
    type(grid_t),             intent(in)    :: gr !< grid
    type(states_elec_t),      intent(inout) :: st !< States
    type(hamiltonian_elec_t), intent(in)    :: hm !< Hamiltonian
    type(space_t),            intent(in)    :: space
    FLOAT,                    intent(out)   :: energy    

    integer :: ist, jst
    FLOAT, allocatable ::  lambda(:,:), fo(:,:)
    type(profile_t), save :: prof_orb_basis
    
    PUSH_SUB(scf_orb)
    call profiling_in(prof_orb_basis, "SCF_ORB_BASIS")

    !matrix of Lagrange Multipliers from  Equation (8), Piris and Ugalde, Vol. 30, No. 13, J. Comput. Chem. 
    SAFE_ALLOCATE(lambda(1:st%nst,1:st%nst)) 
    SAFE_ALLOCATE(fo(1:st%nst, 1:st%nst))    !Generalized Fockian Equation (11) 

    lambda = M_ZERO
    fo = M_ZERO

    call construct_lambda(namespace, hm, st, gr, lambda, rdm)
    
    !Set up fo matrix 
    if (rdm%iter==1) then
      do ist = 1, st%nst
        do jst = 1, ist
          fo(ist, jst) = M_HALF*(lambda(ist, jst) + lambda(jst, ist))
          fo(jst, ist) = fo(ist, jst)
        end do
      end do
    else
      do ist = 1, st%nst
        do jst = 1, ist - 1
          fo(jst, ist) = - ( lambda(jst, ist) - lambda(ist ,jst))
        end do
      end do
      rdm%maxfo = maxval(abs(fo))
      do ist = 1, st%nst
        fo(ist, ist) = rdm%evalues(ist)
        do jst = 1, ist-1
          if(abs(fo(jst, ist)) > rdm%scale_f) then
            fo(jst, ist) = rdm%scale_f*fo(jst,ist)/abs(fo(jst, ist))
          end if
          fo(ist, jst) = fo(jst, ist)
        end do
      end do
    end if
 
    call lalg_eigensolve(st%nst, fo, rdm%evalues)
    call assign_eigfunctions(rdm, st, fo)
    call sum_integrals(rdm) ! to calculate rdm%Coul and rdm%Exch with the new rdm%vecnat 
    call rdm_derivatives(rdm, namespace, hm, st, gr, space)
    call total_energy_rdm(rdm, st%occ(:,1), energy)

    SAFE_DEALLOCATE_A(lambda) 
    SAFE_DEALLOCATE_A(fo) 

    call profiling_out(prof_orb_basis)
    POP_SUB(scf_orb)
  end subroutine scf_orb

  
  !-----------------------------------------------------------------
  ! Minimize the total energy wrt. an orbital by conjugate gradient
  !-----------------------------------------------------------------
  subroutine scf_orb_cg(rdm, namespace, space, gr, ions, st, ks, hm, energy)
    type(rdm_t),              intent(inout) :: rdm
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(grid_t),             intent(in)    :: gr !< grid
    type(ions_t),             intent(in)    :: ions !< geometry
    type(states_elec_t),      intent(inout) :: st !< States
    type(v_ks_t),             intent(inout) :: ks !< Kohn-Sham
    type(hamiltonian_elec_t), intent(inout) :: hm !< Hamiltonian
    FLOAT,                    intent(out)   :: energy

    integer            :: ik, ist
    
    type(profile_t), save :: prof_orb_cg
    
    PUSH_SUB(scf_orb_cg)
    call profiling_in(prof_orb_cg, "CG")
    
    call v_ks_calc(ks, namespace, space, hm, st, ions)
    call hamiltonian_elec_update(hm, gr%mesh, namespace, space)
    
    rdm%eigens%converged = 0
    if(mpi_grp_is_root(mpi_world) .and. .not. debug%info) then
      call loct_progress_bar(-1, st%lnst*st%d%kpt%nlocal)
    end if
    do ik = st%d%kpt%start, st%d%kpt%end
      rdm%eigens%matvec = 0  
      call deigensolver_cg2(namespace, gr%mesh, st, hm, hm%xc, rdm%eigens%pre, rdm%eigens%tolerance, rdm%eigens%es_maxiter, &
        rdm%eigens%converged(ik), ik, rdm%eigens%diff(:, ik), rdm%eigens%orthogonalize_to_all, &
        rdm%eigens%conjugate_direction, rdm%eigens%additional_terms, rdm%eigens%energy_change_threshold)
  
      if(st%calc_eigenval .and. .not. rdm%eigens%folded_spectrum) then
        ! recheck convergence after subspace diagonalization, since states may have reordered (copied from eigensolver_run)
        rdm%eigens%converged(ik) = 0
        do ist = 1, st%nst
          if(rdm%eigens%diff(ist, ik) < rdm%eigens%tolerance) then
            rdm%eigens%converged(ik) = ist
          else
            exit
          end if
        end do
      end if
    end do
    
    if(mpi_grp_is_root(mpi_world) .and. .not. debug%info) then
      write(stdout, '(1x)')
    end if
    
    ! calculate total energy with new states
    call density_calc (st, gr, st%rho)
    call v_ks_calc(ks, namespace, space, hm, st, ions)
    call hamiltonian_elec_update(hm, gr%mesh, namespace, space)
    call rdm_derivatives(rdm, namespace, hm, st, gr, space)
    
    call total_energy_rdm(rdm, st%occ(:,1), energy)

    call profiling_out(prof_orb_cg)
    POP_SUB(scf_orb_cg)
  end subroutine scf_orb_cg


  ! ----------------------------------------
  ! constructs the Lagrange multiplyers needed for the orbital minimization
  subroutine construct_lambda(namespace, hm, st, gr, lambda, rdm)
    type(namespace_t),        intent(in)    :: namespace
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(states_elec_t),      intent(inout) :: st
    type(grid_t),             intent(in)    :: gr
    FLOAT,                    intent(out)   :: lambda(:,:) !< (1:st%nst, 1:st%nst)
    type(rdm_t),              intent(inout) :: rdm

    FLOAT, allocatable :: hpsi(:,:), hpsi1(:,:), dpsi(:,:), dpsi1(:,:) 
    FLOAT, allocatable :: fock(:,:,:), fvec(:)
    integer :: ist, iorb, jorb, jst

    PUSH_SUB(construct_lambda)

    lambda = M_ZERO

    !calculate the Lagrange multiplyer lambda matrix on the grid, Eq. (9), Piris and Ugalde, Vol. 30, No. 13, J. Comput. Chem.
    if (.not. rdm%do_basis) then
      SAFE_ALLOCATE(hpsi(1:gr%mesh%np,1:st%d%dim))
      SAFE_ALLOCATE(hpsi1(1:gr%mesh%np,1:st%d%dim))
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part ,1:st%d%dim))
      SAFE_ALLOCATE(dpsi1(1:gr%mesh%np_part ,1:st%d%dim))

      do iorb = 1, st%nst
        call states_elec_get_state(st, gr%mesh, iorb, 1, dpsi)
        call dhamiltonian_elec_apply_single(hm, namespace, gr%mesh, dpsi, hpsi, iorb, 1)

        do jorb = iorb, st%nst  
          ! calculate <phi_j|H|phi_i> =lam_ji
          call states_elec_get_state(st, gr%mesh, jorb, 1, dpsi1)
          lambda(jorb, iorb) = dmf_dotp(gr%mesh, dpsi1(:,1), hpsi(:,1))
          
          ! calculate <phi_i|H|phi_j>=lam_ij
          if (.not. iorb == jorb ) then
            call dhamiltonian_elec_apply_single(hm, namespace, gr%mesh, dpsi1, hpsi1, jorb, 1)
            lambda(iorb, jorb) = dmf_dotp(gr%mesh, dpsi(:,1), hpsi1(:,1))
          end if
        end do
      end do
 

    else ! calculate the same lambda matrix on the basis
      !call sum_integrals(rdm)
      SAFE_ALLOCATE(fvec(1:st%nst))
      SAFE_ALLOCATE(fock(1:st%nst, 1:st%nst, 1:st%nst))
      fock = M_ZERO
      
      do iorb = 1, st%nst
        do ist = 1, st%nst
          do jst = 1, ist
            fock(ist, jst, iorb) = st%occ(iorb, 1)*rdm%eone_int(ist,jst)
            do jorb = 1, st%nst
              !The coefficient of the Exchange term below is only for the Mueller functional
              fock(ist ,jst, iorb) =  fock(ist, jst, iorb) + st%occ(iorb, 1)*st%occ(jorb, 1)*rdm%Coul(ist, jst, jorb)  &
                                      - sqrt(st%occ(iorb, 1))*sqrt(st%occ(jorb, 1))*rdm%Exch(ist, jst, jorb)
            end do
            fock(jst, ist, iorb) = fock(ist, jst, iorb)
          end do
        end do
      end do

      do jorb = 1, st%nst
        do ist = 1, st%nst
          fvec(ist) = M_ZERO
          do jst = 1, st%nst
            fvec(ist) = fvec(ist) + fock(ist, jst, jorb)*rdm%vecnat(jst, jorb)
          end do
        end do
        do iorb= 1, st%nst
          lambda(iorb, jorb) = M_ZERO
          do ist = 1, st%nst
            lambda(iorb, jorb) = lambda(iorb, jorb) + rdm%vecnat(ist, iorb)*fvec(ist)
          end do
        end do
      end do
    end if


    if (.not. rdm%do_basis) then
      SAFE_DEALLOCATE_A(hpsi)
      SAFE_DEALLOCATE_A(hpsi1)
      SAFE_DEALLOCATE_A(dpsi)
      SAFE_DEALLOCATE_A(dpsi1)
    else
      SAFE_DEALLOCATE_A(fvec) 
      SAFE_DEALLOCATE_A(fock) 
    end if
   
    POP_SUB(construct_lambda)
  end subroutine construct_lambda
   
  ! ----------------------------------------
  
  ! finds the new states after the minimization of the orbitals (Piris method)
  subroutine assign_eigfunctions(rdm, st, lambda)
    type(rdm_t),         intent(inout) :: rdm
    type(states_elec_t), intent(inout) :: st
    FLOAT,               intent(in)    :: lambda(:, :)
    
    integer :: iorb, jorb, ist
    FLOAT, allocatable :: vecnat_new(:,:)

    PUSH_SUB(assign_eigenfunctions)

    SAFE_ALLOCATE(vecnat_new(1:st%nst, 1:st%nst))
    do iorb = 1, st%nst
      do ist = 1, st%nst
        vecnat_new(ist, iorb) = M_ZERO
        do jorb = 1, st%nst
          vecnat_new(ist , iorb) = vecnat_new(ist, iorb) + rdm%vecnat(ist, jorb)*lambda(jorb, iorb)
        end do
      end do
    end do
    
    rdm%vecnat = vecnat_new

    SAFE_DEALLOCATE_A(vecnat_new)

    POP_SUB(assign_eigenfunctions)
  end subroutine assign_eigfunctions
   
  ! --------------------------------------------

  ! calculates the total energy when only the occupation numbers are updated
  subroutine total_energy_rdm(rdm, occ, energy, dE_dn)
    type(rdm_t),          intent(in)  :: rdm
    FLOAT,                intent(in)  :: occ(:)
    FLOAT,                intent(out) :: energy
    FLOAT, optional,      intent(out) :: dE_dn(:) !< (1:st%nst)
     
    integer :: ist, jst
    FLOAT, allocatable :: V_h(:), V_x(:)
     
    PUSH_SUB(total_energy_rdm)
  
    SAFE_ALLOCATE(V_h(1:rdm%nst))
    SAFE_ALLOCATE(V_x(1:rdm%nst))

    energy = M_ZERO
    V_h = M_ZERO
    V_x = M_ZERO

    !Calculate hartree and exchange contribution 
    !This is only for the Mueller functional and has to be changed
    do ist = 1, rdm%nst
      do jst = 1, rdm%nst
        V_h(ist) = V_h(ist) + occ(jst)*rdm%hartree(ist, jst)
        V_x(ist) = V_x(ist) - sqrt(occ(jst))*rdm%exchange(ist, jst) 
      end do
      V_x(ist) = V_x(ist)*M_HALF/max(sqrt(occ(ist)), CNST(1.0e-16))
    end do


    !Calculate the energy derivative with respect to the occupation numbers
    if (present(dE_dn)) then
      dE_dn(:) = rdm%eone(:) + V_h(:) + V_x(:)
    end if

    !Total energy calculation without nuclei interaction  
    do ist = 1, rdm%nst
      energy = energy + occ(ist)*rdm%eone(ist) &
                      + M_HALF*occ(ist)*V_h(ist) & 
                      + occ(ist)*V_x(ist)
    end do

    SAFE_DEALLOCATE_A(V_h)
    SAFE_DEALLOCATE_A(V_x)
    
    POP_SUB(total_energy_rdm)
  end subroutine total_energy_rdm

  ! ----------------------------------------
  ! calculates the derivatives of the energy terms with respect to the occupation numbers
  subroutine rdm_derivatives(rdm, namespace, hm, st, gr, space)
    type(rdm_t),              intent(inout) :: rdm
    type(namespace_t),        intent(in)    :: namespace
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(states_elec_t),      intent(in)    :: st 
    type(grid_t),             intent(in)    :: gr
    type(space_t),            intent(in)    :: space

    
    FLOAT, allocatable  :: hpsi(:,:), rho1(:), rho(:), dpsi(:,:), dpsi2(:,:)
    FLOAT, allocatable  :: v_ij(:,:,:,:,:)
    FLOAT               :: dd
    type(states_elec_t) :: xst
    
    integer :: ist, jst, nspin_, iorb, jorb

    PUSH_SUB(rdm_derivatives) 


    nspin_ = min(st%d%nspin, 2)
   
    if (rdm%do_basis.eqv..false.) then 
      SAFE_ALLOCATE(hpsi(1:gr%mesh%np, 1:st%d%dim))
      SAFE_ALLOCATE(rho1(1:gr%mesh%np))
      SAFE_ALLOCATE(rho(1:gr%mesh%np))
      SAFE_ALLOCATE(dpsi(1:gr%mesh%np_part, 1:st%d%dim))
      SAFE_ALLOCATE(dpsi2(1:gr%mesh%np, 1:st%d%dim))
      SAFE_ALLOCATE(v_ij(1:gr%mesh%np, 1:st%nst, 1:st%nst, 1:st%d%nik, 1:st%d%nik))

      v_ij = M_ZERO
      rdm%eone = M_ZERO
      rdm%hartree = M_ZERO
      rdm%exchange = M_ZERO

      !derivative of one-electron energy with respect to the natural orbitals occupation number
      do ist = 1, st%nst
        call states_elec_get_state(st, gr%mesh, ist, 1, dpsi)

        ! calculate one-body energy
        call dhamiltonian_elec_apply_single(hm, namespace, gr%mesh, dpsi, hpsi, ist, 1, &
                              terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
        rdm%eone(ist) = dmf_dotp(gr%mesh, dpsi(:, 1), hpsi(:, 1))
      end do

      ! integrals used for the hartree and exchange parts of the total energy and their derivatives
      ! maybe better to let that be done from the lower level routines like hamiltonian apply?
      !
      ! only used to calculate total energy
      call xst%nullify()
      call dexchange_operator_compute_potentials(hm%exxop, namespace, space, gr%mesh, st, xst, hm%kpoints, F_out = v_ij)
      call states_elec_end(xst)

      do ist = 1, st%nst
        call states_elec_get_state(st, gr%mesh, ist, 1, dpsi)

        rho1(1:gr%mesh%np) = dpsi(1:gr%mesh%np, 1)**2

        do jst = ist, st%nst
          rdm%hartree(ist, jst) = dmf_dotp(gr%mesh, rho1, v_ij(:,jst, jst, 1, 1))
          rdm%hartree(jst, ist) = rdm%hartree(ist, jst)
          call states_elec_get_state(st, gr%mesh, jst, 1, dpsi2)
          rho(1:gr%mesh%np) = dpsi2(1:gr%mesh%np, 1)*dpsi(1:gr%mesh%np, 1)
          rdm%exchange(ist, jst) = dmf_dotp(gr%mesh, rho, v_ij(:, ist, jst, 1, 1))
          rdm%exchange(jst, ist) = rdm%exchange(ist, jst)
        end do
      end do

    
      SAFE_DEALLOCATE_A(hpsi)
      SAFE_DEALLOCATE_A(rho)
      SAFE_DEALLOCATE_A(rho1)
      SAFE_DEALLOCATE_A(dpsi)
      SAFE_DEALLOCATE_A(dpsi2)
      SAFE_DEALLOCATE_A(v_ij)

    else !if energy derivatives are expanded in a basis set

      do iorb = 1, st%nst
        rdm%eone(iorb) = M_ZERO
        do ist = 1, st%nst
          do jst = 1, st%nst
            dd = rdm%vecnat(ist, iorb)*rdm%vecnat(jst, iorb)  
            rdm%eone(iorb) = rdm%eone(iorb) + dd*rdm%eone_int(ist, jst)
          end do
        end do
      end do

      do iorb = 1, st%nst
        do jorb =1 , iorb
          rdm%hartree(iorb ,jorb) = M_ZERO; rdm%exchange(iorb,jorb) = M_ZERO
          do ist =1, st%nst
            do jst =1, st%nst
              dd = rdm%vecnat(ist, iorb)*rdm%vecnat(jst, iorb)
              rdm%hartree(iorb ,jorb) = rdm%hartree(iorb ,jorb)+rdm%Coul(ist,jst, jorb)*dd 
              rdm%exchange(iorb ,jorb) = rdm%exchange(iorb ,jorb)+rdm%Exch(ist,jst, jorb)*dd 
            end do
          end do
          rdm%hartree(jorb, iorb) = rdm%hartree(iorb, jorb); rdm%exchange(jorb, iorb) = rdm%exchange(iorb, jorb)
        end do
      end do
    end if
     
    POP_SUB(rdm_derivatives)
  end subroutine rdm_derivatives
  
  ! --------------------------------------------
  !calculates the one electron integrals in the basis of the initial orbitals
  subroutine rdm_integrals(rdm, namespace, hm, st, mesh)
    type(rdm_t),              intent(inout) :: rdm
    type(namespace_t),        intent(in)    :: namespace
    type(hamiltonian_elec_t), intent(in)    :: hm
    type(states_elec_t),      intent(in)    :: st 
    type(mesh_t),             intent(in)    :: mesh
    
    FLOAT, allocatable :: hpsi(:,:)
    FLOAT, allocatable :: dpsi(:,:), dpsi2(:,:)
    integer :: ist, jst

    PUSH_SUB(rdm_integrals)
 
    SAFE_ALLOCATE(dpsi(1:mesh%np_part, 1:st%d%dim))
    SAFE_ALLOCATE(dpsi2(1:mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:mesh%np, 1:st%d%dim))

    !calculate integrals of the one-electron energy term with respect to the initial orbital basis
    do ist = 1, st%nst
      call states_elec_get_state(st, mesh, ist, 1, dpsi)
      do jst = ist, st%nst
        call states_elec_get_state(st, mesh, jst, 1, dpsi2)

        ! calculate one-body integrals
        call dhamiltonian_elec_apply_single(hm, namespace, mesh, dpsi, hpsi, ist, 1, &
                              terms = TERM_KINETIC + TERM_LOCAL_EXTERNAL + TERM_NON_LOCAL_POTENTIAL)
        rdm%eone_int(jst, ist) = dmf_dotp(mesh, dpsi2(:, 1), hpsi(:, 1))
        rdm%eone_int(ist, jst) = rdm%eone_int(jst, ist)
      end do
    end do

    SAFE_DEALLOCATE_A(hpsi)    
    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(dpsi2)
    
    POP_SUB(rdm_integrals) 
  end subroutine rdm_integrals

  ! --------------------------------------------
  ! constructs the Hartree and Exchange part of the RDMFT Fock matrix
  subroutine sum_integrals(rdm)
    type(rdm_t),        intent(inout) :: rdm

    integer :: ist, jst, kst, lst, iorb, icount
    logical :: inv_pairs
    FLOAT :: two_int, wij, wik, wil, wjk, wjl, wkl
    FLOAT, allocatable :: dm(:,:,:)

    PUSH_SUB(sum_integrals)

    SAFE_ALLOCATE(dm(1:rdm%nst, 1:rdm%nst, 1:rdm%nst))
     
    rdm%Coul = M_ZERO  
    rdm%Exch = M_ZERO 
    dm = M_ZERO 
  
    do iorb = 1, rdm%nst
      do ist = 1, rdm%nst
        do jst = 1, ist
          dm(ist, jst, iorb) = rdm%vecnat(ist, iorb)*rdm%vecnat(jst, iorb)
          dm(jst, ist, iorb) = dm(ist, jst, iorb)
        end do
      end do
    end do

    do icount = 1, rdm%n_twoint

      ist = rdm%i_index(1,icount) 
      jst = rdm%j_index(1,icount) 
      kst = rdm%k_index(1,icount) 
      lst = rdm%l_index(1,icount) 

      two_int = rdm%twoint(icount)
         
      ! create weights of unique integrals 
      if(ist == jst) then
        wij = M_ONE
      else
        wij = M_TWO
      endif
      if(kst == lst) then
        wkl = M_ONE
      else
        wkl = M_TWO
      endif

      if(ist == kst .and. jst /= lst) then
        wik = M_TWO
      else
        wik = M_ONE
      endif
      if(ist == lst .and. jst /= kst) then
        wil = M_TWO
      else
        wil = M_ONE
      endif
      if(jst == kst .and. ist /= lst) then
        wjk = M_TWO
      else
        wjk = M_ONE
      endif
      if(jst == lst .and. ist /= kst) then
        wjl = M_TWO
      else
        wjl = M_ONE
      endif

      inv_pairs = (ist /= kst .or. jst /= lst)

      do iorb = 1, rdm%nst 

        !the Hartree terms
        rdm%Coul(ist, jst, iorb) = rdm%Coul(ist, jst, iorb) + dm(kst, lst, iorb)*wkl*two_int
        if (inv_pairs) rdm%Coul(kst, lst, iorb) = rdm%Coul(kst, lst, iorb) + dm(ist, jst, iorb)*wij*two_int

        !the exchange terms
        !weights are only included if they can differ from one
        rdm%Exch(ist, kst, iorb) = rdm%Exch(ist, kst, iorb) + two_int*dm(jst, lst, iorb)*wik
        if (kst /= lst) then 
          rdm%Exch(ist, lst, iorb) = rdm%Exch(ist, lst, iorb) + two_int*dm(jst, kst, iorb)*wil
        end if
        if (ist /= jst) then
          if(jst >= kst) then
            rdm%Exch(jst, kst, iorb) = rdm%Exch(jst, kst, iorb) + two_int*dm(ist, lst, iorb)*wjk
          else
            if (inv_pairs) rdm%Exch(kst, jst, iorb) = rdm%Exch(kst, jst, iorb) + two_int*dm(ist, lst, iorb)
          end if
        end if
        if (ist /=jst .and. kst /= lst) then
          if (jst >= lst) then
            rdm%Exch(jst, lst, iorb) = rdm%Exch(jst, lst, iorb) + two_int*dm(ist, kst, iorb)*wjl
          else
            if (inv_pairs) rdm%Exch(lst, jst, iorb) = rdm%Exch(lst, jst, iorb) + two_int*dm(ist, kst, iorb)
          end if
        end if

      end do !iorb
    end do !icount
   
    do iorb =1, rdm%nst
      do ist = 1, rdm%nst
        do jst = 1, ist-1
          rdm%Coul(jst, ist, iorb) = rdm%Coul(ist, jst, iorb)
          rdm%Exch(jst, ist, iorb) = rdm%Exch(ist, jst, iorb)
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(dm)

    POP_SUB(sum_integrals)
  end subroutine sum_integrals

end module rdmft_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:




