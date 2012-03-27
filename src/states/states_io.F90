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

module states_io_m
  use datasets_m
  use geometry_m
  use global_m
  use grid_m
  use io_m
  use kpoints_m
  use lalg_basic_m
  use loct_m
  use math_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use mpi_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2
  use mpi_lib_m
  use ob_green_m
  use ob_interface_m
  use parser_m
  use profiling_m
  use simul_box_m
  use states_m
  use states_dim_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m

  implicit none

  private

  public ::                           &
    states_write_eigenvalues,         &
    states_write_dos,                 &
    states_write_tpa,                 &
    states_write_bands,               &
    states_write_fermi_energy,        &
    states_dump,                      &
    states_init_self_energy,          &
    states_write_proj_lead_wf,        &
    states_read_proj_lead_wf

contains

  ! ---------------------------------------------------------

  subroutine states_write_eigenvalues(iunit, nst, st, sb, error, st_start)
    integer,           intent(in) :: iunit, nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb
    FLOAT, optional,   intent(in) :: error(nst, st%d%nik)
    integer, optional, intent(in) :: st_start

    integer :: ik, ist, ns, is, idir, st_start_
    FLOAT :: occ, kpoint(1:MAX_DIM)
    character(len=80) tmp_str(MAX_DIM), cspin

    PUSH_SUB(states_write_eigenvalues)

    st_start_ = 1
    if(present(st_start)) st_start_ = st_start

    ns = 1
    if(st%d%nspin == 2) ns = 2

    message(1) = 'Eigenvalues [' // trim(units_abbrev(units_out%energy)) // ']'
    call messages_info(1, iunit)
    if (st%d%nik > ns) then
      message(1) = 'k-points [' // trim(units_abbrev(unit_one/units_out%length)) //']'
      call messages_info(1, iunit)
    end if

    if(.not. mpi_grp_is_root(mpi_world)) then
      POP_SUB(states_write_eigenvalues)
      return
    end if

    if(present(error)) then
      if(st%d%ispin .eq. SPINORS) then
        write(message(1), '(a4,1x,a5,1x,a12,1x,a12,2x,a4,4x,a4,4x,a4,5x,a5)')   &
          '#st',' Spin',' Eigenvalue', 'Occupation ', '<Sx>', '<Sy>', '<Sz>', 'Error'
      else
        write(message(1), '(a4,1x,a5,1x,a12,4x,a12,1x,a10)')   &
          '#st',' Spin',' Eigenvalue', 'Occupation ', 'Error'
      end if
    else
      if(st%d%ispin .eq. SPINORS) then
        write(message(1), '(a4,1x,a5,1x,a12,1x,a12,2x,a4,4x,a4,4x,a4)')   &
          '#st',' Spin',' Eigenvalue', 'Occupation ', '<Sx>', '<Sy>', '<Sz>'
      else
        write(message(1), '(a4,1x,a5,1x,a12,4x,a12,1x)')       &
          '#st',' Spin',' Eigenvalue', 'Occupation '
      end if
    end if
    call messages_info(1, iunit)

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then
        kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik))
        kpoint(1:sb%dim) = units_from_atomic(unit_one/units_out%length, kpoint(1:sb%dim))
        write(message(1), '(a,i4,a)') '#k =', ik, ', k = ('
        do idir = 1, sb%dim
          write(tmp_str(1), '(f12.6)') kpoint(idir)
          message(1) = trim(message(1))//trim(tmp_str(1))
          if(idir < sb%dim) message(1) = trim(message(1))//','
        enddo
        message(1) = trim(message(1))//')'
        call messages_info(1, iunit)
      end if

      do ist = st_start_, nst
        do is = 0, ns-1
          if(ist > st%nst) then
            occ = M_ZERO
          else
            occ = st%occ(ist, ik+is)
          end if

          if(is .eq. 0) cspin = 'up'
          if(is .eq. 1) cspin = 'dn'
          if(st%d%ispin .eq. UNPOLARIZED .or. st%d%ispin .eq. SPINORS) cspin = '--'

          write(tmp_str(1), '(i4,3x,a2)') ist, trim(cspin)
          if(simul_box_is_periodic(sb)) then
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,3x,4f5.2)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik)), occ, st%spin(1:3, ist, ik)
              if(present(error)) write(tmp_str(3), '(a7,es8.1,a1)')'      (', error(ist, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik+is)), occ
              if(present(error)) write(tmp_str(3), '(a7,es8.1,a1)')'      (', error(ist, ik), ')'
            end if
          else
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,5x,f5.2,3x,3f8.4)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik)), occ, st%spin(1:3, ist, ik)
              if(present(error)) write(tmp_str(3), '(a3,es8.1,a1)')'  (', error(ist, ik+is), ')'
            else
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik+is)), occ
              if(present(error)) write(tmp_str(3), '(a7,es8.1,a1)')'      (', error(ist, ik), ')'
            end if
          end if
          if(present(error)) then
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))//trim(tmp_str(3))
          else
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))
          end if
          call messages_info(1, iunit)
        end do
      end do
    end do

    POP_SUB(states_write_eigenvalues)
  end subroutine states_write_eigenvalues


  ! ---------------------------------------------------------
  subroutine states_write_bands(dir, nst, st, sb)
    character(len=*),  intent(in) :: dir    
    integer,           intent(in) :: nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb

    integer :: idir, ist, ik, ns, is
    integer, allocatable :: iunit(:)
    FLOAT   :: factor(MAX_DIM), kpoint(1:MAX_DIM)
    logical :: grace_mode, gnuplot_mode
    character(len=80) :: filename    

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(states_write_bands)

    !%Variable OutputBandsGnuplotMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in Gnuplot-friendly format to <tt>bands-gp.dat</tt>
    !% (or <tt>band-gp-is.dat</tt> if spin-polarized).
    !%End
    call parse_logical(datasets_check('OutputBandsGnuplotMode'), .true., gnuplot_mode)

    !%Variable OutputBandsGraceMode
    !%Type logical
    !%Default no
    !%Section Output
    !%Description
    !% The band file will be written in Grace-friendly format to <tt>bands-grace.dat</tt>
    !% (or <tt>bands-grace-is.dat</tt> if spin-polarized).
    !%End
    call parse_logical(datasets_check('OutputBandsGraceMode'), .false., grace_mode)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    SAFE_ALLOCATE(iunit(0:ns-1))

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do idir = 1, sb%dim
      factor(idir) = M_ONE
      if (sb%klattice(idir, idir) /= M_ZERO) factor(idir) = sb%klattice(idir, idir)
    end do

    if (gnuplot_mode) then
      do is = 0, ns-1
        if (ns .gt. 1) then
          write(filename, '(a,i1.1,a)') 'bands-gp-', is+1,'.dat'
        else
          write(filename, '(a)') 'bands-gp.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    

        ! write header
        write(iunit(is),'(a)',advance='no') '# '
        do idir = 1, sb%dim
          write(iunit(is),'(3a)',advance='no') 'k', index2axis(idir), ' '
        enddo
        write(iunit(is),'(a)',advance='no') '(unscaled), '
        do idir = 1, sb%dim
          write(iunit(is),'(3a)',advance='no') 'k', index2axis(idir), ' '
        enddo
        write(iunit(is),'(a, i6)') '(scaled), bands:', nst
      end do

      ! output bands in gnuplot format
      do ist = 1, nst
        do ik = 1, st%d%nik, ns
          do is = 0, ns - 1
            kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik + is))
            write(iunit(is),'(1x)',advance='no')
            do idir = 1, sb%dim
              write(iunit(is),'(f14.8)',advance='no') kpoint(idir)
            enddo
            do idir = 1, sb%dim
              write(iunit(is),'(f14.8)',advance='no') kpoint(idir) / factor(idir)
            enddo
            write(iunit(is),'(3x,f14.8)') units_from_atomic(units_out%energy, st%eigenval(ist, ik + is))
          end do
        end do
        do is = 0, ns-1
          write(iunit(is), '(a)')
        end do
      end do
      do is = 0, ns-1
        call io_close(iunit(is))
      end do
    end if

    if (grace_mode) then
      do is = 0, ns-1
        if (ns .gt. 1) then
          write(filename, '(a,i1.1,a)') 'bands-grace-', is+1,'.dat'
        else
          write(filename, '(a)') 'bands-grace.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    

        ! write header
        write(iunit(is),'(a)',advance='no') '# '
        do idir = 1, sb%dim
          write(iunit(is),'(3a)',advance='no') 'k', index2axis(idir), ' '
        enddo
        write(iunit(is),'(a)',advance='no') '(unscaled), '
        do idir = 1, sb%dim
          write(iunit(is),'(3a)',advance='no') 'k', index2axis(idir), ' '
        enddo
        write(iunit(is),'(a, i6)') '(scaled), bands:', nst
      end do

      ! output bands in xmgrace format, i.e.:
      ! k_x, k_y, k_z, e_1, e_2, ..., e_n
      do ik = 1, st%d%nik, ns
        do is = 0, ns-1
          kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik + is))
          write(iunit(is),'(1x)',advance='no')
          do idir = 1, sb%dim
            write(iunit(is),'(f14.8)',advance='no') kpoint(idir)
          enddo
          do idir = 1, sb%dim
            write(iunit(is),'(f14.8)',advance='no') kpoint(idir) / factor(idir)
          enddo
          write(iunit(is),'(3x)',advance='no')
          do ist = 1, nst
            write(iunit(is),'(f14.8)',advance='no') units_from_atomic(units_out%energy, st%eigenval(ist, ik + is))
          enddo
          write(iunit(is),'(a)')
        end do
      end do
      do is = 0, ns-1
        call io_close(iunit(is))
      end do        
    end if

    SAFE_DEALLOCATE_A(iunit)

    POP_SUB(states_write_bands)
  end subroutine states_write_bands

  ! ---------------------------------------------------------
  subroutine states_write_tpa(dir, gr, st)
    character(len=*), intent(in) :: dir
    type(grid_t),     intent(in) :: gr
    type(states_t),   intent(in) :: st

    type(block_t) :: blk
    integer       :: ncols, icoord, ist, ik, tpa_initialst, tpa_initialk
    integer       :: iunit

    FLOAT, allocatable  :: ff(:)
    CMPLX, allocatable  :: cff(:)
    FLOAT, allocatable  :: osc(:)
    FLOAT               :: transition_energy, osc_strength, dsf

    FLOAT, parameter    :: M_THRESHOLD = CNST(1.0e-6)
    logical             :: use_qvector = .false.
    FLOAT, allocatable  :: qvector(:)

    PUSH_SUB(states_write_tpa)

    ! find the orbital with half-occupation
    tpa_initialst = -1
    do ist = 1, st%nst
      do ik = 1, st%d%nik
        if (abs(st%occ(ist,ik)-0.5) .lt. M_THRESHOLD) then
          tpa_initialst = ist
          tpa_initialk  = ik
        end if
      end do
    end do

    ! make sure that half-occupancy was found
    if(tpa_initialst.eq.-1) then
      if(mpi_grp_is_root(mpi_world)) then

        message(1) = 'No orbital with half-occupancy found. TPA output is not written.'
        call messages_warning(1)

    POP_SUB(states_write_tpa)
return

      end if
    end if

    !%Variable MomentumTransfer
    !%Type block
    !%Section States
    !%Description
    !% Momentum-transfer vector <i>q</i> to be used when calculating matrix elements
    !% &lt;f|exp(iq.r)|i&gt;. This enables the calculation of the dynamical structure factor,
    !% which is closely related to generalized oscillator strengths.
    !% If the vector is not given, but TPA output is requested (<tt>Output = TPA</tt>),
    !% only the oscillator strengths are written in the output file.
    !% For example, to use <i>q</i> = (0.1, 0.2, 0.3), set
    !%
    !% <tt>%MomentumTransfer
    !% <br>&nbsp;&nbsp; 0.1 | 0.2 | 0.3
    !% <br>%</tt>
    !%End
    if(parse_block(datasets_check('MomentumTransfer'),blk)==0) then

      ! check if input makes sense
      ncols = parse_block_cols(blk, 0)

      if(ncols .ne. gr%mesh%sb%dim ) then ! wrong size

        if(mpi_grp_is_root(mpi_world)) then
          message(1) = 'Inconsistent size of momentum-transfer vector. It will not be used in the TPA calculation.'
          call messages_warning(1)
        end if

      else ! correct size

        use_qvector = .true.
        SAFE_ALLOCATE(qvector(1:gr%mesh%sb%dim))

        do icoord = 1,gr%mesh%sb%dim    !for x,y,z
          call parse_block_float(blk, 0, icoord-1, qvector(icoord))
          qvector(icoord) = units_to_atomic(unit_one / units_inp%length, qvector(icoord))
        end do

      end if

    end if

    ! calculate the matrix elements

    SAFE_ALLOCATE(ff(1:gr%mesh%np))
    if(use_qvector) then
      SAFE_ALLOCATE(cff(1:gr%mesh%np))
    end if
    SAFE_ALLOCATE(osc(1:gr%mesh%sb%dim))

    ! root writes output to file

    if(mpi_grp_is_root(mpi_world)) then

      iunit = io_open(trim(dir)//'/'//trim('tpa_xas'), action='write')    

      ! header
      if(use_qvector) then
        write (message(1),'(a1,a30,3(es14.5,1x),a1)') '#', ' momentum-transfer vector : (', &
          (units_from_atomic(unit_one / units_out%length, qvector(icoord)), icoord=1, gr%mesh%sb%dim),')'
        select case(gr%mesh%sb%dim)
          case(1); write(message(2), '(a1,4(a15,1x))') '#', 'E' , '<x>', '<f>', 'S(q,omega)'
          case(2); write(message(2), '(a1,5(a15,1x))') '#', 'E' , '<x>', '<y>', '<f>', 'S(q,omega)'
          case(3); write(message(2), '(a1,6(a15,1x))') '#', 'E' , '<x>', '<y>', '<z>', '<f>', 'S(q,omega)'
        end select
        call messages_info(2,iunit)
      else
        select case(gr%mesh%sb%dim)
          case(1); write(message(1), '(a1,3(a15,1x))') '#', 'E' , '<x>', '<f>'
          case(2); write(message(1), '(a1,4(a15,1x))') '#', 'E' , '<x>', '<y>', '<f>'
          case(3); write(message(1), '(a1,5(a15,1x))') '#', 'E' , '<x>', '<y>', '<z>', '<f>'
        end select
        call messages_info(1,iunit)
      end if

    end if

    ! loop through every state
    do ist = 1,st%nst

      ! final states are the unoccupied ones
      if (abs(st%occ(ist,tpa_initialk)) .lt. M_THRESHOLD) then

        osc_strength=M_ZERO
        transition_energy=st%eigenval(ist,tpa_initialk)-st%eigenval(tpa_initialst,tpa_initialk)

        ! dipole matrix elements <f|x|i> etc. -> oscillator strengths
        do icoord=1,gr%mesh%sb%dim    ! for x,y,z

          ff(1:gr%mesh%np) = st%dpsi(1:gr%mesh%np,1,tpa_initialst,tpa_initialk) * &
                       &  gr%mesh%x(1:gr%mesh%np,icoord)                        * &
                       &  st%dpsi(1:gr%mesh%np,1,ist,tpa_initialk)
          osc(icoord)  = dmf_integrate(gr%mesh, ff)
          osc_strength = osc_strength + 2.0/real(gr%mesh%sb%dim)*transition_energy*abs(osc(icoord))**2.0

        end do

        ! matrix elements <f|exp(iq.r)|i> -> dynamic structure factor
        if (use_qvector) then

          cff(1:gr%mesh%np) = TOCMPLX(st%dpsi(1:gr%mesh%np,1,tpa_initialst,tpa_initialk), M_ZERO) * &
                       &   TOCMPLX(st%dpsi(1:gr%mesh%np,1,ist,tpa_initialk), M_ZERO)
          do icoord=1,gr%mesh%sb%dim    ! for x,y,z
            cff(1:gr%mesh%np) = cff(1:gr%mesh%np) * exp(M_zI*gr%mesh%x(1:gr%mesh%np,icoord)*qvector(icoord))
          end do

          dsf = abs(zmf_integrate(gr%mesh, cff))**2.0
        end if

        ! write oscillator strengths (+ dynamic structure factor if qvector if given) into file
        if(mpi_grp_is_root(mpi_world)) then

          if(use_qvector) then
            write(message(1), '(1x,6(es15.8,1x))') units_from_atomic(units_out%energy, transition_energy), osc(:), osc_strength, &
                                                   units_from_atomic(unit_one/units_out%energy, dsf)
          else
            write(message(1), '(1x,6(es15.8,1x))') units_from_atomic(units_out%energy, transition_energy), osc(:), osc_strength
          endif

          call messages_info(1,iunit)

        end if

      end if

    end do

    ! finally close the file
    if(mpi_grp_is_root(mpi_world)) then
      call io_close(iunit)
    end if

    SAFE_DEALLOCATE_A(ff)
    if(use_qvector) then
      SAFE_DEALLOCATE_A(cff)
    end if
    SAFE_DEALLOCATE_A(osc)
    if (use_qvector) then
      SAFE_DEALLOCATE_A(qvector)
    end if
  
    POP_SUB(states_write_tpa)
 
  end subroutine states_write_tpa

  ! ---------------------------------------------------------
  subroutine states_write_dos(dir, st)
    character(len=*), intent(in) :: dir
    type(states_t),   intent(in) :: st

    integer :: ie, ik, ist, epoints, is, ns
    integer, allocatable :: iunit(:)
    FLOAT   :: emin, emax, de, gamma, energy
    FLOAT   :: evalmax, evalmin, tdos, eextend
    FLOAT, allocatable :: dos(:,:,:)
    character(len=64)  :: filename

    PUSH_SUB(states_write_dos)

    evalmin = minval(st%eigenval)
    evalmax = maxval(st%eigenval)
    ! we extend the energy mesh by this amount
    eextend  = (evalmax - evalmin) / M_FOUR

    !%Variable DOSEnergyMin
    !%Type float
    !%Section Output
    !%Description
    !% Lower bound for the energy mesh of the DOS.
    !% The default is the lowest eigenvalue, minus a quarter of the total range of eigenvalues.
    !%End
    call parse_float(datasets_check('DOSEnergyMin'), units_from_atomic(units_inp%energy, evalmin - eextend), emin)
    emin = units_to_atomic(units_inp%energy, emin)

    !%Variable DOSEnergyMax
    !%Type float
    !%Section Output
    !%Description
    !% Upper bound for the energy mesh of the DOS.
    !% The default is the highest eigenvalue, plus a quarter of the total range of eigenvalues.
    !%End
    call parse_float(datasets_check('DOSEnergyMax'), units_from_atomic(units_inp%energy, evalmax + eextend), emax)
    emax = units_to_atomic(units_inp%energy, emax)

    !%Variable DOSEnergyPoints
    !%Type integer
    !%Default 500
    !%Section Output
    !%Description
    !% Determines how many energy points <tt>Octopus</tt> should use for 
    !% the DOS energy grid.
    !%End
    call parse_integer(datasets_check('DOSEnergyPoints'), 500, epoints)

    !%Variable DOSGamma
    !%Type float
    !%Default 0.008 Ha
    !%Section Output
    !%Description
    !% Determines the width of the Lorentzian which is used for the DOS sum.
    !%End
    call parse_float(datasets_check('DOSGamma'), &
      units_from_atomic(units_inp%energy, CNST(0.008)), gamma)
    gamma = units_to_atomic(units_inp%energy, gamma)

    ! spacing for energy mesh
    de = (emax - emin) / (epoints - 1)

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    ! space for state-dependent DOS
    SAFE_ALLOCATE(dos(1:epoints, 1:st%nst, 0:ns-1))
    SAFE_ALLOCATE(iunit(0:ns-1))    

    ! compute band/spin-resolved density of states
    do ist = 1, st%nst

      do is = 0, ns-1
        if (ns.gt.1) then
          write(filename, '(a,i4.4,a,i1.1,a)') 'dos-', ist, '-', is+1,'.dat'
        else
          write(filename, '(a,i4.4,a)') 'dos-', ist, '.dat'
        end if
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a)') '# energy, band resolved DOS'
      end do

      do ie = 1, epoints
        energy = emin + (ie - 1) * de
        dos(ie, ist, :) = M_ZERO
        ! sum up Lorentzians
        do ik = 1, st%d%nik, ns
          do is = 0, ns-1
            dos(ie, ist, is) = dos(ie, ist, is) + st%d%kweights(ik+is) * M_ONE/M_Pi * &
              gamma / ( (energy - st%eigenval(ist, ik+is))**2 + gamma**2 )
          end do
        end do
        do is = 0, ns-1
          write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                        units_from_atomic(unit_one / units_out%energy, dos(ie, ist, is))
          call messages_info(1, iunit(is))
        end do
      end do

      do is = 0, ns-1
        call io_close(iunit(is))
      end do
    end do

    ! for spin-polarized calculations also output spin-resolved tDOS
    if(st%d%nspin .gt. 1) then    
      do is = 0, ns-1
        write(filename, '(a,i1.1,a)') 'total-dos-', is+1,'.dat'
        iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    
        ! write header
        write(iunit(is), '(a)') '# energy, total DOS (spin-resolved)'

        do ie = 1, epoints
          energy = emin + (ie - 1) * de
          tdos = M_ZERO
          do ist = 1, st%nst
            tdos = tdos + dos(ie, ist, is)
          end do
          write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                        units_from_atomic(unit_one / units_out%energy, tdos)
          call messages_info(1, iunit(is))
        end do

        call io_close(iunit(is))
      end do
    end if


    iunit(0) = io_open(trim(dir)//'/'//'total-dos.dat', action='write')    

    ! compute total density of states
    do ie = 1, epoints
      energy = emin + (ie - 1) * de
      tdos = M_ZERO
      do ist = 1, st%nst
        do is = 0, ns-1
          tdos = tdos + dos(ie, ist, is)
        end do
      end do
      write(message(1), '(2f12.6)') units_from_atomic(units_out%energy, energy), &
                                    units_from_atomic(unit_one / units_out%energy, tdos)
      call messages_info(1, iunit(0))
    end do

    call io_close(iunit(0))

    SAFE_DEALLOCATE_A(iunit)
    SAFE_DEALLOCATE_A(dos)

    POP_SUB(states_write_dos)
  end subroutine states_write_dos


  ! ---------------------------------------------------------
  subroutine states_write_fermi_energy(dir, st, mesh, sb)
    character(len=*),  intent(in) :: dir
    type(states_t), intent(inout) :: st
    type(mesh_t),      intent(in) :: mesh
    type(simul_box_t), intent(in) :: sb

    integer :: iunit, idir
    FLOAT :: maxdos
    FLOAT :: factor(MAX_DIM)
    character(len=100) :: str_tmp

    PUSH_SUB(states_write_fermi_energy)

    call states_fermi(st, mesh)

    iunit = io_open(trim(dir)//'/'//'bands-efermi.dat', action='write')    

    ! define the scaling factor to output k_i/G_i, instead of k_i
    do idir = 1, sb%dim
      factor(idir) = M_ONE
      if (sb%klattice(idir, idir) /= M_ZERO) factor(idir) = sb%klattice(idir, idir)
    end do

    ! write Fermi energy in a format that can be used together 
    ! with bands.dat
    write(message(1), '(a)') '# Fermi energy in a format compatible with bands-gp.dat'

    message(2)=""
    message(3)=""
    message(4)=""
    do idir = 1, sb%dim
      write(str_tmp, '(f12.6)') minval(sb%kpoints%reduced%point(idir, 1:sb%kpoints%reduced%npoints))
      message(2) = trim(message(2)) // trim(str_tmp)
      write(str_tmp, '(f12.6)') M_ZERO     ! Gamma point
      message(3) = trim(message(3)) // trim(str_tmp)
      write(str_tmp, '(f12.6)') maxval(sb%kpoints%reduced%point(idir, 1:sb%kpoints%reduced%npoints))
      message(4) = trim(message(4)) // trim(str_tmp)
    enddo
    do idir = 1, sb%dim
      write(str_tmp, '(f12.6)') minval(sb%kpoints%reduced%point(idir, 1:sb%kpoints%reduced%npoints) / factor(idir))
      message(2) = trim(message(2)) // trim(str_tmp)
      write(str_tmp, '(f12.6)') M_ZERO     ! Gamma point
      message(3) = trim(message(3)) // trim(str_tmp)
      write(str_tmp, '(f12.6)') maxval(sb%kpoints%reduced%point(idir, 1:sb%kpoints%reduced%npoints) / factor(idir))
      message(4) = trim(message(4)) // trim(str_tmp)
    enddo
    write(str_tmp, '(f12.6)') units_from_atomic(units_out%energy, st%smear%e_fermi)
    message(2) = trim(message(2)) // trim(str_tmp)
    message(3) = trim(message(3)) // trim(str_tmp)
    message(4) = trim(message(4)) // trim(str_tmp)

    call messages_info(4, iunit)
    call io_close(iunit)

    ! now we write the same information so that it can be used 
    ! together with total-dos.dat
    iunit = io_open(trim(dir)//'/'//'total-dos-efermi.dat', action='write')    

    write(message(1), '(a)') '# Fermi energy in a format compatible with total-dos.dat'    

    ! this is the maximum that tdos can reach
    maxdos = sum(st%d%kweights) * st%nst

    write(message(2), '(4f12.6)') units_from_atomic(units_out%energy, st%smear%e_fermi), M_ZERO
    write(message(3), '(4f12.6)') units_from_atomic(units_out%energy, st%smear%e_fermi), maxdos

    call messages_info(3, iunit)
    call io_close(iunit)

    POP_SUB(states_write_fermi_energy)
  end subroutine states_write_fermi_energy

  ! ---------------------------------------------------------

  subroutine states_dump(st, iunit)
    type(states_t), intent(in) :: st
     integer,       intent(in) :: iunit

     PUSH_SUB(states_dump)
     
     write(iunit, '(a20,1i10)')  'nst=                ', st%nst
     write(iunit, '(a20,1i10)')  'dim=                ', st%d%dim
     write(iunit, '(a20,1i10)')  'nik=                ', st%d%nik

     POP_SUB(states_dump)
  end subroutine states_dump

  ! ---------------------------------------------------------
  !> initialize the self-energy of the leads
  subroutine states_init_self_energy(st, gr, nspin, d_ispin, lead)
    type(states_t),      intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    integer,             intent(in)    :: nspin
    integer,             intent(in)    :: d_ispin
    type(lead_t),        intent(in)    :: lead(:) ! Diagonal and off-diagonal block of the lead Hamiltonian.

    character(len=2)      :: spin
    character(len=256)    :: fmt, fname_real, fname_imag
    FLOAT                 :: energy
    integer  :: np, ik, ist, il, ispin, s1, s2, k1, k2
    integer  :: green_real, green_imag, irow

    PUSH_SUB(states_init_self_energy)

    ! Calculate self-energy of the leads.
    ! FIXME: For spinors, this calculation is almost certainly wrong.
    ASSERT(st%ob_nst == st%nst)
    ASSERT(st%ob_d%nik == st%d%nik)

    s1 = st%st_start
    s2 = st%st_end
    k1 = st%d%kpt%start
    k2 = st%d%kpt%end

    do il = 1, NLEADS
      np = gr%intf(il)%np_intf
      SAFE_ALLOCATE(st%ob_lead(il)%self_energy(1:np, 1:np, 1:nspin, s1:s2, k1:k2))
    end do
    call messages_print_stress(stdout, "Lead self-energy")
    message(1) = ' st#     k#  Spin      Lead     Energy'
    call messages_info(1)
#ifdef HAVE_MPI 
    ! wait for all processors to finish 
    if(st%d%kpt%parallel) then 
      call MPI_Barrier(st%d%kpt%mpi_grp%comm, mpi_err) 
    end if
#endif
    do ik = k1, k2
      do ist = s1, s2
        energy = st%ob_eigenval(ist, ik)
        do il = 1, NLEADS
          np = gr%intf(il)%np_intf
          do ispin = 1, nspin
            select case(d_ispin)
            case(UNPOLARIZED)
              spin = '--'
            case(SPIN_POLARIZED)
              if(is_spin_up(ik)) then
                spin = 'up'
              else
                spin = 'dn'
              end if
              ! This is nonsense, but at least all indices are present.
            case(SPINORS)
              if(ispin .eq. 1) then
                spin = 'up'
              else
                spin = 'dn'
              end if
            end select
            write(message(1), '(i4,3x,i4,3x,a2,5x,a6,1x,f12.6)') ist, ik, &
              trim(spin), trim(LEAD_NAME(il)), energy
            call messages_info(1)

            call lead_self_energy(energy, lead(il)%h_diag(:, :, ispin), lead(il)%h_offdiag(:, :), &
              gr%intf(il), st%ob_lead(il)%self_energy(:, :, ispin, ist, ik), .true.)

            ! Write the entire self-energy to a file.
            if(in_debug_mode) then
              call io_mkdir('debug/open_boundaries')
              write(fname_real, '(3a,i4.4,a,i3.3,a,i1.1,a)') 'debug/open_boundaries/self-energy-', &
                trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', ispin, '.real'
              write(fname_imag, '(3a,i4.4,a,i3.3,a,i1.1,a)') 'debug/open_boundaries/self-energy-', &
                trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', ispin, '.imag'
              green_real = io_open(fname_real, action='write', grp=st%d%kpt%mpi_grp, is_tmp=.false.)
              green_imag = io_open(fname_imag, action='write', grp=st%d%kpt%mpi_grp, is_tmp=.false.)

              write(fmt, '(a,i6,a)') '(', np, 'e24.16)'
              do irow = 1, np
                write(green_real, fmt) real(st%ob_lead(il)%self_energy(irow, :, ispin, ist, ik))
                write(green_imag, fmt) aimag(st%ob_lead(il)%self_energy(irow, :, ispin, ist, ik))
              end do
              call io_close(green_real)
              call io_close(green_imag)
            end if
          end do
        end do
      end do
    end do
    call messages_print_stress(stdout)

    POP_SUB(states_init_self_energy)
  end subroutine states_init_self_energy


  ! ---------------------------------------------------------
  ! write H_(C,apha)*Psi_(alpha) without using Psi_(alpha)
  subroutine states_write_proj_lead_wf(sb, dir, intf, st)
    type(simul_box_t), intent(in) :: sb
    character(len=*),  intent(in) :: dir   ! directory
    type(interface_t), intent(in) :: intf(:)
    type(states_t),    intent(in) :: st

    integer :: ik, ist, idim, il, np, ip, iunit
    CMPLX, allocatable :: psi(:, :), phi(:, :), hpsi(:, :), self_energy(:, :)
    character(len=256) :: fname
    FLOAT :: kpoint(1:MAX_DIM)

    PUSH_SUB(states_write_proj_lead_wf)

    np = maxval(intf(1:NLEADS)%np_intf)

    SAFE_ALLOCATE(psi(1:np, 1:st%d%dim))
    SAFE_ALLOCATE(phi(1:np, 1:st%d%dim))
    SAFE_ALLOCATE(hpsi(1:np, 1:st%d%dim))
    SAFE_ALLOCATE(self_energy(1:np, 1:np))

    hpsi(:, :) = M_z0

    if(mpi_grp_is_root(mpi_world)) call io_mkdir(dir, is_tmp=.true.)
#ifdef HAVE_MPI
    if(st%parallel_in_states.or.st%d%kpt%parallel) call MPI_Barrier(st%d%kpt%mpi_grp%comm, mpi_err)
#endif

    do ik = st%d%kpt%start, st%d%kpt%end

      kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik))

      do ist = st%st_start, st%st_end
        do il = 1, NLEADS
          np = intf(il)%np_intf
          forall(ip = 1:np)
            psi(ip, :) = st%zpsi(intf(il)%index(ip), :, ist, ik)
            phi(ip, :) = st%zphi(intf(il)%index(ip), :, ist, ik)
          end forall
          do idim = 1, st%d%dim
            self_energy(1:np, 1:np) = st%ob_lead(il)%self_energy(1:np, 1:np, idim, ist, ik)
            call lalg_gemv(np, np, M_z1, self_energy(1:np, 1:np), psi(1:np, idim), M_z0, hpsi(1:np, idim))
            if((il.eq.LEFT).and.(-kpoint(1).gt.M_ZERO) .or. (il.eq.RIGHT).and.(-kpoint(1).lt.M_ZERO)) then
              ! add the reflecting part
              self_energy(1:np, 1:np) = transpose(conjg(self_energy(1:np, 1:np))) - self_energy(1:np, 1:np)
              call lalg_gemv(np, np, M_z1, self_energy(1:np, 1:np), phi(1:np, idim), M_z1, hpsi(1:np, idim))
            end if
            write(fname, '(3a,i4.4,a,i3.3,a,i1.1)') 'src0-', trim(LEAD_NAME(il)), '-', ist, '-', ik, '-', idim
            iunit = io_open(trim(dir)//trim(fname), action='write', form='unformatted', is_tmp=.true.)
            if(iunit.lt.0) then
              message(1) = 'Cannot write term for source term to file.'
              call messages_warning(1)
              call io_close(iunit)
              POP_SUB(states_write_proj_lead_wf)
              return
            end if
            ! Write parameters.
            write(iunit) np
            ! Write matrix.
            write(iunit) hpsi(1:np, idim)
            call io_close(iunit)
          end do
        end do
      end do
    end do

    SAFE_DEALLOCATE_A(psi)
    SAFE_DEALLOCATE_A(phi)
    SAFE_DEALLOCATE_A(hpsi)
    SAFE_DEALLOCATE_A(self_energy)

    POP_SUB(states_write_proj_lead_wf)
  end subroutine states_write_proj_lead_wf

  ! ---------------------------------------------------------
  subroutine states_read_proj_lead_wf(dir, intf, st, src0)
    character(len=*),  intent(in)    :: dir   ! directory
    type(interface_t), intent(in)    :: intf
    type(states_t),    intent(in)    :: st
    CMPLX,             intent(inout) :: src0(:, : ,:, :)

    integer :: ik, ist, idim, np, iunit
    character(len=256) :: fname

    PUSH_SUB(states_read_proj_lead_wf)

    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
          ! Try to open file.
          write(fname, '(3a,i4.4,a,i3.3,a,i1.1)') 'src0-', trim(LEAD_NAME(intf%il)), '-', ist, '-', ik, '-', idim
          iunit = io_open(trim(dir)//trim(fname), action='read', status='old', die=.false., is_tmp=.true., form='unformatted')
          if(iunit.lt.0) then ! no file found
            message(1) = 'Cannot read src(0) from file.'
            call messages_fatal(1)
          end if

          ! Now read the data.
          read(iunit) np

          if(np.ne.size(src0, 1)) then
            message(1) = 'Size mismatch! Cannot read src(0) from file.'
            call messages_fatal(1)
          end if

          ! because we use a sliced array we have to remap the index
          read(iunit) src0(1:np, idim, ist-st%st_start+1, ik-st%d%kpt%start+1)

          call io_close(iunit)
        end do
      end do
    end do

    POP_SUB(states_read_proj_lead_wf)
  end subroutine states_read_proj_lead_wf

end module states_io_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
