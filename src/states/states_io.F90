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

#include "global.h"

module states_io_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use kpoints_oct_m
  use lalg_basic_oct_m
  use loct_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m ! if not before parser_m, ifort 11.072 can`t compile with MPI2
  use mpi_lib_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use sort_oct_m
  use states_oct_m
  use states_dim_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m

  implicit none

  private

  public ::                           &
    states_write_eigenvalues,         &
    states_write_tpa,                 &
    states_write_bandstructure

contains

  ! ---------------------------------------------------------

  subroutine states_write_eigenvalues(iunit, nst, st, sb, error, st_start, compact)
    integer,           intent(in) :: iunit, nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb
    FLOAT, optional,   intent(in) :: error(:,:) !< (nst, st%d%nik)
    integer, optional, intent(in) :: st_start
    logical, optional, intent(in) :: compact
    
    integer :: ik, ikk, ist, ns, is, idir, st_start_, iflat, iqn, homo_index, not_printed
    logical :: print_eigenval
    FLOAT :: kpoint(1:MAX_DIM), max_error
    character(len=120) :: tmp_str(max(MAX_DIM, 3)), cspin

    FLOAT, allocatable :: flat_eigenval(:)
    integer, allocatable :: flat_indices(:, :)
    integer, parameter :: print_range = 8

    PUSH_SUB(states_write_eigenvalues)

    if(.not. st%calc_eigenval) then
      POP_SUB(states_write_eigenvalues)
      return
    end if
    
    st_start_ = 1
    if(present(st_start)) st_start_ = st_start
    ASSERT(nst <= st%nst)


    if(.not. mpi_grp_is_root(mpi_world)) then
      POP_SUB(states_write_eigenvalues)
      return
    end if

    if(.not. optional_default(compact, .false.)) then

      ns = 1
      if(st%d%nspin == 2) ns = 2
      
      message(1) = 'Eigenvalues [' // trim(units_abbrev(units_out%energy)) // ']'
      call messages_info(1, iunit)

      if(st%d%ispin  ==  SPINORS) then
        write(message(1), '(a4,1x,a5,1x,a12,1x,a12,2x,a4,4x,a4,4x,a4)')   &
          '#st',' Spin',' Eigenvalue', 'Occupation ', '<Sx>', '<Sy>', '<Sz>'
      else
        if(st%cmplxscl%space) then
          write(message(1), '(a4,1x,a5,1x,a12,1x,a15,4x,a12)')   &
            '#st',' Spin',' Eigenvalue', ' Im(Eigenvalue)', 'Occupation'
        else
          write(message(1), '(a4,1x,a5,1x,a12,4x,a12)')       &
            '#st',' Spin',' Eigenvalue', 'Occupation'
        end if
      end if
      if(present(error)) &
        write(message(1),'(a,a10)') trim(message(1)), ' Error'
      call messages_info(1, iunit)

      do ik = 1, st%d%nik, ns
        if(simul_box_is_periodic(sb)) then
          ikk = states_dim_get_kpoint_index(st%d, ik)
          kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, ikk, absolute_coordinates = .false.)
          write(message(1), '(a,i4,a)') '#k =', ikk, ', k = ('
          do idir = 1, sb%dim
            write(tmp_str(1), '(f10.6)') kpoint(idir)
            message(1) = trim(message(1))//trim(tmp_str(1))
            if(idir < sb%dim) message(1) = trim(message(1))//','
          end do
          message(1) = trim(message(1))//')'
          call messages_info(1, iunit)
        end if

        do ist = st_start_, nst
          do is = 0, ns-1
            if(is  ==  0) cspin = 'up'
            if(is  ==  1) cspin = 'dn'
            if(st%d%ispin  ==  UNPOLARIZED .or. st%d%ispin  ==  SPINORS) cspin = '--'

            write(tmp_str(1), '(i4,3x,a2)') ist, trim(cspin)
            if(st%d%ispin == SPINORS) then
              write(tmp_str(2), '(1x,f12.6,5x,f5.2,3x,3f8.4)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik)), st%occ(ist, ik), st%spin(1:3, ist, ik)
              if(present(error)) write(tmp_str(3), '(a3,es7.1,a1)')'  (', error(ist, ik), ')'
            else
              if(st%cmplxscl%space) then !cmplxscl
                write(tmp_str(2), '(1x,f12.6,3x,f12.6,3x,f12.6)') &
                  units_from_atomic(units_out%energy, st%zeigenval%Re(ist, ik+is)), &
                  units_from_atomic(units_out%energy, st%zeigenval%Im(ist, ik+is)), st%occ(ist, ik+is)
              else
                write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                  units_from_atomic(units_out%energy, st%eigenval(ist, ik+is)), st%occ(ist, ik+is)
              end if
              if(present(error)) write(tmp_str(3), '(a7,es7.1,a1)')'      (', error(ist, ik+is), ')'
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

    else

      call messages_info(1, iunit)

      SAFE_ALLOCATE(flat_eigenval(1:st%d%nik*nst))
      SAFE_ALLOCATE(flat_indices(1:2, 1:st%d%nik*nst))
      
      iflat = 1
      do iqn = 1, st%d%nik
        do ist = 1, nst
          
          flat_eigenval(iflat) = st%eigenval(ist, iqn)
          flat_indices(1:2, iflat) = (/iqn, ist/)          

          iflat = iflat + 1
        end do
      end do
    
      call sort(flat_eigenval, flat_indices(:, :))

      homo_index = st%d%nik*nst
      do iflat = 1, st%d%nik*nst
        iqn = flat_indices(1, iflat)
        ist = flat_indices(2, iflat)
        if(abs(st%occ(ist, iqn)) < CNST(0.1)) then
          homo_index = iflat - 1
          exit
        end if
      end do
      
      tmp_str(1) = '#  State'

      if(sb%periodic_dim > 0) tmp_str(1) = trim(tmp_str(1))//'  KPoint'

      if(st%d%ispin  ==  SPIN_POLARIZED) tmp_str(1) = trim(tmp_str(1))//'  Spin'

      tmp_str(1) = trim(tmp_str(1))//'  Eigenvalue ['// trim(units_abbrev(units_out%energy)) // ']'

      if(st%cmplxscl%space) then
        tmp_str(1) = trim(tmp_str(1))//'  Im(Eigenvalue)'
      end if
        
      tmp_str(1) = trim(tmp_str(1))//'  Occupation'

      if(st%d%ispin  ==  SPINORS) then
        tmp_str(1) = trim(tmp_str(1))//'      <Sx>     <Sy>     <Sz>'
      end if
      
      if(present(error)) tmp_str(1) = trim(tmp_str(1))//'    Error'

      call messages_write(tmp_str(1))
      call messages_info(iunit = iunit)

      not_printed = 0
      max_error = CNST(0.0)
      do iflat = 1, st%d%nik*nst
        iqn = flat_indices(1, iflat)
        ist = flat_indices(2, iflat)
        ik = states_dim_get_kpoint_index(st%d, iqn)
        is = states_dim_get_spin_index(st%d, iqn)

        print_eigenval = iflat <= print_range
        print_eigenval = print_eigenval .or. st%d%nik*nst - iflat < print_range
        print_eigenval = print_eigenval .or. abs(iflat - homo_index) <= print_range

        if(print_eigenval) then
          
          if(not_printed > 0) then
            call messages_write('')
            call messages_new_line()
            call messages_write('  [output of ')
            call messages_write(not_printed)
            call messages_write(' eigenvalues skipped')
            if(present(error)) then
              call messages_write(': maximum error =')
              call messages_write(max_error, fmt = '(es7.1)', align_left = .true.)
            end if
            call messages_write(']')
            call messages_new_line()
            call messages_write('')
            call messages_info(iunit = iunit)

            not_printed = 0
            max_error = CNST(0.0)

          end if

          write(tmp_str(1), '(i7)') ist

          if(sb%periodic_dim > 0) then
            write(tmp_str(1), '(2a,i7)') trim(tmp_str(1)), ' ', ik
          end if

          if(st%d%ispin  ==  SPIN_POLARIZED) then
            if(is  ==  1) cspin = '   up'
            if(is  ==  2) cspin = '   dn'
            write(tmp_str(1), '(2a,a5)') trim(tmp_str(1)), ' ', cspin
          end if


          if(.not. st%cmplxscl%space) then

            if(len(units_abbrev(units_out%energy)) == 1) then
              write(tmp_str(1), '(2a,f14.6)') trim(tmp_str(1)), ' ', units_from_atomic(units_out%energy, st%eigenval(ist, iqn))
            else
              write(tmp_str(1), '(2a,f15.6)') trim(tmp_str(1)), ' ', units_from_atomic(units_out%energy, st%eigenval(ist, iqn))
            end if

          else

            if(len(units_abbrev(units_out%energy)) == 1) then
              write(tmp_str(1), '(2a,f14.6)') trim(tmp_str(1)), ' ', units_from_atomic(units_out%energy, st%zeigenval%Re(ist, iqn))
            else
              write(tmp_str(1), '(2a,f15.6)') trim(tmp_str(1)), ' ', units_from_atomic(units_out%energy, st%zeigenval%Re(ist, iqn))
            end if

            write(tmp_str(1), '(2a,f15.6)') trim(tmp_str(1)), ' ', units_from_atomic(units_out%energy, st%zeigenval%Im(ist, iqn))

          end if

          write(tmp_str(1), '(2a,f11.6)') trim(tmp_str(1)), ' ', st%occ(ist, iqn)

          if(st%d%ispin  ==  SPINORS) then
            write(tmp_str(1), '(2a,3f9.4)') trim(tmp_str(1)), ' ', st%spin(1:3, ist, iqn)
          end if
          
          if(present(error)) then 
            write(tmp_str(1), '(2a,es7.1,a)') trim(tmp_str(1)), '   (', error(ist, iqn), ')'
          end if

          call messages_write(tmp_str(1))
          call messages_info(iunit = iunit)

        else
          
          not_printed = not_printed + 1

          if(present(error)) then
            max_error = max(max_error, error(ist, iqn))
          end if
          
        end if
        
      end do
      
      SAFE_DEALLOCATE_A(flat_indices)

      if(nst*st%d%nik > 1) call print_dos()

      SAFE_DEALLOCATE_A(flat_eigenval)

    end if

    if(st%smear%method /= SMEAR_SEMICONDUCTOR .and. st%smear%method /= SMEAR_FIXED_OCC) then
      write(message(1), '(a,f12.6,1x,a)') "Fermi energy = ", &
        units_from_atomic(units_out%energy, st%smear%e_fermi), units_abbrev(units_out%energy)
      call messages_info(1, iunit)
    end if

    POP_SUB(states_write_eigenvalues)
    
  contains

    subroutine print_dos()

      integer, parameter :: ndiv = 70, height = 10
      integer :: histogram(1:ndiv), iev, ien, iline, maxhist, ife
      character(len=ndiv) :: line
      FLOAT :: emin, emax, de
      
      PUSH_SUB(states_write_eigenvalues.print_dos)

      emin = flat_eigenval(1)
      emax = flat_eigenval(st%d%nik*nst)
      de = (emax - emin)/(ndiv - 1.0)

      if(de < M_EPSILON) then
        POP_SUB(states_write_eigenvalues.print_dos)
        return
      end if
      
      ife = nint((st%smear%e_fermi - emin)/de) + 1

      histogram = 0
      do iev = 1, st%d%nik*nst
        ien = nint((flat_eigenval(iev) - emin)/de) + 1
        ASSERT(ien >= 1)
        ASSERT(ien <= ndiv)
        histogram(ien) = histogram(ien) + 1
      end do

      !normalize
      if(maxval(histogram) > height) then
        maxhist = maxval(histogram)
        do ien = 1, ndiv
          histogram(ien) = (histogram(ien)*height)/maxhist
        end do
      end if

      call messages_new_line()
      call messages_write('Density of states:')
      call messages_new_line()
      call messages_info()
      
      !print histogram
      do iline = height, 1, -1
        do ien = 1, ndiv
          if(histogram(ien) >= iline) then
            call messages_write('%')
          else
            call messages_write('-')
          end if
        end do
        call messages_info()
      end do

      line(1:ndiv) = ' '
      line(ife:ife) = '^'
      call messages_write(line)
      call messages_new_line()
      call messages_info()

      POP_SUB(states_write_eigenvalues.print_dos)
    end subroutine print_dos
    
  end subroutine states_write_eigenvalues

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
    FLOAT, allocatable  :: qvector(:), psi_initial(:, :), psi_ist(:, :)

    PUSH_SUB(states_write_tpa)

    ! find the orbital with half-occupation
    tpa_initialst = -1
    do ist = 1, st%nst
      do ik = 1, st%d%nik
        if (abs(st%occ(ist,ik)-0.5)  <  M_THRESHOLD) then
          tpa_initialst = ist
          tpa_initialk  = ik
        end if
      end do
    end do

    ! make sure that half-occupancy was found
    if(tpa_initialst == -1) then
      if(mpi_grp_is_root(mpi_world)) then

        call messages_write('No orbital with half-occupancy found. TPA output is not written.')
        call messages_warning()

        POP_SUB(states_write_tpa)
        return
      end if
    end if

    !%Variable MomentumTransfer
    !%Type block
    !%Section Output
    !%Description
    !% Momentum-transfer vector <math>\vec{q}</math> to be used when calculating matrix elements
    !% <math>\left< f \left| e^{i \vec{q} \cdot \vec{r}} \right| i \right></math>.
    !% This enables the calculation of the dynamical structure factor,
    !% which is closely related to generalized oscillator strengths.
    !% If the vector is not given, but TPA output is requested (<tt>Output = TPA</tt>),
    !% only the oscillator strengths are written in the output file.
    !% For example, to use <math>\vec{q}</math> = (0.1, 0.2, 0.3), set
    !%
    !% <tt>%MomentumTransfer
    !% <br>&nbsp;&nbsp; 0.1 | 0.2 | 0.3
    !% <br>%</tt>
    !%End
    if(parse_block('MomentumTransfer', blk) == 0) then

      ! check if input makes sense
      ncols = parse_block_cols(blk, 0)

      if(ncols /= gr%mesh%sb%dim ) then ! wrong size

        if(mpi_grp_is_root(mpi_world)) then
          call messages_write('Inconsistent size of momentum-transfer vector. It will not be used in the TPA calculation.')
          call messages_warning()
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

    SAFE_ALLOCATE(psi_initial(1:gr%mesh%np, 1:st%d%dim))
    SAFE_ALLOCATE(psi_ist(1:gr%mesh%np, 1:st%d%dim))

    call states_get_state(st, gr%mesh, tpa_initialst, tpa_initialk, psi_initial)

    ! loop through every state
    do ist = 1, st%nst

      call states_get_state(st, gr%mesh, ist, tpa_initialk, psi_ist)

      ! final states are the unoccupied ones
      if (abs(st%occ(ist,tpa_initialk))  <  M_THRESHOLD) then

        osc_strength = M_ZERO
        transition_energy = st%eigenval(ist, tpa_initialk) - st%eigenval(tpa_initialst, tpa_initialk)

        ! dipole matrix elements <f|x|i> etc. -> oscillator strengths
        do icoord = 1, gr%mesh%sb%dim    ! for x,y,z

          ff(1:gr%mesh%np) = psi_initial(1:gr%mesh%np, 1)*gr%mesh%x(1:gr%mesh%np, icoord)* &
            psi_ist(1:gr%mesh%np, 1)
          osc(icoord)  = dmf_integrate(gr%mesh, ff)
          osc_strength = osc_strength + 2.0/real(gr%mesh%sb%dim)*transition_energy*abs(osc(icoord))**2.0

        end do

        ! matrix elements <f|exp(iq.r)|i> -> dynamic structure factor
        if (use_qvector) then

          cff(1:gr%mesh%np) = psi_initial(1:gr%mesh%np, 1)*psi_ist(1:gr%mesh%np, 1)

          do icoord = 1, gr%mesh%sb%dim    ! for x,y,z
            cff(1:gr%mesh%np) = cff(1:gr%mesh%np)*exp(M_zI*gr%mesh%x(1:gr%mesh%np, icoord)*qvector(icoord))
          end do

          dsf = abs(zmf_integrate(gr%mesh, cff))**2
        end if

        ! write oscillator strengths (+ dynamic structure factor if qvector is given) into file
        if(mpi_grp_is_root(mpi_world)) then

          if(use_qvector) then
            write(message(1), '(1x,6(es15.8,1x))') units_from_atomic(units_out%energy, transition_energy), osc(:), osc_strength, &
                                                   units_from_atomic(unit_one/units_out%energy, dsf)
          else
            write(message(1), '(1x,6(es15.8,1x))') units_from_atomic(units_out%energy, transition_energy), osc(:), osc_strength
          end if

          call messages_info(1,iunit)

        end if

      end if

    end do

    ! finally close the file
    if(mpi_grp_is_root(mpi_world)) then
      call io_close(iunit)
    end if

    SAFE_DEALLOCATE_A(psi_initial)
    SAFE_DEALLOCATE_A(psi_ist)

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
  subroutine states_write_fermi_for_bands(dir, st, sb)
    character(len=*),  intent(in) :: dir
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb

    integer :: iunit, idir
    character(len=100) :: str_tmp

    PUSH_SUB(states_write_fermi_for_bands)

    iunit = io_open(trim(dir)//'/'//'bands-efermi.dat', action='write')    

    write(message(1), '(3a)') '# Fermi energy [', trim(units_abbrev(units_out%energy)), &
      '] in a format compatible with bands-gp.dat'

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
    end do
    do idir = 1, sb%dim
      write(str_tmp, '(f12.6)') minval(sb%kpoints%reduced%red_point(idir, 1:sb%kpoints%reduced%npoints))
      message(2) = trim(message(2)) // trim(str_tmp)
      write(str_tmp, '(f12.6)') M_ZERO     ! Gamma point
      message(3) = trim(message(3)) // trim(str_tmp)
      write(str_tmp, '(f12.6)') maxval(sb%kpoints%reduced%red_point(idir, 1:sb%kpoints%reduced%npoints))
      message(4) = trim(message(4)) // trim(str_tmp)
    end do
    write(str_tmp, '(f12.6)') units_from_atomic(units_out%energy, st%smear%e_fermi)
    message(2) = trim(message(2)) // trim(str_tmp)
    message(3) = trim(message(3)) // trim(str_tmp)
    message(4) = trim(message(4)) // trim(str_tmp)

    call messages_info(4, iunit)
    call io_close(iunit)

    POP_SUB(states_write_fermi_for_bands)
  end subroutine states_write_fermi_for_bands


    ! ---------------------------------------------------------

  subroutine states_write_bandstructure(dir, nst, st, sb)
    character(len=*),  intent(in) :: dir
    integer,           intent(in) :: nst
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb

    integer :: idir, ist, ik, ns, is,npath
    integer, allocatable :: iunit(:)
    FLOAT   :: red_kpoint(1:MAX_DIM)
    character(len=80) :: filename

    if(.not. mpi_grp_is_root(mpi_world)) return

    PUSH_SUB(states_write_bandstructure)   

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    SAFE_ALLOCATE(iunit(0:ns-1))

    do is = 0, ns-1
      if (ns > 1) then
        write(filename, '(a,i1.1,a)') 'bandstructure-sp', is+1
      else
        write(filename, '(a)') 'bandstructure'
      end if
      iunit(is) = io_open(trim(dir)//'/'//trim(filename), action='write')    

      ! write header
      write(iunit(is),'(a)',advance='no') '# coord. '
      do idir = 1, sb%dim
        write(iunit(is),'(3a)',advance='no') 'k', index2axis(idir), ' '
      end do
      write(iunit(is),'(a,i6,3a)') '(red. coord.), bands:', nst, ' [', trim(units_abbrev(units_out%energy)), ']'
    end do

    npath = SIZE(sb%kpoints%coord_along_path)*ns

    ! output bands
    do ik = st%d%nik-npath+1, st%d%nik, ns
      do is = 0, ns - 1
        red_kpoint(1:sb%dim) = kpoints_get_point(sb%kpoints, states_dim_get_kpoint_index(st%d, ik + is), &
                                                   absolute_coordinates=.false.)
        write(iunit(is),'(1x)',advance='no')
        write(iunit(is),'(f14.8)',advance='no') kpoints_get_path_coord(sb%kpoints, & 
                                                    states_dim_get_kpoint_index(st%d, ik + is)-(st%d%nik -npath)) 
        do idir = 1, sb%dim
          write(iunit(is),'(f14.8)',advance='no') red_kpoint(idir)
        end do
        do ist = 1, nst
          write(iunit(is),'(1x,f14.8)',advance='no') units_from_atomic(units_out%energy, st%eigenval(ist, ik + is))
        end do
      end do
      do is = 0, ns-1
        write(iunit(is), '(a)')
      end do
    end do
    do is = 0, ns-1
      call io_close(iunit(is))
    end do

    SAFE_DEALLOCATE_A(iunit)


  POP_SUB(states_write_bandstructure)
    
  end subroutine states_write_bandstructure

end module states_io_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
