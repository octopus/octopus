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
  use atomic_orbital_oct_m
  use comm_oct_m
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use io_oct_m
  use kpoints_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use orbitalset_oct_m
  use orbitalset_utils_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use smear_oct_m
  use sort_oct_m
  use species_oct_m
  use states_oct_m
  use states_dim_oct_m
  use submesh_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m

  implicit none

  private

  public ::                           &
    states_write_eigenvalues,         &
    states_write_tpa,                 &
    states_write_bandstructure,       &
    states_write_gaps

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
        write(message(1), '(a4,1x,a5,1x,a12,4x,a12)')       &
          '#st',' Spin',' Eigenvalue', 'Occupation'
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
              write(tmp_str(2), '(1x,f12.6,3x,f12.6)') &
                units_from_atomic(units_out%energy, st%eigenval(ist, ik+is)), st%occ(ist, ik+is)
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


          if(len(units_abbrev(units_out%energy)) == 1) then
            write(tmp_str(1), '(2a,f14.6)') trim(tmp_str(1)), ' ', units_from_atomic(units_out%energy, st%eigenval(ist, iqn))
          else
            write(tmp_str(1), '(2a,f15.6)') trim(tmp_str(1)), ' ', units_from_atomic(units_out%energy, st%eigenval(ist, iqn))
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
      
      if(nst*st%d%nik > 1) call print_dos()

      SAFE_DEALLOCATE_A(flat_indices)
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
      FLOAT :: dhistogram(1:ndiv)
      character(len=ndiv) :: line
      FLOAT :: emin, emax, de
      
      PUSH_SUB(states_write_eigenvalues.print_dos)

      emin = flat_eigenval(1)
      emax = flat_eigenval(st%d%nik*nst)
      de = (emax - emin)/(ndiv - M_ONE)

      if(de < M_EPSILON) then
        POP_SUB(states_write_eigenvalues.print_dos)
        return
      end if
      
      ife = nint((st%smear%e_fermi - emin)/de) + 1

      histogram = 0
      dhistogram = M_ZERO
      do iev = 1, st%d%nik*nst
        ien = nint((flat_eigenval(iev) - emin)/de) + 1
        ASSERT(ien >= 1)
        ASSERT(ien <= ndiv)
        dhistogram(ien) = dhistogram(ien) + st%d%kweights(flat_indices(1, iev))*sb%kpoints%full%npoints
      end do

      !normalize
      if(maxval(dhistogram) > real(height)) then
        maxhist = nint(maxval(dhistogram))
        do ien = 1, ndiv
          histogram(ien) = nint((dhistogram(ien)*height)/maxhist)
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
  subroutine states_write_gaps(iunit, st, sb)
    integer,           intent(in) :: iunit
    type(states_t),    intent(in) :: st
    type(simul_box_t), intent(in) :: sb
    
    integer :: ik, ikk, ist

    FLOAT :: homo, lumo, egdir, egindir
    integer :: homok, lumok, egdirk

    PUSH_SUB(states_write_gaps)

    if(.not. st%calc_eigenval .or. .not. mpi_grp_is_root(mpi_world) &
     .or. .not. simul_box_is_periodic(sb) .or. st%smear%method  /=  SMEAR_SEMICONDUCTOR) then 
      POP_SUB(states_write_gaps)
      return
    end if

    homo = -M_HUGE
    homok = -1
    lumo =  M_HUGE 
    lumok = -1
    egdir = M_HUGE
    egdirk = -1
    
    !TODO: This definition depends on the type of smearing
    do ik = 1, st%d%nik
      if(abs(st%d%kweights(ik)) < M_EPSILON) cycle
      do ist = 1,st%nst-1
        if(st%occ(ist,ik) > M_EPSILON .and. st%eigenval(ist,ik) > homo) then
          homo = st%eigenval(ist,ik)
          homok = ik
        end if

        if(st%occ(ist+1,ik) <= M_EPSILON .and. st%eigenval(ist+1,ik) < lumo) then
          lumo = st%eigenval(ist+1,ik)
          lumok = ik
        end if

        if(st%occ(ist,ik) > M_EPSILON .and. st%occ(ist+1,ik) <= M_EPSILON &
          .and. (st%eigenval(ist+1,ik)-st%eigenval(ist,ik))< egdir) then
          egdir = (st%eigenval(ist+1,ik)-st%eigenval(ist,ik))
          egdirk = ik 
        end if
      end do
    end do

    egindir = lumo-homo

    if(lumo == -1 .or. egdir <= M_EPSILON) then
      write(message(1),'(a)') 'The system seems to have no gap.'
      call messages_info(1, iunit)
    else
      write(message(1),'(a,i5,a,f7.4,a)') 'Direct gap at ik=', egdirk, ' of ', &
              units_from_atomic(units_out%energy, egdir),  ' ' // trim(units_abbrev(units_out%energy))
      write(message(2),'(a,i5,a,i5,a,f7.4,a)') 'Indirect gap between ik=', homok, ' and ik=', lumok, &
            ' of ', units_from_atomic(units_out%energy, egindir),  ' ' // trim(units_abbrev(units_out%energy))
      call messages_info(2, iunit)
    end if


    POP_SUB(states_write_gaps)
  end subroutine states_write_gaps


  ! ---------------------------------------------------------
  subroutine states_write_tpa(dir, namespace, gr, st)
    character(len=*),  intent(in) :: dir
    type(namespace_t), intent(in) :: namespace
    type(grid_t),      intent(in) :: gr
    type(states_t),    intent(in) :: st

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
    if(parse_block(namespace, 'MomentumTransfer', blk) == 0) then

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

      iunit = io_open_old(trim(dir)//'/'//trim('tpa_xas'), action='write')    

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
          osc_strength = osc_strength + CNST(2.0)/real(gr%mesh%sb%dim)*transition_energy*abs(osc(icoord))**2

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

  subroutine states_write_bandstructure(dir, namespace, nst, st, sb, geo, mesh, phase, vec_pot, vec_pot_var)
    character(len=*),         intent(in)      :: dir
    type(namespace_t),        intent(in)      :: namespace
    integer,                  intent(in)      :: nst
    type(states_t),           intent(in)      :: st
    type(simul_box_t),        intent(in)      :: sb
    type(geometry_t), target, intent(in)      :: geo
    type(mesh_t),             intent(in)      :: mesh
    CMPLX, pointer                            :: phase(:, :)
    FLOAT, optional, allocatable, intent(in)  :: vec_pot(:) !< (sb%dim)
    FLOAT, optional, allocatable, intent(in)  :: vec_pot_var(:, :) !< (1:sb%dim, 1:ns)

    integer :: idir, ist, ik, ns, is,npath
    integer, allocatable :: iunit(:)
    FLOAT   :: red_kpoint(1:MAX_DIM)
    character(len=80) :: filename

    logical :: projection
    integer :: ii, ll, mm, nn, work, norb, work2
    integer :: ia, iorb, idim, maxnorb
    FLOAT   :: norm
    FLOAT, allocatable :: dpsi(:,:), ddot(:,:)
    CMPLX, allocatable :: zpsi(:,:), zdot(:,:)
    FLOAT, allocatable :: weight(:,:,:,:,:)
    type(orbitalset_t) :: os


    PUSH_SUB(states_write_bandstructure)   

    ! shortcuts
    ns = 1
    if(st%d%nspin == 2) ns = 2

    !%Variable BandStructureComputeProjections
    !%Type logical
    !%Default false
    !%Section Output
    !%Description
    !% Determines if projections of wavefunctions on the atomic orbitals 
    !% are computed or not for obtaining the orbital resolved band-structure.
    !%End
    call parse_variable(namespace, 'BandStructureComputeProjections', .false., projection)


    if(mpi_grp_is_root(mpi_world)) then
      SAFE_ALLOCATE(iunit(0:ns-1))

      do is = 0, ns-1
        if (ns > 1) then
          write(filename, '(a,i1.1,a)') 'bandstructure-sp', is+1
        else
          write(filename, '(a)') 'bandstructure'
        end if
        iunit(is) = io_open_old(trim(dir)//'/'//trim(filename), action='write')    

        ! write header
        write(iunit(is),'(a)',advance='no') '# coord. '
        do idir = 1, sb%dim
          write(iunit(is),'(3a)',advance='no') 'k', index2axis(idir), ' '
        end do
        if(.not.projection) then
          write(iunit(is),'(a,i6,3a)') '(red. coord.), bands:', nst, ' [', trim(units_abbrev(units_out%energy)), ']'
        else 
         write(iunit(is),'(a,i6,3a)',advance='no') '(red. coord.), bands:', nst, ' [', trim(units_abbrev(units_out%energy)), '] '
         do ia = 1, geo%natoms
            work = orbitalset_utils_count(geo, ia)
            do norb = 1, work
             work2 = orbitalset_utils_count(geo, ia, norb)
              write(iunit(is),'(a, i3.3,a,i1.1,a)',advance='no') 'w(at=',ia,',os=',norb,') '
            end do
          end do
          write(iunit(is),'(a)') ''
        end if
      end do
    end if
 
    npath = SIZE(sb%kpoints%coord_along_path)*ns


    !We need to compute the projections of each wavefunctions on the localized basis
    if(projection) then    

      if(states_are_real(st)) then
        SAFE_ALLOCATE(dpsi(1:mesh%np, 1:st%d%dim))
      else
        SAFE_ALLOCATE(zpsi(1:mesh%np, 1:st%d%dim))
      end if

      maxnorb = 0
      do ia = 1, geo%natoms
        maxnorb = max(maxnorb, orbitalset_utils_count(geo, ia))
      end do

      SAFE_ALLOCATE(weight(1:st%d%nik,1:st%nst, 1:maxnorb, 1:MAX_L, 1:geo%natoms))
      weight(1:st%d%nik,1:st%nst, 1:maxnorb, 1:MAX_L, 1:geo%natoms) = M_ZERO
 
      do ia = 1, geo%natoms

        !We first count how many orbital set we have
        work = orbitalset_utils_count(geo, ia)

        !We loop over the orbital sets of the atom ia
        do norb = 1, work
          call orbitalset_nullify(os)

          !We count the orbitals
          work2 = 0
          do iorb = 1, species_niwfs(geo%atom(ia)%species)
           call species_iwf_ilm(geo%atom(ia)%species, iorb, 1, ii, ll, mm)
            call species_iwf_n(geo%atom(ia)%species, iorb, 1, nn )
            if(ii == norb) then
              os%ll = ll
              os%nn = nn
              os%ii = ii
              os%radius = atomic_orbital_get_radius(geo, mesh, ia, iorb, 1, OPTION__AOTRUNCATION__AO_FULL, CNST(0.01))
              work2 = work2 + 1
            end if
          end do
          os%norbs = work2
          os%ndim = 1
          os%submeshforperiodic = .false.
          os%spec => geo%atom(ia)%species
          call submesh_null(os%sphere)

          do iorb = 1, os%norbs
            ! We obtain the orbital
            if(states_are_real(st)) then
              call dget_atomic_orbital(geo, mesh, os%sphere, ia, os%ii, os%ll, os%jj, &
                                                os, iorb, os%radius, os%ndim)
              norm = M_ZERO
              do idim = 1, os%ndim
                norm = norm + dsm_nrm2(os%sphere, os%dorb(1:os%sphere%np,idim,iorb))**2
              end do
              norm = sqrt(norm)
              do idim = 1, os%ndim
                os%dorb(1:os%sphere%np,idim,iorb) =  os%dorb(1:os%sphere%np,idim,iorb)/norm
              end do
            else
              call zget_atomic_orbital(geo, mesh, os%sphere, ia, os%ii, os%ll, os%jj, &
                                                os, iorb, os%radius, os%ndim)
              norm = M_ZERO
              do idim = 1, os%ndim
                norm = norm + zsm_nrm2(os%sphere, os%zorb(1:os%sphere%np,idim,iorb))**2
              end do
              norm = sqrt(norm)
              do idim = 1, os%ndim
                os%zorb(1:os%sphere%np,idim,iorb) =  os%zorb(1:os%sphere%np,idim,iorb)/norm
              end do
            end if
          end do !iorb

          nullify(os%phase)
          if(associated(phase)) then
            ! In case of complex wavefunction, we allocate the array for the phase correction
            SAFE_ALLOCATE(os%phase(1:os%sphere%np, st%d%kpt%start:st%d%kpt%end))
            os%phase(:,:) = M_ZERO
            if(simul_box_is_periodic(mesh%sb) .and. .not. os%submeshforperiodic) then
              SAFE_ALLOCATE(os%eorb_mesh(1:mesh%np, 1:os%norbs, 1:os%ndim, st%d%kpt%start:st%d%kpt%end))
              os%eorb_mesh(:,:,:,:) = M_ZERO
            else
              SAFE_ALLOCATE(os%eorb_submesh(1:os%sphere%np, 1:os%ndim, 1:os%norbs, st%d%kpt%start:st%d%kpt%end))
              os%eorb_submesh(:,:,:,:) = M_ZERO
            end if
            call orbitalset_update_phase(os, sb, st%d%kpt, (st%d%ispin==SPIN_POLARIZED), &
                              vec_pot, vec_pot_var)
          end if

          if(states_are_real(st)) then
            SAFE_ALLOCATE(ddot(1:st%d%dim,1:os%norbs))
          else
            SAFE_ALLOCATE(zdot(1:st%d%dim,1:os%norbs))
          end if

          do ist = st%st_start, st%st_end
           do ik = st%d%kpt%start, st%d%kpt%end
            if(ik < st%d%nik-npath+1 ) cycle ! We only want points inside the k-point path
            if(states_are_real(st)) then
              call states_get_state(st, mesh, ist, ik, dpsi )
              call dorbitalset_get_coefficients(os, st%d%dim, dpsi, ik, .false., .false., ddot(1:st%d%dim,1:os%norbs))
              do iorb = 1, os%norbs
                do idim = 1, st%d%dim
                  weight(ik,ist,iorb,norb,ia) = weight(ik,ist,iorb,norb,ia) + abs(ddot(idim,iorb))**2
                end do
              end do
            else
              call states_get_state(st, mesh, ist, ik, zpsi )
              if(associated(phase)) then
                ! Apply the phase that contains both the k-point and vector-potential terms.
                call states_set_phase(st%d, zpsi, phase(:,ik), mesh%np, .false.)
              end if
              call zorbitalset_get_coefficients(os, st%d%dim, zpsi, ik, associated(phase), .false.,&
                                 zdot(1:st%d%dim,1:os%norbs))
              do iorb = 1, os%norbs
                do idim = 1, st%d%dim
                  weight(ik,ist,iorb,norb,ia) = weight(ik,ist,iorb,norb,ia) + abs(zdot(idim,iorb))**2
                end do
              end do
            end if
          end do
         end do

         SAFE_DEALLOCATE_A(ddot)
         SAFE_DEALLOCATE_A(zdot)

         call orbitalset_end(os)
       end do !norb

       if(st%parallel_in_states .or. st%d%kpt%parallel) then
         call comm_allreduce(st%st_kpt_mpi_grp%comm, weight(1:st%d%nik,1:st%nst, 1:maxnorb, 1:MAX_L,ia))
       end if
     end do !ia

    SAFE_DEALLOCATE_A(dpsi)
    SAFE_DEALLOCATE_A(zpsi)

  end if  !projection

  if(mpi_grp_is_root(mpi_world)) then
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
        if(projection) then
          do ia = 1, geo%natoms
            work = orbitalset_utils_count(geo, ia)
            do norb = 1, work
              work2 = orbitalset_utils_count(geo, ia, norb)              
              do iorb = 1, work2
                do ist = 1, nst
                  write(iunit(is),'(es15.8)',advance='no') weight(ik+is,ist,iorb,norb,ia)
                end do
              end do
            end do
          end do
        end if
      end do

      do is = 0, ns-1
        write(iunit(is), '(a)')
      end do
    end do

    do is = 0, ns-1
      call io_close(iunit(is))
    end do
  end if

  if(projection) then
     SAFE_DEALLOCATE_A(weight)
  end if

  SAFE_DEALLOCATE_A(iunit)


  POP_SUB(states_write_bandstructure)
    
  end subroutine states_write_bandstructure

end module states_io_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
