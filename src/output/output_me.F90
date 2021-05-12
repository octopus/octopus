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

module output_me_oct_m
  use boundaries_oct_m
  use derivatives_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_elec_oct_m
  use io_function_oct_m
  use io_oct_m
  use ions_oct_m
  use kpoints_oct_m
  use loct_math_oct_m
  use lda_u_oct_m
  use mpi_oct_m
  use mpi_lib_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use projector_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use singularity_oct_m
  use space_oct_m
  use states_abst_oct_m
  use states_elec_oct_m
  use states_elec_calc_oct_m
  use states_elec_dim_oct_m
  use scissor_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  use xc_oct_m

  implicit none

  private
  public ::           &
    output_me_t,      &
    output_me_init,   &
    output_me

  type output_me_t
    private
    logical, public :: what(MAX_OUTPUT_TYPES)               !< what to output 
    !> If output_ksdipole, this number sets up which matrix elements will
    !! be printed: e.g. if ksmultipoles = 3, the dipole, quadrupole and 
    !! octopole matrix elements (between Kohn-Sham or single-particle orbitals).
    !! In 2D, only the dipole moments are printed.
    integer :: ks_multipoles      

    integer :: st_start !Start index for the output
    integer :: st_end   !Stop index for the output
    integer :: nst      !Number of states computed
  end type output_me_t

contains
  
  ! ---------------------------------------------------------
  subroutine output_me_init(this, namespace, space, st, nst)
    type(output_me_t),   intent(out) :: this
    type(namespace_t),   intent(in)  :: namespace
    type(space_t),       intent(in)  :: space
    type(states_elec_t), intent(in)  :: st
    integer,             intent(in)  :: nst

    integer(8) :: how(0:MAX_OUTPUT_TYPES)
    integer :: output_interval(0:MAX_OUTPUT_TYPES)

    PUSH_SUB(output_me_init)

    !%Variable OutputMatrixElements
    !%Type block
    !%Default none
    !%Section Output
    !%Description
    !% Specifies what matrix elements to print.
    !% Enabled only if <tt>Output</tt> block includes <tt>matrix_elements</tt>.
    !% The output files go into the <tt>static</tt> directory, except when
    !% running a time-dependent simulation, when the directory <tt>td.XXXXXXX</tt> is used.
    !%
    !% Example:
    !% <br><br><tt>%OutputMatrixElements
    !% <br>&nbsp;&nbsp;momentum
    !% <br>&nbsp;&nbsp;ks_multipoles
    !% <br>%<br></tt>
    !%
    !% It is possible to specify only compute the matrix elements for some of the states
    !% using the variables <tt>OutptMEStart</tt> and <tt>OutputMEEnd</tt>.
    !%Option momentum 1
    !% Momentum. Filename: <tt>ks_me_momentum</tt>.
    !%Option ang_momentum 2
    !% Dimensionless angular momentum <math>\vec{r} \times \vec{k}</math>. Filename: <tt>ks_me_angular_momentum</tt>.
    !%Option one_body 3
    !% <math>\left< i \left| \hat{T} + V_{ext} \right| j \right></math>. Not available with states parallelization.
    !%Option two_body 4
    !% <math>\left< ij \left| \frac{1}{\left|\vec{r}_1-\vec{r}_2\right|} \right| kl \right></math>.
    !% Not available with states parallelization.
    !% Not available with states parallelization. For periodic system, this is not available for k-point parallelization neither.
    !%Option two_body_exc_k 5
    !% <math>\left< n1-k1, n2-k2 \left| \frac{1}{\left|\vec{r}_1-\vec{r}_2\right|} \right| n2-k1 n1-k2 \right></math>.
    !% Not available with states parallelization. For periodic system, this is not available for k-point parallelization neither.
    !%Option ks_multipoles 6
    !% See <tt>OutputMEMultipoles</tt>. Not available with states parallelization.
    !%Option dipole 7
    !% Prints the dipole matrix elements. Not available with states parallelization.
    !% For periodic systems, the intraband terms (dipole matrix elements between degenerated states)
    !% are set to zero, and only the absolute value of the dipole matrix element is printed.
    !% Not yet supported for spinors.
    !%End

    this%what = .false.
    call io_function_read_what_how_when(namespace, space, this%what, how, output_interval, &
    'OutputMatrixElements')

    if (st%parallel_in_states) then
      if (this%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY)) &
        call messages_not_implemented("OutputMatrixElements=two_body is not implemented in states parallelization.", &
        namespace=namespace)
      if (this%what(OPTION__OUTPUTMATRIXELEMENTS__DIPOLE)) &
        call messages_not_implemented("OutputMatrixElements=dipole is not implemented in states parallelization.", &
        namespace=namespace)
    end if

    if (space%dim /= 2 .and. space%dim /= 3) this%what(OPTION__OUTPUTMATRIXELEMENTS__ANG_MOMENTUM) = .false.

    if (this%what(OPTION__OUTPUTMATRIXELEMENTS__KS_MULTIPOLES)) then
      !%Variable OutputMEMultipoles
      !%Type integer
      !%Default 1
      !%Section Output
      !%Description
      !% This variable decides which multipole moments are printed out for
      !% <tt>OutputMatrixElements = ks_multipoles</tt>:
      !%
      !% In 3D, if, for example, <tt>OutputMEMultipoles = 1</tt>, then the program will print three 
      !% files, <tt>ks_me_multipoles.x</tt> (<tt>x</tt>=1,2,3), containing
      !% respectively the (1,-1), (1,0) and (1,1) multipole matrix elements
      !% between Kohn-Sham states.
      !%
      !% In 2D, this variable is ignored: it will always print two files, 
      !% <tt>ks_me_multipoles.i</tt> (<tt>i</tt>=1,2), containing the <math>x</math> and
      !% <math>y</math> dipole matrix elements.
      !%
      !% In 1D, if, for example, <tt>OutputMEMultipoles = 2</tt>, the program will print two files, containing the
      !% <math>x</math> and <math>x^2</math> matrix elements between Kohn-Sham states.
      !%End
      call parse_variable(namespace, 'OutputMEMultipoles', 1, this%ks_multipoles)
    end if

    !%Variable OutputMEStart
    !%Type integer
    !%Default 1
    !%Section Output
    !%Description
    !% Specifies the state/band index for starting to compute the matrix element.
    !% So far, this is only used for dipole matrix elements.
    !%End
    call parse_variable(namespace, 'OutputMEStart', 1, this%st_start)
    ASSERT(this%st_start > 0 .and. this%st_start <= nst)

    !%Variable OutputMEEnd
    !%Type integer
    !%Default 1
    !%Section Output
    !%Description
    !% Specifies the highest state/band index used to compute the matrix element.
    !% So far, this is only used for dipole matrix elements.
    !%End
    call parse_variable(namespace, 'OutputMEEnd', nst, this%st_end)
    ASSERT(this%st_end > 0 .and. this%st_end <= nst)
    ASSERT(this%st_start <= this%st_end)
    this%nst = this%st_end - this%st_start +1

    POP_SUB(output_me_init)
  end subroutine output_me_init


  ! ---------------------------------------------------------
  subroutine output_me(this, namespace, space, dir, st, gr, ions, hm)
    type(output_me_t),        intent(in)    :: this
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=*),         intent(in)    :: dir
    type(states_elec_t),      intent(inout) :: st
    type(grid_t),             intent(in)    :: gr
    type(ions_t),             intent(in)    :: ions
    type(hamiltonian_elec_t), intent(inout) :: hm

    integer :: id, ll, mm, ik, iunit
    character(len=256) :: fname
    FLOAT, allocatable :: doneint(:), dtwoint(:)
    CMPLX, allocatable :: zoneint(:), ztwoint(:)
    integer, allocatable :: iindex(:,:), jindex(:,:), kindex(:,:), lindex(:,:)
    type(singularity_t) :: singul
    
    PUSH_SUB(output_me)

    if (this%what(OPTION__OUTPUTMATRIXELEMENTS__MOMENTUM)) then
      write(fname,'(2a)') trim(dir), '/ks_me_momentum'
      call output_me_out_momentum(fname, st, space, gr, namespace, hm%kpoints)
    end if

    if (this%what(OPTION__OUTPUTMATRIXELEMENTS__ANG_MOMENTUM)) then
      write(fname,'(2a)') trim(dir), '/ks_me_angular_momentum'
      call output_me_out_ang_momentum(fname, st, gr, namespace, hm%kpoints)
    end if

    if (this%what(OPTION__OUTPUTMATRIXELEMENTS__KS_MULTIPOLES)) then
      ! The content of each file should be clear from the header of each file.
      id = 1
      do ik = 1, st%d%nik
        select case(space%dim)
        case(3)
          do ll = 1, this%ks_multipoles
            do mm = -ll, ll
              write(fname,'(i4)') id
              write(fname,'(a)') trim(dir)//'/ks_me_multipoles.'//trim(adjustl(fname))
              if (states_are_real(st)) then
                call doutput_me_ks_multipoles(fname, namespace, st, gr, ll, mm, ik)
              else
                call zoutput_me_ks_multipoles(fname, namespace, st, gr, ll, mm, ik)
              end if

              id = id + 1
            end do
          end do
        case(2)
          do ll = 1, 2
            write(fname,'(i4)') id
            write(fname,'(a)') trim(dir)//'/ks_me_multipoles.'//trim(adjustl(fname))
            if (states_are_real(st)) then
              call doutput_me_ks_multipoles2d(fname, namespace, st, gr%mesh, ll, ik)
            else
              call zoutput_me_ks_multipoles2d(fname, namespace, st, gr%mesh, ll, ik)
            end if

            id = id + 1

          end do
        case(1)
          do ll = 1, this%ks_multipoles
            write(fname,'(i4)') id
            write(fname,'(a)') trim(dir)//'/ks_me_multipoles.'//trim(adjustl(fname))
            if (states_are_real(st)) then
              call doutput_me_ks_multipoles1d(fname, namespace, st, gr%mesh, ll, ik)
            else
              call zoutput_me_ks_multipoles1d(fname, namespace, st, gr%mesh, ll, ik)
            end if

            id = id + 1
          end do
        end select
      end do
    end if

    if (this%what(OPTION__OUTPUTMATRIXELEMENTS__DIPOLE)) then
      ASSERT(.not. st%parallel_in_states)
      ! The content of each file should be clear from the header of each file.
      do ik = st%d%kpt%start, st%d%kpt%end
        write(fname,'(i4)') ik
        write(fname,'(a)') trim(dir)//'/ks_me_dipole.k'//trim(adjustl(fname))//'_'
          if (states_are_real(st)) then
            call doutput_me_dipole(this, fname, namespace, space, st, gr, hm, ions, ik)
          else
            call zoutput_me_dipole(this, fname, namespace, space, st, gr, hm, ions, ik)
          end if
      end do
    end if


    if (this%what(OPTION__OUTPUTMATRIXELEMENTS__ONE_BODY)) then
      message(1) = "Computing one-body matrix elements"
      call messages_info(1)

      if (st%parallel_in_states) call messages_not_implemented("OutputMatrixElements=one_body with states parallelization", &
        namespace=namespace)
      if (st%d%kpt%parallel) call messages_not_implemented("OutputMatrixElements=one_body with k-points parallelization", &
        namespace=namespace)
      if (family_is_mgga_with_exc(hm%xc)) &
      call messages_not_implemented("OutputMatrixElements=one_body with MGGA", namespace=namespace)
      ! how to do this properly? states_elec_matrix
      iunit = io_open(trim(dir)//'/output_me_one_body', namespace, action='write')

      id = st%nst*(st%nst+1)/2

      SAFE_ALLOCATE(iindex(1:id,1:1))
      SAFE_ALLOCATE(jindex(1:id,1:1))

      if (states_are_real(st)) then
        SAFE_ALLOCATE(doneint(1:id))
        call dstates_elec_me_one_body(st, namespace, gr, hm%d%nspin, hm%vhxc, id, iindex(:,1), jindex(:,1), doneint)
        do ll = 1, id
          write(iunit, *) iindex(ll,1), jindex(ll,1), doneint(ll)
        enddo
        SAFE_DEALLOCATE_A(doneint)
      else
        SAFE_ALLOCATE(zoneint(1:id))
        call zstates_elec_me_one_body(st, namespace, gr, hm%d%nspin, hm%vhxc, id, iindex(:,1), jindex(:,1), zoneint)
        do ll = 1, id
          write(iunit, *) iindex(ll,1), jindex(ll,1), zoneint(ll)
        enddo
        SAFE_DEALLOCATE_A(zoneint)
      end if

      SAFE_DEALLOCATE_A(iindex)
      SAFE_DEALLOCATE_A(jindex)
      call io_close(iunit)

    end if

    if (this%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY) .or. this%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY_EXC_K)) then
      message(1) = "Computing two-body matrix elements"
      call messages_info(1)

      ASSERT(.not. st%parallel_in_states)
        if (this%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY)) then
        if(st%parallel_in_states)  call messages_not_implemented("OutputMatrixElements=two_body with states parallelization")
        if(st%d%kpt%parallel) call messages_not_implemented("OutputMatrixElements=two_body with k-points parallelization")
        ! how to do this properly? states_matrix
        iunit = io_open(trim(dir)//'/output_me_two_body', namespace, action='write')
        write(iunit, '(a)') '#(n1,k1) (n2,k2) (n3,k3) (n4,k4) (n1-k1, n2-k2|n3-k3, n4-k4)'
      else
        if(st%parallel_in_states)  call messages_not_implemented("OutputMatrixElements=two_body_exc_k with states parallelization")
        if(st%d%kpt%parallel) call messages_not_implemented("OutputMatrixElements=two_body_exc_k with k-points parallelization")
        ! how to do this properly? states_matrix
        iunit = io_open(trim(dir)//'/output_me_two_body_density', namespace, action='write')
        write(iunit, '(a)') '#(n1,k1) (n2,k2) (n1-k1, n1-k2|n2-k2, n2-k1)'
      end if

      if (this%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY)) then
        if(states_are_real(st)) then
          id = st%d%nik*this%nst*(st%d%nik*this%nst+1)*(st%d%nik**2*this%nst**2+st%d%nik*this%nst+2)/8
        else
          id = (st%d%nik*this%nst)**4
        end if
      else
        id = (st%d%nik*this%nst)**2
      end if

      if(states_are_complex(st)) then
        call singularity_init(singul, namespace, space, st, hm%kpoints)
      end if

      SAFE_ALLOCATE(iindex(1:2, 1:id))
      SAFE_ALLOCATE(jindex(1:2, 1:id))
      SAFE_ALLOCATE(kindex(1:2, 1:id))
      SAFE_ALLOCATE(lindex(1:2, 1:id))

      if(states_are_real(st)) then
        SAFE_ALLOCATE(dtwoint(1:id))
        call dstates_elec_me_two_body(st, namespace, space, gr, hm%kpoints, hm%exxop%psolver, this%st_start, &
                this%st_end, iindex, jindex, kindex, lindex, dtwoint)
        do ll = 1, id
          write(iunit, '(4(i4,i5),e15.6)') iindex(1:2,ll), jindex(1:2,ll), kindex(1:2,ll), lindex(1:2,ll), dtwoint(ll)
        enddo
        SAFE_DEALLOCATE_A(dtwoint)
      else
        SAFE_ALLOCATE(ztwoint(1:id))
        if (allocated(hm%hm_base%phase)) then
          !We cannot pass the phase array like that if kpt%start is not 1.  
          ASSERT(.not.st%d%kpt%parallel) 
          call zstates_elec_me_two_body(st, namespace, space, gr, hm%kpoints, hm%exxop%psolver, this%st_start, this%st_end, &
                     iindex, jindex, kindex, lindex, ztwoint, phase = hm%hm_base%phase, &
                     singularity = singul, exc_k = (this%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY_EXC_K))) 
        else
          call zstates_elec_me_two_body(st, namespace, space, gr, hm%kpoints, hm%exxop%psolver, this%st_start, this%st_end, &
                     iindex, jindex, kindex, lindex, ztwoint, exc_k = (this%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY_EXC_K)))
        end if

        if (this%what(OPTION__OUTPUTMATRIXELEMENTS__TWO_BODY)) then
          do ll = 1, id
            write(iunit, '(4(i4,i5),2e15.6)') iindex(1:2,ll), jindex(1:2,ll), kindex(1:2,ll), lindex(1:2,ll), ztwoint(ll)
          enddo
        else
          do ll = 1, id
            write(iunit, '(2(i4,i5),2e15.6)') iindex(1:2,ll), kindex(1:2,ll), ztwoint(ll)
          enddo
        end if
        SAFE_DEALLOCATE_A(ztwoint)
      end if
      
      SAFE_DEALLOCATE_A(iindex)
      SAFE_DEALLOCATE_A(jindex)
      SAFE_DEALLOCATE_A(kindex)
      SAFE_DEALLOCATE_A(lindex)
      call io_close(iunit)

      call singularity_end(singul)

    end if

    POP_SUB(output_me)
  end subroutine output_me


  ! ---------------------------------------------------------
  subroutine output_me_out_momentum(fname, st, space, gr, namespace, kpoints)
    character(len=*),    intent(in)    :: fname
    type(states_elec_t), intent(inout) :: st
    type(space_t),       intent(in)    :: space
    type(grid_t),        intent(in)    :: gr
    type(namespace_t),   intent(in)    :: namespace
    type(kpoints_t),     intent(in)    :: kpoints

    integer            :: ik, ist, is, ns, iunit, idir
    character(len=80)  :: cspin, str_tmp
    FLOAT              :: kpoint(1:MAX_DIM)
    FLOAT, allocatable :: momentum(:,:,:)

    PUSH_SUB(output_me_out_momentum)

    SAFE_ALLOCATE(momentum(1:space%dim, 1:st%nst, 1:st%d%nik))

    call states_elec_calc_momentum(st, space, gr%der, kpoints, momentum)

    iunit = io_open(fname, namespace, action='write')

    ns = 1
    if(st%d%nspin == 2) ns = 2

    write(message(1),'(a)') 'Momentum of the KS states [a.u.]:'
    call messages_info(1, iunit)      
    if (st%d%nik > ns) then
      message(1) = 'k-points [' // trim(units_abbrev(unit_one/units_out%length)) // ']'
      call messages_info(1, iunit)
    end if

    do ik = 1, st%d%nik, ns
      kpoint = M_ZERO
      kpoint(1:space%dim) = kpoints%get_point(st%d%get_kpoint_index(ik))

      if(st%d%nik > ns) then
        write(message(1), '(a,i4, a)') '#k =', ik, ', k = ('
        do idir = 1, space%dim
          write(str_tmp, '(f12.6, a)') units_from_atomic(unit_one/units_out%length, kpoint(idir)), ','
          message(1) = trim(message(1)) // trim(str_tmp)
          if(idir == space%dim) then
            message(1) = trim(message(1)) // ")"
          else
            message(1) = trim(message(1)) // ","
          end if
        end do
        call messages_info(1, iunit)
      end if

      write(message(1), '(a4,1x,a5)') '#st',' Spin'
      do idir = 1, space%dim
        write(str_tmp, '(a,a1,a)') '        <p', index2axis(idir), '>'
        message(1) = trim(message(1)) // trim(str_tmp)
      end do
      write(str_tmp, '(4x,a12,1x)') 'Occupation '
      message(1) = trim(message(1)) // trim(str_tmp)
      call messages_info(1, iunit)
      
      do ist = 1, st%nst
        do is = 0, ns-1

          if(is  ==  0) cspin = 'up'
          if(is  ==  1) cspin = 'dn'
          if(st%d%ispin  ==  UNPOLARIZED .or. st%d%ispin  ==  SPINORS) cspin = '--'
          
          write(message(1), '(i4,3x,a2,1x)') ist, trim(cspin)
          do idir = 1, space%dim
            write(str_tmp, '(f12.6)') momentum(idir, ist, ik+is)
            message(1) = trim(message(1)) // trim(str_tmp)
          end do
          write(str_tmp, '(3x,f12.6)') st%occ(ist, ik+is)
          message(1) = trim(message(1)) // trim(str_tmp)
          call messages_info(1, iunit)
          
        end do
      end do
      
      write(message(1),'(a)') ''
      call messages_info(1, iunit)      
      
    end do
    
    SAFE_DEALLOCATE_A(momentum)
    call io_close(iunit)

    POP_SUB(output_me_out_momentum)
  end subroutine output_me_out_momentum


  ! ---------------------------------------------------------
  subroutine output_me_out_ang_momentum(fname, st, gr, namespace, kpoints)
    character(len=*),    intent(in)    :: fname
    type(states_elec_t), intent(inout) :: st
    type(grid_t),        intent(in)    :: gr
    type(namespace_t),   intent(in)    :: namespace
    type(kpoints_t),     intent(in)    :: kpoints

    integer            :: iunit, ik, ist, is, ns, idir, kstart, kend
    character(len=80)  :: tmp_str(MAX_DIM), cspin
    FLOAT              :: angular(3), lsquare, kpoint(1:MAX_DIM)
    FLOAT, allocatable :: ang(:, :, :), ang2(:, :)
#if defined(HAVE_MPI)
    integer            :: tmp
    FLOAT, allocatable :: lang(:, :)
    integer            :: kn
#endif

    PUSH_SUB(output_me_out_ang_momentum)

    ns = 1
    if(st%d%nspin == 2) ns = 2
    ASSERT(gr%sb%dim == 3)

    iunit = io_open(fname, namespace, action='write')

    if(mpi_grp_is_root(mpi_world)) then
      write(iunit,'(a)') 'Warning: When non-local pseudopotentials are used '
      write(iunit,'(a)') '         the numbers below may be meaningless.    '
      write(iunit,'(a)') '                                                  '
      write(iunit,'(a)') 'Angular Momentum of the KS states [dimensionless]:'
      ! r x k is dimensionless. we do not include hbar.
      if (st%d%nik > ns) then
        message(1) = 'k-points [' // trim(units_abbrev(unit_one/units_out%length)) // ']'
        call messages_info(1, iunit)
      end if
    end if

    SAFE_ALLOCATE(ang (1:st%nst, 1:st%d%nik, 1:3))
    SAFE_ALLOCATE(ang2(1:st%nst, 1:st%d%nik))

    if (states_are_real(st)) then
      call dstates_elec_angular_momentum(st, gr, ang, ang2)
    else
      call zstates_elec_angular_momentum(st, gr, ang, ang2)
    end if

    kstart = st%d%kpt%start
    kend = st%d%kpt%end
    do idir = 1, 3
      angular(idir) = states_elec_eigenvalues_sum(st, ang(st%st_start:st%st_end, kstart:kend, idir))
    end do
    lsquare = states_elec_eigenvalues_sum(st, ang2(st%st_start:st%st_end, kstart:kend))

#if defined(HAVE_MPI)
    if(st%d%kpt%parallel) then
      kn = st%d%kpt%nlocal
      
      ASSERT(.not. st%parallel_in_states)
      
      ! note: could use lmpi_gen_allgatherv here?
      SAFE_ALLOCATE(lang(1:st%lnst, 1:kn))
      do idir = 1, 3
        lang(1:st%lnst, 1:kn) = ang(st%st_start:st%st_end, kstart:kend, idir)
        call MPI_Allgatherv(lang, st%nst*kn, MPI_FLOAT, &
          ang(:, :, idir), st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
          st%d%kpt%mpi_grp%comm, mpi_err)
      end do
      lang(1:st%lnst, 1:kn) = ang2(st%st_start:st%st_end, kstart:kend)
      call MPI_Allgatherv(lang, st%nst*kn, MPI_FLOAT, &
        ang2, st%d%kpt%num(:)*st%nst, (st%d%kpt%range(1, :) - 1)*st%nst, MPI_FLOAT, &
        st%d%kpt%mpi_grp%comm, mpi_err)
      SAFE_DEALLOCATE_A(lang)
   end if

   if(st%parallel_in_states) then
      SAFE_ALLOCATE(lang(1:st%lnst, 1))
    end if
#endif

    do ik = 1, st%d%nik, ns
      if(st%d%nik > ns) then

        kpoint = M_ZERO
        kpoint(1:gr%sb%dim) = kpoints%get_point(st%d%get_kpoint_index(ik))
        
        write(message(1), '(a,i4, a)') '#k =', ik, ', k = ('
        do idir = 1, gr%sb%dim
          write(tmp_str(1), '(f12.6, a)') units_from_atomic(unit_one/units_out%length, kpoint(idir)), ','
          message(1) = trim(message(1)) // trim(tmp_str(1))
          if(idir == gr%sb%dim) then
            message(1) = trim(message(1)) // ")"
          else
            message(1) = trim(message(1)) // ","
          end if
        end do
        call messages_info(1, iunit)
      end if
      
      ! Exchange ang and ang2.
#if defined(HAVE_MPI)
      if(st%parallel_in_states) then
        ASSERT(.not. st%d%kpt%parallel)

        do is = 1, ns
          do idir = 1, 3
            lang(1:st%lnst, 1) = ang(st%st_start:st%st_end, ik+is-1, idir)
            call lmpi_gen_allgatherv(st%lnst, lang(:, 1), tmp, ang(:, ik+is-1, idir), st%mpi_grp)
          end do
          lang(1:st%lnst, 1) = ang2(st%st_start:st%st_end, ik+is-1)
          call lmpi_gen_allgatherv(st%lnst, lang(:, 1), tmp, ang2(:, ik+is-1), st%mpi_grp)
        end do
      end if
#endif
      write(message(1), '(a4,1x,a5,4a12,4x,a12,1x)')       &
        '#st',' Spin','        <Lx>', '        <Ly>', '        <Lz>', '        <L2>', 'Occupation '
      call messages_info(1, iunit)

      if(mpi_grp_is_root(mpi_world)) then
        do ist = 1, st%nst
          do is = 0, ns-1
            
            if(is  ==  0) cspin = 'up'
            if(is  ==  1) cspin = 'dn'
            if(st%d%ispin  ==  UNPOLARIZED .or. st%d%ispin  ==  SPINORS) cspin = '--'
            
            write(tmp_str(1), '(i4,3x,a2)') ist, trim(cspin)
            write(tmp_str(2), '(1x,4f12.6,3x,f12.6)') &
              (ang(ist, ik+is, idir), idir = 1, 3), ang2(ist, ik+is), st%occ(ist, ik+is)
            message(1) = trim(tmp_str(1))//trim(tmp_str(2))
            call messages_info(1, iunit)
          end do
        end do
      end if
      write(message(1),'(a)') ''
      call messages_info(1, iunit)
      
    end do

    write(message(1),'(a)') 'Total Angular Momentum L [dimensionless]'
    write(message(2),'(10x,4f12.6)') angular(1:3), lsquare
    call messages_info(2, iunit)

    call io_close(iunit)

#if defined(HAVE_MPI)
    if(st%parallel_in_states) then
      SAFE_DEALLOCATE_A(lang)
    end if
#endif
    
    SAFE_DEALLOCATE_A(ang)
    SAFE_DEALLOCATE_A(ang2)
    
    POP_SUB(output_me_out_ang_momentum)
  end subroutine output_me_out_ang_momentum


#include "undef.F90"
#include "real.F90"
#include "output_me_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "output_me_inc.F90"


end module output_me_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
