!! Copyright (C) 2002-2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module mix_oct_m
  use comm_oct_m
  use derivatives_oct_m
  use global_oct_m
  use lalg_adv_oct_m
  use lalg_basic_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use nl_operator_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use stencil_cube_oct_m
  use types_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                     &
    mix_set_mixing,             &
    mix_t,                      &
    mix_init,                   &
    mix_clear,                  &
    mix_end,                    &
    mix_dump,                   &
    mix_load,                   &
    mixing,                     &
    dmixing,                    &
    zmixing,                    &
    mix_coefficient,            &
    mix_scheme,                 &
    mix_d4,                     &
    mix_get_field,              &
    mixfield_t,                 &
    mixfield_nullify,           &
    mixfield_init,              &
    mixfield_clear,             &
    mixfield_end,               &
    mixfield_set_vin,           &
    mixfield_set_vout,          &
    mixfield_get_vnew,          &
    mix_add_auxmixfield

  type mixfield_t
    private
    FLOAT, pointer :: ddf(:, :, :, :)
    FLOAT, pointer :: ddv(:, :, :, :)
    FLOAT, pointer :: df_old(:, :, :)
    FLOAT, pointer :: dvin_old(:, :, :)
    FLOAT, pointer :: dvin(:, :, :)
    FLOAT, pointer :: dvout(:, :, :)
    FLOAT, pointer :: dvnew(:, :, :)

    CMPLX, pointer :: zdf(:, :, :, :)
    CMPLX, pointer :: zdv(:, :, :, :)
    CMPLX, pointer :: zf_old(:, :, :)
    CMPLX, pointer :: zvin_old(:, :, :)
    CMPLX, pointer :: zvin(:, :, :)
    CMPLX, pointer :: zvout(:, :, :)
    CMPLX, pointer :: zvnew(:, :, :)

    type(type_t) :: func_type   !< type of the functions to be mixed
    integer :: d1, d2, d3, d4   !< the dimensions of the arrays that store the information from the previous iterations
  end type mixfield_t

  type mixfield_ptr_t
    private
    type(mixfield_t), pointer :: p
  end type mixfield_ptr_t

  integer, parameter :: MAX_AUXMIXFIELD = 5

  type mix_t
    private
    integer :: scheme           !< the mixing scheme used (linear, broyden, etc)

    FLOAT :: coeff              !< the mixing coefficient (in linear mixing: vnew = (1-coeff)*vin + coeff*vout)

    FLOAT :: dgamma             !< The gamma coefficient of the Broyden mixing scheme
    CMPLX :: zgamma             !< which is stored here for possible auxiliarry mixfields

    integer :: iter             !< number of SCF iterations already done. In case of restart, this number must
                                !< include the iterations done in previous calculations.

    integer :: ns               !< number of steps used to extrapolate the new vector

    integer :: ipos             !< For auxiliary mixing fields
    integer :: last_ipos        !< where is the information about the last iteration stored in arrays df and dv

    integer :: interval

    type(mixfield_t) :: mixfield    !< The field to be mixed

    integer :: nauxmixfield            !< Number of auxiliary mixing fields
    type(mixfield_ptr_t) :: auxmixfield(MAX_AUXMIXFIELD) !< Auxiliary mixing fields

    type(derivatives_t), pointer :: der
    logical                      :: precondition
    type(nl_operator_t)          :: preconditioner

    FLOAT :: residual_coeff
    
  end type mix_t

  interface mixfield_set_vin
    module procedure dmixfield_set_vin2, dmixfield_set_vin3, &
                     zmixfield_set_vin2, zmixfield_set_vin3, &
                     ddmixfield_set_vin2
  end interface mixfield_set_vin
 
  interface mixfield_set_vout
    module procedure dmixfield_set_vout2, dmixfield_set_vout3, &
                     zmixfield_set_vout2, zmixfield_set_vout3, &
                     ddmixfield_set_vout2
  end interface mixfield_set_vout

  interface mixfield_get_vnew
    module procedure dmixfield_get_vnew, ddmixfield_get_vnew, &
                     zmixfield_get_vnew
  end interface mixfield_get_vnew


contains

  ! ---------------------------------------------------------
  subroutine mix_init(smix, namespace, der, d1, d2, d3, def_, func_type_, prefix_)
    type(mix_t),                   intent(out) :: smix
    type(namespace_t),             intent(in)  :: namespace
    type(derivatives_t), target,   intent(in)  :: der
    integer,                       intent(in)  :: d1, d2, d3
    integer,             optional, intent(in)  :: def_
    type(type_t),        optional, intent(in)  :: func_type_
    character(len=*),    optional, intent(in)  :: prefix_

    integer :: def, ii
    character(len=32) :: prefix
    type(type_t) :: func_type

    PUSH_SUB(mix_init)

    smix%der => der

    def = OPTION__MIXINGSCHEME__BROYDEN
    if(present(def_)) def = def_
    if(present(func_type_)) then 
      func_type = func_type_
    else 
      func_type = TYPE_FLOAT
    end if
    prefix = ""
    if(present(prefix_)) prefix = prefix_

    call messages_obsolete_variable(namespace, 'TypeOfMixing', 'MixingScheme')
    
    !%Variable MixingScheme
    !%Type integer
    !%Default broyden
    !%Section SCF::Mixing
    !%Description
    !% The scheme used to produce, at each iteration in the self-consistent cycle
    !% that attempts to solve the Kohn-Sham equations, the input density from the value
    !% of the input and output densities of previous iterations.
    !%Option linear 0
    !% Simple linear mixing.
    !%Option broyden 2
    !% Broyden scheme [C. G Broyden, <i>Math. Comp.</i> <b>19</b>, 577 (1965); 
    !% D. D. Johnson, <i>Phys. Rev. B</i> <b>38</b>, 12807 (1988)].
    !% For complex functions (e.g. Sternheimer with <tt>EMEta</tt> > 0), we use the generalization
    !% with a complex dot product.
    !%Option diis 9
    !% Direct inversion in the iterative subspace (diis)
    !% scheme [P. Pulay, <i>Chem. Phys. Lett.</i>, <b>73</b>, 393
    !% (1980)] as described in [G. Kresse, and J. Hurthmueller,
    !% <i>Phys. Rev. B</i> <b>54</b>, 11169 (1996)].
    !%Option bowler_gillan 1
    !% The Guaranteed-reduction modification of the Pulay scheme by
    !% Bowler and Gillan [D. R. Bowler and M. J. Gillan,
    !% <i>Chem. Phys.  Lett.</i> <b>325</b>, 473 (2000)].
    !%End
    call parse_variable(namespace, trim(prefix)//'MixingScheme', def, smix%scheme)
    if(.not.varinfo_valid_option('MixingScheme', smix%scheme)) call messages_input_error('MixingScheme', 'invalid option')
    call messages_print_var_option(stdout, "MixingScheme", smix%scheme)

    if(smix%scheme == OPTION__MIXINGSCHEME__DIIS) call messages_experimental('MixingScheme = diis')

    !%Variable MixingPreconditioner
    !%Type logical
    !%Default false
    !%Section SCF::Mixing
    !%Description
    !% (Experimental) If set to yes, Octopus will use a preconditioner
    !% for the mixing operator.
    !%End
    call parse_variable(namespace, trim(prefix)+'MixingPreconditioner', .false., smix%precondition)
    if(smix%precondition) call messages_experimental('MixingPreconditioner')
    
    !%Variable Mixing
    !%Type float
    !%Default 0.3
    !%Section SCF::Mixing
    !%Description
    !% The linear, Broyden and DIIS scheme depend on a "mixing parameter", set by this variable. 
    !% Must be 0 < <tt>Mixing</tt> <= 1.
    !%End
    call parse_variable(namespace, trim(prefix)+'Mixing', CNST(0.3), smix%coeff)
    if(smix%coeff <= M_ZERO .or. smix%coeff > M_ONE) then
      call messages_input_error('Mixing', 'Value should be positive and smaller than one.')
    end if
    
    !%Variable MixingResidual
    !%Type float
    !%Default 0.05
    !%Section SCF::Mixing
    !%Description
    !% In the DIIS mixing it is benefitial to include a bit of
    !% residual into the mixing. This parameter controls this amount.
    !%End
    call parse_variable(namespace, trim(prefix)+'MixingResidual', CNST(0.05), smix%residual_coeff)
    if(smix%residual_coeff <= M_ZERO .or. smix%residual_coeff > M_ONE) then
      call messages_input_error('MixingResidual', 'Value should be positive and smaller than one.')
    end if
    
    !%Variable MixNumberSteps
    !%Type integer
    !%Default 3
    !%Section SCF::Mixing
    !%Description
    !% In the Broyden and Bowler_Gillan schemes, the new input density or potential is constructed
    !% from the values of the densities/potentials of a given number of previous iterations.
    !% This number is set by this variable. Must be greater than 1.
    !%End
    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      call parse_variable(namespace, trim(prefix)//'MixNumberSteps', 3, smix%ns)
      if(smix%ns <= 1) call messages_input_error('MixNumberSteps')
    else
      smix%ns = 0
    end if
    
    !%Variable MixInterval
    !%Type integer
    !%Default 1
    !%Section SCF::Mixing
    !%Description
    !% When this variable is set to a value different than 1 (the
    !% defaul) a combined mixing scheme will be used, with MixInterval
    !% - 1 steps of linear mixing followed by 1 step of the selected
    !% mixing. For the moment this variable only works with DIIS mixing.
    !%End
    call parse_variable(namespace, trim(prefix)//'MixInterval', 1, smix%interval)
    if(smix%interval < 1) call messages_input_error('MixInterval', 'MixInterval must be larger or equal than 1')
    
    smix%iter = 0

    smix%nauxmixfield = 0
    do ii=1,MAX_AUXMIXFIELD
      nullify(smix%auxmixfield(ii)%p)   
    end do


    select case (smix%scheme)
    case (OPTION__MIXINGSCHEME__LINEAR, OPTION__MIXINGSCHEME__BROYDEN, OPTION__MIXINGSCHEME__DIIS)
      call mixfield_init(smix, smix%mixfield, d1, d2, d3, smix%ns, func_type)
    case (OPTION__MIXINGSCHEME__BOWLER_GILLAN)
      call mixfield_init(smix, smix%mixfield, d1, d2, d3, smix%ns+1, func_type)
    end select

    call mix_clear(smix)

    if(smix%precondition) call init_preconditioner()

    POP_SUB(mix_init)

  contains

    subroutine init_preconditioner()

      integer :: ns, maxp, ip, is
      FLOAT, parameter :: weight = CNST(50.0)
      
      ! This the mixing preconditioner from GPAW:
      !
      !   https://wiki.fysik.dtu.dk/gpaw/documentation/densitymix/densitymix.html
      !

      ASSERT(.not. der%mesh%use_curvilinear)
      ASSERT(der%mesh%sb%dim == 3)
      
      call nl_operator_init(smix%preconditioner, "Mixing preconditioner")
      call stencil_cube_get_lapl(smix%preconditioner%stencil, der%mesh%sb%dim, 1)
      call nl_operator_build(der%mesh, smix%preconditioner, der%mesh%np, const_w = .not. der%mesh%use_curvilinear)
      
      ns = smix%preconditioner%stencil%size

      if (smix%preconditioner%const_w) then
        maxp = 1
      else
        maxp = der%mesh%np
      end if

      do ip = 1, maxp

        do is = 1, ns
          select case(sum(abs(smix%preconditioner%stencil%points(1:der%mesh%sb%dim, is))))
          case(0)
            smix%preconditioner%w(is, ip) = CNST(1.0) + weight/CNST(8.0)
          case(1)
            smix%preconditioner%w(is, ip) = weight/CNST(16.0)
          case(2)
            smix%preconditioner%w(is, ip) = weight/CNST(32.0)
          case(3)
            smix%preconditioner%w(is, ip) = weight/CNST(64.0)
          case default
            ASSERT(.false.)
          end select

        end do
      end do
      
      call nl_operator_update_weights(smix%preconditioner)

    end subroutine init_preconditioner

  end subroutine mix_init


  ! ---------------------------------------------------------
  subroutine mix_clear(smix)
    type(mix_t),             intent(inout) :: smix
    
    PUSH_SUB(mix_clear)

    call mixfield_clear(smix%scheme, smix%mixfield)

    smix%iter = 0
    smix%last_ipos = 0

    POP_SUB(mix_clear)
  end subroutine mix_clear


  ! ---------------------------------------------------------
  subroutine mix_end(smix)
    type(mix_t), intent(inout) :: smix

    integer :: ii 

    PUSH_SUB(mix_end)

    if(smix%precondition) call nl_operator_end(smix%preconditioner)
   
    call mixfield_end(smix, smix%mixfield)

    smix%nauxmixfield = 0
    do ii=1,MAX_AUXMIXFIELD
      nullify(smix%auxmixfield(ii)%p)
    end do

    POP_SUB(mix_end)
  end subroutine mix_end


  ! ---------------------------------------------------------
  subroutine mix_set_mixing(smix, newmixing)
    type(mix_t), intent(inout) :: smix
    FLOAT, intent(in):: newmixing

    PUSH_SUB(mix_set_mixing)
    
    if(smix%scheme == OPTION__MIXINGSCHEME__LINEAR) then
      smix%coeff = newmixing
    else
    !  message(1) = "Mixing can only be adjusted in linear mixing scheme."
    !  call messages_fatal(1)
    end if
    
    POP_SUB(mix_set_mixing)
  end subroutine mix_set_mixing


  ! ---------------------------------------------------------
  subroutine mix_dump(restart, smix, mesh, ierr)
    type(restart_t), intent(in)  :: restart
    type(mix_t),     intent(in)  :: smix
    type(mesh_t),    intent(in)  :: mesh
    integer,         intent(out) :: ierr

    integer :: iunit, id2, id3, id4, err, err2(4)
    character(len=40) :: lines(8)
    character(len=80) :: filename

    PUSH_SUB(mix_dump)

    ierr = 0

    if (restart_skip(restart)) then
      POP_SUB(mix_dump)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Writing mixing restart."
      call messages_info(1)
    end if

    ! functions to be written need to be compatible with the mesh
    ASSERT(mesh%np == smix%mixfield%d1)

    ! First we write some information about the mixing
    iunit = restart_open(restart, 'mixing')
    write(lines(1), '(a11,i1)')  'scheme=    ', smix%scheme
    ! Number of global mesh points have to be written, not only smix%d1
    write(lines(2), '(a11,i10)') 'd1=        ', mesh%np_global
    write(lines(3), '(a11,i10)') 'd2=        ', smix%mixfield%d2
    write(lines(4), '(a11,i10)') 'd3=        ', smix%mixfield%d3
    write(lines(5), '(a11,i10)') 'd4=        ', smix%mixfield%d4
    write(lines(6), '(a11,i10)') 'iter=      ', smix%iter
    write(lines(7), '(a11,i10)') 'ns=        ', smix%ns
    write(lines(8), '(a11,i10)') 'last_ipos= ', smix%last_ipos
    call restart_write(restart, iunit, lines, 8, err)
    if (err /= 0) ierr = ierr + 1
    call restart_close(restart, iunit)

    ! Now we write the different functions. 
    ! These are not needed when using linear mixing, so we will make sure we skip this step in that case.
    err2 = 0
    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      do id2 = 1, smix%mixfield%d2
        do id3 = 1, smix%mixfield%d3
          do id4 = 1, smix%mixfield%d4

            write(filename,'(a3,i2.2,i2.2,i2.2)') 'df_', id2, id3, id4
            if (smix%mixfield%func_type == TYPE_FLOAT) then
              call drestart_write_mesh_function(restart, filename, mesh, smix%mixfield%ddf(1:mesh%np, id2, id3, id4), err)
            else
              call zrestart_write_mesh_function(restart, filename, mesh, smix%mixfield%zdf(1:mesh%np, id2, id3, id4), err)
            end if
            if (err /= 0) err2(1) = err2(1) + 1

            write(filename,'(a3,i2.2,i2.2,i2.2)') 'dv_', id2, id3, id4
            if (smix%mixfield%func_type == TYPE_FLOAT) then
              call drestart_write_mesh_function(restart, filename, mesh, smix%mixfield%ddv(1:mesh%np, id2, id3, id4), err)
            else
              call zrestart_write_mesh_function(restart, filename, mesh, smix%mixfield%zdv(1:mesh%np, id2, id3, id4), err)
            end if
            if (err /= 0) err2(2) = err2(2) + 1
              
          end do

          write(filename,'(a6,i2.2,i2.2)') 'f_old_', id2, id3
          if (smix%mixfield%func_type == TYPE_FLOAT) then
            call drestart_write_mesh_function(restart, filename, mesh, smix%mixfield%df_old(1:mesh%np, id2, id3), err)
          else
            call zrestart_write_mesh_function(restart, filename, mesh, smix%mixfield%zf_old(1:mesh%np, id2, id3), err)
          end if
          if (err /= 0) err2(3) = err2(3) + 1

          write(filename,'(a8,i2.2,i2.2)') 'vin_old_', id2, id3
          if (smix%mixfield%func_type == TYPE_FLOAT) then
            call drestart_write_mesh_function(restart, filename, mesh, smix%mixfield%dvin_old(1:mesh%np, id2, id3), err)
          else
            call zrestart_write_mesh_function(restart, filename, mesh, smix%mixfield%zvin_old(1:mesh%np, id2, id3), err)
          end if
          if (err /= 0) err2(4) = err2(4) + 1
          
        end do
      end do

      if (err2(1) /= 0) ierr = ierr + 2
      if (err2(2) /= 0) ierr = ierr + 4
      if (err2(3) /= 0) ierr = ierr + 8
      if (err2(4) /= 0) ierr = ierr + 16
    end if

    if (debug%info) then
      message(1) = "Debug: Writing mixing restart done."
      call messages_info(1)
    end if

    POP_SUB(mix_dump)
  end subroutine mix_dump


  !---------------------------------------------------------
  subroutine mix_load(restart, smix, mesh, ierr)
    type(restart_t), intent(in)    :: restart
    type(mix_t),     intent(inout) :: smix
    type(mesh_t),    intent(in)    :: mesh
    integer,         intent(out)   :: ierr

    integer :: iunit, err, err2(4)
    integer :: scheme, d1, d2, d3, d4, ns
    integer :: id2, id3, id4
    character(len=11)  :: str
    character(len=80)  :: filename
    character(len=256) :: lines(8)

    PUSH_SUB(mix_load)

    ierr = 0

    if (restart_skip(restart)) then
      ierr = -1
      POP_SUB(mix_load)
      return
    end if

    if (debug%info) then
      message(1) = "Debug: Reading mixing restart."
      call messages_info(1)
    end if

    ! First we read some information about the mixing
    iunit = restart_open(restart, 'mixing')
    call restart_read(restart, iunit, lines, 8, err)
    if (err /= 0) then
      ierr = ierr + 1
    else
      read(lines(1), *) str, scheme
      read(lines(2), *) str, d1
      read(lines(3), *) str, d2
      read(lines(4), *) str, d3
      read(lines(5), *) str, d4
      read(lines(6), *) str, smix%iter
      read(lines(7), *) str, ns
      read(lines(8), *) str, smix%last_ipos
    end if
    call restart_close(restart, iunit)


    if (ierr == 0) then
      ! We can only use the restart information if the mixing scheme and the number of steps used remained the same
      if (scheme /= smix%scheme .or. ns /= smix%ns) then
        message(1) = "The mixing scheme from the restart data is not the same as the one used in the current calculation."
        call messages_warning(1)
        ierr = ierr + 2
      end if

      ! Check the dimensions of the arrays to be read
      if (mesh%np_global /= d1 .or. mesh%np /= smix%mixfield%d1 .or. d2 /= smix%mixfield%d2 .or. d3 /= smix%mixfield%d3 ) then
        message(1) = "The dimensions of the arrays from the mixing restart data"
        message(2) = "are not the same as the ones used in this calculation."
        call messages_warning(2)
        ierr = ierr + 4
      end if
    end if


    ! Now we read the different functions.
    ! Note that we may have more or less functions than the ones needed (d4 /= smix%d4)
    if (ierr == 0) then
      if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
        err2 = 0
        do id2 = 1, smix%mixfield%d2
          do id3 = 1, smix%mixfield%d3
            do id4 = 1, smix%mixfield%d4

              write(filename,'(a3,i2.2,i2.2,i2.2)') 'df_', id2, id3, id4
              if (smix%mixfield%func_type == TYPE_FLOAT) then
                call drestart_read_mesh_function(restart, trim(filename), mesh, smix%mixfield%ddf(1:mesh%np, id2, id3, id4), err)
              else
                call zrestart_read_mesh_function(restart, trim(filename), mesh, smix%mixfield%zdf(1:mesh%np, id2, id3, id4), err)
              end if
              if (err /= 0) err2(1) = err2(1) + 1
          
              write(filename,'(a3,i2.2,i2.2,i2.2)') 'dv_', id2, id3, id4
              if (smix%mixfield%func_type == TYPE_FLOAT) then
                call drestart_read_mesh_function(restart, trim(filename), mesh, smix%mixfield%ddv(1:mesh%np, id2, id3, id4), err)
              else
                call zrestart_read_mesh_function(restart, trim(filename), mesh, smix%mixfield%zdv(1:mesh%np, id2, id3, id4), err)
              end if
              if (err /= 0) err2(2) = err2(2) + 1

            end do

            write(filename,'(a6,i2.2,i2.2)') 'f_old_', id2, id3
            if (smix%mixfield%func_type == TYPE_FLOAT) then
              call drestart_read_mesh_function(restart, trim(filename), mesh, smix%mixfield%df_old(1:mesh%np, id2, id3), err)
            else
              call zrestart_read_mesh_function(restart, trim(filename), mesh, smix%mixfield%zf_old(1:mesh%np, id2, id3), err)
            end if
            if (err /= 0) err2(3) = err2(3) + 1

            write(filename,'(a8,i2.2,i2.2)') 'vin_old_', id2, id3
            if (smix%mixfield%func_type == TYPE_FLOAT) then
              call drestart_read_mesh_function(restart, trim(filename), mesh, smix%mixfield%dvin_old(1:mesh%np, id2, id3), err)
            else
              call zrestart_read_mesh_function(restart, trim(filename), mesh, smix%mixfield%zvin_old(1:mesh%np, id2, id3), err)
            end if
            if (err /= 0) err2(4) = err2(4) + 1

          end do
        end do

        if (err2(1) /= 0) ierr = ierr + 8
        if (err2(2) /= 0) ierr = ierr + 16
        if (err2(3) /= 0) ierr = ierr + 32
        if (err2(4) /= 0) ierr = ierr + 64
      end if
    end if

    if (ierr /= 0) then
      ! Something went wront, so make sure we start from scratch
      call mix_clear(smix)
    end if

    if (debug%info) then
      message(1) = "Debug: Reading mixing restart done."
      call messages_info(1)
    end if

    POP_SUB(mix_load)
  end subroutine mix_load

  FLOAT pure function mix_coefficient(this) result(coefficient)
    type(mix_t), intent(in) :: this
    
    coefficient = this%coeff
  end function mix_coefficient

  integer pure function mix_scheme(this) result(scheme)
    type(mix_t), intent(in) :: this

    scheme = this%scheme
  end function mix_scheme

  integer pure function mix_d4(this) 
    type(mix_t), intent(in) :: this

    mix_d4 = this%mixfield%d4
  end function mix_d4

  subroutine mix_get_field(this, mixfield)
    type(mix_t), target,  intent(in) :: this
    type(mixfield_t), pointer, intent(out) :: mixfield

    mixfield => this%mixfield
  end subroutine mix_get_field

  ! ---------------------------------------------------------
  subroutine mixing(smix)
    type(mix_t),  intent(inout) :: smix
  
    PUSH_SUB(mixing)

    if(smix%mixfield%func_type == TYPE_FLOAT) then
      call dmixing(smix, smix%mixfield%dvin, smix%mixfield%dvout, smix%mixfield%dvnew)
    else
      call zmixing(smix, smix%mixfield%zvin, smix%mixfield%zvout, smix%mixfield%zvnew)
    end if
  
    POP_SUB(mixing)
  end subroutine mixing

  subroutine mix_add_auxmixfield( smix, mixfield )
    type(mix_t),      intent(inout)      :: smix
    type(mixfield_t), target, intent(in) :: mixfield

    PUSH_SUB(mix_add_auxmixfield)

    smix%nauxmixfield = smix%nauxmixfield + 1
    smix%auxmixfield(smix%nauxmixfield)%p => mixfield 

    if( smix%scheme == OPTION__MIXINGSCHEME__DIIS) &
      call messages_input_error('Mixing scheme  DIIS is not implemented for auxiliary mixing fields')

    if( smix%scheme == OPTION__MIXINGSCHEME__BOWLER_GILLAN) &
      call messages_input_error('Mixing scheme  Bowler Gillan is not implemented for auxiliary mixing fields')

    POP_SUB(mix_add_auxmixfield)
  end subroutine mix_add_auxmixfield

  subroutine mixfield_nullify( mixfield )
    type(mixfield_t), intent(inout) :: mixfield

    PUSH_SUB(mixfield_nullify)

    nullify(mixfield%ddf)
    nullify(mixfield%ddv)
    nullify(mixfield%df_old)
    nullify(mixfield%dvin_old)
    nullify(mixfield%dvin)
    nullify(mixfield%dvout)
    nullify(mixfield%dvnew)

    nullify(mixfield%zdf)
    nullify(mixfield%zdv)
    nullify(mixfield%zf_old)
    nullify(mixfield%zvin_old)
    nullify(mixfield%zvin)
    nullify(mixfield%zvout)
    nullify(mixfield%zvnew)

    POP_SUB(mixfield_nullify)
  end subroutine mixfield_nullify

  subroutine mixfield_init( smix, mixfield, d1, d2, d3, d4, func_type ) 
    type(mix_t),      intent(inout) :: smix
    type(mixfield_t), intent(inout) :: mixfield
    integer,          intent(in)    :: d1, d2, d3, d4
    type(type_t),     intent(in)    :: func_type

    PUSH_SUB(mixfield_init)

    call mixfield_nullify( mixfield )

    mixfield%d1 = d1
    mixfield%d2 = d2
    mixfield%d3 = d3
    mixfield%d4 = d4

    mixfield%func_type = func_type 

    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      if(mixfield%func_type == TYPE_FLOAT) then
        SAFE_ALLOCATE(     mixfield%ddf(1:d1, 1:d2, 1:d3, 1:d4))
        SAFE_ALLOCATE(mixfield%dvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     mixfield%ddv(1:d1, 1:d2, 1:d3, 1:d4))
        SAFE_ALLOCATE(  mixfield%df_old(1:d1, 1:d2, 1:d3))
      else
        SAFE_ALLOCATE(     mixfield%zdf(1:d1, 1:d2, 1:d3, 1:d4))
        SAFE_ALLOCATE(mixfield%zvin_old(1:d1, 1:d2, 1:d3))
        SAFE_ALLOCATE(     mixfield%zdv(1:d1, 1:d2, 1:d3, 1:d4))
        SAFE_ALLOCATE(  mixfield%zf_old(1:d1, 1:d2, 1:d3))
      end if
    end if
  
    if(mixfield%func_type == TYPE_FLOAT) then
      SAFE_ALLOCATE(mixfield%dvin(1:d1, 1:d2, 1:d3))
      SAFE_ALLOCATE(mixfield%dvout(1:d1, 1:d2, 1:d3))
      SAFE_ALLOCATE(mixfield%dvnew(1:d1, 1:d2, 1:d3))
    else
      SAFE_ALLOCATE(mixfield%zvin(1:d1, 1:d2, 1:d3))
      SAFE_ALLOCATE(mixfield%zvout(1:d1, 1:d2, 1:d3))
      SAFE_ALLOCATE(mixfield%zvnew(1:d1, 1:d2, 1:d3))
    end if 

    POP_SUB(mixfield_init)
  end subroutine mixfield_init 

    ! ---------------------------------------------------------
  subroutine mixfield_end(smix, mixfield)
    type(mix_t),      intent(inout) :: smix
    type(mixfield_t), intent(inout) :: mixfield

    PUSH_SUB(mixfield_end)

    ! Arrays got allocated for all mixing schemes, except linear mixing
    if (smix%scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      SAFE_DEALLOCATE_P(mixfield%ddf)
      SAFE_DEALLOCATE_P(mixfield%ddv)
      SAFE_DEALLOCATE_P(mixfield%dvin_old)
      SAFE_DEALLOCATE_P(mixfield%df_old)

      SAFE_DEALLOCATE_P(mixfield%zdf)
      SAFE_DEALLOCATE_P(mixfield%zdv)
      SAFE_DEALLOCATE_P(mixfield%zvin_old)
      SAFE_DEALLOCATE_P(mixfield%zf_old)
    end if

    SAFE_DEALLOCATE_P(mixfield%dvin)
    SAFE_DEALLOCATE_P(mixfield%dvout)
    SAFE_DEALLOCATE_P(mixfield%dvnew)
    SAFE_DEALLOCATE_P(mixfield%zvin)
    SAFE_DEALLOCATE_P(mixfield%zvout)
    SAFE_DEALLOCATE_P(mixfield%zvnew)

    POP_SUB(mixfield_end)
  end subroutine mixfield_end

  ! ---------------------------------------------------------
  subroutine mixfield_clear(scheme, mixfield)
    integer,             intent(in) :: scheme
    type(mixfield_t), intent(inout) :: mixfield

    PUSH_SUB(mixfield_clear)

    if (scheme /= OPTION__MIXINGSCHEME__LINEAR) then
      if(mixfield%func_type == TYPE_FLOAT) then
        mixfield%ddf = M_ZERO
        mixfield%ddv = M_ZERO
        mixfield%dvin_old = M_ZERO
        mixfield%df_old = M_ZERO
      else
        mixfield%zdf = M_z0
        mixfield%zdv = M_z0
        mixfield%zvin_old = M_z0
        mixfield%zf_old = M_z0
      end if
    end if

    if(mixfield%func_type == TYPE_FLOAT) then
      mixfield%dvin  = M_ZERO
      mixfield%dvout = M_ZERO
      mixfield%dvnew = M_ZERO
    else
      mixfield%zvin  = M_z0
      mixfield%zvout = M_z0
      mixfield%zvnew = M_z0
    end if

    POP_SUB(mixfield_clear)
  end subroutine mixfield_clear

  ! --------------------------------------------------------------
  subroutine ddmixfield_set_vin2(mixfield, vin_re, vin_im)
    type(mixfield_t), intent(inout) :: mixfield
    FLOAT,              intent(in)  :: vin_re(:,:), vin_im(:,:)

    PUSH_SUB(ddmixfield_set_vin2)

    mixfield%zvin(1:mixfield%d1, 1, 1:mixfield%d3) = vin_re(1:mixfield%d1, 1:mixfield%d3) &
                                            + M_zI * vin_im(1:mixfield%d1, 1:mixfield%d3)

    POP_SUB(ddmixfield_set_vin2)
  end subroutine ddmixfield_set_vin2

  ! --------------------------------------------------------------
  subroutine ddmixfield_set_vout2(mixfield, vout_re, vout_im)
    type(mixfield_t),  intent(inout) :: mixfield
    FLOAT,               intent(in)  :: vout_re(:,:), vout_im(:,:)

    PUSH_SUB(ddmixfield_set_vout2)

    mixfield%zvout(1:mixfield%d1, 1, 1:mixfield%d3) = vout_re(1:mixfield%d1, 1:mixfield%d3) &
                                            + M_zI * vout_im(1:mixfield%d1, 1:mixfield%d3)

    POP_SUB(ddmixfield_set_vout2)
  end subroutine ddmixfield_set_vout2

  ! --------------------------------------------------------------
  subroutine mixfield_get_dvnew(mixfield, vnew)
    type(mixfield_t),   intent(in) :: mixfield
    FLOAT,          intent(inout)  :: vnew(:,:)

    PUSH_SUB(mixfield_get_dvnew)

    vnew(1:mixfield%d1, 1:mixfield%d3) = mixfield%dvnew(1:mixfield%d1, 1, 1:mixfield%d3)

    POP_SUB(mixfield_get_dvnew)
  end subroutine mixfield_get_dvnew

  ! --------------------------------------------------------------
  subroutine ddmixfield_get_vnew(mixfield, re, im)
    type(mixfield_t),   intent(in) :: mixfield
    FLOAT,          intent(inout)  :: re(:,:), im(:,:)

    PUSH_SUB(mixfield_get_ddvnew)

    re(1:mixfield%d1, 1:mixfield%d3) =  real(mixfield%zvnew(1:mixfield%d1, 1, 1:mixfield%d3), REAL_PRECISION)
    im(1:mixfield%d1, 1:mixfield%d3) = aimag(mixfield%zvnew(1:mixfield%d1, 1, 1:mixfield%d3))

    POP_SUB(mixfield_get_ddvnew)
  end subroutine ddmixfield_get_vnew


#include "undef.F90"
#include "real.F90"

#include "mix_inc.F90"

#include "undef.F90"
#include "complex.F90"

#include "mix_inc.F90"

end module mix_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
