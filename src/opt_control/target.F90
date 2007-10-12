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
!! $Id: opt_control.F90 2870 2007-04-28 06:26:47Z acastro $

#include "global.h"

module opt_control_target_m
  use datasets_m
  use varinfo_m
  use messages_m
  use lib_oct_parser_m
  use lib_oct_m
  use io_m
  use global_m
  use string_m
  use states_m
  use grid_m
  use output_m
  use geometry_m
  use states_output_m
  use mesh_m
  use mesh_function_m
  use restart_m
  use timedep_m
  use opt_control_constants_m
  use opt_control_tdtarget_m

  implicit none

  private
  public :: target_t,       &
            target_init,    &
            target_end,     &
            target_output,  &
            calc_tdfitness, &
            calc_inh

  type target_t
    integer :: totype
    integer :: targetmode
    type(states_t) :: st
    FLOAT, pointer :: rho(:)

    integer :: no_tdtargets
    FLOAT, pointer :: td_fitness(:)
    CMPLX, pointer :: tdtarget(:)
    type(td_target_t), pointer :: tdtg(:)
  end type target_t

  contains

  ! ----------------------------------------------------------------------
  subroutine target_init(gr, geo, stin, td, target)
    type(grid_t), intent(in)      :: gr
    type(geometry_t), intent(in)  :: geo
    type(states_t), intent(in)    :: stin
    type(td_t), intent(in)        :: td
    type(target_t), intent(inout) :: target

    integer           :: ierr, ip, ist, jst
    C_POINTER         :: blk
    FLOAT             :: x(MAX_DIM), r, psi_re, psi_im
    CMPLX, allocatable :: rotation_matrix(:, :)
    type(states_t)    :: tmp_st
    character(len=1024) :: expression

    call push_sub('target.target_init')

    !%Variable OCTTargetMode
    !%Type integer
    !%Section Optimal Control
    !%Description
    !%Option oct_targetmode_static  0
    !% Static or time-independent targets
    !%Option oct_targetmode_td      1
    !% Time-dependent targets, specify block OCTTdTarget
    !%End
    call loct_parse_int(check_inp('OCTTargetMode'), oct_targetmode_static, target%targetmode)
    if(.not.varinfo_valid_option('OCTTargetMode', target%targetmode)) &
      call input_error('OCTTargetMode')

    !%Variable OCTTargetOperator
    !%Type integer
    !%Section Optimal Control
    !%Default 2
    !%Description
    !% The string OCTTargetOperator describes the initial state of the quantum system
    !% Possible arguments are:
    !%Option oct_tg_groundstate 1 
    !% Targetoperator is a projection operator on the ground state
    !%Option oct_tg_excited 2
    !% Currently not in use.
    !%Option oct_tg_gstransformation 3
    !% Targetoperator is a projection operator on a transformation of the ground state 
    !% orbitals defined by the block OCTTargetTransformStates
    !%Option oct_tg_userdefined 4
    !% Targetoperator is a projection operator on a user defined state
    !%Option oct_tg_density 5
    !% Targetoperator is a given density.
    !%Option oct_tg_local 6
    !% Target operator is a local operator.
    !%End
    call loct_parse_int(check_inp('OCTTargetOperator'),oct_tg_excited, target%totype)
    if(.not.varinfo_valid_option('OCTTargetOperator', target%totype)) call input_error('OCTTargetOperator')    

    call states_copy(target%st, stin)

    select case(target%totype)
    case(oct_tg_groundstate)
      message(1) =  'Info: Using Ground State for TargetOperator'
      call write_info(1)
      call restart_read(trim(tmpdir)//'gs', target%st, gr, geo, ierr)
      if(ierr.ne.0) then
        write(message(1),'(a)') 'Could not read ground-state wavefunctions from '//trim(tmpdir)//'gs.'
        call write_fatal(1)
      end if
      
    case(oct_tg_excited) 
      message(1) = 'Error: using an excited state as the target state for an '
      message(2) = 'optimal control run is not possible yet.'
      message(3) = 'Try using "OCTInitialState = oct_is_transformation" instead.'
      call write_fatal(3)

    case(oct_tg_gstransformation)  

      message(1) =  'Info: Using Superposition of States for TargetOperator'
      call write_info(1)

      !%Variable OCTTargetTransformStates
      !%Type block
      !%Default no
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_gstransformation, you must specify one
      !% OCTTargetTransformStates block, in order to specify which linear
      !% combination of the states present in "restart/gs" is used to
      !% create the target state.
      !% 
      !% The syntax is equivalent to the one used for the TransformStates
      !% block.
      !%End
      if(loct_parse_isdef(check_inp('OCTTargetTransformStates')).ne.0) then
        if(loct_parse_block(check_inp('OCTTargetTransformStates'), blk) == 0) then
          tmp_st = target%st
          deallocate(tmp_st%zpsi)
          call restart_look_and_read("tmp", tmp_st, gr, geo, ierr)
          ALLOCATE(rotation_matrix(target%st%nst, tmp_st%nst), target%st%nst*tmp_st%nst)
          rotation_matrix = M_z0
          do ist = 1, target%st%nst
            do jst = 1, loct_parse_block_cols(blk, ist-1)
              call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
            end do
          end do
          call rotate_states(gr%m, target%st, tmp_st, rotation_matrix)
          deallocate(rotation_matrix)
          call states_end(tmp_st)
        else
          message(1) = '"OCTTargetTransformStates" has to be specified as block.'
          call write_info(1)
          call input_error('OCTTargetTransformStates')
        end if
      else
        message(1) = 'Error: if "OCTTargetOperator = oct_tg_superposition", then you must'
        message(2) = 'supply one "OCTTargetTransformStates" block to create the superposition.'
        call write_info(2)
        call input_error('OCTTargetTransformStates')
      end if

    case(oct_tg_userdefined) 
      message(1) =  'Error: Option oct_tg_userdefined is disabled in this version'
      call write_fatal(1)

    case(oct_tg_density) 

      message(1) =  'Info: Target is a density.'
      call write_info(1)


      !%Variable OCTTargetDensity
      !%Type string
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_local, then one must supply the target density
      !% that should be searched for. This one can do by supplying a string through
      !% the variable OCTLocalTarget.
      !%End


      !%Variable OCTTargetDensityFromState
      !%Type block
      !%Default no
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_local, and OCTLocalTarget = "OCTTargetDensityFromState",
      !% you must specify one OCTTargetDensityState block, in order to specify which linear
      !% combination of the states present in "restart/gs" is used to
      !% create the target density.
      !%
      !% The syntax is equivalent to the one used for the TransformStates
      !% block.
      !%End

      if(loct_parse_isdef('OCTTargetDensity').ne.0) then
        ALLOCATE(target%rho(NP), NP)
        target%rho = M_ZERO
        call loct_parse_string('OCTTargetDensity', "0", expression)


        if(trim(expression).eq.'OCTTargetDensityFromState') then

          if(loct_parse_block(check_inp('OCTTargetDensityFromState'), blk) == 0) then
            tmp_st = target%st
            deallocate(tmp_st%zpsi)
            call restart_look_and_read("tmp", tmp_st, gr, geo, ierr)
            ALLOCATE(rotation_matrix(target%st%nst, tmp_st%nst), target%st%nst*tmp_st%nst)
            rotation_matrix = M_z0
            do ist = 1, target%st%nst
              do jst = 1, loct_parse_block_cols(blk, ist-1)
                call loct_parse_block_cmplx(blk, ist-1, jst-1, rotation_matrix(ist, jst))
              end do
            end do
            call rotate_states(gr%m, target%st, tmp_st, rotation_matrix)
            deallocate(rotation_matrix)
            call states_calc_dens(target%st, NP_PART, target%st%rho)
            do ip = 1, NP
              target%rho(ip) = sum(target%st%rho(ip, 1:target%st%d%spin_channels))
            end do
            call states_end(tmp_st)
          else
            message(1) = '"OCTTargetDensityState" has to be specified as block.'
            call write_info(1)
            call input_error('OCTTargetDensity')
          end if

        else

          call conv_to_C_string(expression)
          do ip = 1, NP
            call mesh_r(gr%m, ip, r, x = x)
            ! parse user defined expression
            call loct_parse_expression(psi_re, psi_im, &
              x(1), x(2), x(3), r, M_ZERO, expression)
            target%rho(ip) = psi_re
          end do
          ! Normalize
          r = dmf_integrate(gr%m, target%rho)
          target%rho = (-target%st%val_charge) * target%rho/r
        end if

      else
        message(1) = 'If OCTTargetOperator = oct_tg_density, then you must give the shape'
        message(2) = 'of this target in variable "OCTTargetDensity"'
        call write_fatal(2)
      end if

    case(oct_tg_local)
      !%Variable OCTTargetLocal
      !%Type string
      !%Section Optimal Control
      !%Description
      !% If OCTTargetOperator = oct_tg_local, then one must supply the target density
      !% that should be searched for. This one can do by supplying a string through
      !% the variable OCTLocalTarget.
      !%End

      if(loct_parse_isdef('OCTTargetLocal').ne.0) then
        ALLOCATE(target%rho(NP), NP)
        target%rho = M_ZERO
        call loct_parse_string('OCTTargetLocal', "0", expression)
        call conv_to_C_string(expression)
        do ip = 1, NP
          call mesh_r(gr%m, ip, r, x = x)
          ! parse user defined expression
          call loct_parse_expression(psi_re, psi_im, &
            x(1), x(2), x(3), r, M_ZERO, expression)
          target%rho(ip) = psi_re
        end do
      else
        message(1) = 'If OCTTargetOperator = oct_tg_local, then you must give the shape'
        message(2) = 'of this target in variable "OCTTargetLocal"'
        call write_fatal(2)
      end if

    case default
      write(message(1),'(a)') "Target Operator not properly defined."
      call write_fatal(1)
    end select

    call tdtargetset_init(target, gr, td)

    if(target%no_tdtargets.eq.0) then
      nullify(target%td_fitness)
    else
      ALLOCATE(target%td_fitness(0:td%max_iter), td%max_iter+1)
    end if
    
    call pop_sub()
  end subroutine target_init


  subroutine target_end(target)
    type(target_t), intent(inout) :: target
    integer :: i

    call push_sub('target.target_end')

    call states_end(target%st)
    if(target%totype .eq. oct_tg_local .or. target%totype .eq. oct_tg_density) then
      deallocate(target%rho)
      nullify(target%rho)
    end if
    if(associated(target%td_fitness)) then
      deallocate(target%td_fitness); nullify(target%td_fitness)
    end if

    do i = 1, target%no_tdtargets
      call tdtarget_end(target%tdtg(i))
    end do
    if(target%no_tdtargets > 0) then 
      deallocate(target%tdtg); nullify(target%tdtg)
      deallocate(target%tdtarget); nullify(target%tdtarget)
    end if
    target%no_tdtargets = 0

    call pop_sub()
  end subroutine target_end


  ! ----------------------------------------------------------------------
  subroutine target_output(gr, dir, outp, target)
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(output_t),   intent(in)  :: outp
    type(target_t), intent(inout) :: target

    integer :: ierr
    call push_sub('target.target_output')

    select case(target%totype)
    case(oct_tg_local)
      call doutput_function(outp%how, trim(dir), 'local_target', gr%m, gr%sb, &
        target%rho, M_ONE, ierr, is_tmp = .false.)
    case(oct_tg_density)
      call doutput_function(outp%how, trim(dir), 'density_target', gr%m, gr%sb, &
        target%rho, M_ONE, ierr, is_tmp = .false.)
    case default
      call states_output(target%st, gr, trim(dir), outp)
    end select

    call pop_sub()
  end subroutine target_output



  ! ---------------------------------------------------------
  subroutine calc_tdfitness(target, gr, psi, i)
    type(target_t), intent(inout) :: target
    type(grid_t),      intent(in) :: gr
    type(states_t),    intent(in) :: psi
    integer, intent(in) :: i

    integer             :: jj, ik, p, dim

    call push_sub('target.calc_tdfitness')

    target%td_fitness(i) = M_ZERO
    do jj = 1, target%no_tdtargets
      if(target%tdtg(jj)%type.eq.oct_tgtype_local) then
        do ik = 1, psi%d%nik
          do p  = psi%st_start, psi%st_end
            do dim = 1, psi%d%dim
             target%td_fitness(i) = target%td_fitness(i) + &
                zmf_integrate(gr%m, target%tdtarget(:)* &
                abs(psi%zpsi(:,dim,ik,p))**2)
            end do
          end do
        end do
      else
        do ik = 1, psi%d%nik
          do p  = psi%st_start, psi%st_end
            do dim = 1, psi%d%dim
              target%td_fitness(i) = target%td_fitness(i) + &
                abs(zmf_integrate(gr%m, target%tdtarget(:)* &
                conjg(psi%zpsi(:,dim,ik,p))))**2

            end do
          end do
        end do
      end if
    end do

    call pop_sub()
  end subroutine calc_tdfitness


  ! ---------------------------------------------------------------
  subroutine calc_inh(psi_n, gr, target, iter, max_iter, dt, chi_n)
    type(states_t),    intent(in)        :: psi_n
    type(grid_t),      intent(in)        :: gr
    type(target_t),    intent(inout)     :: target
    integer,           intent(in)        :: iter
    integer,           intent(in)        :: max_iter
    FLOAT,             intent(in)        :: dt
    type(states_t),    intent(inout)     :: chi_n
    
    CMPLX               :: tgt(NP_PART)
    integer             :: jj, dim
    CMPLX               :: olap

    call push_sub('target.calc_inh')
    
    ! TODO: change to right dimensions
    !       build target operator
    
    ! FIXME: more flexibility when defining multiple target operators
    ! usually one has only one specie of operators, i.e., only local or only non-local, 
    ! when mixing these, this routine has to be improved.
    if(target%tdtg(1)%type.eq.oct_tgtype_local) then
      do jj = 1, target%no_tdtargets
        call build_tdtarget(target%tdtg(jj), gr, tgt, iter) ! tdtarget is build
        target%tdtarget(:) = target%tdtarget(:) + tgt(:) 
      end do
      do dim=1, chi_n%d%dim
        chi_n%zpsi(:,dim,1,1) = chi_n%zpsi(:,dim,1,1) &
          - sign(M_ONE,dt)/real(max_iter)*target%tdtarget(:)* &
          psi_n%zpsi(:,dim,1,1)
      enddo
    else 
      do jj = 1, target%no_tdtargets
        call build_tdtarget(target%tdtg(jj), gr, tgt, iter)
        target%tdtarget = tgt
        olap= m_z0
        do dim=1, psi_n%d%dim
          olap = zmf_integrate(gr%m, psi_n%zpsi(:,dim,1,1) * conjg(target%tdtarget(:)))
          chi_n%zpsi(:,dim,1,1) = chi_n%zpsi(:,dim,1,1) &
            - sign(M_ONE,dt)/real(max_iter)*olap*target%tdtarget(:)
        end do
      end do
    end if
    
    call pop_sub()
  end subroutine calc_inh


  !----------------------------------------------------------
  ! 
  !----------------------------------------------------------
  subroutine tdtargetset_init(target, gr, td)
    type(target_t), intent(inout) :: target
    type(grid_t),      intent(in) :: gr
    type(td_t),        intent(in) :: td

    integer             :: i, jj, iunit
    FLOAT               :: mxloc(MAX_DIM)
    integer             :: tt, pos(1)
    FLOAT               :: t
    CMPLX               :: tgt(NP_PART)
    character(len=80)  :: filename

    integer :: no_c, no_tds, pol
    C_POINTER         :: blk

    call push_sub('opt_control_tdtarget.tdtargetset_init')

    !%Variable OCTTdTarget
    !%Type block
    !%Section Optimal Control
    !%Description
    !% octopus features also time-dependent targets, i.e., one wants to optimize a laser 
    !% field that achieves a predefined time-dependent target. An example, could be the 
    !% evolution of occupation numbers in time. A time-dependent target consists of two 
    !% parts, i.e., the operator itself (a projection or a local operator) and its 
    !% time-dependence. Both are specified in one row of the block. You may enter as many 
    !% rows as you like. For time-dependent occupation targets OCTOPUS takes care of the 
    !% normalization.
    !% 
    !% The structure of the block is as follows:
    !%
    !% <tt>%OCTTdTarget
    !% <br>&nbsp;&nbsp;type | ftype | width | weight | g_x(t) | g_y(t) | g_z(t) 
    !% <br>%</tt>
    !%  
    !% Type:
    !% Choose betweem local and state target. 
    !% *oct_tgtype_local*
    !% *oct_tgtype_state*
    !%
    !% Ftype: 
    !% Spatial function of target operator.
    !% *oct_ftype_gauss*
    !%
    !% width:
    !% width of the Gaussian 
    !%
    !% weight:
    !% If multiple operators ae defined weight the importance between them
    !% g_x(t), g_y(t), g_z(t):
    !% describe the time-dependence
    !% 
    !% Example: 
    !% R0=1.0
    !% <tt>%OCTTdTarget
    !% <br>&nbsp;&nbsp;oct_tgtype_local | oct_ftype_gauss | 20 | 1 | "R0*sin(1.0*t)" | "R0*cos(1.0*t)" | "0"
    !% <br>%</tt>
    !% 
    !% In this example we define a narrow Gaussian which circles with radius R0 around the center. 
    !% 
    !%End
    no_tds           = 0
    target%no_tdtargets = no_tds
    if((target%targetmode==oct_targetmode_td) &
        .AND.(loct_parse_block(check_inp('OCTTdTarget'),blk)==0)) then

      no_tds = loct_parse_block_n(blk)
      if(no_tds > 0) then
        target%no_tdtargets = no_tds
        ALLOCATE(target%tdtg(no_tds), no_tds)
        do i=1, no_tds
          do pol=1, MAX_DIM
            target%tdtg(i)%expression(pol) = " "
          end do
          ! The structure of the block is:
          ! domain | function_type | center | width | weight 
          no_c = loct_parse_block_cols(blk, i-1)
          !td_tg(i)%type = oct_tgtype_local
          call loct_parse_block_int(blk, i-1, 0, target%tdtg(i)%type)
          call loct_parse_block_int(blk, i-1, 1, target%tdtg(i)%ftype)
          call loct_parse_block_float(blk, i-1, 2, target%tdtg(i)%width)
          call loct_parse_block_float(blk, i-1, 3, target%tdtg(i)%weight)

          do pol=1, NDIM
            call loct_parse_block_string(blk, i-1, 3+pol, target%tdtg(i)%expression(pol))
          end do
          !
          ALLOCATE(target%tdtg(i)%tdshape(NDIM,0:td%max_iter), NDIM*(td%max_iter+1))
          call build_tdshape(target%tdtg(i), gr, td%max_iter, td%dt)
        end do

      end if
      
    end if
    ! calc norm, give warning when to small

    if(no_tds.eq.0) then
      nullify(target%tdtarget)
      nullify(target%tdtg)
      call pop_sub(); return
    end if

    ! td target fitness
    ! to avoid double parsing and build up tdtarget array
    ALLOCATE(target%tdtarget(1:NP_PART), NP_PART)
    target%tdtarget   = M_z0

    ! generate some output to analyse the defined target
    do jj=1, target%no_tdtargets
       write(message(1),'(a,i2.2)') 'Info: Dumping td target ', jj
       call write_info(1)
       write(filename,'(a,i2.2)') 'opt-control/td_target_',jj
       iunit = io_open(filename, action='write')
       ! build target_function and calc inhomogeneity
       do tt=0, td%max_iter, 20 
          t = real(tt,REAL_PRECISION)*td%dt
          call build_tdtarget(target%tdtg(jj), gr, tgt, tt)
          target%tdtarget = tgt
          pos = maxloc(abs(target%tdtarget))
          mxloc(:) = gr%m%x(pos(1),:) 
       end do
       close(iunit)
    end do
  
    call pop_sub()
  end subroutine tdtargetset_init


end module opt_control_target_m
