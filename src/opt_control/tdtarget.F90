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

module opt_control_tdtarget_m
  use datasets_m
  use varinfo_m
  use messages_m
  use lib_oct_parser_m
  use lib_oct_m
  use io_m
  use global_m
  use grid_m
  use mesh_function_m
  use states_m
  use hamiltonian_m
  use timedep_m
  use td_write_m
  use opt_control_constants_m

  implicit none

  private
  public :: td_target_t,      &
            td_target_set_t,  &
            tdtargetset_init, &
            tdtargetset_end,  &
            calc_inh

  type td_target_t
    integer :: type       ! local or state 
    FLOAT   :: par 
    FLOAT   :: width      ! width parameter
    FLOAT   :: weight     ! relative weight of filter
    character(len=1024) :: expression(MAX_DIM) ! expression of user def func
    integer             :: ftype               ! which function: gaussian etc
    CMPLX, pointer :: tdshape(:,:) ! evaluate user defined functions
  end type td_target_t

  type td_target_set_t
    integer :: no_tdtargets
    type(td_target_t), pointer :: tdtg(:)
    FLOAT, pointer :: td_fitness(:)
    CMPLX, pointer :: tdtarget(:)
  end type td_target_set_t

  contains

  !----------------------------------------------------------
  ! 
  !----------------------------------------------------------
  subroutine tdtargetset_init(targetmode, gr, td, tdt)
    integer, intent(in) :: targetmode
    type(grid_t),      intent(in) :: gr
    type(td_t),        intent(in) :: td
    type(td_target_set_t), intent(inout) :: tdt

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
    no_tds = 0
    if((targetmode==oct_targetmode_td) &
        .AND.(loct_parse_block(check_inp('OCTTdTarget'),blk)==0)) then

      no_tds = loct_parse_block_n(blk)
      if(no_tds > 0) then
        tdt%no_tdtargets = no_tds
        ALLOCATE(tdt%tdtg(no_tds), no_tds)
        do i=1, no_tds
          do pol=1, MAX_DIM
            tdt%tdtg(i)%expression(pol) = " "
          end do
          ! The structure of the block is:
          ! domain | function_type | center | width | weight 
          no_c = loct_parse_block_cols(blk, i-1)
          !td_tg(i)%type = oct_tgtype_local
          call loct_parse_block_int(blk, i-1, 0, tdt%tdtg(i)%type)
          call loct_parse_block_int(blk, i-1, 1, tdt%tdtg(i)%ftype)
          call loct_parse_block_float(blk, i-1, 2, tdt%tdtg(i)%width)
          call loct_parse_block_float(blk, i-1, 3, tdt%tdtg(i)%weight)
          do pol=1, NDIM
            call loct_parse_block_string(blk, i-1, 3+pol, tdt%tdtg(i)%expression(pol))
          end do
          !
          ALLOCATE(tdt%tdtg(i)%tdshape(NDIM,0:td%max_iter), NDIM*(td%max_iter+1))
          call build_tdshape(tdt%tdtg(i), gr, td%max_iter, td%dt)
        end do

      end if
      
    end if
    ! calc norm, give warning when to small

    if(no_tds.eq.0) then
      nullify(tdt%td_fitness); nullify(tdt%tdtarget); nullify(tdt%tdtg)
      call pop_sub(); return
    end if

    ! td target fitness
    ALLOCATE(tdt%td_fitness(0:td%max_iter), td%max_iter+1)
    tdt%td_fitness = M_ZERO
    ! to avoid double parsing and build up tdtarget array
    ALLOCATE(tdt%tdtarget(1:NP_PART), NP_PART)
    tdt%tdtarget   = M_z0

    ! generate some output to analyse the defined target
    do jj=1, tdt%no_tdtargets
       write(message(1),'(a,i2.2)') 'Info: Dumping td target ', jj
       call write_info(1)
       write(filename,'(a,i2.2)') 'opt-control/td_target_',jj
       iunit = io_open(filename, action='write')
       ! build target_function and calc inhomogeneity
       do tt=0, td%max_iter, 20 
          t = real(tt,REAL_PRECISION)*td%dt
          !write(6,*) tt, trim(td_tg(jj)%expression(1))
          !write(6,*) tt, trim(td_tg(jj)%expression(2))
          !write(6,*) tt, trim(td_tg(jj)%expression(3))
          call build_tdtarget(tdt%tdtg(jj), gr, tgt, tt)
          tdt%tdtarget = tgt
          pos = maxloc(abs(tdt%tdtarget))
          !write(6,*) pos
          mxloc(:) = gr%m%x(pos(1),:) 
       end do
       close(iunit)
    end do
  
    call pop_sub()
  end subroutine tdtargetset_init


  !----------------------------------------------------------
  subroutine tdtargetset_end(tdt)
    type(td_target_set_t), intent(inout) :: tdt
    integer :: i
    call push_sub('opt_control_tdtarget.tdtargetset_end')

    do i = 1, tdt%no_tdtargets
      call tdtarget_end(tdt%tdtg(i))
    end do
    if(tdt%no_tdtargets > 0) then 
      deallocate(tdt%tdtg); nullify(tdt%tdtg)
      deallocate(tdt%tdtarget); nullify(tdt%tdtarget)
      deallocate(tdt%td_fitness); nullify(tdt%td_fitness)
    end if

  end subroutine tdtargetset_end


  !----------------------------------------------------------
  subroutine tdtarget_end(tdtg)
    type(td_target_t), intent(inout) :: tdtg
    call push_sub('opt_control_tdtarget.tdtarget_end')
    deallocate(tdtg%tdshape); nullify(tdtg%tdshape)
    call pop_sub()
  end subroutine tdtarget_end
  



  !----------------------------------------------------------
  ! Parses tdtd%expression to build tgtd%tdshape 
  !----------------------------------------------------------
  subroutine build_tdshape(tdtg, gr, steps, dt)
    type(td_target_t), intent(inout) :: tdtg
    type(grid_t),      intent(in)    :: gr
    integer,           intent(in)    :: steps
    FLOAT,             intent(in)    :: dt

    FLOAT   :: f_re, f_im, t
    integer :: kk, pol
    call push_sub('opt_control_tdtarget.build_tdshape')

    do kk=0, steps
      t = (kk-1)*steps
      do pol=1, NDIM
        f_re = M_ZERO
        f_im = M_ZERO
        call loct_parse_expression(f_re, f_im, & 	
          M_ZERO, M_ZERO, M_ZERO, M_ZERO, t, &
          tdtg%expression(pol)) 
        tdtg%tdshape(pol,kk) = f_re + M_zI*f_im
      end do
    end do

    call pop_sub()
  end subroutine build_tdshape


  ! ---------------------------------------------------------
  subroutine build_tdtarget(tdtg, gr, tgt, iter)
    type(td_target_t), intent(in)   :: tdtg
    type(grid_t),      intent(in)   :: gr
    integer,           intent(in)   :: iter
    CMPLX  ,           intent(out)  :: tgt(:)

    integer :: pol
    call push_sub('opt_control_tdtarget.build_tdtarget')
    
    select case(tdtg%ftype)
    case(1)
      ! gauss
      ! Product x = gr%m%x(ip, :) 
      tgt = M_ONE
      do pol=1, NDIM
        tgt = tgt * exp( - (gr%m%x(:, pol) - &
          real(tdtg%tdshape(pol,iter),REAL_PRECISION ))**2/(M_TWO*tdtg%width**2))
      enddo
      
      !case(2)
      ! tgt = M-zero
      !   ! step
      !   Nwidth = tg(jj)%width
      !   Nwmin  = tg(jj)%numerical
      !   Nwmax  = 
      !   do kk=1,  gr%m%np  tg(jj)%width  tg(jj)%tdshape
      !   
      !   end do
    case default
      message(1)="Unknown TD Target Operator: Choose Gaussian(1) or Step(2)"
      call write_fatal(1)
    end select
    
    call pop_sub()
  end subroutine build_tdtarget


  ! ---------------------------------------------------------------
  subroutine calc_inh(psi_n, gr, tdt, iter, max_iter, dt, chi_n)
    type(states_t),    intent(in)        :: psi_n
    type(grid_t),      intent(in)        :: gr
    type(td_target_set_t), intent(inout) :: tdt
    integer,           intent(in)        :: iter
    integer,           intent(in)        :: max_iter
    FLOAT,             intent(in)        :: dt
    type(states_t),    intent(inout)     :: chi_n
    
    CMPLX               :: tgt(NP_PART)
    integer             :: jj, dim
    CMPLX               :: olap

    call push_sub('opt_control_tdtarget.calc_inh')
    
    ! TODO: change to right dimensions
    !       build target operator
    
    ! FIXME: more flexibility when defining multiple target operators
    ! usually one has only one specie of operators, i.e., only local or only non-local, 
    ! when mixing these, this routine has to be improved.
    if(tdt%tdtg(1)%type.eq.oct_tgtype_local) then
      do jj = 1, tdt%no_tdtargets
        call build_tdtarget(tdt%tdtg(jj), gr, tgt, iter) ! tdtarget is build
        tdt%tdtarget(:) = tdt%tdtarget(:) + tgt(:) 
      end do
      do dim=1, chi_n%d%dim
        chi_n%zpsi(:,dim,1,1) = chi_n%zpsi(:,dim,1,1) &
          - sign(M_ONE,dt)/real(max_iter)*tdt%tdtarget(:)* &
          psi_n%zpsi(:,dim,1,1)
      enddo
    else 
      do jj = 1, tdt%no_tdtargets
        call build_tdtarget(tdt%tdtg(jj), gr, tgt, iter)
        tdt%tdtarget = tgt
        olap= m_z0
        do dim=1, psi_n%d%dim
          olap = zmf_integrate(gr%m,psi_n%zpsi(:,dim,1,1)*conjg(tdt%tdtarget(:)))
          chi_n%zpsi(:,dim,1,1) = chi_n%zpsi(:,dim,1,1) &
            - sign(M_ONE,dt)/real(max_iter)*olap*tdt%tdtarget(:)
        end do
      end do
    end if
    
    call pop_sub()
  end subroutine calc_inh


end module opt_control_tdtarget_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

