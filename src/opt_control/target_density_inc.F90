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
!! $Id: target_density_inc.F90 $

  ! ----------------------------------------------------------------------
  !> 
  subroutine target_init_density(gr, tg, stin, td)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg
    type(states_t),   intent(inout) :: stin 
    type(td_t),       intent(in)    :: td

    integer             :: ip, ist, jst, cstr_dim(MAX_DIM), ib, idim, jj, no_constraint, no_ptpair
    type(block_t)       :: blk
    FLOAT               :: psi_re, psi_im, xx(MAX_DIM), rr, fact, xend, xstart
    FLOAT, allocatable  :: xp(:), tmp_box(:, :)
    CMPLX, allocatable  :: rotation_matrix(:, :)
    type(states_t)      :: tmp_st
    character(len=1024) :: expression

    PUSH_SUB(target_init_density)

    message(1) =  'Info: Target is a combination of a density and/or a current functional.'
    call messages_info(1)

    tg%curr_functional = oct_no_curr
    tg%curr_weight = M_ZERO
    tg%strt_iter_curr_tg = 0

    !%Variable OCTTargetDensity
    !%Type string
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If <tt>OCTTargetOperator = oct_tg_density</tt>, then one must supply the target density
    !% that should be searched for. This one can do by supplying a string through
    !% the variable <tt>OCTTargetDensity</tt>.
    !%End

    !%Variable OCTTargetDensityFromState
    !%Type block
    !%Default no
    !%Section Calculation Modes::Optimal Control
    !%Description
    !% If <tt>OCTTargetOperator = oct_tg_density</tt>, and <tt>OCTLocalTarget = "OCTTargetDensityFromState"</tt>,
    !% you must specify a <tt>OCTTargetDensityState</tt> block, in order to specify which linear
    !% combination of the states present in <tt>restart/gs</tt> is used to
    !% create the target density.
    !%
    !% The syntax is the same as the <tt>TransformStates</tt> block.
    !%End

    if(parse_isdef('OCTTargetDensity') /= 0) then
      tg%density_weight = M_ONE
      SAFE_ALLOCATE(tg%rho(1:gr%mesh%np))
      tg%rho = M_ZERO
      call parse_string(datasets_check('OCTTargetDensity'), "0", expression)

      if(trim(expression)  ==  'OCTTargetDensityFromState') then

        if(parse_block(datasets_check('OCTTargetDensityFromState'), blk) == 0) then
          call states_copy(tmp_st, tg%st)
          SAFE_DEALLOCATE_P(tmp_st%zpsi)
          call restart_look_and_read(tmp_st, gr)
          SAFE_ALLOCATE(rotation_matrix(1:tg%st%nst, 1:tmp_st%nst))
          rotation_matrix = M_z0
          do ist = 1, tg%st%nst
            do jst = 1, parse_block_cols(blk, ist - 1)
              call parse_block_cmplx(blk, ist - 1, jst - 1, rotation_matrix(ist, jst))
            end do
          end do
          call states_rotate(gr%mesh, tg%st, tmp_st, rotation_matrix)
          SAFE_DEALLOCATE_A(rotation_matrix)
          call density_calc(tg%st, gr, tg%st%rho)
          do ip = 1, gr%mesh%np
            tg%rho(ip) = sum(tg%st%rho(ip, 1:tg%st%d%spin_channels))
          end do
          call states_end(tmp_st)
          call parse_block_end(blk)
        else
          message(1) = '"OCTTargetDensityState" has to be specified as block.'
          call messages_info(1)
          call input_error('OCTTargetDensity')
        end if

      else

        call conv_to_C_string(expression)
        do ip = 1, gr%mesh%np
          call mesh_r(gr%mesh, ip, rr, coords = xx)
          ! parse user-defined expression
          call parse_expression(psi_re, psi_im, gr%sb%dim, xx, rr, M_ZERO, expression)
          tg%rho(ip) = psi_re
        end do
        ! Normalize
        rr = dmf_integrate(gr%mesh, tg%rho)
        tg%rho = (-tg%st%val_charge) * tg%rho/rr
      end if

    else
      tg%density_weight = M_ZERO
    end if

      !%Variable OCTCurrentFunctional
      !%Type integer
      !%Section Calculation Modes::Optimal Control
      !%Default oct_no_curr
      !%Description
      !% The variable <tt>OCTCurrentFunctional</tt> describes which kind of current target functional J1_c[j] is
      !% to be used. EXPERIMENTAL!
      !%Option oct_no_curr 0
      !% No current functional is used, no current calculated.
      !%Option oct_curr_square 1
      !% Calculates the square of current j: J1_c[j] = <tt>OCTCurrentWeight</tt>*\int|j(r)|^2 dr. 
      !% For <tt>OCTCurrentWeight</tt> .LT. 0 the current will be minimized (useful in combination with 
      !% target density in order to obtain stable final target density), while for 
      !% OCTCurrentWeight</tt> .GT. 0 it will be maximized (useful in combination with a target density 
      !% in order to obtain a high-velocity impact, for instance). It is a static target, to be reached at
      !% total time. 
      !%Option oct_max_curr_ring 2
      !% Maximizes the current of a quantum ring in one direction. The functional maximizes the z projection of the 
      !% outer product between the position \vec{r} and the current \vec{j}: 
      !% J1[j] = <tt>OCTCurrentWeight</tt>*\int (\vec{r} \times \vec{j}) \hat{z} dr. For <tt>OCTCurrentWeight</tt> .GT. 0. the
      !% current flows in counter-clockwise direction, while for <tt>OCTCurrentWeight</tt> .LT. 0 the current is clockwise.
      !%Option oct_curr_square_td 3
      !% The time dependent version of <tt>oct_curr_square</tt>. In fact, calculates the 
      !% square of current in time interval [<tt>OCTStartTimeCurrTg</tt>, 
      !% total time = <tt>TDMaximumIter</tt> * <tt>TDTimeStep</tt>]. 
      !% Set <tt>TDPropagator</tt> = <tt>crank_nicolson</tt>.
      !%End 

      !%Variable OCTCurrentWeight
      !%Type float
      !%Section Calculation Modes::Optimal Control
      !%Default 0.0
      !%Description
      !% In the case of simultaneous optimization of density n and current j, one can tune the importance
      !% of the current functional J1_c[j], as the respective functionals might not provide results on the
      !% same scale of magnitude. J1[n,j]= J1_d[n]+ <tt>OCTCurrentWeight</tt> * J1_c[j]. Be aware that its
      !% sign is crucial for the chosen <tt>OCTCurrentFunctional</tt> as explained there.
      !%End

      !%Variable OCTStartIterCurrTg
      !%Type integer
      !%Section Calculation Modes::Optimal Control
      !%Default 0
      !%Description
      !% Allows for a time dependent target for the current without defining it for the total 
      !% time-interval of the simulation.
      !% Thus it can be switched on at the iteration desired, <tt>OCTStartIterCurrTg</tt> >= 0
      !% and  <tt>OCTStartIterCurrTg</tt>  <  <tt>TDMaximumIter</tt>. 
      !% Tip: If you would like to specify a real time for switching
      !% the functional on rather than the number of steps, just use something
      !% like:
      !% <tt>OCTStartIterCurrTg</tt> = 100.0 / <tt>TDTimeStep</tt>
      !%End
 
      !%Variable OCTSpatialCurrWeight
      !%Type block
      !%Section Calculation Modes::Optimal Control
      !%Description
      !% Can be seen as a position dependent <tt>OCTCurrentWeight</tt>. Consequently, it
      !% weights contribution of current j to its functional J1_c[j] according to the position in space. 
      !% For example, <tt>oct_curr_square</tt> thus becomes 
      !% J1_c[j] = <tt>OCTCurrentWeight</tt> \int |j(r)|^2 <tt>OCTSpatialCurrWeight</tt>(r) dr.
      !%
      !% It is defined <tt>OCTSpatialCurrWeight</tt>(r) = g(x)*g(y)*g(z), where   
      !% g(x) = \sum_{i} 1/(1+exp( -fact*(x-startpoint_i) )) - 1/(1+exp( -fact*(x-endpoint_i) )).
      !% If not specified, g(x) = 1. 
      !% 
      !% Each g(x) is represented by one line of the block that has the following form
      !%
      !% <tt>%OCTSpatialCurrWeight
      !% <br>&nbsp;&nbsp;  dimension  |  fact |  startpoint_1  | endpoint_1  | startpoint_2 | endpoint_2 |...
      !% <br>%</tt>
      !%
      !% There are no restrictions on the number of lines, nor on the number of pairs of start- and endpoints. 
      !% Attention: <tt>startpoint</tt> and <tt>endpoint</tt> have to be supplied pairwise 
      !% with <tt>startpoint  <  endpoint</tt>. <tt>dimension > 0</tt> is integer, <tt>fact</tt> is float.
      !%End

    call parse_integer(datasets_check('OCTCurrentFunctional'), oct_no_curr, tg%curr_functional)
    select case(tg%curr_functional)
    case(oct_no_curr)
    case(oct_curr_square, oct_max_curr_ring, oct_curr_square_td)
      SAFE_ALLOCATE(stin%current( 1:gr%mesh%np_part, 1:gr%mesh%sb%dim, 1:stin%d%nspin ) )
      stin%current= M_ZERO
    end select

    call parse_float(datasets_check('OCTCurrentWeight'), M_ZERO, tg%curr_weight)
    write(message(1), '(a,i3)')   'Info: OCTCurrentFunctional = ', tg%curr_functional
    write(message(2), '(a,f8.3)') 'Info: OCTCurrentWeight = ',  tg%curr_weight
    call messages_info(2)

    if (target_mode(tg)  ==  oct_targetmode_td) then
      call parse_integer(datasets_check('OCTStartIterCurrTg'), 0, tg%strt_iter_curr_tg)
      if (tg%strt_iter_curr_tg  <  0) then
        message(1) = 'OCTStartIterCurrTg must be positive.'
        call messages_fatal(1)
      elseif (tg%strt_iter_curr_tg >= td%max_iter) then
        message(1) = 'OCTStartIterCurrTg has to be  <  TDMaximumIter.'
        call messages_fatal(1)
      end if
      write(message(1), '(a,i3)')   'Info: TargetMode = ', target_mode(tg)
      write(message(2), '(a,i8)') 'Info: OCTStartIterCurrTg = ',  tg%strt_iter_curr_tg
      call messages_info(2)
      tg%dt = td%dt
      SAFE_ALLOCATE(tg%td_fitness(0:td%max_iter))
      tg%td_fitness = M_ZERO
    else
      tg%strt_iter_curr_tg = 0
    end if

    if(parse_isdef(datasets_check('OCTSpatialCurrWeight')) /= 0) then
      if(parse_block(datasets_check('OCTSpatialCurrWeight'), blk) == 0) then
        SAFE_ALLOCATE(tg%spatial_curr_wgt(1:gr%mesh%np_part))
        SAFE_ALLOCATE(xp(1:gr%mesh%np_part))
        SAFE_ALLOCATE(tmp_box(1:gr%mesh%np_part, 1:gr%mesh%sb%dim))
          
        no_constraint = parse_block_n(blk)
        tmp_box = M_ZERO
        cstr_dim = 0
        do ib = 1, no_constraint
          call parse_block_integer(blk, ib - 1, 0, idim)
          if( idim  <=  0 .or. idim > gr%mesh%sb%dim) then
            write(message(1), '(a,i3)') 'Error in "OCTSpatialCurrWeight" block, line:', ib
            write(message(2), '(a)'   ) '"dimension" has to be positive'
            write(message(3), '(a)'   ) 'and must not exceed dimensions of the system.'
            call messages_fatal(3)
          end if
          cstr_dim(idim) = 1
          xp(1:gr%mesh%np_part) = gr%mesh%x(1:gr%mesh%np_part, idim)

          call parse_block_float(blk, ib - 1, 1, fact)

          no_ptpair = parse_block_cols(blk, ib-1) - 2
          if (mod(no_ptpair,2) /= 0) then
            write(message(1), '(a,i3)') 'Error in "OCTSpatialCurrWeight" block, line:', ib
            write(message(2), '(a)'   ) 'Each interval needs start and end point!'
            call messages_fatal(2)
          end if
              
          do jj= 2, no_ptpair, 2
            call parse_block_float(blk, ib - 1, jj, xstart)
            call parse_block_float(blk, ib - 1, jj+1, xend)
                           
            if (xstart >= xend) then
              write(message(1), '(a,i3)') 'Error in "OCTSpatialCurrWeight" block, line:', ib
              write(message(2), '(a)'   ) 'Set "startpoint"  <  "endpoint" ' 
              call messages_fatal(2)
            end if

            do ip = 1, gr%mesh%np_part
              tmp_box(ip,idim) = tmp_box(ip,idim) + M_ONE/(M_ONE+exp(-fact*(xp(ip)-xstart) )) -  &
                                                    M_ONE/(M_ONE+exp(-fact*(xp(ip)-xend) ))
            end do
          end do
            
        end do
          
        do idim = 1, gr%mesh%sb%dim
          if(cstr_dim(idim) == 0) tmp_box(:,idim) = M_ONE
        end do
        tg%spatial_curr_wgt(1:gr%mesh%np_part) = product(tmp_box(1:gr%mesh%np_part, 1:gr%mesh%sb%dim),2) 
        SAFE_DEALLOCATE_A(xp)
        SAFE_DEALLOCATE_A(tmp_box)
                             
        call parse_block_end(blk)     
      else
        message(1) = '"OCTSpatialCurrWeight" has to be specified as a block.'
        call messages_info(1)
        call input_error('OCTEvalBoxCurrTg')
      end if
    end if

    POP_SUB(target_init_density)
  end subroutine target_init_density


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_end_density(tg)
    type(target_t),   intent(inout) :: tg
    PUSH_SUB(target_end_density)
    SAFE_DEALLOCATE_P(tg%rho)
    SAFE_DEALLOCATE_P(tg%spatial_curr_wgt)
    select case(tg%curr_functional)
    case(oct_curr_square_td) 
      SAFE_DEALLOCATE_P(tg%td_fitness)
    end select
    POP_SUB(target_end_density)
  end subroutine target_end_density


  ! ----------------------------------------------------------------------
  subroutine target_output_density(tg, gr, dir, geo, outp)
    type(target_t), intent(inout) :: tg
    type(grid_t), intent(inout)   :: gr
    character(len=*), intent(in)  :: dir
    type(geometry_t),       intent(in)  :: geo
    type(output_t),         intent(in)  :: outp

    integer :: ierr
    PUSH_SUB(target_output_density)
    
    call loct_mkdir(trim(dir))
    if(outp%how /= 0) then
      if(tg%density_weight > M_ZERO) then
        call dio_function_output(outp%how, trim(dir), 'density_target', gr%mesh, &
          tg%rho, units_out%length**(-gr%sb%dim), ierr, geo = geo)
      end if
    end if

    POP_SUB(target_output_density)
  end subroutine target_output_density
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> 
  FLOAT function target_j1_density(gr, tg, psi) result(j1)
    type(grid_t),     intent(in)    :: gr
    type(target_t),   intent(inout) :: tg
    type(states_t), intent(inout)   :: psi

    integer :: ip, maxiter
    FLOAT :: currfunc_tmp
    FLOAT, allocatable :: local_function(:)

    PUSH_SUB(target_j1_density)

    if (tg%density_weight > M_ZERO) then
      SAFE_ALLOCATE(local_function(1:gr%mesh%np))
      do ip = 1, gr%mesh%np
        local_function(ip) = - ( sqrt(psi%rho(ip, 1)) - sqrt(tg%rho(ip)) )**2
      end do
      j1 = tg%density_weight * dmf_integrate(gr%mesh, local_function)
      SAFE_DEALLOCATE_A(local_function)
    else
      j1 = M_ZERO
    end if

    ! current functionals are conveniently combined with others
    if (tg%curr_functional /= oct_no_curr) then
      select case(target_mode(tg))
      case(oct_targetmode_static)
        currfunc_tmp = jcurr_functional(tg, gr, psi )
      case(oct_targetmode_td)
        maxiter = size(tg%td_fitness) - 1
        currfunc_tmp = M_HALF * tg%dt * tg%td_fitness(tg%strt_iter_curr_tg) + & 
                       M_HALF * tg%dt * tg%td_fitness(maxiter) + & 
                       tg%dt * sum(tg%td_fitness(tg%strt_iter_curr_tg+1:maxiter-1))  
      end select
      if(conf%devel_version) then
        write(message(1), '(6x,a,f12.5)')    " => Other functional   = ", j1
        write(message(2), '(6x,a,f12.5)')    " => Current functional = ", currfunc_tmp
        call messages_info(2)
      end if
      ! accumulating functional values
      j1 = j1 + currfunc_tmp
    end if 

    POP_SUB(target_j1_density)
  end function target_j1_density


  ! ----------------------------------------------------------------------
  !> 
  subroutine target_chi_density(tg, gr, psi_in, chi_out)
    type(target_t),    intent(inout) :: tg
    type(grid_t),      intent(inout) :: gr
    type(states_t),    intent(inout) :: psi_in
    type(states_t),    intent(inout) :: chi_out

    integer :: no_electrons, ip, ist

    PUSH_SUB(target_chi_density)

    no_electrons = -nint(psi_in%val_charge)

    chi_out%zpsi = M_ZERO

    if( tg%density_weight > M_ZERO ) then

    select case(psi_in%d%ispin)
    case(UNPOLARIZED)
      ASSERT(psi_in%d%nik  ==  1)

      if(no_electrons  ==  1) then
        do ip = 1, gr%mesh%np
          chi_out%zpsi(ip, 1, 1, 1) = sqrt(tg%rho(ip)) * &
            exp(M_zI * atan2(aimag(psi_in%zpsi(ip, 1, 1, 1)), &
                              real(psi_in%zpsi(ip, 1, 1, 1) )) )
        end do
      else
        do ist = psi_in%st_start, psi_in%st_end
          do ip = 1, gr%mesh%np
            if(psi_in%rho(ip, 1) > CNST(1.0e-8)) then
              chi_out%zpsi(ip, 1, ist, 1) = psi_in%occ(ist, 1) * sqrt(tg%rho(ip) / psi_in%rho(ip, 1)) * &
                psi_in%zpsi(ip, 1, ist, 1)
            else
              chi_out%zpsi(ip, 1, ist, 1) = M_ZERO !sqrt(tg%rho(ip))
            end if
          end do
        end do
      end if

    case(SPIN_POLARIZED)
       message(1) = 'Error in target.target_chi_density: spin_polarized.'
       call messages_fatal(1)
    case(SPINORS)
       message(1) = 'Error in target.target_chi_density: spinors.'
       call messages_fatal(1)
    end select

    end if

    if(tg%curr_functional /= oct_no_curr) then
      if (target_mode(tg)  ==  oct_targetmode_static ) then
        chi_out%zpsi = chi_out%zpsi + chi_current(tg, gr, psi_in)
      end if
    end if 

    POP_SUB(target_chi_density)
  end subroutine target_chi_density


  ! ---------------------------------------------------------
  !> 
  !!
  subroutine target_tdcalc_density(tg, gr, psi, time)
    type(target_t),      intent(inout) :: tg
    type(grid_t),        intent(inout) :: gr
    type(states_t),      intent(inout) :: psi
    integer,             intent(in)    :: time

    PUSH_SUB(target_tdcalc_density)

    if (time >= tg%strt_iter_curr_tg) then
      tg%td_fitness(time) = jcurr_functional(tg, gr, psi)
    end if 


    POP_SUB(target_tdcalc_density)
  end subroutine target_tdcalc_density
  ! ----------------------------------------------------------------------


  ! ----------------------------------------------------------------------
  !> Calculates a current functional that may be combined with
  !! other functionals found in function target_j1.
  FLOAT function jcurr_functional(tg, gr, psi) result(jcurr)
    type(target_t), intent(in)    :: tg
    type(grid_t),   intent(in)    :: gr
    type(states_t), intent(inout) :: psi

    integer :: ip
    FLOAT, allocatable :: semilocal_function(:)
 
    PUSH_SUB(jcurr_functional)
    
    jcurr = M_ZERO    
    ASSERT(psi%d%nik  ==  1)
    SAFE_ALLOCATE(semilocal_function(1:gr%mesh%np))
    semilocal_function = M_ZERO

    select case(tg%curr_functional)
    case(oct_no_curr)
      semilocal_function = M_ZERO

    case(oct_curr_square,oct_curr_square_td)
      call states_calc_quantities(gr%der, psi, paramagnetic_current=psi%current) 
      do ip = 1, gr%mesh%np
        semilocal_function(ip) =  sum(psi%current(ip, 1:gr%sb%dim, 1)**2)  
      end do
      
    case(oct_max_curr_ring)
      call states_calc_quantities(gr%der, psi, paramagnetic_current=psi%current) 

      if(gr%sb%dim /= M_TWO) then
        call messages_not_implemented('Target for dimension != 2')
      end if

      do ip = 1, gr%mesh%np
        ! func = j_y * x - j_x * y 
        semilocal_function (ip) = psi%current(ip, 2, 1) *  gr%mesh%x(ip,1) -  &
                                    psi%current(ip, 1, 1) * gr%mesh%x(ip,2)
      end do
    case default
      message(1) = 'Error in target.jcurr_functional: chosen target does not exist'
      call messages_fatal(1)
    end select

 
    if( is_spatial_curr_wgt(tg) ) then
      do ip = 1, gr%mesh%np
        semilocal_function(ip) = semilocal_function(ip) * tg%spatial_curr_wgt(ip) 
      end do
    end if

    jcurr = tg%curr_weight * dmf_integrate(gr%mesh, semilocal_function)
        
    SAFE_DEALLOCATE_A(semilocal_function)

    POP_SUB(jcurr_functional)
  end function jcurr_functional
  !-----------------------------------------------------------------------

  ! ----------------------------------------------------------------------
  ! Calculates current-specific boundary condition
  !-----------------------------------------------------------------------
  function chi_current(tg, gr, psi_in) result(chi)
    type(target_t),    intent(in)    :: tg
    type(grid_t),      intent(in)    :: gr
    type(states_t),    intent(inout) :: psi_in
    CMPLX                            :: chi( size(psi_in%zpsi,1), 1, size(psi_in%zpsi,3), 1)
        
    CMPLX, allocatable :: grad_psi_in(:,:,:)
    FLOAT, allocatable :: div_curr_psi_in(:,:)   
    integer :: ip, ist, idim

    PUSH_SUB(chi_current)

    chi = M_ZERO    

    SAFE_ALLOCATE(grad_psi_in(1:gr%der%mesh%np_part, 1:gr%der%mesh%sb%dim, 1))

    if(target_mode(tg)  ==  oct_targetmode_td ) then 
      call states_calc_quantities(gr%der, psi_in, paramagnetic_current=psi_in%current) 
    end if

    select case(tg%curr_functional)
    case(oct_no_curr)
      chi(ip, 1, ist, 1) = M_ZERO

    case(oct_curr_square,oct_curr_square_td)
       ! components current weighted by its position in the mesh, np_part included,
       ! since needed for the divergence of current.
       if( is_spatial_curr_wgt(tg) ) then
        do idim = 1, gr%sb%dim
          do ip = 1, gr%mesh%np_part 
            psi_in%current(ip, idim, 1) = psi_in%current(ip, idim, 1) * tg%spatial_curr_wgt(ip) 
          end do
        end do
      end if
      
      SAFE_ALLOCATE(div_curr_psi_in(1:gr%der%mesh%np_part,1))
      call dderivatives_div(gr%der, psi_in%current(1:gr%der%mesh%np_part, 1:gr%mesh%sb%dim, 1), &
                                    div_curr_psi_in(1:gr%der%mesh%np_part,1)) 
      
      ! the boundary condition  
      do ist = psi_in%st_start, psi_in%st_end
        call zderivatives_grad(gr%der, psi_in%zpsi(1:gr%der%mesh%np_part,1, ist, 1), & 
                               grad_psi_in(1:gr%der%mesh%np_part, 1:gr%der%mesh%sb%dim,1))
            
        do ip = 1, gr%mesh%np 
          chi(ip, 1, ist, 1) =  -M_zI * tg%curr_weight  * &
               ( M_TWO * sum(psi_in%current(ip, 1:gr%sb%dim, 1) * grad_psi_in(ip, 1:gr%sb%dim, 1))+ &
               div_curr_psi_in(ip,1) * psi_in%zpsi(ip, 1, ist, 1) )
        end do
      end do
      SAFE_DEALLOCATE_A(div_curr_psi_in)

    case(oct_max_curr_ring)

      if( is_spatial_curr_wgt(tg) ) then

        do ist = psi_in%st_start, psi_in%st_end
          call zderivatives_grad(gr%der, psi_in%zpsi(:,1, ist, 1), grad_psi_in(:,:,1))
          do ip = 1, gr%mesh%np 
            chi(ip, 1, ist, 1) =  M_zI * tg%curr_weight * tg%spatial_curr_wgt(ip) * &
                 ( grad_psi_in(ip, 1, 1)  * gr%mesh%x(ip,2) - &
                 grad_psi_in(ip, 2, 1)  * gr%mesh%x(ip,1)   ) 
          end do
        end do

      else

        do ist = psi_in%st_start, psi_in%st_end
          call zderivatives_grad(gr%der, psi_in%zpsi(:,1, ist, 1), grad_psi_in(:,:,1))
          do ip = 1, gr%mesh%np 
            chi(ip, 1, ist, 1) =  M_zI * tg%curr_weight * &
                 ( grad_psi_in(ip, 1, 1)  * gr%mesh%x(ip,2) - &
                 grad_psi_in(ip, 2, 1)  * gr%mesh%x(ip,1)   ) 
          end do
        end do
      end if

    case default
      message(1) = 'Error in target.chi_current: chosen target does not exist'
      call messages_fatal(1)
    end select  

    SAFE_DEALLOCATE_A(grad_psi_in)       
    POP_SUB(chi_current)
  end function chi_current
  ! ----------------------------------------------------------------------



!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
