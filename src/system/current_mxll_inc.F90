


    !----------------------------------------------------------
  subroutine get_rs_density_ext(st, mesh, time, rs_current_density_ext)
    type(states_mxll_t),  intent(inout) :: st
    type(mesh_t),    intent(in)    :: mesh
    FLOAT,           intent(in)    :: time
    CMPLX, optional, intent(inout) :: rs_current_density_ext(:,:)

    integer :: ip, idim
    FLOAT, allocatable :: current(:,:,:)

    PUSH_SUB(get_rs_density_ext)

    SAFE_ALLOCATE(current(1:mesh%np_part, 1:mesh%sb%dim,1))  !< The 1 in the last column is a dummy to use batch routines

    call external_current_calculation(st, mesh, time, current(:,:,1))
    call build_rs_current_state(current(:,:,1), rs_current_density_ext(:,:), st%ep(:), mesh%np_part)
    rs_current_density_ext = - rs_current_density_ext

    SAFE_DEALLOCATE_A(current)

    POP_SUB(get_rs_density_ext)
  end subroutine get_rs_density_ext


    !----------------------------------------------------------
  subroutine external_current_init(st, namespace, mesh)
    type(states_mxll_t), intent(inout) :: st
    type(mesh_t),   intent(in)    :: mesh
    type(namespace_t), intent(in)    :: namespace

    type(block_t)        :: blk
    integer              :: ip, il, nlines, ncols, idim, ierr
    FLOAT                :: j_vector(MAX_DIM), dummy(MAX_DIM), xx(MAX_DIM), rr, omega
    character(len=1024)  :: tdf_expression, phase_expression

    PUSH_SUB(external_current_init)

    !%Variable UserDefinedMaxwellExternalCurrent
    !%Type block
    !%Section MaxwellStates
    !%Description
    !%
    !% Example:
    !%
    !% <tt>%UserDefinedMaxwellExternalCurrent
    !% <br>&nbsp;&nbsp; external_current_parser      | "expression_x_dir1" | "expression_y_dir1" | "expression_z_dir1" 
    !% <br>&nbsp;&nbsp; external_current_parser      | "expression_x_dir2" | "expression_y_dir2" | "expression_z_dir2" 
    !% <br>&nbsp;&nbsp; external_current_td_function | "amplitude_j0_x"    | "amplitude_j0_y"    | "amplitude_j0_z"    | omega   | envelope_td_function_name | phase
    !% <br>%</tt>
    !%
    !% Description about UserDefinedMaxwellExternalCurrent follows
    !%
    !%Option external_current_parser 0
    !% description follows
    !%Option external_current_td_function 1
    !% description follows
    !%End

    if(parse_block(namespace, 'UserDefinedMaxwellExternalCurrent', blk) == 0) then

      ! find out how many lines (i.e. states) the block has
      nlines = parse_block_n(blk)

      st%external_current_number = nlines
      SAFE_ALLOCATE(st%external_current_modus(nlines))
      SAFE_ALLOCATE(st%external_current_string(MAX_DIM,nlines))
      SAFE_ALLOCATE(st%external_current_amplitude(1:mesh%np,MAX_DIM,nlines))
      SAFE_ALLOCATE(st%external_current_td_function(nlines))
      SAFE_ALLOCATE(st%external_current_omega(nlines))
      SAFE_ALLOCATE(st%external_current_td_phase(nlines))

      ! read all lines
      do il = 1, nlines
        ! Check that number of columns is four or seven.
        ncols = parse_block_cols(blk, il - 1)
        if((ncols  /=  4) .and. (ncols /= 6) .and. (ncols /= 7) .and. (ncols /= 5)) then
      !  if (ncols /= 3) then
          message(1) = 'Each line in the MaxwellExternalCurrent block must have'
          message(2) = 'four or seven columns.'
          call messages_fatal(2)
        end if

        call parse_block_integer(blk, il - 1, 0, st%external_current_modus(il))

        if (st%external_current_modus(il) &
            == OPTION__USERDEFINEDMAXWELLEXTERNALCURRENT__EXTERNAL_CURRENT_PARSER) then
          ! parse formula string
          do idim=1, st%d%dim
            call parse_block_string( blk, il - 1, idim, st%external_current_string(idim,il))
            call conv_to_C_string(st%external_current_string(idim,il))
          end do
        else if (st%external_current_modus(il) &
                 == OPTION__USERDEFINEDMAXWELLEXTERNALCURRENT__EXTERNAL_CURRENT_TD_FUNCTION)then
          do ip=1, mesh%np
            xx(:) = mesh%x(ip,:)
            rr    = sqrt(sum(xx(:)**2))
            do idim=1, st%d%dim
              call parse_block_string( blk, il - 1, idim, st%external_current_string(idim,il))
              call conv_to_C_string(st%external_current_string(idim,il))
              call parse_expression(j_vector(idim), dummy(idim), st%d%dim, xx, rr, M_ZERO, &
                                    st%external_current_string(idim,il))
              j_vector(idim) = units_to_atomic(units_inp%energy/(units_inp%length**2),j_vector(idim))
            end do
            st%external_current_amplitude(ip,:,il) = j_vector(:)
          end do
          call parse_block_float(blk, il-1, 4, omega, unit_one/units_inp%time)
          st%external_current_omega(il) = omega
          call parse_block_string(blk, il-1, 5, tdf_expression)
          call tdf_read(st%external_current_td_function(il), trim(tdf_expression), ierr)
          if(parse_block_cols(blk, il-1) > 6) then
            call parse_block_string(blk, il-1, 6, phase_expression)
            call tdf_read(st%external_current_td_phase(il), trim(phase_expression), ierr)
            if (ierr /= 0) then            
              write(message(1),'(3A)') 'Error in the "', trim(tdf_expression), '" field defined in the TDExternalFields block:'
              write(message(2),'(3A)') 'Time-dependent phase function "', trim(phase_expression), '" not found.'
              call messages_warning(2)
            end if
          else
            call tdf_init(st%external_current_td_phase(il))
          end if
        end if
      end do
    end if
    POP_SUB(external_current_init)
  end subroutine external_current_init

  !----------------------------------------------------------
  subroutine external_current_calculation(st, mesh, time, current)
    type(states_mxll_t), intent(inout) :: st
    type(mesh_t),   intent(in)    :: mesh
    FLOAT,          intent(in)    :: time
    FLOAT,          intent(inout) :: current(:,:)

    integer :: ip, jn, idim, il
    FLOAT   :: xx(MAX_DIM), rr, tt, j_vector(MAX_DIM), dummy(MAX_DIM), omega, shift, width, amp(MAX_DIM)
    CMPLX   :: exp_arg

    PUSH_SUB(external_current_calculation)

    current(:,:) = M_ZERO
    do jn=1, st%external_current_number
       if (st%external_current_modus(jn) &
            &== OPTION__USERDEFINEDMAXWELLEXTERNALCURRENT__EXTERNAL_CURRENT_PARSER) then
        do ip=1, mesh%np
          do idim=1, st%d%dim
            xx(:) = mesh%x(ip,:)
            rr    = sqrt(sum(xx(:)**2))
            tt    = time
            call parse_expression(j_vector(idim), dummy(idim), st%d%dim, xx, rr, tt, &
                 & trim(st%external_current_string(idim,jn)))
            j_vector(idim) = units_to_atomic(units_inp%energy/(units_inp%length**2),j_vector(idim))
          end do
          current(ip,:) = current(ip,:) + j_vector(:)
        end do
     else if(st%external_current_modus(jn)&
          & == OPTION__USERDEFINEDMAXWELLEXTERNALCURRENT__EXTERNAL_CURRENT_TD_FUNCTION) then
        do ip=1, mesh%np
           exp_arg = st%external_current_omega(jn) * time &
                 & + tdf(st%external_current_td_phase(jn),time)
           amp(:)  = st%external_current_amplitude(ip,:,jn) &
                 & * tdf(st%external_current_td_function(jn),time) 
          j_vector(:) = real(amp(:) * exp(-M_zI *(exp_arg)))
          current(ip,:) = current(ip,:) + j_vector(:)
        end do
      end if
    end do

    POP_SUB(external_current_calculation)
  end subroutine external_current_calculation
