subroutine td_init(td, sys, m, st)
  type(td_type), intent(out) :: td
  type(system_type), intent(IN) :: sys
  type(mesh_type), intent(inout) :: m
  type(states_type), intent(inout) :: st

  integer :: iunit

  sub_name = 'td_init'; call push_sub()

  td%iter = 0

  call oct_parse_int(C_string("TDMaximumIter"), 1500, td%max_iter)
  if(td%max_iter < 1) then
    write(message(1), '(a,i6,a)') "Input: '", td%max_iter, "' is not a valid TDMaximumIter"
    message(2) = '(1 <= TDMaximumIter)'
    call write_fatal(2)
  end if
  
  call oct_parse_int(C_string("TDSaveIter"), 100, td%save_iter)
  if (td%save_iter < 0 .or. td%save_iter>td%max_iter) then
    write(message(1), '(a,i6,a)') "Input: '", td%save_iter, "' is not a valid TDSaveIter"
    message(2) = '(0 <= TDSaveIter <= TDMaximumIter)'
    call write_fatal(2)
  end if

  call oct_parse_double(C_string("TDTimeStep"), 0.07_r8/units_inp%time%factor, td%dt)
  td%dt = td%dt * units_inp%time%factor
  if (td%dt <= 0._r8) then
    write(message(1),'(a,f14.6,a)') "Input: '", td%dt, "' is not a valid TDTimeStep"
    message(2) = '(0 < TDTimeStep)'
    call write_fatal(2)
  end if
  
  call oct_parse_int(C_string("TDEvolutionMethod"), 2, td%evolution_method)
  if (td%evolution_method<1 .or. td%evolution_method>3) then
    write(message(1), '(a,i6,a)') "Input: '", td%evolution_method, "' is not a valid TDEvolutionMethod"
    message(2) = '(1 <= TDEvolutionMethod <= 3)'
    call write_fatal(2)
  end if
  if(td%evolution_method == 3) then ! split operator method
    call td_init_evolution_splitop()
  end if

  call oct_parse_int(C_string("TDDipoleLmax"), 1, td%lmax)
  if (td%lmax < 0 .or. td%lmax > 4) then
    write(message(1), '(a,i6,a)') "Input: '", td%lmax, "' is not a valid TDDipoleLmax"
    message(2) = '(0 <= TDDipoleLmax <= 4 )'
    call write_fatal(2)
  end if

  ! delta impulse used to calculate optical spectrum
  ! units are 1/length
  call oct_parse_double(C_string("TDDeltaStrength"), 0._r8, td%delta_strength)
  td%delta_strength = td%delta_strength / units_inp%length%factor
  if(oct_parse_isdef(C_string('TDPolarization')) .ne. 0) then
    call oct_parse_block_double(C_string('TDPolarization'), 0, 0, td%pol(1))
    call oct_parse_block_double(C_string('TDPolarization'), 0, 1, td%pol(2))
    call oct_parse_block_double(C_string('TDPolarization'), 0, 2, td%pol(3))
  else  !default along the z-direction
    td%pol(1:2) = 0._r8
    td%pol(3)   = 1._r8
  endif

  call oct_parse_int(C_string("TDGauge"), 1, td%gauge)
  if (td%gauge < 1 .or. td%gauge > 2) then
    write(message(1), '(a,i6,a)') "Input: '", td%gauge, "' is not a valid TDGauge"
    message(2) = 'Accepted values are:'
    message(3) = '   1 = length gauge'
    message(4) = '   2 = velocity gauge'
    call write_fatal(4)
  end if  

  ! now the lasers
  call td_init_lasers()

  ! now come the absorbing boundaries
  call td_init_ab()

  ! now the photoelectron stuff
  call PES_init(td%PESv, m, sys%st, td%ab, td%save_iter)

  ! occupational analysis stuff
  call oct_parse_logical(C_string("TDOccupationalAnalysis"), .false., td%occ_analysis)

  ! should we move the ions during the simulation?
  call oct_parse_int(C_string("MoveIons"), 0, td%move_ions)
  if(td%move_ions.ne.0 .and. td%move_ions<3 .and. td%move_ions>4) then
    write(message(1),'(a,i4)') "Input: '", td%move_ions, &
         "' is not a valid MoveIons"
    message(2) = '  MoveIons = 0 <= do not move'
    message(3) = '  MoveIons = 3 <= verlet'
    message(4) = '  MoveIons = 4 <= velocity verlet'
    call write_fatal(4)
  endif
  
  call td_init_states()

  ! get name of continuation file
  write(td%filename, '(a,a,i3.3,a)') trim(sys%sysname), &
       '.', mpiv%node, '.cont'

contains
  subroutine td_init_evolution_splitop()
    integer :: ix, iy, iz, ixx, iyy, izz
    real(r8) :: temp(3), vec

    allocate(td%kin_2(m%fft_n(1), m%fft_n(2), m%fft_n(3)))
    temp(1:3) = 2.0_r8*M_PI/(m%fft_n(1:3)*m%h(1:3))    
    do ix = 1, m%fft_n(1)
      ixx = pad_feq(ix, m%fft_n(1), .true.)
      do iy = 1, m%fft_n(2)
        iyy = pad_feq(iy, m%fft_n(2), .true.)
        do iz = 1, m%fft_n(3)
          izz = pad_feq(iz, m%fft_n(3), .true.)
          vec = real(temp(1)**2*ixx**2 + temp(2)**2*iyy**2 + temp(3)**2*izz**2, r8)/2._r8
          
          td%kin_2(ix, iy, iz) = exp(- M_zI*td%dt/2._r8* vec)
        end do
      end do
    end do
  end subroutine td_init_evolution_splitop

  subroutine td_init_lasers()
    call laser_init(m, td%no_lasers, td%lasers)
    td%output_laser = .false.
    if(td%no_lasers>0 ) then
      message(1) = 'Info: Lasers'
      call write_info(1)
      if(conf%verbose > 20) call laser_write_info(td%no_lasers, td%lasers, stdout)
      
      td%delta_strength = 0._r8 ! no delta impulse if lasers exist
      call oct_parse_logical(C_string("TDOutputLaser"), .false., td%output_laser)
    end if

  end subroutine td_init_lasers

  subroutine td_init_ab()
    integer :: i, dummy
    real(r8) :: d, r, x(3)

    call oct_parse_int(C_string("TDAbsorbingBoundaries"), 0, dummy)
    if(dummy .eq. 1 .or. dummy .eq. 2) then
      td%ab = dummy
      call oct_parse_double(C_string("TDABWidth"), 4._r8/units_inp%length%factor, td%ab_width)
      td%ab_width  = td%ab_width * units_inp%length%factor
      if(td%ab == 1) then
        call oct_parse_double(C_string("TDABHeight"), -0.2_r8/units_inp%energy%factor, td%ab_height)
        td%ab_height = td%ab_height * units_inp%energy%factor
      else
        td%ab_height = 1._r8
      end if

      ! generate boundary potential...
      allocate(td%ab_pot(m%np))
      td%ab_pot = 0._r8
      pot: do i = 1, m%np
        call mesh_r(m, i, r, x=x)
        
        select case(m%box_shape)
        case(1)
          d = r - (m%rsize - td%ab_width)
          if(d.gt.0._r8) then
            td%ab_pot(i) = td%ab_height * sin(d*M_PI/(2._r8*td%ab_width))**2
            !td%ab_pot(i) = 1._r8
          end if
        
        case(2)
          d = sqrt(x(1)**2 + x(2)**2) - (m%rsize - td%ab_width)
          if(d.gt.0._r8)  &
               td%ab_pot(i) = td%ab_height * sin(d*M_PI/(2._r8*td%ab_width))**2
          d = abs(x(3)) - (m%zsize - td%ab_width)
          if(d.gt.0._r8)  &
               td%ab_pot(i) = td%ab_pot(i) + td%ab_height * sin(d*M_PI/(2._r8*td%ab_width))**2
          if(abs(td%ab_pot(i)) > abs(td%ab_height)) td%ab_pot(i) = td%ab_height
        
        case default
          message(1) = "Absorbing boundaries are not implemented for"
          message(2) = "Box_shape = 3"
          call write_warning(2)
          exit pot
        end select
      end do pot
      else
        ! some compilers have problems passing unallocated
        allocate(td%ab_pot(1))
    end if
  end subroutine td_init_ab
  
  subroutine td_init_states()
    integer :: i, ix

    ! MPI stuff
#if defined(HAVE_MPI) && defined(MPI_TD)
    if(st%nst < mpiv%numprocs) then
      message(1) = "Have more processors than necessary"
      write(message(2),'(i4,a,i4,a)') mpiv%numprocs, " processors and ", st%nst, " states."
      call write_fatal(2)
    end if

    i = st%nst / mpiv%numprocs
    ix = st%nst - i*mpiv%numprocs
    if(ix > 0 .and. mpiv%node < ix) then
      i = i + 1
      st%st_start = mpiv%node*i + 1
      st%st_end = st%st_start + i - 1
    else
      st%st_end = st%nst - (mpiv%numprocs - mpiv%node - 1)*i
      st%st_start = st%st_end - i + 1
    end if
    write(stdout, '(a,i4,a,i4,a,i4)') "Info: Node ", mpiv%node, " will propagate state ", &
         st%st_start, " - ", st%st_end

    ! syncronize to get the output properly
    call MPI_Barrier(MPI_COMM_WORLD, i)
  
#else
    st%st_start = 1
    st%st_end   = st%nst
#endif

    ! allocate memory
    allocate(td%v_old1(m%np, st%nspin), td%v_old2(m%np, st%nspin))
    allocate(st%zpsi(0:m%np, st%dim, st%st_start:st%st_end, st%nik))
  end subroutine td_init_states

end subroutine td_init

subroutine td_end(td)
  type(td_type), intent(inout) :: td

  sub_name = 'td_end'; call push_sub()

#ifndef NO_PES
  call PES_end(td%PESv)
!  if(td%calc_PES_mask) call PES_mask_end(td%PES_mask)
#endif

  if(td%evolution_method == 3) then ! split operator method
    deallocate(td%kin_2); nullify(td%kin_2)
  end if

  if(associated(td%ab_pot)) then
    deallocate(td%ab_pot); nullify(td%ab_pot)
  end if

  if(associated(td%v_old1)) then
    deallocate(td%v_old1); nullify(td%v_old1)
    deallocate(td%v_old2); nullify(td%v_old2)
  end if

  call laser_end(td%no_lasers, td%lasers)

  call pop_sub()
end subroutine td_end
