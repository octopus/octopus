! Orthonormalizes nst orbital in mesh m
subroutine R_FUNC(states_gram_schmidt)(nst, m, dim, psi, start)
  integer, intent(in) :: nst, dim
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(inout) :: psi(0:m%np, dim, nst) 
  integer, intent(in), optional :: start

  integer :: p, q, id, stst
  real(r8) :: nrm2
  R_TYPE :: ss

  if(present(start)) then
    stst = start
  else
    stst = 1
  end if

  do p = stst, nst
    do q = 1, p - 1
      ss = R_FUNC(states_dotp)(m, dim, psi(1:m%np, :, q), psi(1:m%np, :, p))
      do id = 1, dim
        call R_FUNC(axpy) (m%np, -ss, psi(1:m%np, id, q), 1, psi(1:m%np, id, p), 1)
      end do
    enddo
    nrm2 = R_FUNC(states_nrm2)(m, dim, psi(1:m%np, :, p))
    ss = R_TOTYPE(1.0_r8/nrm2)
    do id = 1, dim
      call R_FUNC(scal) (m%np, ss, psi(1:m%np, id, p), 1)
    end do
  end do

  return
end subroutine R_FUNC(states_gram_schmidt)

R_TYPE function R_FUNC(states_dotp)(m, dim, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  R_TYPE, intent(IN) :: f1(1:m%np, dim), f2(1:m%np, dim)
  R_TYPE :: R_FUNC(states_ddot)

  integer :: i

  dotp = R_TOTYPE(0._r8)
  do i = 1, dim
    dotp = dotp + R_FUNC(mesh_dotp)(m, f1, f2)
  end do

end function R_FUNC(states_dotp)

real(r8) function R_FUNC(states_nrm2)(m, dim, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  R_TYPE, intent(IN) :: f(1:m%np, dim)

  integer :: i

  nrm2 = 0._r8
  do i = 1, dim
    nrm2 = nrm2 + R_FUNC(mesh_nrm2)(m, f)**2
  end do

  nrm2 = sqrt(nrm2)
end function R_FUNC(states_nrm2)

! TODO use netcdf
subroutine R_FUNC(states_write_restart)(filename, m, st, iter, v1, v2)
  character(len=*), intent(in) :: filename
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  integer, intent(in), optional :: iter ! used in TD
  real(r8), intent(in), optional :: v1(m%np, st%nspin), v2(m%np, st%nspin)

  integer :: iunit, ik, ist, id
  integer(i4) :: mode

  sub_name = 'states_write_restart'; call push_sub()

  call io_assign(iunit)
  open(iunit, status='unknown', file=trim(filename), form='unformatted')
    
#ifdef R_TREAL
  mode = 0_i4
#else
  mode = 1_i4
#endif

  write(iunit) int(m%box_shape, i4), m%h, m%rsize, m%zsize
  write(iunit) int(m%np, i4), int(st%dim, i4), int(st%st_start, i4), &
       int(st%st_end, i4), int(st%nik, i4), int(st%ispin, i4), mode

  ! psi has to be written in parts, or segmentation fault in some machines
  do ik = 1, st%nik
    do ist = st%st_start, st%st_end
      do id = 1, st%dim
        write(iunit) st%R_FUNC(psi)(1:m%np, id, ist, ik)
      end do
    end do
  end do

  ! eigenvalues are also needed ;)
  do ik = 1, st%nik
    do ist = st%st_start, st%st_end
      write(iunit) st%eigenval(ist, ik)
    end do
  end do

  if(present(iter)) then
    write(iunit) int(iter, i4), v1, v2
  end if
  call io_close(iunit)
  
  call pop_sub()
end subroutine R_FUNC(states_write_restart)

logical function R_FUNC(states_load_restart)(filename, m, st, iter, v1, v2) result(ok)
  character(len=*), intent(in) :: filename
  type(mesh_type), intent(in) :: m
  type(states_type), intent(inout) :: st
  integer, intent(out), optional :: iter ! used in TD
  real(r8), intent(out), optional :: v1(m%np, st%nspin), v2(m%np, st%nspin)

  integer :: iunit, ik, ist, id
  integer(i4) :: ii, old_np, old_dim, old_start, old_end, old_nik, old_ispin, mode
  real(r8) :: e
  real(r8), allocatable :: dpsi(:)
  complex(r8), allocatable :: zpsi(:)

  sub_name = 'systm_load_psi'; call push_sub()
  ok = .true.

  if(conf%verbose > 20) then
    write(stdout, '(3a)')"Info: Reading wavefunctions from file '", &
         trim(filename), "'"
  end if

  call io_assign(iunit)
  open(iunit, status='old', file=trim(filename), form='unformatted', err=999)

  read(iunit, err=999) ! mesh stuff is now skipped
  read(iunit, err=100) old_np, old_dim, old_start, old_end, old_nik, old_ispin, mode
  go to 200

  ! ugly thing to handle old restarts
100 continue
  backspace(iunit)
  read(iunit, err=999) old_np, old_dim, old_start, old_end, old_nik, old_ispin
  mode = -1 ! old restart file

200 continue
  if(old_np.ne.m%np .or. old_dim.ne.st%dim .or. & ! different mesh, cannot read
       old_start.gt.st%st_start .or. old_end.lt.st%st_end .or. old_nik.ne.st%nik) then
    message(1) = 'Restart file has a different mesh!'
    write(message(2), '(a,i6,a,i6,a)') '  m%np        = ', m%np,        ' != ', old_np, ' or'
    write(message(3), '(a,i6,a,i6,a)') '  st%dim      = ', st%dim,      ' != ', old_dim, ' or'
    write(message(4), '(a,i6,a,i6,a)') '  st%st_start = ', st%st_start, ' <  ', old_start, ' or'
    write(message(5), '(a,i6,a,i6,a)') '  st%st_end   = ', st%st_end,   ' >  ', old_end, ' or'
    write(message(6), '(a,i6,a,i6)')   '  st%nik      = ', st%nik,      ' != ', old_nik

    call write_warning(6)
    go to 999 ! one go to does not harm :)
  else
    st%R_FUNC(psi) = R_TOTYPE(0._r8)

    if(mode == 0 .or. mode == -1) then
      allocate(dpsi(1:m%np))
    else
      allocate(zpsi(1:m%np))
    end if

    do ik = 1, st%nik
      do ist = old_start, old_end
        do id = 1, st%dim
          imode: if(mode == 0 .or. mode == -1) then
            read(iunit, err=999) dpsi(:)
            if(ist >= st%st_start .and. ist <= st%st_end) then
#             ifdef R_TREAL
                st%dpsi(1:m%np, id, ist, ik) = dpsi(:)
#             else
                st%zpsi(1:m%np, id, ist, ik) = cmplx(dpsi(:), 0._r8, r8)
#             endif
            else
              st%R_FUNC(psi) (1:m%np, id, ist, ik) = R_TOTYPE(1._r8) ! stupid
            end if
          else
            read(iunit, err=999) zpsi(:)
            if(ist >= st%st_start .and. ist <= st%st_end) then
#             ifdef R_TREAL
                st%dpsi(1:m%np, id, ist, ik) = real(zpsi(:), r8)
#             else
                st%zpsi(1:m%np, id, ist, ik) = zpsi(:)
#             endif
            else
              st%R_FUNC(psi) (1:m%np, id, ist, ik) = R_TOTYPE(1._r8) ! stupid
            end if
          end if imode

        end do
      end do
    end do

    if(mode == 0 .or. mode == -1) then
      deallocate(dpsi)
    else
      deallocate(zpsi)
    end if

    if(mode == -1) then
      read(iunit, err=999) st%eigenval(st%st_start:st%st_end, 1:st%nik)
    else
      do ik = 1, st%nik
        do ist = old_start, old_end
          read(iunit, err=999) e
          if(ist >= st%st_start .and. ist <= st%st_end) then
            st%eigenval(ist, ik) = e
          end if
        end do
      end do
    end if

    if(present(iter)) then ! read the time-dependent stuff
      read(iunit, err=999) ii, v1, v2
      iter = int(ii)
    end if
  end if
  
  call io_close(iunit)
  return

999 continue
  message(1) = 'Error reading from file '//trim(filename)//"'"
  call write_warning(1)
  ok = .false.
  return

end function R_FUNC(states_load_restart)
