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
    ss = R_TOTYPE(1.0_r8/sqrt(nrm2))
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

  call io_assign(iunit)
  open(iunit, status='unknown', file=trim(filename), form='unformatted')
    
  write(iunit) m%box_shape, m%h, m%rsize, m%zsize
  write(iunit) m%np, st%dim, 1, st%nst, st%nik, st%ispin
  write(iunit) st%R_FUNC(psi)(1:m%np, 1:st%dim, st%st_start:st%st_end, 1:st%nik)
  ! eigenvalues are also needed ;)
  write(iunit) st%eigenval(st%st_start:st%st_end, 1:st%nik)

  if(present(iter)) then
    write(iunit) iter, v1, v2
  end if
  call io_close(iunit)

end subroutine R_FUNC(states_write_restart)

logical function R_FUNC(states_load_restart)(filename, m, st, iter, v1, v2) result(ok)
  character(len=*), intent(in) :: filename
  type(mesh_type), intent(in) :: m
  type(states_type), intent(inout) :: st
  integer, intent(out), optional :: iter ! used in TD
  real(r8), intent(out), optional :: v1(m%np, st%nspin), v2(m%np, st%nspin)

  integer :: iunit, ik, ist, id, old_np, old_dim, old_start, old_end, old_nik

  sub_name = 'systm_load_psi'; call push_sub()
  ok = .true.

  if(conf%verbose > 20) then
    write(stdout, '(3a)')"Info: Reading wavefunctions from file '", &
         trim(filename), "'"
  end if

  call io_assign(iunit)
  open(iunit, status='old', file=trim(filename), form='unformatted', err=999)

  read(iunit, err=999) ! mesh stuff is now skipped
  read(iunit, err=999) old_np, old_dim, old_start, old_end, old_nik
  
  if(old_np.ne.m%np .or. old_dim.ne.st%dim .or. & ! different mesh, cannot read
       old_start.ne.st%st_start .or. old_end.ne.st%st_end .or. old_nik.ne.st%nik) then
    message(1) = 'Restart file has a different mesh!'
    write(message(2), '(a,i6,a,i6,a)') '  m%np        = ', m%np,        ' != ', old_np, ' or'
    write(message(3), '(a,i6,a,i6,a)') '  st%dim      = ', st%dim,      ' != ', old_dim, ' or'
    write(message(4), '(a,i6,a,i6,a)') '  st%st_start = ', st%st_start, ' != ', old_start, ' or'
    write(message(5), '(a,i6,a,i6,a)') '  st%st_end   = ', st%st_end,   ' != ', old_end, ' or'
    write(message(6), '(a,i6,a,i6)')   '  st%nik      = ', st%nik,      ' != ', old_nik

    call write_warning(6)
    go to 999 ! one go to does not harm :)
  else
    st%R_FUNC(psi) = R_TOTYPE(0._r8)
    read(iunit, err=999) st%R_FUNC(psi)(1:m%np, 1:st%dim, st%st_start:st%st_end, 1:st%nik)
    read(iunit, err=999) st%eigenval(st%st_start:st%st_end, 1:st%nik)

    if(present(iter)) then ! read the time-dependent stuff
      read(iunit, err=999) iter, v1, v2
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
