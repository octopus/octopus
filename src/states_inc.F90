!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

! Orthonormalizes nst orbital in mesh m
subroutine R_FUNC(states_gram_schmidt)(nst, m, dim, psi, start)
  integer, intent(in) :: nst, dim
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(inout) :: psi(0:m%np, 1:dim, 1:nst)
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

subroutine R_FUNC(states_gram_schmidtq)(nst, m, dim, psi, start)
  integer, intent(in) :: nst, dim
  type(mesh_type), intent(IN) :: m
  complex(r8), intent(inout), dimension(*) :: psi
  integer, intent(in), optional :: start

  integer :: p, q, id, stst, n
  real(r8) :: nrm2
  complex(r8) :: ss

  n = m%R_FUNC(npw)*dim

  if(present(start)) then
    stst = start
  else
    stst = 1
  end if

  do p = stst, nst
    do q = 1, p - 1
      ss = R_FUNC(states_dotpq)(m, dim, psi((q-1)*n+1), psi((p-1)*n+1))
      call R_FUNC(axpy) (n, -ss, psi((q-1)*n+1), 1, psi((p-1)*n+1), 1)
    enddo
    nrm2 = R_FUNC(states_nrm2q)(m, dim, psi((p-1)*n+1))
    ss = cmplx(1.0_r8/nrm2)
    call R_FUNC(scal) (n, ss, psi((p-1)*n+1), 1)
  end do

  return
end subroutine R_FUNC(states_gram_schmidtq)

subroutine R_FUNC(states_rs2fs)(m, dim, psirs, psifs)
  type(mesh_type), intent(in)            :: m
  integer, intent(in)                    :: dim
  R_TYPE, intent(in), dimension(*)       :: psirs
  complex(r8), intent(out), dimension(*) :: psifs
  integer :: i

  do i = 1, dim
     call R_FUNC(mesh_rs2fs)(m, psirs(m%np*(i-1)+1), psifs(m%R_FUNC(npw)*(i-1)+1))
  enddo

end subroutine R_FUNC(states_rs2fs)

subroutine R_FUNC(states_fs2rs)(m, dim, psifs, psirs)
  type(mesh_type), intent(in)            :: m
  integer, intent(in)                    :: dim
  complex(r8), intent(inout), dimension(*)  :: psifs
  R_TYPE, intent(out), dimension(*)         :: psirs
  integer :: i
  do i = 1, dim
     call R_FUNC(mesh_fs2rs)(m, psifs(m%R_FUNC(npw)*(i-1)+1:m%R_FUNC(npw)*i), &
                                psirs(m%np*(i-1)+1:m%np*i))
  enddo

end subroutine R_FUNC(states_fs2rs)

R_TYPE function R_FUNC(states_dotp)(m, dim, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  R_TYPE, intent(IN) :: f1(1:m%np, dim), f2(1:m%np, dim)
  R_TYPE :: R_FUNC(states_ddot)
  R_TYPE, external :: R_DOT

  dotp = R_DOT(m%np*dim, f1, 1, f2, 1)*m%vol_pp
end function R_FUNC(states_dotp)

complex(r8) function R_FUNC(states_dotpq)(m, dim, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  complex(r8), intent(IN), dimension(*) :: f1, f2
  complex(r8) :: R_FUNC(states_ddot)
  complex(r8), external :: zdotc

  integer :: i

  dotp = R_TOTYPE(M_ZERO)
  do i = 1, dim
     dotp = dotp + R_FUNC(mesh_dotpq)(m, f1(m%R_FUNC(npw)*(i-1)+1), f2(m%R_FUNC(npw)*(i-1)+1))
  enddo
end function R_FUNC(states_dotpq)

real(r8) function R_FUNC(states_nrm2)(m, dim, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  R_TYPE, intent(IN) :: f(1:m%np, dim)
  real(r8), external :: R_NRM2

  nrm2 = R_NRM2(m%np*dim, f, 1)*sqrt(m%vol_pp)
end function R_FUNC(states_nrm2)

real(r8) function R_FUNC(states_nrm2q)(m, dim, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  complex(r8), intent(IN), dimension(*) :: f
  real(r8), external :: dznrm2

  integer :: i

  nrm2 = M_ZERO
  do i = 1, dim
     nrm2 = nrm2 + R_FUNC(mesh_nrm2q)(m, f(m%R_FUNC(npw)*(i-1)+1))
  enddo
end function R_FUNC(states_nrm2q)

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

  write(iunit) int(m%box_shape, i4), m%h, m%rsize, m%xsize
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

  sub_name = 'states_load_restart'; call push_sub()
  ok = .true.

  if(conf%verbose > 20) then
    write(stdout, '(3a)')"Info: Reading wavefunctions from file '", &
         trim(filename), "'"
  end if

  call io_assign(iunit)
  open(iunit, status='old', file=trim(filename), form='unformatted', err=998)
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
    st%R_FUNC(psi) = R_TOTYPE(M_ZERO)

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
                st%zpsi(1:m%np, id, ist, ik) = cmplx(dpsi(:), M_ZERO, r8)
#             endif
            end if
          else
            read(iunit, err=999) zpsi(:)
            if(ist >= st%st_start .and. ist <= st%st_end) then
#             ifdef R_TREAL
                st%dpsi(1:m%np, id, ist, ik) = real(zpsi(:), r8)
#             else
                st%zpsi(1:m%np, id, ist, ik) = zpsi(:)
#             endif
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
  call pop_sub()
  return

999 call io_close(iunit)
998 continue
  message(1) = 'Error reading from file '//trim(filename)//"'"
  call write_warning(1)
  ok = .false.
  call pop_sub()
  return

end function R_FUNC(states_load_restart)

subroutine R_FUNC(states_output) (st, m, dir, outp)
  type(states_type), intent(IN) :: st
  type(mesh_type), intent(IN) :: m
  character(len=*), intent(IN) :: dir
  type(output_type), intent(IN) :: outp

  integer :: ik, ist, idim, is
  character(len=80) :: fname
  real(r8) :: u

  u = 1._r8/units_out%length%factor**conf%dim

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif
    if(outp%what(output_density)) then
      do is = 1, st%nspin
        write(fname, '(a,i1)') 'density-', is
        call doutput_function(outp, dir, fname, m, st%rho(:, is), u)
      end do
    end if
#ifdef HAVE_MPI
  end if
#endif

  if(outp%what(output_wfs)) then
    do ist = 1, st%nst
      is = outp%wfs((ist-1)/32 + 1)
      if(iand(is, 2**(modulo(ist-1, 32))).ne.0) then
        do ik = 1, st%nik
          do idim = 1, st%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'wf-', ik, '-', ist, '-', idim
            call R_FUNC(output_function) (outp, dir, fname, m, &
                 st%R_FUNC(psi) (1:, idim, ist, ik), sqrt(u))
          end do
        end do
      end if
    end do
  end if

  if(conf%dim==3.and.outp%what(output_elf)) then
    call elf()
  end if

contains
  ! WARNING some constants are probably wrong for 1 and 2D
  subroutine elf()
    real(r8) :: f, d, s
    real(r8), allocatable :: c(:), j(:,:,:), r(:), gr(:,:)
    R_TYPE, allocatable :: gpsi(:,:)
    integer :: i, is, ik

    ! single or double occupancy
    if(st%nspin == 1) then
      s = 2._r8
    else
      s = 1._r8
    end if

#if defined(R_TCOMPLEX)
    allocate(j(3, m%np, st%nspin))
    call calc_current(m, st, j)
#endif

    allocate(c(m%np))
    do_is: do is = 1, st%nspin
      ! first term
      allocate(r(0:m%np), gr(3, m%np))
      r(1:m%np) = st%rho(1:m%np, is)/s
      r(0) = 0._r8
      call dmesh_derivatives(m, r(0:m%np), grad=gr)
      do i = 1, m%np
        if(r(i) >= 1d-10) then
          c(i) = -0.25_r8*sum(gr(1:conf%dim, i)**2)/r(i)
#if defined(R_TCOMPLEX)
          c(i) = c(i) - sum(j(1:conf%dim, i, is)**2)/(s*s*r(i))
#endif
        end if
      end do
      deallocate(gr)

      ! now the second term
      allocate(gpsi(3, m%np))
      do ik = is, st%nik, st%nspin
        do ist = 1, st%nst
          do idim = 1, st%dim
            call R_FUNC(mesh_derivatives) (m, st%R_FUNC(psi)(0:m%np, idim, ist, ik), grad=gpsi(:,:))
            do i = 1, m%np
              if(r(i) >= 1d-10) then
                c(i) = c(i) + st%occ(ist, ik)/s*sum(gpsi(1:conf%dim, i)*R_CONJ(gpsi(1:conf%dim, i)))
              end if
            end do
          end do
        end do
      end do
      deallocate(gpsi)
      
      f = 3._r8/5._r8*(6._r8*M_PI**2)**(2._r8/3._r8)
      do i = 1, m%np
        if(r(i) >= 1d-10) then
          d    = f*r(i)**(5._r8/3._r8)
          c(i) = 1._r8/(1._r8 + (c(i)/d)**2)
        else
          c(i) = M_ZERO
        end if
      end do
      
      deallocate(r)
      
      write(fname, '(a,i1)') 'elf-', is
      call doutput_function(outp, dir, fname, m, c, 1._r8)
      
    end do do_is
#if defined(R_TCOMPLEX)
    deallocate(j)
#endif
    deallocate(c)

  end subroutine elf
end subroutine R_FUNC(states_output)
