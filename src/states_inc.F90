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

! Calculates the new density out the wavefunctions and occupations...
subroutine R_FUNC(calcdens)(st, np, rho, reduce)
  type(states_type), intent(in) :: st
  integer, intent(in) :: np
  real(r8), intent(out) :: rho(np, st%nspin)
  logical, intent(in), optional :: reduce

  integer :: i, ik, p, sp
#ifdef HAVE_MPI
  real(r8), allocatable :: reduce_rho(:,:)
  R_TYPE, allocatable :: reduce_rho_off(:)
  integer :: ierr
#endif

  call push_sub('calc_dens')

  if(st%ispin == SPIN_POLARIZED) then
    sp = 2
  else
    sp = 1
  end if

  rho = M_ZERO
  do ik = 1, st%nik, sp
    do p  = st%st_start, st%st_end
      do i = 1, np
           rho(i, 1) = rho(i, 1) + st%kweights(ik)  *st%occ(p, ik)  *R_ABS(st%R_FUNC(psi)(i, 1, p, ik))**2
         select case(st%ispin)
         case(SPIN_POLARIZED)
           rho(i, 2) = rho(i, 2) + st%kweights(ik+1)*st%occ(p, ik+1)*R_ABS(st%R_FUNC(psi)(i, 1, p, ik+1))**2
         case(SPINORS)
           rho(i, 2) = rho(i, 2) + st%kweights(ik)  *st%occ(p, ik)  *R_ABS(st%R_FUNC(psi)(i, 2, p, ik))**2
           rho(i, 3) = rho(i, 3) + st%kweights(ik)*st%occ(p, ik)  * &
                       R_REAL (st%R_FUNC(psi)(i, 1, p, ik) * R_CONJ(st%R_FUNC(psi)(i, 2, p, ik)))
           rho(i, 4) = rho(i, 4) + st%kweights(ik)*st%occ(p, ik)  * &
                       R_AIMAG(st%R_FUNC(psi)(i, 1, p, ik) * R_CONJ(st%R_FUNC(psi)(i, 2, p, ik)))

         end select
      end do
    end do
  end do

#if defined(HAVE_MPI) && defined(MPI_TD)
  ! reduce density (assumes memory is contiguous)
  if(present(reduce)) then
  if(reduce) then
    allocate(reduce_rho(1:np, st%nspin))
    call MPI_ALLREDUCE(rho(1, 1), reduce_rho(1, 1), np*st%nspin, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    rho = reduce_rho
    deallocate(reduce_rho)
  end if
  end if
#endif

  call pop_sub()
  return
end subroutine R_FUNC(calcdens)

!!! This subroutine generates a gaussian wave-function in a random
!!! position in space
subroutine X(states_random)(m, f)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(out) :: f(1:m%np)

  integer, save :: iseed = 123
  integer :: i
  real(r8) :: a(3), rnd, r

  call push_sub('mesh_random')

  call quickrnd(iseed, rnd)
  a(1) = 2.0_r8*(2*rnd - 1)
  call quickrnd(iseed, rnd)
  a(2) = 2.0_r8*(2*rnd - 1)
  call quickrnd(iseed, rnd)
  a(3) = 2.0_r8*(2*rnd - 1)

  do i = 1, m%np
     call mesh_r(m, i, r, a=a)
     f(i) = exp(-M_HALF*r*r)
  end do

  r = X(mf_nrm2)(m, f)
  call X(scal)(m%np, R_TOTYPE(M_ONE/r), f, 1) 

  call pop_sub()
end subroutine X(states_random)


! Orthonormalizes nst orbital in mesh m
subroutine R_FUNC(states_gram_schmidt)(nst, m, dim, psi, start)
  integer, intent(in) :: nst, dim
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(inout) :: psi(:, :, :)
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
  R_TYPE, intent(IN), dimension(*) :: f1, f2
  R_TYPE :: R_FUNC(states_ddot)
  R_TYPE, external :: R_DOT

  dotp = R_DOT(m%np*dim, f1, 1, f2, 1)*m%vol_pp
end function R_FUNC(states_dotp)

real(r8) function R_FUNC(states_nrm2)(m, dim, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  R_TYPE, intent(IN), dimension(*) :: f
  real(r8), external :: R_NRM2

  nrm2 = R_NRM2(m%np*dim, f, 1)*sqrt(m%vol_pp)
end function R_FUNC(states_nrm2)

real(r8) function R_FUNC(states_residue)(m, dim, hf, e, f, res) result(r)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  R_TYPE, intent(IN), dimension(:, :) :: hf, f
  real(r8), intent(in)  :: e
  R_TYPE, intent(out), optional :: res(:, :)

  if(present(res)) then
    res = hf - e*f
    r = R_FUNC(states_nrm2)(m, dim, res)
  else
    r = R_FUNC(states_nrm2)(m, dim, hf - e*f)
  endif
end function R_FUNC(states_residue)

! TODO use netcdf
subroutine R_FUNC(states_write_restart)(filename, m, st, iter, v1, v2)
  character(len=*), intent(in) :: filename
  type(mesh_type), intent(in) :: m
  type(states_type), intent(in) :: st
  integer, intent(in), optional :: iter ! used in TD
  real(r8), intent(in), optional :: v1(m%np, st%nspin), v2(m%np, st%nspin)

  integer :: iunit, ik, ist, id
  integer(i4) :: mode

  call push_sub('states_write_restart')

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

  call push_sub('states_load_restart')
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
  message(1) = 'Error reading from file '+trim(filename)+"'"
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

  if(outp%what(output_wfs_sqmod)) then
    do ist = 1, st%nst
      is = outp%wfs((ist-1)/32 + 1)
      if(iand(is, 2**(modulo(ist-1, 32))).ne.0) then
        do ik = 1, st%nik
          do idim = 1, st%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'sqm-wf-', ik, '-', ist, '-', idim
            call doutput_function (outp, dir, fname, m, &
                 abs(st%R_FUNC(psi) (1:, idim, ist, ik))**2, sqrt(u))
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

    real(r8), parameter :: dmin = 1e-10_r8
    
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
      allocate(r(m%np), gr(3, m%np))
      r(1:m%np) = st%rho(1:m%np, is)/s
      call df_gradient(m, r, gr)
      do i = 1, m%np
        if(r(i) >= dmin) then
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
            call R_FUNC(f_gradient) (m, st%R_FUNC(psi)(:, idim, ist, ik), gpsi)
            do i = 1, m%np
              if(r(i) >= dmin) then
                c(i) = c(i) + st%occ(ist, ik)/s*sum(gpsi(1:conf%dim, i)*R_CONJ(gpsi(1:conf%dim, i)))
              end if
            end do
          end do
        end do
      end do
      deallocate(gpsi)
      
      f = 3._r8/5._r8*(6._r8*M_PI**2)**(2._r8/3._r8)
      do i = 1, m%np
        if(r(i) >= dmin) then
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the dot product of two many-body states st1 and st2, for a given
! irreducible subspace (k-point) ik.
! Warning: it does not do any check on the "conformance" of the two states.
!  (maybe a subroutine should be written for that purpose...)
! Warning: it disregards completely the occupation number problem.
!  This is not important, unless the occupations are different in the
!  two states (this is not the case for the moment), or if any of the
!  occupation is null (this can be problem, and will have to be cared about.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
R_TYPE function R_FUNC(states_mpdotp)(m, ik, st1, st2) result(dotp)
  implicit none
  type(mesh_type), intent(in) :: m
  integer, intent(in) :: ik
  type(states_type), intent(in) :: st1, st2

  R_TYPE, allocatable :: a(:, :)
  R_TYPE :: det

  call push_sub('states_mpdotp')

  allocate(a(st1%nst, st1%nst))

  call R_FUNC(calculate_matrix)(m, ik, st1, st2, st1%nst, a)
  dotp = R_FUNC(det)(a, st1%nst)

  select case(st1%ispin)
   case(UNPOLARIZED)
     dotp = dotp**2
   case(SPIN_POLARIZED)
     ! We assume that ik is an odd number (maybe this should be checked. So one
     ! spin component is ik and another spin component is ik + 1.
     call R_FUNC(calculate_matrix)(m, ik+1, st1, st2, st1%nst, a)
     dotp = dotp*R_FUNC(det)(a, st1%nst)
  end select

  deallocate(a)
  call pop_sub()
  contains

    subroutine R_FUNC(calculate_matrix)(m, ik, st1, st2, n, a)
      implicit none
      type(mesh_type), intent(in) :: m
      integer, intent(in) :: ik
      type(states_type), intent(in) :: st1, st2
      integer, intent(in) :: n
      R_TYPE :: a(n, n)

      R_TYPE, allocatable :: phi2(:, :)
      integer :: i, j, k, l, r, sn, sn1, ierr, dim
#if defined(HAVE_MPI) && defined(MPI_TD)
      integer :: status(MPI_STATUS_SIZE)
      integer :: request, node(st1%nst)
#endif

      dim = st1%dim
#if defined(HAVE_MPI) && defined(MPI_TD)

      ! This code just builds an integer array, node, that assigns to each
      ! state the process that holds it. If such a thing is ever needed for any other
      ! purpose, maybe it should go to other place, more general (such as "st%node")
      sn  = st1%nst/mpiv%numprocs
      sn1 = sn + 1
      r  = mod(st1%nst, mpiv%numprocs)
      do j = 1, r
         node((j-1)*sn1+1:j*sn1) = j - 1
      enddo
      k = sn1*r
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      do j = 1, mpiv%numprocs - r
         node(k+(j-1)*sn+1:k+j*sn) = r + j - 1
      enddo

      ! Each process sends the states in st2 to the rest of the processes.
      do i = st1%st_start, st1%st_end
         do j = 0, mpiv%numprocs-1
            if(mpiv%node.ne.j) then
               call MPI_ISEND(st2%R_FUNC(psi)(1, 1, i, ik), st1%dim*m%np, MPI_DOUBLE_COMPLEX, &
                              j, i, MPI_COMM_WORLD, request, ierr)
            endif
         enddo
      enddo

      ! Processes are received, and then the matrix elements are calculated.
      allocate(phi2(m%np, st1%dim))
      do j = 1, n
         l = node(j)
         if(l.ne.mpiv%node) then
            call MPI_RECV(phi2(1, 1), st1%dim*m%np, MPI_DOUBLE_COMPLEX, &
                          l, j, MPI_COMM_WORLD, status, ierr)
         else
            phi2(:, :) = st2%R_FUNC(psi)(:, :, j, ik)
         endif
         do i = st1%st_start, st1%st_end
            a(i, j) = R_FUNC(states_dotp)(m, dim, st1%R_FUNC(psi)(1:, :, i, ik), &
                                                phi2(1:, :))
         enddo
      enddo
      deallocate(phi2)

      ! Each process holds some lines of the matrix. So it is broadcasted (All processes
      ! should get the whole matrix)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      do i = 1, n
         k = node(i)
         do j = 1, n
            call MPI_BCAST(a(i, j), 1, MPI_DOUBLE_COMPLEX, k, MPI_COMM_WORLD, ierr)
         enddo
      enddo
#else
      do i = 1, n
         do j = 1, n
            a(i, j) = R_FUNC(states_dotp)(m, dim, st1%R_FUNC(psi)(1:, :, i, ik), &
                                                  st2%R_FUNC(psi)(1:, :, j, ik))
         enddo
      enddo
#endif
    end subroutine R_FUNC(calculate_matrix)

end function R_FUNC(states_mpdotp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine R_FUNC(states_calculate_angular)(m, st, angular)
  type(mesh_type), intent(IN) :: m
  type(states_type), intent(IN) :: st
  real(r8), intent(out) :: angular(3)

  R_TYPE, allocatable :: lpsi(:, :)
  integer :: idim, ik, j, i

  angular = M_ZERO
  allocate(lpsi(conf%dim, m%np))
  do ik = 1, st%nik
     do j = 1, st%nst
        do idim = 1, st%dim
           call R_FUNC(mesh_angular_momentum)(m, st%R_FUNC(psi)(:, idim, j, ik), lpsi)
           angular(1) = angular(1) + st%occ(j, ik)*R_FUNC(mf_dotp)(m, st%R_FUNC(psi)(:, idim, j, ik), lpsi(1, :))
           angular(2) = angular(2) + st%occ(j, ik)*R_FUNC(mf_dotp)(m, st%R_FUNC(psi)(:, idim, j, ik), lpsi(2, :))
           angular(3) = angular(3) + st%occ(j, ik)*R_FUNC(mf_dotp)(m, st%R_FUNC(psi)(:, idim, j, ik), lpsi(3, :))
        enddo
     enddo
  enddo
  deallocate(lpsi)
end subroutine R_FUNC(states_calculate_angular)

subroutine R_FUNC(mesh_angular_momentum)(m, f, lf)
  type(mesh_type), intent(in) :: m
  R_TYPE, intent(in)  :: f(m%np)
  R_TYPE, intent(out) :: lf(3, m%np)

  R_TYPE, allocatable :: gf(:, :)
  real(r8) :: x(3)
  integer :: i

  allocate(gf(3, m%np))
  call R_FUNC(f_gradient)(m, f, grad = gf)

  do i = 1, m%np
     call mesh_xyz(m, i, x)
     lf(1, i) = (x(2)*gf(3, i)-x(3)*gf(2, i))
     lf(2, i) = (x(3)*gf(1, i)-x(1)*gf(3, i))
     lf(3, i) = (x(1)*gf(2, i)-x(2)*gf(1 ,i))
  enddo
#if defined(R_TCOMPLEX)
  call zscal(3*m%np, -M_zI, lf, 1)
#endif
  deallocate(gf)
end subroutine R_FUNC(mesh_angular_momentum)
