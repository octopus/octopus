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
subroutine X(calcdens)(st, np, rho, reduce)
  type(states_type), intent(IN)  :: st
  integer,           intent(in)  :: np
  FLOAT,             intent(out) :: rho(np, st%d%nspin)
  logical,           intent(in), optional :: reduce

  integer :: i, ik, p, sp
#ifdef HAVE_MPI
  FLOAT,  allocatable :: reduce_rho(:,:)
  R_TYPE, allocatable :: reduce_rho_off(:)
  integer :: ierr
#endif

  call push_sub('calc_dens')

  if(st%d%ispin == SPIN_POLARIZED) then
    sp = 2
  else
    sp = 1
  end if

  rho = M_ZERO
  do ik = 1, st%d%nik, sp
    do p  = st%st_start, st%st_end
      do i = 1, np
        rho(i, 1) = rho(i, 1) + st%d%kweights(ik)  *st%occ(p, ik)  *R_ABS(st%X(psi)(i, 1, p, ik))**2
        select case(st%d%ispin)
        case(SPIN_POLARIZED)
          rho(i, 2) = rho(i, 2) + st%d%kweights(ik+1)*st%occ(p, ik+1)*R_ABS(st%X(psi)(i, 1, p, ik+1))**2
        case(SPINORS)
          rho(i, 2) = rho(i, 2) + st%d%kweights(ik)  *st%occ(p, ik)  *R_ABS(st%X(psi)(i, 2, p, ik))**2
          rho(i, 3) = rho(i, 3) + st%d%kweights(ik)*st%occ(p, ik)  * &
                      R_REAL (st%X(psi)(i, 1, p, ik) * R_CONJ(st%X(psi)(i, 2, p, ik)))
          rho(i, 4) = rho(i, 4) + st%d%kweights(ik)*st%occ(p, ik)  * &
                      R_AIMAG(st%X(psi)(i, 1, p, ik) * R_CONJ(st%X(psi)(i, 2, p, ik)))
        end select
      end do
    end do
  end do

#if defined(HAVE_MPI) && defined(MPI_TD)
  ! reduce density (assumes memory is contiguous)
  if(present(reduce)) then
    if(reduce) then
      allocate(reduce_rho(1:np, st%d%nspin))
      call MPI_ALLREDUCE(rho(1, 1), reduce_rho(1, 1), np*st%d%nspin, &
           MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
      rho = reduce_rho
      deallocate(reduce_rho)
    end if
  end if
#endif

  call pop_sub()
  return
end subroutine X(calcdens)

!!! This subroutine generates a gaussian wave-function in a random
!!! position in space
subroutine X(states_random)(m, f)
  type(mesh_type), intent(IN)  :: m
  R_TYPE,          intent(out) :: f(1:m%np)

  integer, save :: iseed = 123
  integer :: i
  FLOAT :: a(3), rnd, r

  call push_sub('states_random')

  call quickrnd(iseed, rnd)
  a(1) = M_TWO*(2*rnd - 1)
  call quickrnd(iseed, rnd)
  a(2) = M_TWO*(2*rnd - 1)
  call quickrnd(iseed, rnd)
  a(3) = M_TWO*(2*rnd - 1)

  do i = 1, m%np
     call mesh_r(m, i, r, a=a)
     f(i) = exp(-M_HALF*r*r)
  end do

  r = X(mf_nrm2)(m, f)
  call X(lalg_scal)(m%np, R_TOTYPE(M_ONE/r), f(1)) 

  call pop_sub()
end subroutine X(states_random)


! Orthonormalizes nst orbital in mesh m
subroutine X(states_gram_schmidt)(nst, m, dim, psi, start)
  integer,         intent(in) :: nst, dim
  type(mesh_type), intent(IN) :: m
  R_TYPE,          intent(inout) :: psi(m%np, dim, nst)
  integer,         intent(in), optional :: start

  integer :: p, q, id, stst
  FLOAT :: nrm2
  R_TYPE :: ss

  if(present(start)) then
    stst = start
  else
    stst = 1
  end if

  do p = stst, nst
    do q = 1, p - 1
      ss = X(states_dotp)(m, dim, psi(1:m%np, :, q), psi(1:m%np, :, p))
      do id = 1, dim
        call X(lalg_axpy)(m%np, -ss, psi(1, id, q), psi(1, id, p))
      end do
    enddo
    nrm2 = X(states_nrm2)(m, dim, psi(1:m%np, :, p))
    ss = R_TOTYPE(M_ONE/nrm2)
    do id = 1, dim
      call X(lalg_scal)(m%np, ss, psi(1, id, p))
    end do
  end do

  return
end subroutine X(states_gram_schmidt)

R_TYPE function X(states_dotp)(m, dim, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(IN) :: f1(*), f2(*)

  dotp = X(lalg_dot)(m%np*dim, f1(1), f2(1))*m%vol_pp

end function X(states_dotp)

FLOAT function X(states_nrm2)(m, dim, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(IN) :: f(*)

  nrm2 = X(lalg_nrm2)(m%np*dim, f(1))*sqrt(m%vol_pp)

end function X(states_nrm2)

FLOAT function X(states_residue)(m, dim, hf, e, f) result(r)
  type(mesh_type),   intent(IN)  :: m
  integer,           intent(in)  :: dim
  R_TYPE,            intent(IN)  :: hf(:,:), f(:,:)
  FLOAT,             intent(in)  :: e

  R_TYPE, allocatable :: res(:,:)

  allocate(res(m%np, dim))
  res = hf - e*f
  r = X(states_nrm2)(m, dim, res)
  deallocate(res)

end function X(states_residue)

! TODO use netcdf
subroutine X(states_write_restart)(filename, m, st, iter, v1, v2)
  character(len=*),  intent(in) :: filename
  type(mesh_type),   intent(IN) :: m
  type(states_type), intent(IN) :: st
  integer, optional, intent(in) :: iter ! used in TD
  FLOAT,   optional, intent(IN) :: v1(m%np, st%d%nspin), v2(m%np, st%d%nspin)

  integer :: iunit, ik, ist, id
  integer(i4) :: mode

  call push_sub('states_write_restart')

  call io_assign(iunit)
  if(st%restart_format == 1) then
    open(iunit, status='unknown', file=trim(filename), form='unformatted')
  else
    open(iunit, status='unknown', file=trim(filename) // '.ascii', form='formatted')
  end if

#ifdef R_TREAL
  mode = 0_i4
#else
  mode = 1_i4
#endif

  if(st%restart_format == 1) then
    write(iunit) int(m%box_shape, i4), m%h, m%rsize, m%xsize
    write(iunit) int(m%np, i4), int(st%d%dim, i4), int(st%st_start, i4), &
        int(st%st_end, i4), int(st%d%nik, i4), int(st%d%ispin, i4), mode
  else
    write(iunit,'(i12,7E23.14)') int(m%box_shape, i4), m%h, m%rsize, m%xsize
    write(iunit,'(7i12)') int(m%np, i4), int(st%d%dim, i4), int(st%st_start, i4), &
        int(st%st_end, i4), int(st%d%nik, i4), int(st%d%ispin, i4), mode
  end if

  ! psi has to be written in parts, or segmentation fault in some machines
  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end
      do id = 1, st%d%dim
        if(st%restart_format == 1) then
          write(iunit) st%X(psi)(1:m%np, id, ist, ik)
        else
          write(iunit,*) st%X(psi)(1:m%np, id, ist, ik)
        end if
      end do
    end do
  end do

  ! eigenvalues are also needed ;)
  do ik = 1, st%d%nik
    do ist = st%st_start, st%st_end
      if(st%restart_format == 1) then
        write(iunit) st%eigenval(ist, ik)
      else
        write(iunit,*) st%eigenval(ist, ik)
      end if
    end do
  end do

  if(present(iter)) then
    if(st%restart_format == 1) then
      write(iunit) int(iter, i4), v1, v2
    else
      write(iunit,*) int(iter, i4), v1, v2
    end if
  end if
  call io_close(iunit)
  
  call pop_sub()
end subroutine X(states_write_restart)

logical function X(states_load_restart)(filename, m, st, iter, v1, v2) result(ok)
  character(len=*),  intent(in)    :: filename
  type(mesh_type),   intent(IN)    :: m
  type(states_type), intent(inout) :: st
  integer, optional, intent(out)   :: iter ! used in TD
  FLOAT,   optional, intent(out)   :: v1(m%np, st%d%nspin), v2(m%np, st%d%nspin)

  integer :: iunit, ik, ist, id, restart_format
  integer(i4) :: ii, old_np, old_dim, old_start, old_end, old_nik, old_ispin, mode
  logical :: found
  character(len=256) :: fname
  FLOAT :: e
  FLOAT, allocatable :: dpsi(:)
  CMPLX, allocatable :: zpsi(:)

  call push_sub('states_load_restart')
  ok = .true.

  call io_assign(iunit)

  ! check what kind of restart file we have
  inquire(file=filename, exist=found)
  if(found) then ! unformatted file
    fname = trim(filename)
    restart_format = 1
    
    open(iunit, status='old', file=trim(fname), form='unformatted', err=998)
    read(iunit, err=999) ! mesh stuff is now skipped
    read(iunit, err=100) old_np, old_dim, old_start, old_end, old_nik, old_ispin, mode
    go to 200

    ! ugly thing to handle old restarts
100 continue
    backspace(iunit)
    read(iunit, err=999) old_np, old_dim, old_start, old_end, old_nik, old_ispin
    mode = -1 ! old restart file
200 continue

  else ! unformatted
    inquire(file=trim(filename) // '.ascii', exist=found)
    if(found) then ! unformatted file
      fname = trim(filename) // '.ascii'
      restart_format = 2

      open(iunit, status='old', file=trim(fname), form='formatted', err=998)
      read(iunit, *, err=999) ! mesh stuff is now skipped
      read(iunit, *, err=999) old_np, old_dim, old_start, old_end, old_nik, old_ispin, mode
    else ! error
      go to 999 ! one go to does not harm :)
    end if
  end if

  if(conf%verbose > 20) then
    write(stdout, '(3a)') "Info: Reading wavefunctions from file '", &
        trim(fname), "'"
  end if

  if(old_np.ne.m%np .or. old_dim.ne.st%dim .or. & ! different mesh, cannot read
       old_start.gt.st%st_start .or. old_end.lt.st%st_end .or. old_nik.ne.st%nik) then
    message(1) = 'Restart file has a different mesh!'
    write(message(2), '(a,i6,a,i6,a)') '  m%np        = ', m%np,        ' != ', old_np, ' or'
    write(message(3), '(a,i6,a,i6,a)') '  st%d%dim    = ', st%d%dim,    ' != ', old_dim, ' or'
    write(message(4), '(a,i6,a,i6,a)') '  st%st_start = ', st%st_start, ' <  ', old_start, ' or'
    write(message(5), '(a,i6,a,i6,a)') '  st%st_end   = ', st%st_end,   ' >  ', old_end, ' or'
    write(message(6), '(a,i6,a,i6)')   '  st%d%nik    = ', st%d%nik,    ' != ', old_nik

    call write_warning(6)
    go to 999 ! one go to does not harm :)
  end if

  st%X(psi) = R_TOTYPE(M_ZERO)

  if(mode == 0 .or. mode == -1) then
    allocate(dpsi(1:m%np))
  else
    allocate(zpsi(1:m%np))
  end if
  
  do ik = 1, st%d%nik
    do ist = old_start, old_end
      do id = 1, st%dim
        imode: if(mode == 0 .or. mode == -1) then
          if(restart_format == 1) then
            read(iunit, err=999) dpsi(:)
          else
            read(iunit, *, err=999) dpsi(:)
          end if
          if(ist >= st%st_start .and. ist <= st%st_end) then
#           ifdef R_TREAL
              st%dpsi(1:m%np, id, ist, ik) = dpsi(:)
#           else
              st%zpsi(1:m%np, id, ist, ik) = cmplx(dpsi(:), M_ZERO, PRECISION)
#           endif
          end if
        else
          read(iunit, err=999) zpsi(:)
          if(ist >= st%st_start .and. ist <= st%st_end) then
#           ifdef R_TREAL
              st%dpsi(1:m%np, id, ist, ik) = real(zpsi(:), PRECISION)
#           else
              st%zpsi(1:m%np, id, ist, ik) = zpsi(:)
#           endif
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
    read(iunit, err=999) st%eigenval(st%st_start:st%st_end, 1:st%d%nik)
  else
    do ik = 1, st%d%nik
      do ist = old_start, old_end
        if(restart_format == 1) then
          read(iunit, err=999) e
        else
          read(iunit, *, err=999) e
        end if
        if(ist >= st%st_start .and. ist <= st%st_end) then
          st%eigenval(ist, ik) = e
        end if
      end do
    end do
  end if

  if(present(iter)) then ! read the time-dependent stuff
    if(st%restart_format == 1) then
      read(iunit, err=999) ii, v1, v2
    else
      read(iunit, *, err=999) ii, v1, v2
    end if
    iter = int(ii)
  end if
  
  call io_close(iunit)
  call pop_sub()
  return

999 call io_close(iunit)
998 continue
  message(1) = 'Error reading from file ' // trim(filename) // "'"
  call write_warning(1)
  ok = .false.
  call pop_sub()
  return

end function X(states_load_restart)

subroutine X(states_output) (st, m, dir, outp)
  type(states_type), intent(IN) :: st
  type(mesh_type),   intent(IN) :: m
  character(len=*),  intent(in) :: dir
  type(output_type), intent(IN) :: outp

  integer :: ik, ist, idim, is
  character(len=80) :: fname
  FLOAT :: u
  FLOAT, allocatable :: dtmp(:)

  u = M_ONE/units_out%length%factor**conf%dim

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif
    if(outp%what(output_density)) then
      do is = 1, st%d%nspin
        write(fname, '(a,i1)') 'density-', is
        call doutput_function(outp%how, dir, fname, m, st%rho(:, is), u)
      end do
    end if
#ifdef HAVE_MPI
  end if
#endif

  if(outp%what(output_wfs)) then
    do ist = st%st_start, st%st_end
      is = outp%wfs((ist-1)/32 + 1)
      if(iand(is, 2**(modulo(ist-1, 32))).ne.0) then
        do ik = 1, st%d%nik
          do idim = 1, st%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'wf-', ik, '-', ist, '-', idim
            call X(output_function) (outp%how, dir, fname, m, &
                 st%X(psi) (1:, idim, ist, ik), sqrt(u))
          end do
        end do
      end if
    end do
  end if

  if(outp%what(output_wfs_sqmod)) then
    allocate(dtmp(m%np))
    do ist = 1, st%nst
      is = outp%wfs((ist-1)/32 + 1)
      if(iand(is, 2**(modulo(ist-1, 32))).ne.0) then
        do ik = 1, st%d%nik
          do idim = 1, st%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'sqm-wf-', ik, '-', ist, '-', idim
            dtmp = abs(st%X(psi) (:, idim, ist, ik))**2
            call doutput_function (outp%how, dir, fname, m, dtmp, u)
          end do
        end do
      end if
    end do
    deallocate(dtmp)
  end if
  
  if(conf%dim==3) then
    if(outp%what(output_elf))    call elf_FS(.true.,  'elf_rs')
    if(outp%what(output_elf_FS)) call elf_FS(.false., 'elf_fs')
  end if

contains
  subroutine mf2mf_RS2FS(m, fin, fout, c)
    type(mesh_type), intent(IN) :: m
    R_TYPE, intent(IN) :: fin(m%np)
    CMPLX, intent(out) :: fout(m%np)
    type(X(cf)), intent(inout) :: c

    call X(cf_alloc_RS) (c)
    call X(cf_alloc_FS) (c)
    call X(mf2cf) (m, fin, c)
    call X(cf_RS2FS) (c)
    call X(cf_FS2mf) (m, c, fout)
    call X(cf_free_RS) (c)
    call X(cf_free_FS) (c)
  end subroutine mf2mf_RS2FS

  subroutine elf_FS(rs, filename)
    logical, intent(in) :: rs
    character(len=*), intent(in) :: filename

    FLOAT :: f, d, s
    integer :: i, is, ik
    CMPLX, allocatable :: psi_fs(:), gpsi(:,:)
    FLOAT, allocatable :: c(:), r(:), gr(:,:), j(:,:)
    type(X(cf)) :: cf_tmp

    FLOAT, parameter :: dmin = CNST(1e-10)
    
    ! single or double occupancy
    if(st%d%nspin == 1) then
      s = M_TWO
    else
      s = M_ONE
    end if

    if(.not.rs) then
      call X(cf_new)(m%l, cf_tmp)
      call X(cf_fft_init)(cf_tmp)
    end if

    allocate(c(m%np))

    do_is: do is = 1, st%d%nspin
      allocate(r(m%np), gr(3, m%np), j(3, m%np))
      r = M_ZERO; gr = M_ZERO; j  = M_ZERO
      c = M_ZERO

      allocate(psi_fs(m%np), gpsi(3, m%np))
      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          do idim = 1, st%dim
            
            if(rs) then
              psi_fs(:) = cmplx(st%X(psi)(:, idim, ist, ik), KIND=PRECISION)
            else
              call mf2mf_RS2FS(m, st%X(psi)(:, idim, ist, ik), psi_fs(:), cf_tmp)
            end if
            call zf_gradient(m, psi_fs(:), gpsi)

            if(.not.rs) then
              do i = 1, conf%dim
                gpsi(i,:) = gpsi(i,:) * m%h(i)**2 * real(cf_tmp%n(i), PRECISION) / (M_TWO*M_PI)
              end do
            end if

            r(:) = r(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do i = 1, conf%dim
              gr(i,:) = gr(i,:) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   M_TWO * real(conjg(psi_fs(:))*gpsi(i,:))
              j(i,:)  =  j(i,:) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                   aimag(conjg(psi_fs(:))*gpsi(i,:))
            end do

            do i = 1, m%np
              if(r(i) >= dmin) then
                c(i) = c(i) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                     sum(abs(gpsi(1:conf%dim, i))**2)
              end if
            end do
          end do
        end do
      end do
      deallocate(psi_fs, gpsi)

      do i = 1, m%np
        if(r(i) >= dmin) then
          c(i) = c(i) - (M_FOURTH*sum(gr(1:conf%dim, i)**2) + sum(j(1:conf%dim, i)**2))/(s*r(i))
        end if
      end do

      ! normalization
      f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
      do i = 1, m%np
        if(abs(r(i)) >= dmin) then
          d    = f*(r(i)/s)**(M_FIVE/M_THREE)
          c(i) = M_ONE/(M_ONE + (c(i)/d)**2)
        else
          c(i) = M_ZERO
        end if
      end do
      
      deallocate(r, gr, j)
      write(fname, '(a,a,i1)') trim(filename), '-', is
      call doutput_function(outp%how, dir, trim(fname), m, c, M_ONE)
      
    end do do_is
    
    if(.not.rs) call X(cf_free)(cf_tmp)
    deallocate(c)
    
  end subroutine elf_FS

end subroutine X(states_output)

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
R_TYPE function X(states_mpdotp)(m, ik, st1, st2) result(dotp)
  type(mesh_type),   intent(IN) :: m
  integer,           intent(in) :: ik
  type(states_type), intent(IN) :: st1, st2

  R_TYPE, allocatable :: a(:, :)

  call push_sub('states_mpdotp')

  allocate(a(st1%nst, st1%nst))

  call X(calculate_matrix)(m, ik, st1, st2, st1%nst, a)
  dotp = lalg_det(a, st1%nst)

  select case(st1%d%ispin)
   case(UNPOLARIZED)
     dotp = dotp**2
   case(SPIN_POLARIZED)
     ! We assume that ik is an odd number (maybe this should be checked. So one
     ! spin component is ik and another spin component is ik + 1.
     call X(calculate_matrix)(m, ik+1, st1, st2, st1%nst, a)
     dotp = dotp*lalg_det(a, st1%nst)
  end select

  deallocate(a)
  call pop_sub()
  contains

    subroutine X(calculate_matrix)(m, ik, st1, st2, n, a)
      type(mesh_type),   intent(IN)  :: m
      integer,           intent(in)  :: ik
      type(states_type), intent(IN)  :: st1, st2
      integer,           intent(in)  :: n
      R_TYPE,            intent(out) :: a(n, n)

      integer :: i, j, dim
#if defined(HAVE_MPI) && defined(MPI_TD)
      R_TYPE, allocatable :: phi2(:, :)
      integer :: k, l, r, sn, sn1, ierr
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
               call MPI_ISEND(st2%X(psi)(1, 1, i, ik), st1%dim*m%np, MPI_DOUBLE_COMPLEX, &
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
            phi2(:, :) = st2%X(psi)(:, :, j, ik)
         endif
         do i = st1%st_start, st1%st_end
            a(i, j) = X(states_dotp)(m, dim, st1%X(psi)(1:, :, i, ik), &
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
            a(i, j) = X(states_dotp)(m, dim, st1%X(psi)(1:, :, i, ik), &
                                                  st2%X(psi)(1:, :, j, ik))
         enddo
      enddo
#endif
    end subroutine X(calculate_matrix)

end function X(states_mpdotp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine X(states_calculate_angular)(m, st, angular)
  type(mesh_type),   intent(IN)  :: m
  type(states_type), intent(IN)  :: st
  FLOAT,             intent(out) :: angular(3)

  FLOAT :: temp(3)
  R_TYPE, allocatable :: lpsi(:, :)
  integer :: idim, ik, j
#if defined(HAVE_MPI)
  integer :: i, ierr
#endif

  temp = M_ZERO
  allocate(lpsi(conf%dim, m%np))
  do ik = 1, st%d%nik
     do j  = st%st_start, st%st_end
        do idim = 1, st%dim
           call X(mesh_angular_momentum)(m, st%X(psi)(:, idim, j, ik), lpsi)
           temp(1) = temp(1) + st%occ(j, ik)*X(mf_dotp)(m, st%X(psi)(:, idim, j, ik), lpsi(1, :))
           temp(2) = temp(2) + st%occ(j, ik)*X(mf_dotp)(m, st%X(psi)(:, idim, j, ik), lpsi(2, :))
           temp(3) = temp(3) + st%occ(j, ik)*X(mf_dotp)(m, st%X(psi)(:, idim, j, ik), lpsi(3, :))
        enddo
     enddo
  enddo

#if defined(HAVE_MPI) && defined(MPI_TD)
    call MPI_ALLREDUCE(temp, angular, 3, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  angular = temp
#endif

  deallocate(lpsi)
end subroutine X(states_calculate_angular)

subroutine X(mesh_angular_momentum)(m, f, lf)
  type(mesh_type), intent(IN)  :: m
  R_TYPE,          intent(IN)  :: f(m%np)
  R_TYPE,          intent(out) :: lf(3, m%np)

  R_TYPE, allocatable :: gf(:, :)
  FLOAT :: x(3)
  integer :: i

  allocate(gf(3, m%np))
  call X(f_gradient)(m, f, grad = gf)

  do i = 1, m%np
     call mesh_xyz(m, i, x)
     lf(1, i) = (x(2)*gf(3, i)-x(3)*gf(2, i))
     lf(2, i) = (x(3)*gf(1, i)-x(1)*gf(3, i))
     lf(3, i) = (x(1)*gf(2, i)-x(2)*gf(1 ,i))
  enddo
#if defined(R_TCOMPLEX)
  call X(lalg_scal)(3*m%np, -M_zI, lf(1, 1))
#endif
  deallocate(gf)
end subroutine X(mesh_angular_momentum)
