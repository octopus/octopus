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
  CMPLX   :: c

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

          c = st%d%kweights(ik)*st%occ(p, ik) * &
             st%X(psi)(i, 1, p, ik) * R_CONJ(st%X(psi)(i, 2, p, ik))
          rho(i, 3) = rho(i, 3) + R_REAL(c)
          rho(i, 4) = rho(i, 4) + R_AIMAG(c)
        end select
      end do
    end do
  end do

#if defined(HAVE_MPI)
  ! reduce density (assumes memory is contiguous)
  if(present(reduce)) then
    if(reduce) then
      allocate(reduce_rho(1:np, st%d%nspin))
      call MPI_ALLREDUCE(rho(1, 1), reduce_rho(1, 1), np*st%d%nspin, &
           MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
      rho = reduce_rho
      deallocate(reduce_rho)
    end if
  end if
#endif

  call pop_sub()
  return
end subroutine X(calcdens)

! Orthonormalizes nst orbital in mesh m
subroutine X(states_gram_schmidt)(nst, m, dim, psi, start)
  integer,           intent(in)    :: nst, dim
  type(mesh_type),   intent(in)    :: m
  R_TYPE,            intent(inout) :: psi(:,:,:)   ! psi(m%np, dim, nst)
  integer, optional, intent(in)    :: start

  integer :: p, q, stst
  FLOAT :: nrm2
  R_TYPE :: ss

  call push_sub('states_gram_schmidt')

  if(present(start)) then
    stst = start
  else
    stst = 1
  end if

  do p = stst, nst
    do q = 1, p - 1
      ss = X(states_dotp)(m, dim, psi(:,:, q), psi(:,:, p))
      call lalg_axpy(m%np, dim, -ss, psi(:,:, q), psi(:,:, p))
    end do

    nrm2 = X(states_nrm2)(m, dim, psi(:,:, p))
    ss = R_TOTYPE(M_ONE/nrm2)
    call lalg_scal(m%np, dim, ss, psi(:, :, p))
  end do

  call pop_sub()
end subroutine X(states_gram_schmidt)

R_TYPE function X(states_dotp)(m, dim, f1, f2) result(dotp)
  type(mesh_type), intent(IN) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(IN) :: f1(:,:), f2(:,:)

  integer :: i

  dotp = R_TOTYPE(M_ZERO)
  do i = 1, dim
    dotp = dotp + X(mf_dotp)(m, f1(:,i), f2(:,i))
  end do

end function X(states_dotp)

FLOAT function X(states_nrm2)(m, dim, f) result(nrm2)
  type(mesh_type), intent(IN) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(IN) :: f(:,:)

  integer :: i

  nrm2 = M_ZERO
  do i = 1, dim
    nrm2 = nrm2 + X(mf_nrm2)(m, f(:,i))**2
  end do
  nrm2 = sqrt(nrm2)

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

subroutine X(states_output) (st, m, f_der, dir, outp)
  type(states_type),   intent(in)    :: st
  type(mesh_type),     intent(in)    :: m
  type(f_der_type),    intent(inout) :: f_der
  character(len=*),    intent(in)    :: dir
  type(output_type),   intent(in)    :: outp

  integer :: ik, ist, idim, is, ierr
  character(len=80) :: fname
  FLOAT :: u
  FLOAT, allocatable :: dtmp(:)

  call push_sub('states_output')

  u = M_ONE/units_out%length%factor**conf%dim

#ifdef HAVE_MPI
  if(mpiv%node == 0) then
#endif
    if(outp%what(output_density)) then
      do is = 1, st%d%nspin
        write(fname, '(a,i1)') 'density-', is
        call doutput_function(outp%how, dir, fname, m, st%rho(:, is), u, ierr)
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
          do idim = 1, st%d%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'wf-', ik, '-', ist, '-', idim
            call X(output_function) (outp%how, dir, fname, m, &
               st%X(psi) (1:, idim, ist, ik), sqrt(u), ierr)
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
          do idim = 1, st%d%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'sqm-wf-', ik, '-', ist, '-', idim
            dtmp = abs(st%X(psi) (:, idim, ist, ik))**2
            call doutput_function (outp%how, dir, fname, m, dtmp, u, ierr)
          end do
        end do
      end if
    end do
    deallocate(dtmp)
  end if
  
  if(conf%dim==3) then
    if(outp%what(output_elf))    call elf(.true.,  'elf_rs')
    if(outp%what(output_elf_FS)) call elf(.false., 'elf_fs')
  end if

  call pop_sub()
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

  subroutine elf(rs, filename)
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
      allocate(r(m%np), gr(m%np, conf%dim), j(m%np, conf%dim))
      r = M_ZERO; gr = M_ZERO; j  = M_ZERO
      c = M_ZERO

      allocate(psi_fs(m%np), gpsi(m%np, conf%dim))
      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          do idim = 1, st%d%dim
            
            if(rs) then
              psi_fs(:) = cmplx(st%X(psi)(:, idim, ist, ik), KIND=PRECISION)
            else
              call mf2mf_RS2FS(m, st%X(psi)(:, idim, ist, ik), psi_fs(:), cf_tmp)
            end if
            call zf_gradient(f_der, psi_fs(:), gpsi)

            if(.not.rs) then
              do i = 1, conf%dim
                gpsi(:,i) = gpsi(:,i) * m%h(i)**2 * real(cf_tmp%n(i), PRECISION) / (M_TWO*M_PI)
              end do
            end if

            r(:) = r(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do i = 1, conf%dim
              gr(:,i) = gr(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                 M_TWO * real(conjg(psi_fs(:))*gpsi(:,i))
              j (:,i) =  j(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                 aimag(conjg(psi_fs(:))*gpsi(:,i))
            end do

            do i = 1, m%np
              if(r(i) >= dmin) then
                c(i) = c(i) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                   sum(abs(gpsi(i, 1:conf%dim))**2)
              end if
            end do
          end do
        end do
      end do
      deallocate(psi_fs, gpsi)

      do i = 1, m%np
        if(r(i) >= dmin) then
          c(i) = c(i) - (M_FOURTH*sum(gr(i, 1:conf%dim)**2) + sum(j(i, 1:conf%dim)**2))/(s*r(i))
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
      call doutput_function(outp%how, dir, trim(fname), m, c, M_ONE, ierr)
      
    end do do_is
    
    if(.not.rs) call X(cf_free)(cf_tmp)
    deallocate(c)
    
  end subroutine elf

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
#if defined(HAVE_MPI)
    R_TYPE, allocatable :: phi2(:, :)
    integer :: k, l, r, sn, sn1, ierr
    integer :: status(MPI_STATUS_SIZE)
    integer :: request, node(st1%nst)
#endif
    
    dim = st1%d%dim
#if defined(HAVE_MPI)
    
    ! This code just builds an integer array, node, that assigns to each
    ! state the process that holds it. If such a thing is ever needed for any other
    ! purpose, maybe it should go to other place, more general (such as "st%node")
    sn  = st1%nst/mpiv%numprocs
    sn1 = sn + 1
    r  = mod(st1%nst, mpiv%numprocs)
    do j = 1, r
      node((j-1)*sn1+1:j*sn1) = j - 1
    end do
    k = sn1*r
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do j = 1, mpiv%numprocs - r
      node(k+(j-1)*sn+1:k+j*sn) = r + j - 1
    enddo

    ! Each process sends the states in st2 to the rest of the processes.
    do i = st1%st_start, st1%st_end
      do j = 0, mpiv%numprocs-1
        if(mpiv%node.ne.j) then
          call MPI_ISEND(st2%X(psi)(1, 1, i, ik), st1%d%dim*m%np, R_MPITYPE, &
             j, i, MPI_COMM_WORLD, request, ierr)
        end if
      end do
    end do

    ! Processes are received, and then the matrix elements are calculated.
    allocate(phi2(m%np, st1%d%dim))
    do j = 1, n
      l = node(j)
      if(l.ne.mpiv%node) then
        call MPI_RECV(phi2(1, 1), st1%d%dim*m%np, R_MPITYPE, &
           l, j, MPI_COMM_WORLD, status, ierr)
      else
        phi2(:, :) = st2%X(psi)(:, :, j, ik)
      end if
      do i = st1%st_start, st1%st_end
        a(i, j) = X(states_dotp)(m, dim, st1%X(psi)(1:, :, i, ik), phi2(1:, :))
      end do
    end do
    deallocate(phi2)

    ! Each process holds some lines of the matrix. So it is broadcasted (All processes
    ! should get the whole matrix)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    do i = 1, n
      k = node(i)
      do j = 1, n
        call MPI_BCAST(a(i, j), 1, R_MPITYPE, k, MPI_COMM_WORLD, ierr)
      enddo
    enddo
#else
    do i = 1, n
      do j = 1, n
        a(i, j) = X(states_dotp)(m, dim, st1%X(psi)(1:, :, i, ik), &
           st2%X(psi)(1:, :, j, ik))
      end do
    end do
#endif
  end subroutine X(calculate_matrix)

end function X(states_mpdotp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine X(states_calculate_magnetization)(m, st, mag)
  type(mesh_type),   intent(in)  :: m
  type(states_type), intent(in)  :: st
  FLOAT,             intent(out) :: mag(3)

  integer :: ik, ist, i
  FLOAT :: sign, temp(3)
  R_TYPE :: c
#if defined(HAVE_MPI)
  integer :: ierr
#endif

  call push_sub('states_calculate_magnetization')

  select case(st%d%ispin)
  case(SPIN_POLARIZED) ! collinear spin
    sign = M_ONE
    temp = M_ZERO
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        temp(3) = temp(3) + sign*st%d%kweights(ik)*st%occ(ist, ik)
      end do
      sign = -sign
    end do

  case(SPINORS) ! non-collinear
    temp = M_ZERO
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        do i = 1, m%np
          c = st%X(psi) (i, 1, ist, ik) * R_CONJ(st%X(psi) (i, 2, ist, ik))* m%vol_pp(i)

          temp(1) = temp(1) + st%d%kweights(ik)*st%occ(ist, ik)* M_TWO*R_REAL(c)
          temp(2) = temp(2) - st%d%kweights(ik)*st%occ(ist, ik)* M_TWO*R_AIMAG(c)
          c = R_ABS(st%X(psi) (i, 1, ist, ik))**2 - R_ABS(st%X(psi) (i, 2, ist, ik))**2
          temp(3) = temp(3) + st%d%kweights(ik)*st%occ(ist, ik)* R_REAL(c)
        end do
      end do
    end do
    temp = temp

  end select

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE(temp, mag, 3, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  mag = temp
#endif

  call pop_sub()
end subroutine X(states_calculate_magnetization)

subroutine X(states_calculate_angular)(m, f_der, st, angular, l2)
  type(mesh_type),   intent(IN)  :: m
  type(f_der_type), intent(inout) :: f_der
  type(states_type), intent(IN)  :: st
  FLOAT,             intent(out) :: angular(3)
  FLOAT, optional,   intent(out) :: l2

  FLOAT :: temp(3), ltemp
  R_TYPE, allocatable :: lpsi(:, :)
  integer :: idim, ik, j
#if defined(HAVE_MPI)
  integer :: i, ierr
#endif

  call push_sub('states_calculate_angular')


  temp = M_ZERO; ltemp = M_ZERO
  allocate(lpsi(m%np, conf%dim))
  do ik = 1, st%d%nik
    do j  = st%st_start, st%st_end
      do idim = 1, st%d%dim
#if defined(R_TREAL)
        temp = M_ZERO ! The expectation value of L of *any* real function is null
#else
        call X(f_angular_momentum)(f_der, st%X(psi)(:, idim, j, ik), lpsi)
        temp(1) = temp(1) + st%occ(j, ik)*X(mf_dotp)(m, st%X(psi)(:, idim, j, ik), lpsi(:, 1))
        temp(2) = temp(2) + st%occ(j, ik)*X(mf_dotp)(m, st%X(psi)(:, idim, j, ik), lpsi(:, 2))
        temp(3) = temp(3) + st%occ(j, ik)*X(mf_dotp)(m, st%X(psi)(:, idim, j, ik), lpsi(:, 3))
#endif
        if(present(l2)) then
          call X(f_l2)(f_der, st%X(psi)(:, idim, j, ik), lpsi(:, 1))
          ltemp = ltemp + st%occ(j, ik)*X(mf_dotp)(m, st%X(psi)(:, idim, j, ik), lpsi(:, 1))
        end if
      end do
    end do
  end do

#if defined(HAVE_MPI)
  call MPI_ALLREDUCE(temp, angular, 3, &
     MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
  if(present(l2)) call MPI_ALLREDUCE(ltemp, l2, 1, &
     MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
  angular = temp
  if(present(l2)) l2 = ltemp
#endif

  deallocate(lpsi)
  call pop_sub()

end subroutine X(states_calculate_angular)
