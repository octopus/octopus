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
!!
!! $Id$

! ---------------------------------------------------------
! Allocates the KS wavefunctions defined within an states_t
! structure.
subroutine X(states_allocate_wfns)(st, m)
  type(states_t), intent(inout) :: st
  type(mesh_t),    intent(in)    :: m
  integer :: n

  n = m%np_part * st%d%dim * (st%st_end-st%st_start+1) * st%d%nik
  ALLOCATE(st%X(psi)(m%np_part, st%d%dim, st%st_start:st%st_end, st%d%nik), n)
  st%X(psi) = R_TOTYPE(M_ZERO)

end subroutine X(states_allocate_wfns)

! ---------------------------------------------------------
! Calculates the new density out the wavefunctions and
! occupations...
subroutine X(states_calc_dens)(st, np, rho)
  type(states_t), intent(in)  :: st
  integer,        intent(in)  :: np
  FLOAT,          intent(out) :: rho(:,:)

  integer :: i, ik, p, sp
  CMPLX   :: c

#ifdef HAVE_MPI
  FLOAT,  allocatable :: reduce_rho(:)
  integer :: ierr
#endif

  call push_sub('states_inc.states_calc_dens')

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
  ! reduce density
  if(st%parallel_in_states) then
    ALLOCATE(reduce_rho(1:np), np)
    do i = 1, st%d%nspin
      call MPI_ALLREDUCE(rho(1, i), reduce_rho(1), np, &
         MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, ierr)
      rho(1:np, i) = reduce_rho(1:np)
    end do
    deallocate(reduce_rho)
  end if
#endif

  call pop_sub()
  return
end subroutine X(states_calc_dens)


! ---------------------------------------------------------
! Orthonormalizes nst orbital in mesh m
subroutine X(states_gram_schmidt1)(nst, m, dim, psi, start)
  integer,           intent(in)    :: nst, dim
  type(mesh_t),      intent(in)    :: m
  R_TYPE,            intent(inout) :: psi(:,:,:)   ! psi(m%np_part, dim, nst)
  integer, optional, intent(in)    :: start

  integer :: p, q, stst, idim
  FLOAT   :: nrm2
  R_TYPE  :: ss

  call profiling_in(C_PROFILING_GRAM_SCHMIDT)
  call push_sub('states_inc.states_gram_schmidt')

  if(present(start)) then
    stst = start
  else
    stst = 1
  end if

  do p = stst, nst
    do q = 1, p - 1
      ss = X(states_dotp)(m, dim, psi(:,:, q), psi(:,:, p))
      do idim = 1, dim
        call lalg_axpy(m%np, -ss, psi(:, idim, q), psi(:, idim, p))
      end do
    end do

    nrm2 = X(states_nrm2)(m, dim, psi(:,:, p))
    ss = R_TOTYPE(M_ONE/nrm2)
    do idim = 1, dim
      call lalg_scal(m%np, ss, psi(:, idim, p))
    end do
  end do

  call pop_sub()
  call profiling_out(C_PROFILING_GRAM_SCHMIDT)
end subroutine X(states_gram_schmidt1)

! ---------------------------------------------------------
! Orthonormalizes phi to the nst orbitals psi.
! It also permits to do only the orthogonalization (no normalization).
! And one can pass an extra optional argument, mask, which:
!  - on input, if mask(p) = .true., the p-orbital is not used.
!  - on output, mask(p) = .true. if p was already orthogonal (to within 1e-12).
subroutine X(states_gram_schmidt2)(nst, m, dim, psi, phi, normalize, mask)
  integer,           intent(in)    :: nst, dim
  type(mesh_t),      intent(in)    :: m
  R_TYPE,            intent(inout) :: psi(:,:,:)   ! psi(m%np_part, dim, nst)
  R_TYPE,            intent(inout) :: phi(:,:)     ! phi(m%np_part, dim)
  logical, optional, intent(in)    :: normalize
  logical, optional, intent(inout) :: mask(:)      ! nst

  logical :: normalize_
  integer :: q, idim
  FLOAT   :: nrm2
  R_TYPE  :: ss

  call profiling_in(C_PROFILING_GRAM_SCHMIDT)
  call push_sub('states_inc.states_gram_schmidt')

  do q = 1, nst
    if(present(mask)) then
      if(mask(q)) cycle
    end if
    ss = X(states_dotp)(m, dim, psi(:,:, q), phi)
    if(abs(ss) > CNST(1.0e-13)) then
      do idim = 1, dim
        call lalg_axpy(m%np, -ss, psi(:, idim, q), phi(:, idim))
      end do
    else
      if(present(mask)) mask(q) = .true.
    end if
  end do

  normalize_ = .false.
  if(present(normalize)) normalize_ = normalize
  if(normalize) then
    nrm2 = X(states_nrm2)(m, dim, phi)
    ss = R_TOTYPE(M_ONE/nrm2)
    do idim = 1, dim
      call lalg_scal(m%np, ss, phi(:, idim))
    end do
  end if

  call pop_sub()
  call profiling_out(C_PROFILING_GRAM_SCHMIDT)
end subroutine X(states_gram_schmidt2)


! ---------------------------------------------------------
R_TYPE function X(states_dotp)(m, dim, f1, f2) result(dotp)
  type(mesh_t),    intent(in) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(in) :: f1(:,:), f2(:,:)

  integer :: idim

  call push_sub('states_inc.Xstates_dotp')

  dotp = R_TOTYPE(M_ZERO)
  do idim = 1, dim
    dotp = dotp + X(mf_dotp)(m, f1(:, idim), f2(:, idim))
  end do

  call pop_sub()

end function X(states_dotp)


! ---------------------------------------------------------
subroutine X(states_normalize_orbital)(m, dim, psi)
  type(mesh_t),    intent(in)  :: m
  integer,         intent(in)  :: dim
  R_TYPE,          intent(out) :: psi(:,:)

  FLOAT   :: norm
  integer :: idim

  norm = X(states_nrm2) (m, dim, psi)
  norm = sqrt(norm)

  do idim = 1, dim
    psi(1:m%np, idim) = psi(1:m%np, idim)/norm
  end do

end subroutine X(states_normalize_orbital)


! ---------------------------------------------------------
FLOAT function X(states_nrm2)(m, dim, f) result(nrm2)
  type(mesh_t),    intent(in) :: m
  integer,         intent(in) :: dim
  R_TYPE,          intent(in) :: f(:,:)

  integer :: idim

  call push_sub('states_inc.Xstates_nrm2')

  nrm2 = M_ZERO
  do idim = 1, dim
    nrm2 = nrm2 + X(mf_nrm2)(m, f(:, idim))**2
  end do
  nrm2 = sqrt(nrm2)

  call pop_sub()

end function X(states_nrm2)


! ---------------------------------------------------------
FLOAT function X(states_residue)(m, dim, hf, e, f) result(r)
  type(mesh_t),      intent(in)  :: m
  integer,           intent(in)  :: dim
  R_TYPE,            intent(in)  :: hf(:,:), f(:,:)
  FLOAT,             intent(in)  :: e

  R_TYPE, allocatable :: res(:,:)

  call push_sub('states_inc.Xstates_nrm2')

  ALLOCATE(res(m%np_part, dim), m%np_part*dim)

  res(1:m%np, :) = hf(1:m%np, :) - e*f(1:m%np, :)

  r = X(states_nrm2)(m, dim, res)
  deallocate(res)

  call pop_sub()

end function X(states_residue)


! ---------------------------------------------------------
subroutine X(states_output) (st, gr, dir, outp)
  type(states_t),   intent(inout) :: st
  type(grid_t),     intent(inout) :: gr
  character(len=*), intent(in)    :: dir
  type(output_t),   intent(in)    :: outp

  integer :: ik, ist, idim, is, id, ierr
  character(len=80) :: fname
  FLOAT :: u
  FLOAT, allocatable :: dtmp(:)

  call push_sub('states_inc.states_output')

  u = M_ONE/units_out%length%factor**NDIM

  if(iand(outp%what, output_density).ne.0) then
    do is = 1, st%d%nspin
      write(fname, '(a,i1)') 'density-', is
      call doutput_function(outp%how, dir, fname, gr%m, gr%sb, st%rho(:, is), u, ierr)
    end do
  end if

#ifdef COMPLEX_WFNS
  if(iand(outp%what, output_current).ne.0) then
    ! calculate current first
    call calc_paramagnetic_current(gr, st, st%j)
    do is = 1, st%d%nspin
      do id = 1, NDIM
        write(fname, '(a,i1,a,a)') 'current-', is, '-', index2axis(id)
        call doutput_function(outp%how, dir, fname, gr%m, gr%sb, st%j(:, id, is), u, ierr)
      end do
    end do
  end if
#endif

  if(iand(outp%what, output_wfs).ne.0) then
    do ist = st%st_start, st%st_end
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = 1, st%d%nik
          do idim = 1, st%d%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'wf-', ik, '-', ist, '-', idim
            call X(output_function) (outp%how, dir, fname, gr%m, gr%sb, &
              st%X(psi) (1:, idim, ist, ik), sqrt(u), ierr)
          end do
        end do
      end if
    end do
  end if

  if(iand(outp%what, output_wfs_sqmod).ne.0) then
    ALLOCATE(dtmp(NP_PART), NP_PART)
    do ist = 1, st%nst
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = 1, st%d%nik
          do idim = 1, st%d%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'sqm-wf-', ik, '-', ist, '-', idim
            dtmp = abs(st%X(psi) (:, idim, ist, ik))**2
            call doutput_function (outp%how, dir, fname, gr%m, gr%sb, dtmp, u, ierr)
          end do
        end do
      end if
    end do
    deallocate(dtmp)
  end if

  if(NDIM==3) then
    if(iand(outp%what, output_elf).ne.0)    call elf(.true.,  'elf_rs')
    if(iand(outp%what, output_elf_FS).ne.0) call elf(.false., 'elf_fs')
  end if

  call pop_sub()
contains

  ! ---------------------------------------------------------
  subroutine mf2mf_RS2FS(m, fin, fout, c)
    type(mesh_t),  intent(in)    :: m
    R_TYPE,        intent(in)    :: fin(:)
    CMPLX,         intent(out)   :: fout(:)
    type(X(cf_t)), intent(inout) :: c

    call X(cf_alloc_RS) (c)
    call X(cf_alloc_FS) (c)
    call X(mf2cf) (m, fin, c)
    call X(cf_RS2FS) (c)
    call X(cf_FS2mf) (m, c, fout)
    call X(cf_free_RS) (c)
    call X(cf_free_FS) (c)
  end subroutine mf2mf_RS2FS


  ! ---------------------------------------------------------
  subroutine elf(rs, filename)
    logical, intent(in) :: rs
    character(len=*), intent(in) :: filename

    FLOAT :: f, d, s
    integer :: i, is, ik
    CMPLX, allocatable :: psi_fs(:), gpsi(:,:)
    FLOAT, allocatable :: c(:), r(:), gradr(:,:), j(:,:)
    type(X(cf_t)) :: cf_tmp

    FLOAT, parameter :: dmin = CNST(1e-10)

    ! single or double occupancy
    if(st%d%nspin == 1) then
      s = M_TWO
    else
      s = M_ONE
    end if

    if(.not.rs) then
      call X(cf_new)(gr%m%l, cf_tmp)
      call X(cf_fft_init)(cf_tmp, gr%sb)
    end if

    ALLOCATE(c(NP), NP)

    do_is: do is = 1, st%d%nspin
      ALLOCATE(    r(NP),       NP)
      ALLOCATE(gradr(NP, NDIM), NP*NDIM)
      ALLOCATE(    j(NP, NDIM), NP*NDIM)
      r = M_ZERO; gradr = M_ZERO; j  = M_ZERO
      c = M_ZERO

      ALLOCATE(psi_fs(NP_PART),  NP_PART)
      ALLOCATE(gpsi  (NP, NDIM), NP*NDIM)
      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          do idim = 1, st%d%dim

            if(rs) then
              psi_fs(:) = cmplx(st%X(psi)(:, idim, ist, ik), KIND=PRECISION)
            else
              call mf2mf_RS2FS(gr%m, st%X(psi)(:, idim, ist, ik), psi_fs(:), cf_tmp)
            end if
            call zf_gradient(gr%sb, gr%f_der, psi_fs(:), gpsi)

            if(.not.rs) then
              do i = 1, NDIM
                gpsi(:,i) = gpsi(:,i) * gr%m%h(i)**2 * real(cf_tmp%n(i), PRECISION) / (M_TWO*M_PI)
              end do
            end if

            r(:) = r(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do i = 1, NDIM
              gradr(:,i) = gradr(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                M_TWO * real(conjg(psi_fs(:))*gpsi(:,i))
              j (:,i) =  j(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                aimag(conjg(psi_fs(:))*gpsi(:,i))
            end do

            do i = 1, NP
              if(r(i) >= dmin) then
                c(i) = c(i) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                  sum(abs(gpsi(i, 1:NDIM))**2)
              end if
            end do
          end do
        end do
      end do
      deallocate(psi_fs, gpsi)

      do i = 1, NP
        if(r(i) >= dmin) then
          c(i) = c(i) - (M_FOURTH*sum(gradr(i, 1:NDIM)**2) + sum(j(i, 1:NDIM)**2))/(s*r(i))
        end if
      end do

      ! normalization
      f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
      do i = 1, NP
        if(abs(r(i)) >= dmin) then
          d    = f*(r(i)/s)**(M_FIVE/M_THREE)
          c(i) = M_ONE/(M_ONE + (c(i)/d)**2)
        else
          c(i) = M_ZERO
        end if
      end do

      deallocate(r, gradr, j)
      write(fname, '(a,a,i1)') trim(filename), '-', is
      call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, c, M_ONE, ierr)

    end do do_is

    if(.not.rs) call X(cf_free)(cf_tmp)
    deallocate(c)

  end subroutine elf

end subroutine X(states_output)


! -------------------------------------------------------------
! Returns the dot product of two many-body states st1 and
! st2, for a given irreducible subspace (k-point) ik.
! Warning: it does not do any check on the "conformance" of
!  the two states. (maybe a subroutine should be written
!  for that purpose...)
! Warning: it disregards completely the occupation number
!  problem. This is not important, unless the occupations
!  are different in the two states (this is not the case
!  for the moment), or if any of the occupation is null
!  (this can be problem, and will have to be cared about.
! -------------------------------------------------------------
R_TYPE function X(states_mpdotp)(m, ik, st1, st2) result(dotp)
  type(mesh_t),   intent(in) :: m
  integer,        intent(in) :: ik
  type(states_t), intent(in) :: st1, st2

  R_TYPE, allocatable :: a(:, :)

  call push_sub('states_inc.states_mpdotp')

  ALLOCATE(a(st1%nst, st1%nst), st1%nst*st1%nst)

  call X(calculate_matrix)(m, ik, st1, st2, st1%nst, a)
  dotp = lalg_determinant(st1%nst, a, invert = .false.)

  select case(st1%d%ispin)
  case(UNPOLARIZED)
    dotp = dotp**2
  case(SPIN_POLARIZED)
    ! We assume that ik is an odd number (maybe this should be checked. So one
    ! spin component is ik and another spin component is ik + 1.
    call X(calculate_matrix)(m, ik+1, st1, st2, st1%nst, a)
    dotp = dotp*lalg_determinant(st1%nst, a, invert = .false.)
  end select

  deallocate(a)
  call pop_sub()

contains

  ! ---------------------------------------------------------
  subroutine X(calculate_matrix)(m, ik, st1, st2, n, a)
    type(mesh_t),   intent(in)  :: m
    integer,        intent(in)  :: ik
    type(states_t), intent(in)  :: st1, st2
    integer,        intent(in)  :: n
    R_TYPE,         intent(out) :: a(n, n)

    integer :: i, j, dim
#if defined(HAVE_MPI)
    R_TYPE, allocatable :: phi2(:, :)
    integer :: k, l, ierr
    integer :: status(MPI_STATUS_SIZE)
    integer :: request
#endif

    dim = st1%d%dim
#if defined(HAVE_MPI)
    call mpi_barrier(st1%mpi_grp%comm, ierr)
    ! Each process sends the states in st2 to the rest of the processes.
    do i = st1%st_start, st1%st_end
      do j = 0, st1%mpi_grp%size - 1
        if(st1%mpi_grp%rank.ne.j) then
          call mpi_isend(st2%X(psi)(1, 1, i, ik), st1%d%dim*m%np, R_MPITYPE, &
            j, i, st1%mpi_grp%comm, request, ierr)
        end if
      end do
    end do

    ! Processes are received, and then the matrix elements are calculated.
    ALLOCATE(phi2(m%np, st1%d%dim), m%np*st1%d%dim)
    do j = 1, n
      l = st1%node(j)
      if(l.ne.st1%mpi_grp%rank) then
        call mpi_irecv(phi2(1, 1), st1%d%dim*m%np, R_MPITYPE, l, j, st1%mpi_grp%comm, request, ierr)
        call mpi_wait(request, status, ierr)
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
    call MPI_BARRIER(st1%mpi_grp%comm, ierr)
    do i = 1, n
      k = st1%node(i)
      do j = 1, n
        call MPI_BCAST(a(i, j), 1, R_MPITYPE, k, st1%mpi_grp%comm, ierr)
      end do
    end do
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


! ---------------------------------------------------------
subroutine X(states_calc_angular)(gr, st, angular, l2)
  type(grid_t),    intent(inout) :: gr
  type(states_t),  intent(inout) :: st
  FLOAT,           intent(out)   :: angular(MAX_DIM)
  FLOAT, optional, intent(out)   :: l2

  FLOAT :: temp(MAX_DIM), ltemp
  R_TYPE, allocatable :: lpsi(:, :)
  integer :: idim, ik, j
#if defined(HAVE_MPI)
  integer :: mpi_err
#endif

  call push_sub('states_inc.states_calc_angular')

  temp = M_ZERO; ltemp = M_ZERO
  ALLOCATE(lpsi(NP, NDIM), NP*NDIM)

  do ik = 1, st%d%nik
    do j  = st%st_start, st%st_end
      do idim = 1, st%d%dim
#if defined(R_TREAL)
        temp = M_ZERO ! The expectation value of L of *any* real function is null
#else
        call X(f_angular_momentum)(gr%sb, gr%f_der, st%X(psi)(:, idim, j, ik), lpsi)

        temp(1) = temp(1) + st%occ(j, ik)*X(mf_dotp)(gr%m, st%X(psi)(:, idim, j, ik), lpsi(:, 1))
        temp(2) = temp(2) + st%occ(j, ik)*X(mf_dotp)(gr%m, st%X(psi)(:, idim, j, ik), lpsi(:, 2))
        temp(3) = temp(3) + st%occ(j, ik)*X(mf_dotp)(gr%m, st%X(psi)(:, idim, j, ik), lpsi(:, 3))
#endif
        if(present(l2)) then
          call X(f_l2)(gr%sb, gr%f_der, st%X(psi)(:, idim, j, ik), lpsi(:, 1))
          ltemp = ltemp + st%occ(j, ik)*X(mf_dotp)(gr%m, st%X(psi)(:, idim, j, ik), lpsi(:, 1))
        end if
      end do
    end do
  end do

  if(st%parallel_in_states) then
#if defined(HAVE_MPI)
    call MPI_ALLREDUCE(temp, angular, 3, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
    if(present(l2)) call MPI_ALLREDUCE(ltemp, l2, 1, MPI_FLOAT, MPI_SUM, st%mpi_grp%comm, mpi_err)
#endif
  else
    angular = temp
    if(present(l2)) l2 = ltemp
  end if

  deallocate(lpsi)

  call pop_sub()
end subroutine X(states_calc_angular)
