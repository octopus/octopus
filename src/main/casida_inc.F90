!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2011-2013 D. Strubbe
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
!! $Id$

subroutine X(oscillator_strengths)(cas, mesh, st)
  type(casida_t), intent(inout) :: cas
  type(mesh_t), intent(in) :: mesh
  type(states_t), intent(in) :: st

  FLOAT, allocatable :: deltav(:)
  R_TYPE, allocatable :: xx(:), psi_a(:)
  CMPLX, allocatable :: zf(:), zx(:)
  FLOAT, allocatable :: gaus_leg_points(:), gaus_leg_weights(:)
  integer :: ii, jj, ia, ip, idir
  FLOAT :: theta, phi, qlen
  FLOAT :: qvect(MAX_DIM)

  PUSH_SUB(X(oscillator_strengths))

  if(cas%triplet) then
    cas%X(tm)(:,:) = M_ZERO
    cas%f(:) = M_ZERO

    if(cas%qcalc) then
      cas%qf(:) = M_ZERO
      if(cas%avg_order > 0) cas%qf_avg(:) = M_ZERO
    endif

    POP_SUB(X(oscillator_strengths))
    return
  endif

  call profiling_in(prof, "CASIDA_OSCILLATOR_STRENGTHS")

  if(cas%qcalc) then
    SAFE_ALLOCATE(zf(1:mesh%np))
    SAFE_ALLOCATE(psi_a(1:mesh%np))
    SAFE_ALLOCATE(zx(1:cas%n_pairs))

    ! matrix element
    do ia = 1, cas%n_pairs

      call states_get_state(st, mesh, 1, cas%pair(ia)%i, cas%pair(ia)%sigma, zf)
      call states_get_state(st, mesh, 1, cas%pair(ia)%a, cas%pair(ia)%sigma, psi_a)

      do ip = 1, mesh%np
        zf(ip) = exp(M_zI*dot_product(cas%qvector(1:mesh%sb%dim), mesh%x(ip, 1:mesh%sb%dim)))*aimag(zf(ip))*psi_a(ip)
      end do

      zx(ia) = zmf_integrate(mesh, zf)
    end do

    ! intensity
    do ia = 1, cas%n_pairs
      cas%qf(ia) = abs(ztransition_matrix_element(cas, ia, zx))**2
    end do

    ! do we calculate the average
    if(cas%avg_order > 0) then

      ! use Gauss-Legendre quadrature scheme
      SAFE_ALLOCATE(gaus_leg_points (1:cas%avg_order))
      SAFE_ALLOCATE(gaus_leg_weights(1:cas%avg_order))
      call gauss_legendre_points(cas%avg_order, gaus_leg_points, gaus_leg_weights)

      qlen = sqrt(dot_product(cas%qvector, cas%qvector))
      do ii = 1, cas%avg_order
        do jj = 1, 2 * cas%avg_order

          ! construct the q-vector
          phi   = acos(gaus_leg_points(ii))
          theta = M_PI * jj / cas%avg_order
          qvect(1) = qlen * cos(theta) * sin(phi)
          qvect(2) = qlen * sin(theta) * sin(phi)
          qvect(3) = qlen * cos(phi)

          ! matrix elements
          zx(:) = M_ZERO
          zf(:) = M_ZERO

          call states_get_state(st, mesh, 1, cas%pair(ia)%i, cas%pair(ia)%sigma, zf)
          call states_get_state(st, mesh, 1, cas%pair(ia)%a, cas%pair(ia)%sigma, psi_a)

          do ia = 1, cas%n_pairs
            forall(ip = 1:mesh%np)
              zf(ip) = exp(M_zI*dot_product(qvect(1:mesh%sb%dim), mesh%x(ip, 1:mesh%sb%dim)))*aimag(zf(ip))*psi_a(ip)
            end forall
            zx(ia) = zmf_integrate(mesh, zf)
          end do

          ! intensities
          do ia = 1, cas%n_pairs
            cas%qf_avg(ia) = cas%qf_avg(ia) + &
              gaus_leg_weights(ii)*abs(ztransition_matrix_element(cas, ia, zx))**2
          end do

        end do ! jj (thetas)
      end do ! ii (phis)

      ! normalize: for integral over sphere one would multiply by pi/N, but since
      !            we want the average, the integral must be divided by 4*pi
      forall(ia = 1:cas%n_pairs) cas%qf_avg(ia) = cas%qf_avg(ia) / (4*cas%avg_order)

      ! and finalize
      SAFE_DEALLOCATE_A(gaus_leg_points)
      SAFE_DEALLOCATE_A(gaus_leg_weights)

    end if ! averaging

    SAFE_DEALLOCATE_A(zf)
    SAFE_DEALLOCATE_A(psi_a)
    SAFE_DEALLOCATE_A(zx)

  end if

  SAFE_ALLOCATE(xx(1:cas%n_pairs))
  SAFE_ALLOCATE(deltav(1:mesh%np))

  do idir = 1, mesh%sb%dim
    deltav(1:mesh%np) = mesh%x(1:mesh%np, idir)
    ! let us get now the x vector.
    xx = X(ks_matrix_elements)(cas, st, mesh, deltav)
    ! And now we are able to get the transition matrix elements between many-electron states.
    do ia = 1, cas%n_pairs
      cas%X(tm)(ia, idir) = X(transition_matrix_element)(cas, ia, xx)
    end do
  end do
  SAFE_DEALLOCATE_A(xx)

  ! And the oscillator strengths.
  do ia = 1, cas%n_pairs
    cas%f(ia) = (M_TWO / mesh%sb%dim) * cas%w(ia) * sum( (abs(cas%X(tm)(ia, :)))**2 )
  end do

  call profiling_out(prof)

  POP_SUB(X(oscillator_strengths))
end subroutine X(oscillator_strengths)

! ---------------------------------------------------------
function X(ks_matrix_elements) (cas, st, mesh, dv) result(xx)
  type(casida_t), intent(in) :: cas
  type(states_t), intent(in) :: st
  type(mesh_t),   intent(in) :: mesh
  FLOAT,          intent(in) :: dv(:)
  R_TYPE :: xx(cas%n_pairs)

  R_TYPE, allocatable :: ff(:)
  R_TYPE, allocatable :: psii(:, :), psia(:, :)
  integer :: ip, ia, idim
  type(profile_t), save :: prof

  PUSH_SUB(X(ks_matrix_elements))
  call profiling_in(prof, 'CASIDA_KS')

  SAFE_ALLOCATE(ff(1:mesh%np))
  SAFE_ALLOCATE(psii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psia(1:mesh%np, 1:st%d%dim))

  do ia = 1, cas%n_pairs
    call states_get_state(st, mesh, cas%pair(ia)%i, cas%pair(ia)%sigma, psii)
    call states_get_state(st, mesh, cas%pair(ia)%a, cas%pair(ia)%sigma, psia)

    ! FIXME: parallelize in states
    ! use forall here
    do ip = 1, mesh%np
      ff(ip) = M_ZERO
      do idim = 1, st%d%dim
        ff(ip) = ff(ip) + dv(ip)*R_CONJ(psii(ip, idim))*psia(ip, idim)
      end do
    end do

    xx(ia) = X(mf_integrate)(mesh, ff)
  end do

  SAFE_DEALLOCATE_A(ff)

  call profiling_out(prof)
  POP_SUB(X(ks_matrix_elements))
end function X(ks_matrix_elements)

! ---------------------------------------------------------
R_TYPE function X(transition_matrix_element) (cas, ia, xx) result(zz)
  type(casida_t), intent(in) :: cas
  integer,        intent(in) :: ia
  R_TYPE,         intent(in) :: xx(:) !< these are KS matrix elements

  integer :: jb

  PUSH_SUB(X(transition_matrix_element))

  zz = R_TOTYPE(M_ZERO)
  if(cas%w(ia) > M_ZERO) then
    if(cas%type == CASIDA_EPS_DIFF) then
      zz = sqrt(TOFLOAT(cas%el_per_state)) * xx(ia)
    else if(cas%type == CASIDA_CASIDA) then
      do jb = 1, cas%n_pairs
        zz = zz + xx(jb) * (M_ONE/sqrt(cas%s(jb))) * cas%X(mat)(jb, ia)
      end do
      zz = (M_ONE/sqrt(cas%w(ia))) * zz
    else ! TAMM_DANCOFF, VARIATIONAL, PETERSILKA
      do jb = 1, cas%n_pairs
        zz = zz + xx(jb) * cas%X(mat)(jb, ia)
      end do
      if(cas%nik == 1) zz = sqrt(M_TWO) * zz
    endif
  end if

  POP_SUB(X(transition_matrix_element))
end function X(transition_matrix_element)

! ---------------------------------------------------------
subroutine X(transition_density) (cas, st, mesh, ia, n0I)
  type(casida_t), intent(in)  :: cas
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: mesh
  integer,        intent(in)  :: ia
  R_TYPE,         intent(out) :: n0I(:)

  integer :: ip, jb, idim, block_size, sp, ep
  R_TYPE, allocatable :: xx(:), psi(:, :, :, :)

  PUSH_SUB(X(transition_density))

  SAFE_ALLOCATE(xx(1:cas%n_pairs))

  block_size = 1000

  SAFE_ALLOCATE(psi(st%st_start:st%st_end, 1:st%d%dim, 1:block_size, st%d%kpt%start:st%d%kpt%end))

  ! We do this by blocks of points to avoid allocating an array of the
  ! size of the full states.
  do sp = 1, mesh%np, block_size
    ep = min(sp + block_size - 1, mesh%np)

    call states_get_points(st, sp, ep, psi)

    do ip = sp, ep
      do jb = 1, cas%n_pairs
        do idim = 1, st%d%dim
          xx(jb) = R_CONJ(psi(cas%pair(jb)%i, idim, ip - sp + 1, cas%pair(jb)%sigma))* &
            psi(cas%pair(jb)%a, idim, ip - sp + 1, cas%pair(jb)%sigma)
        end do
      end do
      n0I(ip) = X(transition_matrix_element)(cas, ia, xx)
    end do

  end do

  SAFE_DEALLOCATE_A(psi)
  SAFE_DEALLOCATE_A(xx)
  POP_SUB(X(transition_density))
end subroutine X(transition_density)

! ---------------------------------------------------------
subroutine X(get_transition_densities) (cas, sys)
  type(casida_t),    intent(in) :: cas
  type(system_t),    intent(in) :: sys

  integer :: ia, ierr
  character(len=5) :: intstr
  character(len=130) :: filename
  R_TYPE, allocatable :: n0I(:)
  type(unit_t) :: fn_unit

  if(.not. mpi_grp_is_root(mpi_world)) return

  PUSH_SUB(X(get_transition_densities))

  SAFE_ALLOCATE(n0I(1:sys%gr%mesh%np))
  n0I = M_ZERO
  fn_unit = units_out%length**(-sys%gr%sb%dim)

  do ia = 1, cas%n_pairs
    if(loct_isinstringlist(ia, cas%trandens)) then
      call X(transition_density) (cas, sys%st, sys%gr%mesh, ia, n0I)
      write(intstr,'(i5)') ia
      write(intstr,'(i1)') len(trim(adjustl(intstr)))
      write(filename,'(a,a,i'//trim(intstr)//')') trim(theory_name(cas)), '_rho_n0',ia
      call X(io_function_output)(sys%outp%how, CASIDA_DIR, trim(filename), &
        sys%gr%mesh, n0I, fn_unit, ierr, geo = sys%geo)
    end if
  end do

  SAFE_DEALLOCATE_A(n0I)
  POP_SUB(X(get_transition_densities))
end subroutine X(get_transition_densities)

! -----------------------------------------------------------------------------

subroutine X(casida_get_rho)(st, mesh, ii, ia, sigma, rho) 
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: mesh
  integer,        intent(in)  :: ii
  integer,        intent(in)  :: ia
  integer,        intent(in)  :: sigma
  R_TYPE,         intent(out) :: rho(:)

  R_TYPE, pointer :: psi_i(:), psi_a(:)
  integer :: ip, idim, iblock, ablock, ilin, alin
  type(profile_t), save :: prof

  PUSH_SUB(X(casida_get_rho))
  call profiling_in(prof, 'CASIDA_GET_RHO')

  ! For performance reasons we don`t use states_get_states, but we access the states directly

  iblock = st%group%iblock(ii, sigma)
  ablock = st%group%iblock(ia, sigma)

  ! FIXME: need to take into account spinor dimension here, not just 1
  idim = 1
  ilin = batch_inv_index(st%group%psib(iblock, sigma), (/ii, idim/))
  alin = batch_inv_index(st%group%psib(ablock, sigma), (/ia, idim/))

  ASSERT(.not. batch_is_packed(st%group%psib(iblock, sigma)))
  ASSERT(.not. batch_is_packed(st%group%psib(ablock, sigma)))

  psi_i => st%group%psib(iblock, sigma)%states_linear(ilin)%X(psi)
  psi_a => st%group%psib(ablock, sigma)%states_linear(alin)%X(psi)

  forall(ip = 1:mesh%np) rho(ip) = R_CONJ(psi_i(ip))*psi_a(ip)

  call profiling_out(prof)
  POP_SUB(X(casida_get_rho))
end subroutine X(casida_get_rho)

! -----------------------------------------------------------------------------

!> one-particle matrix elements of perturbation
subroutine X(casida_calc_lr_hmat1)(sys, hm, pert, hvar, lr_hmat1, is_saved, st_start, st_end, ik)
  type(system_t),      intent(in)    :: sys
  type(hamiltonian_t), intent(inout) :: hm
  type(pert_t),        intent(in)    :: pert
  FLOAT,               intent(in)    :: hvar(:,:,:)
  R_TYPE,              intent(out)   :: lr_hmat1(:,:,:)
  logical,             intent(in)    :: is_saved(:,:,:)
  integer,             intent(in)    :: st_start
  integer,             intent(in)    :: st_end
  integer,             intent(in)    :: ik

  integer :: ist, jst, ispin, idim
  R_TYPE, allocatable :: psi(:,:,:), pert_psi(:,:)

  PUSH_SUB(X(casida_calc_lr_hmat1))

  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim, st_start:st_end))
  SAFE_ALLOCATE(pert_psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))

  ! could use batches?
  ! FIXME: parallelize in states

  ispin = states_dim_get_spin_index(sys%st%d, ik)

  do ist = st_start, st_end
    call states_get_state(sys%st, sys%gr%mesh, ist, ik, psi(:, :, ist))
  enddo

  do ist = st_start, st_end
    if(all(is_saved(ist, ist:st_end, ik))) cycle
    call X(pert_apply)(pert, sys%gr, sys%geo, hm, ik, psi(:, :, ist), pert_psi(:, :))
    do idim = 1, sys%st%d%dim
      pert_psi(:, idim) = pert_psi(:, idim) + hvar(:, ispin, 1) * psi(:, idim, ist)
    enddo

    do jst = ist, st_end
      if(.not. is_saved(ist, jst, ik)) then
        lr_hmat1(jst, ist, ik) = X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, psi(:, :, jst), pert_psi(:, :))
        if(jst /= ist) lr_hmat1(ist, jst, ik) = R_CONJ(lr_hmat1(jst, ist, ik)) ! Hermiticity
      endif
    enddo
  enddo

  POP_SUB(X(casida_calc_lr_hmat1))
end subroutine X(casida_calc_lr_hmat1)

! -----------------------------------------------------------------------------

!> two-particle matrix elements of perturbation
subroutine X(casida_lr_hmat2)(cas, st, lr_hmat1, ik)
  type(casida_t), intent(inout) :: cas
  type(states_t), intent(in)    :: st
  R_TYPE,         intent(in)    :: lr_hmat1(:,:,:)
  integer,        intent(in)    :: ik

  integer :: ia, jb

  PUSH_SUB(X(casida_lr_hmat2))

  do ia = 1, cas%n_pairs
    do jb = ia, cas%n_pairs
      ! only matrix elements between degenerate states matter for degenerate perturbation theory
      if((cas%type == CASIDA_PETERSILKA .or. cas%type == CASIDA_EPS_DIFF) &
        .and. isnt_degenerate(cas, st, ia, jb)) cycle

      ! if occ states the same, apply unocc matrix elements
      if(cas%pair(ia)%i == cas%pair(jb)%i) then
        cas%X(lr_hmat2)(ia, jb) = cas%X(lr_hmat2)(ia, jb) + lr_hmat1(cas%pair(ia)%a, cas%pair(jb)%a, ik)
        if(ia /= jb) cas%X(lr_hmat2)(jb, ia) = cas%X(lr_hmat2)(jb, ia) + lr_hmat1(cas%pair(jb)%a, cas%pair(ia)%a, ik)
      endif

      ! if unocc states the same, apply occ matrix elements
      if(cas%pair(ia)%a == cas%pair(jb)%a) then
        cas%X(lr_hmat2)(ia, jb) = cas%X(lr_hmat2)(ia, jb) - lr_hmat1(cas%pair(ia)%i, cas%pair(jb)%i, ik)
        if(ia /= jb) cas%X(lr_hmat2)(jb, ia) = cas%X(lr_hmat2)(jb, ia) - lr_hmat1(cas%pair(jb)%i, cas%pair(ia)%i, ik)
      endif
    enddo
  enddo

  POP_SUB(X(casida_lr_hmat2))

end subroutine X(casida_lr_hmat2)

! -----------------------------------------------------------------------------

subroutine X(casida_get_matrix)(cas, hm, st, mesh, matrix, xc, restart_file, is_forces)
  type(casida_t),      intent(inout) :: cas
  type(hamiltonian_t), intent(in)    :: hm
  type(states_t),      intent(in)    :: st
  type(mesh_t),        intent(in)    :: mesh
  R_TYPE,              intent(out)   :: matrix(:,:)
  FLOAT,               intent(in)    :: xc(:,:,:)
  character(len=*),    intent(in)    :: restart_file
  logical,   optional, intent(in)    :: is_forces

  integer :: ia, jb, iunit, ia_iter, ia_length, jb_tmp
  integer :: maxcount, actual, counter
  R_TYPE :: mtxel_vh, mtxel_xc
  logical, allocatable :: is_saved(:, :), is_calcd(:, :)
  logical :: is_forces_
  type(casida_save_pot_t) :: saved_pot

  PUSH_SUB(X(casida_get_matrix))

  mtxel_vh = M_ZERO
  mtxel_xc = M_ZERO
  is_forces_ = optional_default(is_forces, .false.)

  ! load saved matrix elements
  SAFE_ALLOCATE(is_saved(1:cas%n_pairs, 1:cas%n_pairs))
  call load_saved(matrix, is_saved, restart_file)

  SAFE_ALLOCATE(is_calcd(1:cas%n_pairs, 1:cas%n_pairs))
  is_calcd = .true.
  ! purge saved non-degenerate offdiagonals, mark which are being calculated
  if(cas%type == CASIDA_PETERSILKA .and. mpi_grp_is_root(mpi_world)) then
    do ia = 1, cas%n_pairs
      do jb = ia, cas%n_pairs
        if(isnt_degenerate(cas, st, ia, jb)) then
          matrix(ia, jb) = M_ZERO
          matrix(jb, ia) = M_ZERO
          is_calcd(ia, jb) = .false.
          is_calcd(jb, ia) = .false.
        endif
      enddo
    enddo
  endif

  if(cas%type == CASIDA_PETERSILKA) then
    maxcount = cas%n_pairs
  else
    maxcount = ceiling((cas%n_pairs*(M_ONE + cas%n_pairs)/M_TWO)/cas%mpi_grp%size)
  endif
  counter = 0
  actual = 0
  if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, maxcount)

  ! only root retains the saved values
  if(.not. mpi_grp_is_root(mpi_world)) matrix = M_ZERO

  call X(casida_save_pot_init)(saved_pot, mesh)

  ! calculate the matrix elements of (v + fxc)
  do jb = 1, cas%n_pairs
    actual = actual + 1
    if(mod(actual, cas%mpi_grp%size) /= cas%mpi_grp%rank) cycle

    ! we only count diagonals for Petersilka
    if(cas%type == CASIDA_PETERSILKA) counter = counter + 1

    ! note: the ordering of jb, ia loops are crucial to minimize number of Poisson solves required.
    ia_length = (cas%n_pairs - 1) / 2
    if(mod(cas%n_pairs, 2) == 0) then ! even
      if(jb > cas%n_pairs / 2) then
        jb_tmp = cas%n_pairs - jb + 1
      else
        jb_tmp = jb
      endif
      ia_length = ia_length + mod(jb_tmp, 2)
    endif

    do ia_iter = jb, jb + ia_length

      ! make ia in range [1, cas%n_pairs]
      ia = mod(ia_iter, cas%n_pairs)
      if(ia == 0) ia = cas%n_pairs

      if(cas%type == CASIDA_PETERSILKA) then
        ! only calculate off-diagonals in degenerate subspace
        if(isnt_degenerate(cas, st, ia, jb)) cycle
      else
        counter = counter + 1
      endif

      ! if not loaded, then calculate matrix element
      if(.not. is_saved(ia, jb)) then
        if(is_forces_) then
          call X(K_term)(cas%pair(ia), cas%pair(jb), saved_pot, mtxel_xc = mtxel_xc)
        else
          call X(K_term)(cas%pair(ia), cas%pair(jb), saved_pot, mtxel_vh = mtxel_vh, mtxel_xc = mtxel_xc)
        endif
        matrix(ia, jb) = mtxel_vh + mtxel_xc
      end if
      if(jb /= ia) matrix(jb, ia) = R_CONJ(matrix(ia, jb))
    end do
    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(counter, maxcount)
  end do

  call X(casida_save_pot_end)(saved_pot)

  if(mpi_grp_is_root(mpi_world)) then
    call loct_progress_bar(maxcount, maxcount)
    ! complete progress bar
    write(stdout, '(1x)')
  endif

  ! sum all matrix elements
  if(cas%parallel_in_eh_pairs) then
    call comm_allreduce(cas%mpi_grp%comm, matrix)
  end if

  if(mpi_grp_is_root(mpi_world)) then
    ! output the restart file
    iunit = io_open(trim(restart_file), action='write', &
      position='append', is_tmp=.true.)

    do ia = 1, cas%n_pairs
      do jb = ia, cas%n_pairs
        if(.not. is_saved(ia, jb) .and. is_calcd(ia, jb)) &
          call X(write_K_term)(cas, matrix(ia, jb), iunit, ia, jb)
      enddo
    enddo

    call io_close(iunit)
  endif
  SAFE_DEALLOCATE_A(is_saved)

  POP_SUB(X(casida_get_matrix))

contains

  ! ---------------------------------------------------------
  !> calculates the matrix elements <i(p),a(p)|v|j(q),b(q)> and/or <i(p),a(p)|xc|j(q),b(q)>
  subroutine X(K_term)(pp, qq, saved, mtxel_vh, mtxel_xc)
    type(states_pair_t),               intent(in)    :: pp
    type(states_pair_t),               intent(in)    :: qq
    type(casida_save_pot_t),           intent(inout) :: saved
    R_TYPE,                  optional, intent(out)   :: mtxel_vh
    R_TYPE,                  optional, intent(out)   :: mtxel_xc

    integer :: pi, qi, sigma, pa, qa, mu
    R_TYPE, allocatable :: rho_i(:), rho_j(:), integrand(:)
    FLOAT :: coeff_vh
    type(profile_t), save :: prof

    PUSH_SUB(X(casida_get_matrix).X(K_term))
    call profiling_in(prof, 'CASIDA_K')
    
    if(cas%herm_conj) then
      pi = qq%i
      pa = qq%a
      sigma = qq%sigma

      qi = pp%i
      qa = pp%a
      mu = pp%sigma
    else
      pi = pp%i
      pa = pp%a
      sigma = pp%sigma

      qi = qq%i
      qa = qq%a
      mu = qq%sigma
    endif

    SAFE_ALLOCATE(rho_i(1:mesh%np))
    SAFE_ALLOCATE(rho_j(1:mesh%np))
    SAFE_ALLOCATE(integrand(1:mesh%np))

    call X(casida_get_rho)(st, mesh, pa, pi, sigma, rho_i)
    call X(casida_get_rho)(st, mesh, qi, qa, mu,    rho_j)

    !  first the Hartree part
    if(present(mtxel_vh)) then
      coeff_vh = - cas%kernel_lrc_alpha / (M_FOUR * M_PI)
      if(.not. cas%triplet) coeff_vh = coeff_vh + M_ONE
      if(abs(coeff_vh) > M_EPSILON) then
        if(qi /= saved%qi  .or.   qa /= saved%qa .or.  mu /= saved%mu) then
          saved%X(pot)(1:mesh%np) = M_ZERO
          if(hm%theory_level /= INDEPENDENT_PARTICLES) call X(poisson_solve)(psolver, saved%X(pot), rho_j, all_nodes=.false.)

          saved%qi = qi
          saved%qa = qa
          saved%mu = mu
      else
        endif
        ! value of pot is retained between calls
        mtxel_vh = coeff_vh * X(mf_dotp)(mesh, rho_i(:), saved%X(pot)(:))

      else
        mtxel_vh = M_ZERO
      endif
    end if

    if(present(mtxel_xc)) then
      integrand(1:mesh%np) = rho_i(1:mesh%np)*rho_j(1:mesh%np)*xc(1:mesh%np, sigma, mu)
      mtxel_xc = X(mf_integrate)(mesh, integrand)
    endif

    if(cas%herm_conj) then
      if(present(mtxel_vh)) mtxel_vh = R_CONJ(mtxel_vh)
      if(present(mtxel_xc)) mtxel_vh = R_CONJ(mtxel_xc)
    endif

    SAFE_DEALLOCATE_A(rho_i)
    SAFE_DEALLOCATE_A(rho_j)
    SAFE_DEALLOCATE_A(integrand)

    call profiling_out(prof)
    POP_SUB(X(casida_get_matrix).X(K_term))
  end subroutine X(K_term)

  ! ---------------------------------------------------------
  subroutine load_saved(matrix, is_saved, restart_file)
    R_TYPE,           intent(out) :: matrix(:,:)
    logical,          intent(out) :: is_saved(:,:)
    character(len=*), intent(in)  :: restart_file

    integer :: iunit, err
    integer :: ia, jb, ii, aa, ik, jj, bb, jk
    R_TYPE  :: val

    PUSH_SUB(X(casida_get_matrix).load_saved)

    is_saved = .false.
    matrix = M_ZERO

    ! if fromScratch, we already deleted the restart files
    if(mpi_grp_is_root(mpi_world)) then
      iunit = io_open(trim(restart_file), action='read', &
        status='old', die=.false., is_tmp=.true.)

      if( iunit > 0) then
        do
          read(iunit, fmt=*, iostat=err) ii, aa, ik, jj, bb, jk, val
          if(err /= 0) exit

          if(ii < 1 .or. aa < 1 .or. ik < 1) then
            message(1) = "Illegal indices in '" // trim(restart_file) // "': working from scratch."
            call messages_warning(1)
            call loct_rm(trim(restart_file))
            ! if file is corrupt, do not trust anything that was read
            is_saved = .false.
            exit
          endif

          ia = cas%index(ii, aa, ik)
          jb = cas%index(jj, bb, jk)

          if(ia > 0 .and. jb > 0) then
            matrix(ia, jb) = val
            is_saved(ia, jb) = .true.
            matrix(jb, ia) = R_CONJ(val)
            is_saved(jb, ia) = .true.
          endif
        end do

        call io_close(iunit)
      else if(.not. cas%fromScratch) then
        message(1) = "Could not find restart file '" // trim(restart_file) // "'. Starting from scratch."
        call messages_warning(1)
      endif
    endif

    ! if no file found, root has no new information to offer the others
#ifdef HAVE_MPI
    call MPI_Bcast(is_saved(1, 1), cas%n_pairs**2, MPI_LOGICAL, 0, mpi_world, mpi_err)
    ! No need to bcast these, since they will be obtained from a reduction
    !      call MPI_Bcast(cas%X(mat)(1, 1), cas%n_pairs**2, R_MPITYPE,   0, mpi_world, mpi_err)
#endif

    POP_SUB(X(casida_get_matrix).load_saved)
  end subroutine load_saved

end subroutine X(casida_get_matrix)

! ---------------------------------------------------------
!> write matrix element to casida_restart file
subroutine X(write_K_term)(cas, mat_val, iunit, ia, jb)
  type(casida_t), intent(in) :: cas
  R_TYPE,         intent(in) :: mat_val
  integer,        intent(in) :: iunit
  integer,        intent(in) :: ia
  integer,        intent(in) :: jb

  PUSH_SUB(X(write_K_term))
  
  write(iunit,*) cas%pair(ia)%i, cas%pair(ia)%a, cas%pair(ia)%sigma, &
    cas%pair(jb)%i, cas%pair(jb)%a, cas%pair(jb)%sigma, mat_val
  
  POP_SUB(X(write_K_term))
end subroutine X(write_K_term)

! ---------------------------------------------------------
subroutine X(casida_forces)(cas, sys, mesh, st, hm)
  type(casida_t), intent(inout) :: cas
  type(system_t), intent(inout) :: sys
  type(mesh_t), intent(in) :: mesh
  type(states_t), intent(inout) :: st
  type(hamiltonian_t), intent(inout) :: hm
  
  integer :: ip, iatom, idir, is1, is2, ierr, ik, ia
  FLOAT, allocatable :: dl_rho(:,:), kxc(:,:,:,:)
  FLOAT, target, allocatable :: lr_fxc(:,:,:)
  R_TYPE, allocatable :: lr_hmat1(:,:,:)
  FLOAT :: factor = CNST(1e6) ! FIXME: allow user to set
  character(len=100) :: restart_filename

  PUSH_SUB(X(casida_forces))
  
  if(cas%type == CASIDA_CASIDA) then
    message(1) = "Forces for Casida theory level not implemented"
    call messages_warning(1)
    POP_SUB(X(casida_forces))
    return
  endif
  
  if(cas%type /= CASIDA_EPS_DIFF) then
    SAFE_ALLOCATE(kxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin, 1:st%d%nspin))
    kxc = M_ZERO
    ! not spin polarized so far
    call xc_get_kxc(sys%ks%xc, mesh, cas%rho, st%d%ispin, kxc(:, :, :, :))
  endif
  
  message(1) = "Reading vib_modes density for calculating excited-state forces."
  call messages_info(1)
  
  SAFE_ALLOCATE(dl_rho(1:mesh%np, 1:st%d%nspin))
  if (cas%type /= CASIDA_EPS_DIFF) then
    SAFE_ALLOCATE(lr_fxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
  endif
  
  if(cas%type == CASIDA_EPS_DIFF) then
    cas%X(mat_save) = M_ZERO
    do ia = 1, cas%n_pairs
      cas%X(mat_save)(ia, ia) = cas%w(ia)
    enddo
  endif
  
  SAFE_ALLOCATE(lr_hmat1(cas%nst, cas%nst, cas%nik))
  SAFE_ALLOCATE(cas%X(lr_hmat2)(cas%n_pairs, cas%n_pairs))
  SAFE_ALLOCATE(cas%X(mat2)(cas%n_pairs, cas%n_pairs))
  SAFE_ALLOCATE(cas%X(w2)(cas%n_pairs))

  do iatom = 1, sys%geo%natoms
    do idir = 1, mesh%sb%dim
      
      call drestart_read_lr_rho(dl_rho, sys%gr, st%d%nspin, &
        VIB_MODES_DIR, phn_rho_tag(iatom, idir), ierr)
      
      if(ierr /= 0) then
        message(1) = "Could not load vib_modes density; previous vib_modes calculation required."
        call messages_fatal(1)
      end if

      call X(casida_get_lr_hmat1)(cas, sys, hm, iatom, idir, dl_rho, lr_hmat1)
      
      cas%X(lr_hmat2) = M_ZERO
      ! use them to make two-particle matrix elements (as for eigenvalues)
      do ik = 1, cas%nik
        call X(casida_lr_hmat2)(cas, st, lr_hmat1, ik)
      enddo
      
      if (cas%type /= CASIDA_EPS_DIFF .and. cas%calc_forces_kernel) then
        forall(ip = 1:mesh%np, is1 = 1:st%d%nspin, is2 = 1:st%d%nspin)
          lr_fxc(ip, is1, is2) = sum(kxc(ip, is1, is2, :) * dl_rho(ip, :))
        end forall
        
        write(restart_filename,'(a,a,i6.6,a,i1)') trim(cas%restart_dir), '/lr_kernel_', iatom, '_', idir
        if(cas%triplet) restart_filename = trim(restart_filename)//'_triplet'
        
        call X(casida_get_matrix)(cas, hm, st, mesh, cas%X(mat2), lr_fxc, restart_filename, is_forces = .true.)
        cas%X(mat2) = cas%X(mat2) * casida_matrix_factor(cas, sys)
      else
        cas%X(mat2) = M_ZERO
      endif
      
      cas%X(mat2) = cas%X(mat_save) * factor + cas%X(lr_hmat2) + cas%X(mat2)
      call lalg_eigensolve(cas%n_pairs, cas%X(mat2), cas%X(w2))
      do ia = 1, cas%n_pairs
        cas%forces(iatom, idir, cas%ind(ia)) = factor * cas%w(cas%ind(ia)) - R_REAL(cas%X(w2)(ia))
      enddo
    enddo
  enddo
  
  if(cas%type /= CASIDA_EPS_DIFF) then
    SAFE_DEALLOCATE_A(kxc)
    SAFE_DEALLOCATE_A(lr_fxc)
  endif
  SAFE_DEALLOCATE_A(dl_rho)

  SAFE_DEALLOCATE_A(lr_hmat1)
  SAFE_DEALLOCATE_P(cas%X(mat2))
  SAFE_DEALLOCATE_P(cas%X(w2))
  
  if(cas%calc_forces_scf) then
    call forces_calculate(sys%gr, sys%geo, hm, st)
    do ia = 1, cas%n_pairs
      do iatom = 1, sys%geo%natoms
        do idir = 1, sys%gr%sb%dim
          cas%forces(iatom, idir, ia) = cas%forces(iatom, idir, ia) + sys%geo%atom(iatom)%f(idir)
        enddo
      enddo
    enddo
  endif
  
  POP_SUB(X(casida_forces))

end subroutine X(casida_forces)

! ---------------------------------------------------------
subroutine X(casida_get_lr_hmat1)(cas, sys, hm, iatom, idir, dl_rho, lr_hmat1)
  type(casida_t),      intent(in)     :: cas
  type(system_t),      intent(inout)  :: sys
  type(hamiltonian_t), intent(inout)  :: hm
  integer,             intent(in)     :: iatom
  integer,             intent(in)     :: idir
  FLOAT,               intent(in)     :: dl_rho(:,:)
  R_TYPE,              intent(out)    :: lr_hmat1(:,:,:)

  FLOAT, allocatable :: hvar(:,:,:)
  integer :: ik, ist, jst, iunit, err, ii, aa, num_saved
  type(pert_t) :: ionic_pert
  character(len=100) :: restart_filename
  R_TYPE :: val
  logical :: all_done
  logical, allocatable :: is_saved(:,:,:)

  PUSH_SUB(X(casida_get_lr_hmat1))

  lr_hmat1 = M_ZERO
  SAFE_ALLOCATE(is_saved(cas%nst, cas%nst, cas%nik))
  is_saved = .false.
  num_saved = 0

  ! if fromScratch, we already deleted the restart files
  if(mpi_grp_is_root(mpi_world)) then
    write(restart_filename,'(a,a,i6.6,a,i1)') trim(cas%restart_dir), '/lr_hmat1_', iatom, '_', idir
    iunit = io_open(restart_filename, action = 'read', status = 'old', die = .false., is_tmp = .true.)

    if(iunit > 0) then
      do
        read(iunit, fmt=*, iostat=err) ii, aa, ik, val
        if(err /= 0) exit

        if(ii < 1 .or. aa < 1 .or. ik < 1) then
          message(1) = "Illegal indices in '" // trim(restart_filename) // "': working from scratch."
          call messages_warning(1)
          call loct_rm(trim(restart_filename))
          ! if file is corrupt, do not trust anything that was read
          is_saved = .false.
          exit
        endif

        ! FIXME: what about elements which are always zero?
        if(ii <= cas%nst .and. aa <= cas%nst .and. ik <= cas%nik) then
          lr_hmat1(ii, aa, ik) = val
          is_saved(ii, aa, ik) = .true.
          lr_hmat1(aa, ii, ik) = R_CONJ(val)
          is_saved(aa, ii, ik) = .true.
          num_saved = num_saved + 1
        endif
      enddo
      write(6,'(a,i8,a,a)') 'Read ', num_saved, ' saved elements from ', trim(restart_filename)

      call io_close(iunit)
    else if(.not. cas%fromScratch) then
      message(1) = "Could not find restart file '" // trim(restart_filename) // "'. Starting from scratch."
      call messages_warning(1)
    endif
  endif

#ifdef HAVE_MPI
    call MPI_Bcast(is_saved(1, 1, 1), cas%nst**2, MPI_LOGICAL, 0, mpi_world, mpi_err)
#endif

  all_done = .true.
  do ik = 1, cas%nik
    all_done = all_done .and. all(is_saved(1:cas%n_occ(ik), 1:cas%n_occ(ik), ik)) &
      .and. all(is_saved(cas%n_occ(ik) + 1:cas%nst, cas%n_occ(ik) + 1:cas%nst, ik))
  enddo

  if(all_done) then
    SAFE_DEALLOCATE_A(is_saved)
    POP_SUB(X(casida_get_lr_hmat1))
    return
  endif

  call pert_init(ionic_pert, PERTURBATION_IONIC, sys%gr, sys%geo)
  call pert_setup_atom(ionic_pert, iatom)
  call pert_setup_dir(ionic_pert, idir)

  SAFE_ALLOCATE(hvar(1:sys%gr%mesh%np, 1:sys%st%d%nspin, 1:1))
  call dcalc_hvar(.true., sys, dl_rho, 1, hvar, fxc = cas%fxc)

  ! FIXME: do this only for states called for in CasidaKohnShamStates
  do ik = 1, cas%nik
    ! occ-occ matrix elements
    call X(casida_calc_lr_hmat1)(sys, hm, ionic_pert, hvar, lr_hmat1, is_saved, 1, cas%n_occ(ik), ik)
    ! unocc-unocc matrix elements
    call X(casida_calc_lr_hmat1)(sys, hm, ionic_pert, hvar, lr_hmat1, is_saved, cas%n_occ(ik) + 1, cas%nst, ik)
  enddo

  SAFE_DEALLOCATE_A(hvar)
  call pert_end(ionic_pert)
  
  if(mpi_grp_is_root(mpi_world)) then
    iunit = io_open(restart_filename, action = 'write', position = 'append', is_tmp = .true.)
    do ik = 1, cas%nik
      do ist = 1, cas%nst
        do jst = ist, cas%nst
          if(.not. is_saved(ist, jst, ik)) write(iunit, *) ist, jst, ik, lr_hmat1(ist, jst, ik)
        enddo
      enddo
    enddo
    call io_close(iunit)
  endif

  SAFE_DEALLOCATE_A(is_saved)

  POP_SUB(X(casida_get_lr_hmat1))

end subroutine X(casida_get_lr_hmat1)

! -----------------------------------------------------

subroutine X(casida_save_pot_init)(this, mesh)
  type(casida_save_pot_t), intent(out)   :: this
  type(mesh_t),            intent(in)    :: mesh
  
  SAFE_ALLOCATE(this%X(pot)(1:mesh%np))
  this%qi = -1
  this%qa = -1
  this%mu = -1

end subroutine X(casida_save_pot_init)
    
! -----------------------------------------------------

subroutine X(casida_save_pot_end)(this)
  type(casida_save_pot_t), intent(inout)   :: this

  SAFE_DEALLOCATE_P(this%X(pot))

end subroutine X(casida_save_pot_end)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
