
!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2011-2013 D. Strubbe
!! Copyright (C) 2017-2018 J. Flick, S. Ohlmann
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

subroutine X(oscillator_strengths)(cas, mesh, st)
  type(casida_t),   intent(inout) :: cas
  type(mesh_t),        intent(in) :: mesh
  type(states_elec_t), intent(in) :: st

  FLOAT, allocatable :: deltav(:)
  R_TYPE, allocatable :: xx(:), psi_a(:)
  CMPLX, allocatable :: zf(:), zx(:)
  FLOAT, allocatable :: gaus_leg_points(:), gaus_leg_weights(:)
  integer :: ii, jj, ia, ip, idir
  FLOAT :: theta, phi, qlen
  FLOAT :: qvect(MAX_DIM)
  type(profile_t), save :: prof

  PUSH_SUB(X(oscillator_strengths))

  call profiling_in(prof, TOSTRING(X(CASIDA_OSCILLATOR_STRENGTHS)))

  if(cas%qcalc) then
    SAFE_ALLOCATE(zf(1:mesh%np))
    SAFE_ALLOCATE(psi_a(1:mesh%np))
    SAFE_ALLOCATE(zx(1:cas%n_pairs))

    ! matrix element
    do ia = 1, cas%n_pairs

      call states_elec_get_state(st, mesh, 1, cas%pair(ia)%i, cas%pair(ia)%kk, zf)
      call states_elec_get_state(st, mesh, 1, cas%pair(ia)%a, cas%pair(ia)%kk, psi_a)

      do ip = 1, mesh%np
        zf(ip) = exp(M_zI*dot_product(cas%qvector(1:mesh%sb%dim), mesh%x(ip, 1:mesh%sb%dim)))*aimag(zf(ip))*psi_a(ip)
      end do

      zx(ia) = zmf_integrate(mesh, zf, reduce = .false.)
    end do

    if(mesh%parallel_in_domains) then
      call mesh%allreduce(zx)
    end if

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

          call states_elec_get_state(st, mesh, 1, cas%pair(ia)%i, cas%pair(ia)%kk, zf)
          call states_elec_get_state(st, mesh, 1, cas%pair(ia)%a, cas%pair(ia)%kk, psi_a)

          do ia = 1, cas%n_pairs
            do ip = 1, mesh%np
              zf(ip) = exp(M_zI*dot_product(qvect(1:mesh%sb%dim), mesh%x(ip, 1:mesh%sb%dim)))*aimag(zf(ip))*psi_a(ip)
            end do
            zx(ia) = zmf_integrate(mesh, zf, reduce = .false.)
          end do

          if(mesh%parallel_in_domains) then
            call mesh%allreduce(zx)
          end if
          

          ! intensities
          do ia = 1, cas%n_pairs
            cas%qf_avg(ia) = cas%qf_avg(ia) + &
              gaus_leg_weights(ii)*abs(ztransition_matrix_element(cas, ia, zx))**2
          end do

        end do ! jj (thetas)
      end do ! ii (phis)

      ! normalize: for integral over sphere one would multiply by pi/N, but since
      !            we want the average, the integral must be divided by 4*pi
      do ia = 1, cas%n_pairs
        cas%qf_avg(ia) = cas%qf_avg(ia) / (4*cas%avg_order)
      end do

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
    do ia = 1, cas%n
      cas%X(tm)(ia, idir) = X(transition_matrix_element)(cas, ia, xx)
    end do
  end do
  SAFE_DEALLOCATE_A(xx)
  SAFE_DEALLOCATE_A(deltav)

  ! And the oscillator strengths.
  do ia = 1, cas%n
    cas%f(ia) = (M_TWO / mesh%sb%dim) * cas%w(ia) * sum( (abs(cas%X(tm)(ia, :)))**2 )
  end do

  call profiling_out(prof)

  POP_SUB(X(oscillator_strengths))
end subroutine X(oscillator_strengths)

! ---------------------------------------------------------
function X(ks_matrix_elements) (cas, st, mesh, dv) result(xx)
  type(casida_t),      intent(in) :: cas
  type(states_elec_t), intent(in) :: st
  type(mesh_t),        intent(in) :: mesh
  FLOAT,               intent(in) :: dv(:)
  R_TYPE                          :: xx(cas%n_pairs)

  R_TYPE, allocatable :: ff(:)
  R_TYPE, allocatable :: psii(:, :), psia(:, :)
  integer :: ip, ia
  type(profile_t), save :: prof

  PUSH_SUB(X(ks_matrix_elements))
  call profiling_in(prof, TOSTRING(X(CASIDA_KS)))

  SAFE_ALLOCATE(ff(1:mesh%np))
  SAFE_ALLOCATE(psii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psia(1:mesh%np, 1:st%d%dim))

  do ia = 1, cas%n_pairs
    call states_elec_get_state(st, mesh, cas%pair(ia)%i, cas%pair(ia)%kk, psii)
    call states_elec_get_state(st, mesh, cas%pair(ia)%a, cas%pair(ia)%kk, psia)

    do ip = 1, mesh%np
      ff(ip) = dv(ip)*sum(R_CONJ(psii(ip, 1:st%d%dim))*psia(ip, 1:st%d%dim))
    end do

    xx(ia) = X(mf_integrate)(mesh, ff, reduce = .false.)
  end do

  if(mesh%parallel_in_domains) then
    call mesh%allreduce(xx)
  end if

  SAFE_DEALLOCATE_A(ff)

  call profiling_out(prof)
  POP_SUB(X(ks_matrix_elements))
end function X(ks_matrix_elements)

! ---------------------------------------------------------
!> Casida: \vec{d}_k = \sum_{cv} \vec{d}_{cv} x_{cv} \sqrt{\frac{\epsilon_c - \epsilon_v}{\omega_k}}
!! others: \vec{d}_k = \sum_{cv} \vec{d}_{cv} x_{cv}
R_TYPE function X(transition_matrix_element) (cas, ia, xx) result(zz)
  type(casida_t), intent(in) :: cas
  integer,        intent(in) :: ia
  R_TYPE,         intent(in) :: xx(:) !< these are KS matrix elements
  integer :: jb, jb_local, ia_local
  logical :: on_this_processor

  PUSH_SUB(X(transition_matrix_element))

  zz = R_TOTYPE(M_ZERO)
  if(cas%w(ia) > M_ZERO) then
    if(cas%type == CASIDA_EPS_DIFF) then
      zz = xx(ia)
    else if(cas%type == CASIDA_CASIDA) then
      do jb = 1, cas%n_pairs
        call local_indices(cas, ia, jb, on_this_processor, ia_local, jb_local)
        if(on_this_processor) then
          zz = zz + xx(jb) * cas%X(mat)(jb_local, ia_local) / sqrt(cas%s(jb))
        end if
      end do
      zz = X(allreduce_sum)(cas, zz)
      zz = zz / sqrt(cas%w(ia))
    else ! TAMM_DANCOFF, VARIATIONAL, PETERSILKA
      do jb = 1, cas%n_pairs
        call local_indices(cas, ia, jb, on_this_processor, ia_local, jb_local)
        if(on_this_processor) then
          zz = zz + xx(jb) * cas%X(mat)(jb_local, ia_local)
        end if
      end do
      zz = X(allreduce_sum)(cas, zz)
    end if
    zz = sqrt(TOFLOAT(cas%el_per_state)) * zz
  end if

  POP_SUB(X(transition_matrix_element))
end function X(transition_matrix_element)

! ---------------------------------------------------------
subroutine X(transition_density) (cas, st, mesh, ia, n0I)
  type(casida_t),      intent(in)  :: cas
  type(states_elec_t), intent(in)  :: st
  type(mesh_t),        intent(in)  :: mesh
  integer,             intent(in)  :: ia
  R_TYPE,              intent(out) :: n0I(:)

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

    call states_elec_get_points(st, sp, ep, psi)

    do ip = sp, ep
      do jb = 1, cas%n_pairs
        do idim = 1, st%d%dim
          xx(jb) = R_CONJ(psi(cas%pair(jb)%i, idim, ip - sp + 1, cas%pair(jb)%kk))* &
            psi(cas%pair(jb)%a, idim, ip - sp + 1, cas%pair(jb)%kk)
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
  type(casida_t),      intent(in) :: cas
  type(electrons_t),   intent(in) :: sys

  integer :: ia, ierr
  character(len=5) :: intstr
  character(len=130) :: filename
  R_TYPE, allocatable :: n0I(:)
  type(unit_t) :: fn_unit

  PUSH_SUB(X(get_transition_densities))

  SAFE_ALLOCATE(n0I(1:sys%gr%mesh%np))
  n0I = M_ZERO
  fn_unit = units_out%length**(-sys%space%dim)

  do ia = 1, cas%n_pairs
    if(loct_isinstringlist(ia, cas%trandens)) then
      call X(transition_density) (cas, sys%st, sys%gr%mesh, ia, n0I)
      write(intstr,'(i5)') ia
      write(intstr,'(i1)') len(trim(adjustl(intstr)))
      write(filename,'(a,a,i'//trim(intstr)//')') trim(theory_name(cas)), '_rho_n0',ia
      call X(io_function_output)(sys%outp%how(0), CASIDA_DIR, trim(filename), &
        sys%namespace, sys%space, sys%gr%mesh, n0I, fn_unit, ierr, ions = sys%ions)
    end if
  end do

  SAFE_DEALLOCATE_A(n0I)
  POP_SUB(X(get_transition_densities))
end subroutine X(get_transition_densities)

! -----------------------------------------------------------------------------

subroutine X(casida_get_rho)(st, mesh, ii, ia, kk, rho) 
  type(states_elec_t), intent(in)  :: st
  type(mesh_t),        intent(in)  :: mesh
  integer,             intent(in)  :: ii
  integer,             intent(in)  :: ia
  integer,             intent(in)  :: kk
  R_TYPE,              intent(out) :: rho(:)

  integer :: ip, idim, iblock, ablock, ilin, alin
  type(profile_t), save :: prof

  PUSH_SUB(X(casida_get_rho))
  call profiling_in(prof, TOSTRING(X(CASIDA_GET_RHO)))

  ! For performance reasons we don`t use states_elec_get_states, but we access the states directly

  iblock = st%group%iblock(ii, kk)
  ablock = st%group%iblock(ia, kk)

  ! FIXME: need to take into account spinor dimension here, not just 1
  idim = 1
  ilin = st%group%psib(iblock, kk)%inv_index((/ii, idim/))
  alin = st%group%psib(ablock, kk)%inv_index((/ia, idim/))

  ASSERT(.not. st%group%psib(iblock, kk)%is_packed())
  ASSERT(.not. st%group%psib(ablock, kk)%is_packed())

  do ip = 1, mesh%np
    rho(ip) = R_CONJ(st%group%psib(iblock, kk)%X(ff_linear)(ip, ilin))*st%group%psib(ablock, kk)%X(ff_linear)(ip, alin)
  end do

  call profiling_out(prof)
  POP_SUB(X(casida_get_rho))
end subroutine X(casida_get_rho)

! -----------------------------------------------------------------------------

!> one-particle matrix elements of perturbation
subroutine X(casida_calc_lr_hmat1)(sys, pert, hvar, lr_hmat1, is_saved, st_start, st_end, ik)
  type(electrons_t),   intent(inout) :: sys
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

  ispin = sys%st%d%get_spin_index(ik)

  do ist = st_start, st_end
    call states_elec_get_state(sys%st, sys%gr%mesh, ist, ik, psi(:, :, ist))
  end do

  do ist = st_start, st_end
    if(all(is_saved(ist, ist:st_end, ik))) cycle
    call X(pert_apply)(pert, sys%namespace, sys%gr, sys%ions, sys%hm, ik, psi(:, :, ist), pert_psi(:, :))
    do idim = 1, sys%st%d%dim
      pert_psi(:, idim) = pert_psi(:, idim) + hvar(:, ispin, 1) * psi(:, idim, ist)
    end do

    do jst = ist, st_end
      if(.not. is_saved(ist, jst, ik)) then
        lr_hmat1(jst, ist, ik) = X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, psi(:, :, jst), pert_psi(:, :))
        if(jst /= ist) lr_hmat1(ist, jst, ik) = R_CONJ(lr_hmat1(jst, ist, ik)) ! Hermiticity
      end if
    end do
  end do

  POP_SUB(X(casida_calc_lr_hmat1))
end subroutine X(casida_calc_lr_hmat1)

! -----------------------------------------------------------------------------

!> two-particle matrix elements of perturbation
subroutine X(casida_lr_hmat2)(cas, st, lr_hmat1, ik)
  type(casida_t),      intent(inout) :: cas
  type(states_elec_t), intent(in)    :: st
  R_TYPE,              intent(in)    :: lr_hmat1(:,:,:)
  integer,             intent(in)    :: ik

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
      end if

      ! if unocc states the same, apply occ matrix elements
      if(cas%pair(ia)%a == cas%pair(jb)%a) then
        cas%X(lr_hmat2)(ia, jb) = cas%X(lr_hmat2)(ia, jb) - lr_hmat1(cas%pair(ia)%i, cas%pair(jb)%i, ik)
        if(ia /= jb) cas%X(lr_hmat2)(jb, ia) = cas%X(lr_hmat2)(jb, ia) - lr_hmat1(cas%pair(jb)%i, cas%pair(ia)%i, ik)
      end if
    end do
  end do

  POP_SUB(X(casida_lr_hmat2))

end subroutine X(casida_lr_hmat2)

! -----------------------------------------------------------------------------

subroutine X(casida_get_matrix)(cas, hm, st, ks, mesh, matrix, xc, restart_file, is_forces)
  type(casida_t),           intent(inout) :: cas
  type(hamiltonian_elec_t), intent(in)    :: hm
  type(states_elec_t),      intent(in)    :: st
  type(v_ks_t),             intent(in)    :: ks
  type(mesh_t),             intent(in)    :: mesh
  R_TYPE,                   intent(out)   :: matrix(:,:)
  FLOAT,                    intent(in)    :: xc(:,:,:)
  character(len=*),         intent(in)    :: restart_file
  logical, optional,        intent(in)    :: is_forces

  integer :: ia, jb, iunit, ia_length, ii
  integer :: ia_local, jb_local
  integer :: maxcount, actual, counter
  R_TYPE :: mtxel_vh, mtxel_xc, mtxel_vm, save_value
  R_TYPE :: rhobufferx, rhobuffery, rhobufferz
  logical, allocatable :: is_saved(:, :), is_calcd(:, :)
  logical :: is_forces_, write_value, on_this_processor
  R_TYPE, allocatable :: xx(:)
  R_TYPE, allocatable :: X(pot)(:)
  R_TYPE, allocatable :: buffer(:)
  R_TYPE, allocatable :: bufferx(:)
  R_TYPE, allocatable :: buffery(:)
  R_TYPE, allocatable :: bufferz(:)
#ifdef HAVE_SCALAPACK
  integer :: mpi_status(mpi_status_size), src
  R_TYPE, allocatable :: buffer_transpose(:,:)
#endif
  integer, allocatable :: rank_of_element(:, :)

  type(profile_t), save :: prof

  R_TYPE, allocatable :: rho_i(:), rho_j(:), integrand_xc(:,:)
  FLOAT :: coeff_vh
  integer :: iimode
  FLOAT, allocatable  :: deltav(:)

  PUSH_SUB(X(casida_get_matrix))
  call profiling_in(prof, 'CASIDA_GET_MATRIX')

  SAFE_ALLOCATE(xx(1:cas%n_pairs))

  mtxel_vh = M_ZERO
  mtxel_xc = M_ZERO
  mtxel_vm = M_ZERO
  is_forces_ = optional_default(is_forces, .false.)

  ! load saved matrix elements
  SAFE_ALLOCATE(is_saved(1:cas%n_pairs, 1:cas%n_pairs))
  if(.not. cas%has_photons) then
    call load_saved(matrix, is_saved, restart_file)
  else
    is_saved = .false.
  end if

  SAFE_ALLOCATE(is_calcd(1:cas%n_pairs, 1:cas%n_pairs))
  is_calcd = .true.
  ! purge saved non-degenerate offdiagonals, mark which are being calculated
  if(cas%type == CASIDA_PETERSILKA) then
    do jb_local = 1, cas%nb_rows
      jb = get_global_row(cas, jb_local)
      if(jb > cas%n_pairs) cycle
      do ia_local = 1, cas%nb_cols
        ia = get_global_col(cas, ia_local)
        if(ia > cas%n_pairs) cycle
        if(isnt_degenerate(cas, st, ia, jb)) then
          matrix(jb_local, ia_local) = M_ZERO
          is_calcd(ia, jb) = .false.
          if(.not. cas%distributed_matrix) then
            matrix(ia_local, jb_local) = M_ZERO
            is_calcd(jb, ia) = .false.
          end if
        end if
      end do
    end do
  end if

  if(cas%distributed_matrix) then
    SAFE_ALLOCATE(rank_of_element(1:cas%n_pairs, 1:cas%n_pairs))
    rank_of_element = 0
  end if

  if(cas%type == CASIDA_PETERSILKA) then
    maxcount = cas%n_pairs
  else
    maxcount = ceiling((cas%n_pairs*(M_ONE + cas%n_pairs)/M_TWO + cas%n_pairs*cas%pt_nmodes) &
      * cas%nb_rows/cas%n * cas%nb_cols/cas%n)
  end if
  counter = 0
  actual = 0
  if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(-1, maxcount)

  if(cas%herm_conj) then
    message(1) = "herm_conj currently not implemented due to performance reasons."
    ! look at older versions in git for reimplementing this...
    call messages_fatal(1)
  end if


  if(.not. cas%distributed_matrix) then
    ! only root retains the saved values
    if(.not. mpi_grp_is_root(mpi_world)) matrix = M_ZERO
  end if

  SAFE_ALLOCATE(rho_i(1:mesh%np))
  SAFE_ALLOCATE(rho_j(1:mesh%np))
  SAFE_ALLOCATE(integrand_xc(1:mesh%np, 1:cas%nik))
  SAFE_ALLOCATE(X(pot)(1:mesh%np))
  SAFE_ALLOCATE(buffer(1:mesh%np))
  SAFE_ALLOCATE(bufferx(1:mesh%np))
  SAFE_ALLOCATE(buffery(1:mesh%np))
  SAFE_ALLOCATE(bufferz(1:mesh%np))
  SAFE_ALLOCATE(deltav(1:mesh%np))

  ! coefficients for Hartree potential
  coeff_vh = - cas%kernel_lrc_alpha / (M_FOUR * M_PI)
  if(.not. cas%triplet) coeff_vh = coeff_vh + M_ONE
  if (ks%sic_type == SIC_ADSIC) coeff_vh = coeff_vh*(M_ONE - M_ONE/st%qtot)

  ! precompute buffer once for photon terms
  if ((cas%has_photons).and.(cas%type == CASIDA_CASIDA)) then
    bufferx(1:mesh%np) = R_TOTYPE(M_ZERO)
    buffery(1:mesh%np) = R_TOTYPE(M_ZERO)
    bufferz(1:mesh%np) = R_TOTYPE(M_ZERO)
    do ii = 1, cas%pt_nmodes
        buffer(1:mesh%np) = (cas%pt%lambda(ii)* &
             ( cas%pt%pol(ii,1)*mesh%x(1:mesh%np, 1)  &
             + cas%pt%pol(ii,2)*mesh%x(1:mesh%np, 2)  &
             + cas%pt%pol(ii,3)*mesh%x(1:mesh%np, 3)))
        bufferx(1:mesh%np) = bufferx(1:mesh%np) + &
            buffer(1:mesh%np)*cas%pt%lambda(ii)*cas%pt%pol(ii,1)
        buffery(1:mesh%np) = buffery(1:mesh%np) + &
            buffer(1:mesh%np)*cas%pt%lambda(ii)*cas%pt%pol(ii,2)
        bufferz(1:mesh%np) = bufferz(1:mesh%np) + &
            buffer(1:mesh%np)*cas%pt%lambda(ii)*cas%pt%pol(ii,3)
    end do
  end if

  ! calculate the matrix elements of (v + fxc)
  do jb_local = 1, cas%nb_rows
    if(.not. cas%distributed_matrix) then
      actual = actual + 1
      if(mod(actual, cas%mpi_grp%size) /= cas%mpi_grp%rank) cycle
    end if
    jb = get_global_row(cas, jb_local)

    ! first electron-hole contributions
    if(jb <= cas%n_pairs) then

      ! we only count diagonals for Petersilka
      if(cas%type == CASIDA_PETERSILKA) counter = counter + 1

      ! compute rho (order of indices important!) and potential (depends only on jb)
      call X(casida_get_rho)(st, mesh, cas%pair(jb)%i, cas%pair(jb)%a, cas%pair(jb)%kk, rho_j)

      if(.not. is_forces_ .and. abs(coeff_vh) > M_EPSILON) then
        X(pot)(1:mesh%np) = M_ZERO
        if(hm%theory_level /= INDEPENDENT_PARTICLES) then
          call X(poisson_solve)(hm%psolver, X(pot), rho_j, all_nodes=.false.)
        end if
      end if

      ! compute part of fxc, spin-resolved
      do ii = 1, cas%nik
        integrand_xc(1:mesh%np, ii) = rho_j(1:mesh%np)*xc(1:mesh%np, ii, cas%pair(jb)%kk)
      end do

      if ((cas%has_photons).and.(cas%type == CASIDA_CASIDA)) then
        rhobufferx = X(mf_integrate)(mesh, rho_j(1:mesh%np)*mesh%x(1:mesh%np, 1))
        rhobuffery = X(mf_integrate)(mesh, rho_j(1:mesh%np)*mesh%x(1:mesh%np, 2))
        rhobufferz = X(mf_integrate)(mesh, rho_j(1:mesh%np)*mesh%x(1:mesh%np, 3))
      end if

      ! take care of not computing elements twice for the symmetric matrix
      ia_length = cas%n_pairs / 2
      if(mod(cas%n_pairs, 2) == 0 .and. jb > cas%n_pairs/2) then
        ia_length = ia_length - 1
      end if

      do ia_local = 1, cas%nb_cols
        ia = get_global_col(cas, ia_local)

        if(jb+ia_length <= cas%n_pairs) then
          ! from diagonal to the right
          if(.not.(ia >= jb .and. ia <= jb+ia_length)) cycle
        else
          ! wrap around the end of the electron-hole part
          if(.not.((ia >= jb .and. ia <= cas%n_pairs) .or. &
                    ia <= jb+ia_length-cas%n_pairs)) cycle
        end if
        if(cas%distributed_matrix) then
          ! save rank number for later communcation (+1 is corrected later again)
          rank_of_element(ia, jb) = cas%mpi_grp%rank + 1
        end if

        ! only calculate off-diagonals in degenerate subspace
        if(cas%type == CASIDA_PETERSILKA .and. isnt_degenerate(cas, st, ia, jb)) cycle
        counter = counter + 1


        if(.not. is_saved(ia, jb)) then
          ! ---------------------------------------------------------
          !> calculates the matrix elements <i(p),a(p)|v|j(q),b(q)> and/or <i(p),a(p)|xc|j(q),b(q)>
          if (jb /= ia) then
            call X(casida_get_rho)(st, mesh, cas%pair(ia)%a, cas%pair(ia)%i, cas%pair(ia)%kk, rho_i)
          else
            !> no need to get rho twice
            rho_i(1:mesh%np) = rho_j(1:mesh%np)
          end if

          !  first the Hartree part
          if(.not. is_forces_ .and. abs(coeff_vh) > M_EPSILON) then
            ! value of pot is retained between calls
            mtxel_vh = coeff_vh * X(mf_dotp)(mesh, rho_i(:), X(pot)(:))
          else
            mtxel_vh = M_ZERO
          end if

          ! now the exchange part
          mtxel_xc = X(mf_dotp)(mesh, rho_i(:), integrand_xc(:, cas%pair(ia)%kk))

          mtxel_vm = M_ZERO
          if ((cas%has_photons).and.(cas%type == CASIDA_CASIDA)) then
            ! buffer is precomputed once before double loop
            mtxel_vm = mtxel_vm + &
            rhobufferx * X(mf_dotp)(mesh, rho_i(:), bufferx(:)) + &
            rhobuffery * X(mf_dotp)(mesh, rho_i(:), buffery(:)) + &
            rhobufferz * X(mf_dotp)(mesh, rho_i(:), bufferz(:))
          end if

          matrix(jb_local, ia_local) = mtxel_vh + mtxel_xc + mtxel_vm

        end if
        if(.not. cas%distributed_matrix) then
          if(jb /= ia) matrix(ia, jb) = R_CONJ(matrix(jb, ia))
        end if
        if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(counter, maxcount)
      end do

      ! compute photon part, reuse rho_j to avoid excessive memory access
      if ((cas%has_photons).and.(cas%type == CASIDA_CASIDA)) then
        do ia_local = 1, cas%nb_cols
          ia = get_global_col(cas, ia_local)
          if(ia <= cas%n_pairs) cycle
          iimode = ia - cas%n_pairs
          deltav(1:mesh%np) = cas%pt%lambda(iimode)*&
                       ( cas%pt%pol(iimode,1)*mesh%x(1:mesh%np, 1)  &
                       + cas%pt%pol(iimode,2)*mesh%x(1:mesh%np, 2)  &
                       + cas%pt%pol(iimode,3)*mesh%x(1:mesh%np, 3))


          xx(jb) = X(mf_integrate)(mesh, deltav(1:mesh%np)*rho_j(1:mesh%np))

          matrix(jb_local, ia_local) = sqrt(M_HALF*cas%pt%omega(ia-cas%n_pairs))*xx(jb)
          if(.not. cas%distributed_matrix) then
            matrix(ia, jb) = matrix(jb, ia)
          end if
        end do
      end if
    end if
    if(mpi_grp_is_root(mpi_world)) call loct_progress_bar(counter, maxcount)
  end do

  SAFE_DEALLOCATE_A(rho_i)
  SAFE_DEALLOCATE_A(rho_j)
  SAFE_DEALLOCATE_A(integrand_xc)
  SAFE_DEALLOCATE_A(X(pot))
  SAFE_DEALLOCATE_A(buffer)
  SAFE_DEALLOCATE_A(bufferx)
  SAFE_DEALLOCATE_A(buffery)
  SAFE_DEALLOCATE_A(bufferz)
  SAFE_DEALLOCATE_A(deltav)

  if(mpi_grp_is_root(mpi_world)) then
    call loct_progress_bar(maxcount, maxcount)
    ! complete progress bar
    write(stdout, '(1x)')
  end if

  ! sum all matrix elements
  if(cas%parallel_in_eh_pairs .and. .not. cas%distributed_matrix) then
    call comm_allreduce(cas%mpi_grp, matrix)
  end if
  if(cas%distributed_matrix) then
#ifdef HAVE_SCALAPACK
    ! add transpose of matrix to get full matrix
    SAFE_ALLOCATE(buffer_transpose(1:cas%nb_rows,1:cas%nb_cols))
    buffer_transpose(1:cas%nb_rows,1:cas%nb_cols) = matrix(1:cas%nb_rows,1:cas%nb_cols)
    ! set diagonal to zero and add transpose of buffer to matrix to get full matrix
    do jb_local = 1, cas%nb_rows
      jb = get_global_row(cas, jb_local)
      do ia_local = 1, cas%nb_cols
        ia = get_global_col(cas, ia_local)
        if(ia == jb) buffer_transpose(jb_local,ia_local) = M_ZERO
      end do
    end do
    ! this call adds transpose/hermitian conjugate and saves it to matrix
    call pblas_tran(cas%n, cas%n, R_TOTYPE(M_ONE), buffer_transpose(1,1), 1, 1, cas%desc(1), &
      R_TOTYPE(M_ONE), matrix(1,1), 1, 1, cas%desc(1))
    SAFE_DEALLOCATE_A(buffer_transpose)
#endif
  end if

  if(.not. cas%has_photons) then
    if(cas%distributed_matrix) then
      call comm_allreduce(cas%mpi_grp, rank_of_element)
      ! we subtract 1 again here; thus elements not associated to any rank
      ! are now negative
      rank_of_element = rank_of_element - 1
    end if

    ! now write out the restart files
    iunit = restart_open(cas%restart_dump, restart_file, position='append')
    if(mpi_grp_is_root(mesh%mpi_grp)) then
      do ia = 1, cas%n_pairs
        do jb = 1, cas%n_pairs
          if(cas%distributed_matrix) then
            if(rank_of_element(ia, jb) < 0) cycle
          else
            if(jb < ia) cycle
          end if
          call local_indices(cas, ia, jb, on_this_processor, ia_local, jb_local)
          if(on_this_processor) then
            save_value = matrix(jb_local, ia_local)
            write_value = .not. is_saved(ia, jb) .and. is_calcd(ia, jb)
          end if
          if(mpi_grp_is_root(mpi_world)) then
            if(.not.on_this_processor .and. cas%distributed_matrix) then
#ifdef HAVE_SCALAPACK
              src = rank_of_element(ia, jb)
              call MPI_Recv(write_value, 1, MPI_LOGICAL, src, 0, cas%mpi_grp%comm, mpi_status, mpi_err)
              call MPI_Recv(save_value, 1, R_MPITYPE, src, 0, cas%mpi_grp%comm, mpi_status, mpi_err)
#else
              ! if scalapack is not used, this branch is never taken
              ASSERT(.false.)
#endif
            end if
            if (write_value) call X(write_K_term)(cas, save_value, iunit, ia, jb)
          else
            if(on_this_processor .and. cas%distributed_matrix) then
#ifdef HAVE_SCALAPACK
              call MPI_Send(write_value, 1, MPI_LOGICAL, 0, 0, cas%mpi_grp%comm, mpi_err)
              call MPI_Send(save_value, 1, R_MPITYPE, 0, 0, cas%mpi_grp%comm, mpi_err)
#else
              ! if scalapack is not used, this branch is never taken
              ASSERT(.false.)
#endif
            end if
          end if
        end do
      end do
    end if
    if (iunit > 0) call restart_close(cas%restart_dump, iunit)
  end if

  SAFE_DEALLOCATE_A(is_saved)
  SAFE_DEALLOCATE_A(is_calcd)
  if(cas%distributed_matrix) then
    SAFE_DEALLOCATE_A(rank_of_element)
  end if

  SAFE_DEALLOCATE_A(xx)

  call profiling_out(prof)
  POP_SUB(X(casida_get_matrix))

contains
  ! ---------------------------------------------------------
  subroutine load_saved(matrix, is_saved, restart_file)
    R_TYPE,           intent(out) :: matrix(:,:)
    logical,          intent(out) :: is_saved(:,:)
    character(len=*), intent(in)  :: restart_file

    integer :: iunit, err, num_saved
    integer :: ia, jb, ii, aa, ik, jj, bb, jk, ia_local, jb_local
    R_TYPE  :: val
    logical :: on_this_processor

    PUSH_SUB(X(casida_get_matrix).load_saved)

    is_saved = .false.
    matrix = M_ZERO
    num_saved = 0

    ! if fromScratch, we already deleted the restart files
    iunit = restart_open(cas%restart_load, restart_file, silent = .true.)
    if(mpi_grp_is_root(mpi_world) .or. cas%distributed_matrix) then
      if (iunit > 0) then
        do
          ! for scalapack layout: read on rank 0 and broadcast to others
          if(mpi_grp_is_root(mpi_world)) then
            read(iunit, fmt=*, iostat=err) ii, aa, ik, jj, bb, jk, val
          end if
#ifdef HAVE_MPI
          ! broadcast error first
          if(cas%distributed_matrix) then
            call MPI_Bcast(err, 1, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
          end if
#endif
          ! exit if there was an error or if the end of the file is reached
          if(err /= 0) exit
#ifdef HAVE_MPI
          ! now broadcast values read from file
          if(cas%distributed_matrix) then
            call MPI_Bcast(ii, 1, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
            call MPI_Bcast(aa, 1, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
            call MPI_Bcast(ik, 1, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
            call MPI_Bcast(jj, 1, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
            call MPI_Bcast(bb, 1, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
            call MPI_Bcast(jk, 1, MPI_INTEGER, 0, mpi_world%comm, mpi_err)
            call MPI_Bcast(val, 1, R_MPITYPE, 0, mpi_world%comm, mpi_err)
          end if
#endif

          if(ii < 1 .or. aa < 1 .or. ik < 1) then
            message(1) = "Illegal indices in '" // trim(restart_file) // "': working from scratch."
            call messages_warning(1)
            call restart_rm(cas%restart_load, trim(restart_file))
            ! if file is corrupt, do not trust anything that was read
            is_saved = .false.
            exit
          end if

          ia = cas%index(ii, aa, ik)
          jb = cas%index(jj, bb, jk)

          if(ia > 0 .and. jb > 0) then
            call local_indices(cas, ia, jb, on_this_processor, ia_local, jb_local)
            if(on_this_processor) then
              matrix(jb_local, ia_local) = val
              is_saved(ia, jb) = .true.
              if(.not. cas%distributed_matrix) then
                matrix(ia_local, jb_local) = R_CONJ(val)
                is_saved(jb, ia) = .true.
              end if
              num_saved = num_saved + 1
            end if
          end if
        end do

        write(6,'(a,i5,a,i8,a,a)') 'Rank ', mpi_world%rank, ' read ', num_saved, &
          ' saved elements from ', trim(restart_file)
      else if(.not. cas%fromScratch) then
        message(1) = "Could not find restart file '" // trim(restart_file) // "'. Starting from scratch."
        call messages_warning(1)
      end if
    end if
    if (iunit > 0) call restart_close(cas%restart_load, iunit)

    ! if no file found, root has no new information to offer the others
#ifdef HAVE_MPI
    if(.not. cas%distributed_matrix) then
      call MPI_Bcast(is_saved(1, 1), cas%n_pairs**2, MPI_LOGICAL, 0, mpi_world%comm, mpi_err)
    end if
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
  
  write(iunit,*) cas%pair(ia)%i, cas%pair(ia)%a, cas%pair(ia)%kk, &
    cas%pair(jb)%i, cas%pair(jb)%a, cas%pair(jb)%kk, mat_val
  
  POP_SUB(X(write_K_term))
end subroutine X(write_K_term)

! ---------------------------------------------------------
subroutine X(casida_forces)(cas, sys, mesh, st)
  type(casida_t),      intent(inout) :: cas
  type(electrons_t),   intent(inout) :: sys
  type(mesh_t),        intent(in) :: mesh
  type(states_elec_t), intent(inout) :: st
  
  integer :: ip, iatom, idir, is1, is2, ierr, ik, ia
  FLOAT, allocatable :: ddl_rho(:,:), kxc(:,:,:,:)
  FLOAT, allocatable :: lr_fxc(:,:,:)
  R_TYPE, allocatable :: lr_hmat1(:,:,:)
  CMPLX, allocatable :: zdl_rho(:,:)
  FLOAT :: factor = CNST(1e6) ! FIXME: allow user to set
  character(len=MAX_PATH_LEN) :: restart_filename
  type(restart_t) :: restart_vib

  PUSH_SUB(X(casida_forces))
  
  if(cas%type == CASIDA_CASIDA) then
    message(1) = "Forces for Casida theory level not implemented"
    call messages_warning(1)
    POP_SUB(X(casida_forces))
    return
  end if
  
  if(cas%type /= CASIDA_EPS_DIFF) then
    SAFE_ALLOCATE(kxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin, 1:st%d%nspin))
    kxc = M_ZERO
    ! not spin polarized so far
    call xc_get_kxc(sys%ks%xc, mesh, sys%namespace, cas%rho, st%d%ispin, kxc(:, :, :, :))
  end if
  
  message(1) = "Reading vib_modes density for calculating excited-state forces."
  call messages_info(1)
  
  ! vib_modes for complex wavefunctions will write complex rho, however it is actually real
  ! FIXME: why?
#ifdef R_TCOMPLEX
  SAFE_ALLOCATE(zdl_rho(1:mesh%np, 1:st%d%nspin))
#endif
  SAFE_ALLOCATE(ddl_rho(1:mesh%np, 1:st%d%nspin))

  if (cas%type /= CASIDA_EPS_DIFF) then
    SAFE_ALLOCATE(lr_fxc(1:mesh%np, 1:st%d%nspin, 1:st%d%nspin))
  end if
  
  if(cas%type == CASIDA_EPS_DIFF) then
    cas%X(mat_save) = M_ZERO
    do ia = 1, cas%n_pairs
      cas%X(mat_save)(ia, ia) = cas%w(ia)
    end do
  end if
  
  SAFE_ALLOCATE(lr_hmat1(1:cas%nst, 1:cas%nst, 1:cas%nik))
  SAFE_ALLOCATE(cas%X(lr_hmat2)(1:cas%n_pairs, 1:cas%n_pairs))
  SAFE_ALLOCATE(cas%X(mat2)(1:cas%n_pairs, 1:cas%n_pairs))
  SAFE_ALLOCATE(cas%X(w2)(1:cas%n_pairs))

  call restart_init(restart_vib, sys%namespace, RESTART_VIB_MODES, RESTART_TYPE_LOAD, sys%mc, ierr, mesh = sys%gr%mesh)

  do iatom = 1, sys%ions%natoms
    do idir = 1, mesh%sb%dim
      
      if(ierr == 0) &
        call X(lr_load_rho)(X(dl_rho), sys%space, sys%gr%mesh, st%d%nspin, restart_vib, phn_rho_tag(iatom, idir), ierr)
      if(ierr /= 0) then
        message(1) = "Could not read vib_modes density; previous vib_modes calculation required."
        call messages_fatal(1)
      end if

      if(allocated(zdl_rho)) then
        if(any(abs(aimag(zdl_rho)) > M_EPSILON)) then
          message(1) = "Vib_modes density has an imaginary part."
          call messages_warning(1)
        end if
        ddl_rho = TOFLOAT(zdl_rho)
      end if
      call X(casida_get_lr_hmat1)(cas, sys, iatom, idir, ddl_rho, lr_hmat1)
      
      cas%X(lr_hmat2) = M_ZERO
      ! use them to make two-particle matrix elements (as for eigenvalues)
      do ik = 1, cas%nik
        call X(casida_lr_hmat2)(cas, st, lr_hmat1, ik)
      end do

      if (cas%type /= CASIDA_EPS_DIFF .and. cas%calc_forces_kernel) then
        do is2 = 1, st%d%nspin
          do is1 = 1, st%d%nspin
            do ip = 1, mesh%np
              lr_fxc(ip, is1, is2) = sum(kxc(ip, is1, is2, :) * ddl_rho(ip, :))
            end do
          end do
        end do

        write(restart_filename,'(a,i6.6,a,i1)') 'lr_kernel_', iatom, '_', idir
        if(cas%triplet) restart_filename = trim(restart_filename)//'_triplet'

        call X(casida_get_matrix)(cas, sys%hm, st, sys%ks, mesh, cas%X(mat2), lr_fxc, restart_filename, &
          is_forces = .true.)
        cas%X(mat2) = cas%X(mat2) * casida_matrix_factor(cas, sys)
      else
        cas%X(mat2) = M_ZERO
      end if
      
      cas%X(mat2) = cas%X(mat_save) * factor + cas%X(lr_hmat2) + cas%X(mat2)
      ! use parallel eigensolver here; falls back to serial solver if ScaLAPACK
      ! not available
      call lalg_eigensolve_parallel(cas%n_pairs, cas%X(mat2), cas%X(w2))
      do ia = 1, cas%n_pairs
        cas%forces(idir, iatom, cas%ind(ia)) = factor * cas%w(cas%ind(ia)) - R_REAL(cas%X(w2)(ia))
      end do
    end do
  end do

  call restart_end(restart_vib)
  
  if(cas%type /= CASIDA_EPS_DIFF) then
    SAFE_DEALLOCATE_A(kxc)
    SAFE_DEALLOCATE_A(lr_fxc)
  end if
  SAFE_DEALLOCATE_A(ddl_rho)
  SAFE_DEALLOCATE_A(zdl_rho)

  SAFE_DEALLOCATE_A(lr_hmat1)
  SAFE_DEALLOCATE_A(cas%X(lr_hmat2))
  SAFE_DEALLOCATE_A(cas%X(mat2))
  SAFE_DEALLOCATE_A(cas%X(w2))
  
  if(cas%calc_forces_scf) then
    call forces_calculate(sys%gr, sys%namespace, sys%ions, sys%hm, st, sys%ks)
    do ia = 1, cas%n_pairs
      cas%forces(:, :, ia) = cas%forces(:, :, ia) + sys%ions%tot_force
    end do
  end if
  
  POP_SUB(X(casida_forces))

end subroutine X(casida_forces)

! ---------------------------------------------------------
subroutine X(casida_get_lr_hmat1)(cas, sys, iatom, idir, dl_rho, lr_hmat1)
  type(casida_t),      intent(in)     :: cas
  type(electrons_t),   intent(inout)  :: sys
  integer,             intent(in)     :: iatom
  integer,             intent(in)     :: idir
  FLOAT,               intent(in)     :: dl_rho(:,:)
  R_TYPE,              intent(out)    :: lr_hmat1(:,:,:)

  FLOAT, allocatable :: hvar(:,:,:)
  integer :: ik, ist, jst, iunit, err, ii, aa, num_saved
  type(pert_t) :: ionic_pert
  character(len=MAX_PATH_LEN) :: restart_filename
  R_TYPE :: val
  logical :: all_done
  logical, allocatable :: is_saved(:,:,:)

  PUSH_SUB(X(casida_get_lr_hmat1))

  lr_hmat1 = M_ZERO
  SAFE_ALLOCATE(is_saved(1:cas%nst, 1:cas%nst, 1:cas%nik))
  is_saved = .false.
  num_saved = 0

  ! if fromScratch, we already deleted the restart files
  write(restart_filename,'(a,i6.6,a,i1)') 'lr_hmat1_', iatom, '_', idir
  iunit = restart_open(cas%restart_load, restart_filename, silent = .true.)
  ! NOTE: this restart file is different from 'kernel' because it is one-state not two-state matrix elements
  if(mpi_grp_is_root(mpi_world)) then
    if(iunit > 0) then
      do
        read(iunit, fmt=*, iostat=err) ii, aa, ik, val
        if(err /= 0) exit

        if(ii < 1 .or. aa < 1 .or. ik < 1) then
          message(1) = "Illegal indices in '" // trim(restart_filename) // "': working from scratch."
          call messages_warning(1)
          call restart_rm(cas%restart_load, trim(restart_filename))
          ! if file is corrupt, do not trust anything that was read
          is_saved = .false.
          exit
        end if

        ! FIXME: what about elements which are always zero?
        if(ii <= cas%nst .and. aa <= cas%nst .and. ik <= cas%nik) then
          lr_hmat1(ii, aa, ik) = val
          is_saved(ii, aa, ik) = .true.
          lr_hmat1(aa, ii, ik) = R_CONJ(val)
          is_saved(aa, ii, ik) = .true.
          num_saved = num_saved + 1
        end if
      end do
      write(6,'(a,i8,a,a)') 'Read ', num_saved, ' saved elements from ', trim(restart_filename)

    else if(.not. cas%fromScratch) then
      message(1) = "Could not find restart file '" // trim(restart_filename) // "'. Starting from scratch."
      call messages_warning(1)
    end if
  end if
  if (iunit > 0) call restart_close(cas%restart_load, iunit)

#ifdef HAVE_MPI
    call MPI_Bcast(is_saved(1, 1, 1), cas%nst**2, MPI_LOGICAL, 0, mpi_world%comm, mpi_err)
#endif

  all_done = .true.
  do ik = 1, cas%nik
    all_done = all_done .and. all(is_saved(1:cas%n_occ(ik), 1:cas%n_occ(ik), ik)) &
      .and. all(is_saved(cas%n_occ(ik) + 1:cas%nst, cas%n_occ(ik) + 1:cas%nst, ik))
  end do

  if(all_done) then
    SAFE_DEALLOCATE_A(is_saved)
    POP_SUB(X(casida_get_lr_hmat1))
    return
  end if

  call pert_init(ionic_pert, sys%namespace, PERTURBATION_IONIC, sys%gr, sys%ions)
  call pert_setup_atom(ionic_pert, iatom)
  call pert_setup_dir(ionic_pert, idir)

  SAFE_ALLOCATE(hvar(1:sys%gr%mesh%np, 1:sys%st%d%nspin, 1:1))
  call dcalc_hvar(.true., sys%gr%mesh, sys%st, sys%hm, sys%ks%xc, dl_rho, 1, hvar, fxc = cas%fxc)

  ! FIXME: do this only for states called for in CasidaKohnShamStates
  do ik = 1, cas%nik
    ! occ-occ matrix elements
    call X(casida_calc_lr_hmat1)(sys, ionic_pert, hvar, lr_hmat1, is_saved, 1, cas%n_occ(ik), ik)
    ! unocc-unocc matrix elements
    call X(casida_calc_lr_hmat1)(sys, ionic_pert, hvar, lr_hmat1, is_saved, cas%n_occ(ik) + 1, cas%nst, ik)
  end do

  SAFE_DEALLOCATE_A(hvar)
  call pert_end(ionic_pert)

  iunit = restart_open(cas%restart_dump, restart_filename, position = 'append')
  if (mpi_grp_is_root(mpi_world)) then
    do ik = 1, cas%nik
      do ist = 1, cas%nst
        do jst = ist, cas%nst
          if(.not. is_saved(ist, jst, ik)) write(iunit, *) ist, jst, ik, lr_hmat1(ist, jst, ik)
        end do
      end do
    end do
  end if
  call restart_close(cas%restart_dump, iunit)

  SAFE_DEALLOCATE_A(is_saved)

  POP_SUB(X(casida_get_lr_hmat1))

end subroutine X(casida_get_lr_hmat1)


! ---------------------------------------------------------
subroutine X(casida_solve)(cas, sys)
  type(casida_t),      intent(inout) :: cas
  type(electrons_t),   intent(in)    :: sys

  FLOAT :: eig_diff
  FLOAT, allocatable :: occ_diffs(:)
  integer :: ia, jb, ia_local, jb_local
#ifdef HAVE_MPI
  integer :: info
#endif
  type(profile_t), save :: prof
#ifdef HAVE_SCALAPACK
  R_TYPE, allocatable :: eigenvectors(:,:)
  R_TYPE, allocatable :: work(:)
  R_TYPE :: worksize
#ifdef R_TCOMPLEX
  R_TYPE, allocatable :: rwork(:)
  R_TYPE :: rworksize
#endif
#ifdef HAVE_ELPA
  class(elpa_t), pointer :: elpa
#endif
#endif

  PUSH_SUB(X(casida_solve))

  cas%X(mat) = cas%X(mat) * casida_matrix_factor(cas, sys)

  ! Note: this is not just a matter of implementation. CASIDA_CASIDA assumes real wfns in the derivation.
  if(states_are_complex(sys%st) .and. (cas%type == CASIDA_VARIATIONAL .or. cas%type == CASIDA_CASIDA)) then
    message(1) = "Variational and full Casida theory levels cannot be used with complex wavefunctions."
    call messages_fatal(1, only_root_writes = .true.)
    ! see section II.D of CV(2) paper regarding this assumption. Would be Eq. 30 with complex wfns.
  end if

  if(cas%type == CASIDA_CASIDA) then
    do ia = 1, cas%n_pairs
      cas%s(ia) = cas%el_per_state / ( &
        ( sys%st%eigenval(cas%pair(ia)%a, cas%pair(ia)%kk) - sys%st%eigenval(cas%pair(ia)%i, cas%pair(ia)%kk) ) &
        * ( sys%st%occ(cas%pair(ia)%i, cas%pair(ia)%kk) - sys%st%occ(cas%pair(ia)%a, cas%pair(ia)%kk) ) )
    end do
  end if

  if(cas%type /= CASIDA_CASIDA) then
    SAFE_ALLOCATE(occ_diffs(1:cas%n_pairs))
    do ia = 1, cas%n_pairs
      occ_diffs(ia) = (sys%st%occ(cas%pair(ia)%i, cas%pair(ia)%kk) - sys%st%occ(cas%pair(ia)%a, cas%pair(ia)%kk)) &
        / cas%el_per_state
    end do
  end if

  ! complete the matrix
  do jb_local = 1, cas%nb_rows
    jb = get_global_row(cas, jb_local)
    if(jb > cas%n_pairs) cycle
      eig_diff = sys%st%eigenval(cas%pair(jb)%a, cas%pair(jb)%kk) - sys%st%eigenval(cas%pair(jb)%i, cas%pair(jb)%kk)

    do ia_local = 1, cas%nb_cols
      ia = get_global_col(cas, ia_local)
      if(ia > cas%n_pairs) cycle
      ! FIXME: need the equivalent of this stuff for forces too.
      if(cas%type == CASIDA_CASIDA) then
        cas%X(mat)(jb_local, ia_local) = M_TWO * cas%X(mat(jb_local, ia_local)) &
          / sqrt(cas%s(jb) * cas%s(ia))
      else
        cas%X(mat)(jb_local, ia_local) = cas%X(mat(jb_local, ia_local)) * sqrt(occ_diffs(jb) * occ_diffs(ia))
      end if
      if(ia == jb) then
        if(cas%type == CASIDA_CASIDA) then
          cas%X(mat)(jb_local, ia_local) = eig_diff**2 + cas%X(mat)(jb_local, ia_local)
        else
          cas%X(mat)(jb_local, ia_local) = eig_diff + cas%X(mat)(jb_local, ia_local)
        end if
      end if
    end do
  end do

  if(cas%type /= CASIDA_CASIDA) then
    SAFE_DEALLOCATE_A(occ_diffs)
  end if

  if ((cas%has_photons).and.(cas%type == CASIDA_CASIDA)) then
    ! first diagonal for photon-photon interaction
    do jb_local = 1, cas%nb_rows
      jb = get_global_row(cas, jb_local)
      if(jb <= cas%n_pairs) cycle
      do ia_local = 1, cas%nb_cols
        ia = get_global_col(cas, ia_local)
        if(ia <= cas%n_pairs) cycle
        if(ia == jb) cas%X(mat)(jb_local,ia_local) = (cas%pt%omega(ia-cas%n_pairs))**2
      end do
    end do
    ! now other photon-electron elements
    do jb_local = 1, cas%nb_rows
      jb = get_global_row(cas, jb_local)
      do ia_local = 1, cas%nb_cols
        ia = get_global_col(cas, ia_local)
        ! lower left part
        if(jb > cas%n_pairs .and. ia <= cas%n_pairs) then
          cas%X(mat)(jb_local,ia_local) = M_TWO*cas%X(mat)(jb_local,ia_local) &
             *sqrt(cas%pt%omega(jb-cas%n_pairs))/sqrt(cas%s(ia))/sqrt(M_TWO)
        end if
        ! upper right part
        if(jb <= cas%n_pairs .and. ia > cas%n_pairs) then
          cas%X(mat)(jb_local,ia_local) = M_TWO*cas%X(mat)(jb_local,ia_local) &
             *sqrt(cas%pt%omega(ia-cas%n_pairs))/sqrt(cas%s(jb))/sqrt(M_TWO)
        end if
      end do
    end do
  end if

  if(cas%parallel_in_eh_pairs) then
#ifdef HAVE_MPI
    ! wait here instead of in diagonalization for better profiling
    call MPI_Barrier(cas%mpi_grp%comm, info)
#endif
  end if
  message(1) = "Info: Diagonalizing matrix for resonance energies."
  call messages_info(1)

  ! now we diagonalize the matrix
  call profiling_in(prof, "CASIDA_DIAGONALIZATION")
  if(cas%calc_forces) cas%X(mat_save) = cas%X(mat) ! save before gets turned into eigenvectors

  if(.not. cas%distributed_matrix) then
    call lalg_eigensolve_parallel(cas%n, cas%X(mat), cas%w)
  else
#ifdef HAVE_SCALAPACK
    SAFE_ALLOCATE(eigenvectors(cas%nb_rows, cas%nb_cols))

    if(cas%parallel_solver == SOLVER_ELPA) then
#ifdef HAVE_ELPA
      ! eigensolver settings (allocate workspace)
      if (elpa_init(20170403) /= elpa_ok) then
        write(message(1),'(a)') "ELPA API version not supported"
        call messages_fatal(1)
      endif
      elpa => elpa_allocate()

      ! set parameters describing the matrix
      call elpa%set("na", cas%n, info)
      call elpa%set("nev", cas%n, info)
      call elpa%set("local_nrows", cas%nb_rows, info)
      call elpa%set("local_ncols", cas%nb_cols, info)
      call elpa%set("nblk", cas%block_size, info)
      call elpa%set("mpi_comm_parent", cas%mpi_grp%comm, info)
      call elpa%set("process_row", cas%proc_grid%myrow, info)
      call elpa%set("process_col", cas%proc_grid%mycol, info)

      info = elpa%setup()

      call elpa%set("solver", elpa_solver_2stage, info)

      ! call eigensolver
      call elpa%eigenvectors(cas%X(mat), cas%w, eigenvectors, info)

      ! error handling
      if (info /= elpa_ok) then
        write(message(1),'(a,i6,a,a)') "Error in ELPA, code: ", info, ", message: ", &
          elpa_strerr(info)
        call messages_fatal(1)
      end if

      call elpa_deallocate(elpa)
      call elpa_uninit()
#endif
    else
#ifdef HAVE_SCALAPACK
      ! use ScaLAPACK solver if ELPA not available
      ! eigensolver settings (allocate workspace)
      ! workspace query
#ifdef R_TREAL
      call pdsyev(jobz='V', uplo='L', n=cas%n, &
        a=cas%X(mat)(1, 1), ia=1, ja=1, desca=cas%desc(1), w=cas%w(1), &
        z=eigenvectors(1, 1), iz=1, jz=1, descz=cas%desc(1), &
        work=worksize, lwork=-1, info=info)
#else
      call pzheev(jobz='V', uplo='L', n=cas%n, &
        a=cas%X(mat)(1, 1), ia=1, ja=1, desca=cas%desc(1), w=cas%w(1), &
        z=eigenvectors(1, 1), iz=1, jz=1, descz=cas%desc(1), &
        work=worksize, lwork=-1, rwork=rworksize, lrwork=-1, info=info)
#endif

      if(info /= 0) then
        write(message(1),'(a,i6)') "ScaLAPACK workspace query failure, error code=", info
        call messages_fatal(1)
      end if

      SAFE_ALLOCATE(work(1:int(worksize)))
#ifdef R_TCOMPLEX
      SAFE_ALLOCATE(rwork(1:int(rworksize)))
#endif

      ! call eigensolver
#ifdef R_TREAL
      call pdsyev(jobz='V', uplo='L', n=cas%n, &
        a=cas%X(mat)(1, 1) , ia=1, ja=1, desca=cas%desc(1), w=cas%w(1), &
        z=eigenvectors(1, 1), iz=1, jz=1, descz=cas%desc(1), &
        work=work(1), lwork=int(worksize), info=info)
#else
      call pzheev(jobz='V', uplo='L', n=cas%n, &
        a=cas%X(mat)(1, 1), ia=1, ja=1, desca=cas%desc(1), w=cas%w(1), &
        z=eigenvectors(1, 1), iz=1, jz=1, descz=cas%desc(1), &
        work=work(1), lwork=int(worksize), &
        rwork=rwork(1), lrwork=int(rworksize), info=info)
#endif

      SAFE_DEALLOCATE_A(work)
#ifdef R_TCOMPLEX
      SAFE_DEALLOCATE_A(rwork)
#endif

      ! error handling
      if(info /= 0) then
#ifdef R_TCOMPLEX
        write(message(1),'(3a,i5)') trim(message(1)), TOSTRING(pX(heev)), &
          ' returned error message ', info
#else
        write(message(1),'(3a,i5)') trim(message(1)), TOSTRING(pX(syev)), &
          ' returned error message ', info
#endif
!*  INFO    (global output) INTEGER
!*          = 0:  successful exit
!*          < 0:  If the i-th argument is an array and the j-entry had
!*                an illegal value, then INFO = -(i*100+j), if the i-th
!*                argument is a scalar and had an illegal value, then
!*                INFO = -i.
!*          > 0:  If INFO = 1 through N, the i(th) eigenvalue did not
!*                converge in DSTEQR2 after a total of 30*N iterations.
!*                If INFO = N+1, then PDSYEV has detected heterogeneity
!*                by finding that eigenvalues were not identical across
!*                the process grid.  In this case, the accuracy of
!*                the results from PDSYEV cannot be guaranteed.
        if(info < 0) then
          write(message(2), '(a,i5,a)') 'Argument number ', -info, ' had an illegal value.'
        else if(info == cas%n+1) then
          write(message(2), '(a)') 'Eigenvalues were not identical over the process grid.'
        else
          write(message(2), '(i5,a)') info, 'th eigenvalue did not converge.'
        end if
        call messages_fatal(2)
      end if
#endif
    end if
    ! save eigenvectors to X(mat)
    cas%X(mat)(1:cas%nb_rows,1:cas%nb_cols) = eigenvectors(1:cas%nb_rows,1:cas%nb_cols)
    SAFE_DEALLOCATE_A(eigenvectors)
#endif
  end if
  call profiling_out(prof)

  do ia = 1, cas%n
    if(cas%type == CASIDA_CASIDA) then
      if(cas%w(ia) < -M_EPSILON) then
        write(message(1),'(a,i4,a)') 'Casida excitation energy', ia, ' is imaginary.'
        call messages_warning(1)
        cas%w(ia) = -sqrt(-cas%w(ia))
      else
        cas%w(ia) = sqrt(cas%w(ia))
      end if
    else
      if(cas%w(ia) < -M_EPSILON) then
        write(message(1),'(a,i4,a)') 'For whatever reason, excitation energy', ia, ' is negative.'
        write(message(2),'(a)')      'This should not happen.'
        call messages_warning(2)
      end if
    end if

    cas%ind(ia) = ia ! diagonalization returns eigenvalues in order.
  end do

  POP_SUB(X(casida_solve))
end subroutine X(casida_solve)

! ---------------------------------------------------------
subroutine X(casida_write)(cas, sys)
  type(casida_t),      intent(in) :: cas
  type(electrons_t),   intent(in) :: sys
  
  character(len=5) :: str
  character(len=50) :: dir_name
  integer :: iunit, ia, jb, idim, ia_local, jb_local
  logical :: full_printing, on_this_processor
  FLOAT   :: norm, norm_e, norm_p
  
  PUSH_SUB(X(casida_write))

  full_printing = .false.
  if (cas%print_exst == 'all') full_printing = .true.

  if(mpi_grp_is_root(mpi_world)) then
    ! output excitation energies and oscillator strengths
    call io_mkdir(CASIDA_DIR, sys%namespace)
    iunit = io_open(CASIDA_DIR//trim(theory_name(cas)), sys%namespace, action='write')

    if(cas%type == CASIDA_EPS_DIFF) then
      write(iunit, '(2a4)', advance='no') 'From', '  To'
      if(cas%nik > 1) then
        write(iunit, '(a7)', advance='no') 'Spin/k'
      end if
    else
      write(iunit, '(6x)', advance='no')
    end if
    
    write(iunit, '(1x,a15)', advance='no') 'E [' // trim(units_abbrev(units_out%energy)) // ']' 
    do idim = 1, cas%space_dim
      write(iunit, '(1x,a15)', advance='no') '<' // index2axis(idim) // '> [' // trim(units_abbrev(units_out%length)) // ']' 
#ifdef R_TCOMPLEX
      write(iunit, '(16x)', advance='no')
#endif
    end do
    write(iunit, '(1x,a15)') '<f>'
    
    do ia = 1, cas%n
      if((cas%type == CASIDA_EPS_DIFF)) then
        write(iunit, '(2i4)', advance='no') cas%pair(cas%ind(ia))%i, cas%pair(cas%ind(ia))%a
        if(cas%nik > 1) then
          write(iunit, '(i7)', advance='no') cas%pair(cas%ind(ia))%kk
        end if
      else
        write(iunit, '(i6)', advance='no') cas%ind(ia)
      end if
      write(iunit, '(99(1x,es15.8))') units_from_atomic(units_out%energy, cas%w(cas%ind(ia))), &
        (units_from_atomic(units_out%length, cas%X(tm)(cas%ind(ia), idim)), idim=1,cas%space_dim), cas%f(cas%ind(ia))
    end do
    call io_close(iunit)
  
    if(cas%qcalc) call qcasida_write(cas, sys%namespace)
  
    if (.not.(cas%print_exst == "0" .or. cas%print_exst == "none")) then
      if(cas%type /= CASIDA_EPS_DIFF .or. cas%calc_forces) then
        dir_name = CASIDA_DIR//trim(theory_name(cas))//'_excitations'
        call io_mkdir(trim(dir_name), sys%namespace)
      end if
      
      do ia = 1, cas%n
        if(loct_isinstringlist(ia, cas%print_exst) .or. full_printing) then 
          write(str,'(i5.5)') ia
          
          ! output eigenvectors
          if(cas%type /= CASIDA_EPS_DIFF) then
            iunit = io_open(trim(dir_name)//'/'//trim(str), sys%namespace, action='write')
            ! First, a little header
            write(iunit,'(a,es14.5)') '# Energy ['// trim(units_abbrev(units_out%energy)) // '] = ', &
              units_from_atomic(units_out%energy, cas%w(cas%ind(ia)))
            do idim = 1, cas%space_dim
              write(iunit,'(a,2es14.5)') '# <' // index2axis(idim) // '> ['//trim(units_abbrev(units_out%length))// '] = ', &
                units_from_atomic(units_out%length, cas%X(tm)(cas%ind(ia), idim))
            end do

            ! this stuff should go BEFORE calculation of transition matrix elements!
            if(cas%type == CASIDA_TAMM_DANCOFF .or. cas%type == CASIDA_VARIATIONAL .or. cas%type == CASIDA_PETERSILKA) then
              call X(write_implied_occupations)(cas, iunit, cas%ind(ia))
            end if
            call io_close(iunit)
          end if
          
          if(cas%calc_forces .and. cas%type /= CASIDA_CASIDA) then
            iunit = io_open(trim(dir_name)//'/forces_'//trim(str)//'.xsf', sys%namespace, action='write')
            call write_xsf_geometry(iunit, sys%ions, sys%gr%mesh, forces = cas%forces(:, :, cas%ind(ia)))
            call io_close(iunit)
          end if
        end if
      end do
    end if
  end if

  if(cas%write_matrix .and. mpi_grp_is_root(sys%gr%mesh%mpi_grp)) then
      call X(write_distributed_matrix)(cas, cas%X(mat), &
        CASIDA_DIR//trim(theory_name(cas))//"_matrix")
  end if

  ! Output the norm files always in photon mode
  if(cas%has_photons) then
    ! compute and write norms
    if(mpi_grp_is_root(mpi_world)) then
      iunit = io_open(CASIDA_DIR//trim(theory_name(cas))//"_norms", sys%namespace, action='write')
      ! first, a header line
      write(iunit, '(6x)', advance='no')
      write(iunit, '(1x,a15)', advance='no') 'E [' // trim(units_abbrev(units_out%energy)) // ']'
      write(iunit, '(1x,a15)', advance='no') 'Electron part'
      write(iunit, '(1x,a15)') 'Photon part'
    end if
    do ia = 1, cas%n
      norm = M_ZERO
      norm_e = M_ZERO
      norm_p = M_ZERO
      do jb = 1, cas%n
        call local_indices(cas, cas%ind(ia), jb, on_this_processor, ia_local, jb_local)
        if(on_this_processor) then
          ! sum up norm, electronic norm and photonic norm on each processor
          norm = norm + abs(cas%X(mat)(jb_local, ia_local))**2
          if(jb <= cas%n_pairs) then
            norm_e = norm_e + abs(cas%X(mat)(jb_local, ia_local))**2
          else
            norm_p = norm_p + abs(cas%X(mat)(jb_local, ia_local))**2
          end if
        end if
      end do
      norm = dallreduce_sum(cas, norm)
      norm_e = dallreduce_sum(cas, norm_e)
      norm_p = dallreduce_sum(cas, norm_p)
      if(mpi_grp_is_root(mpi_world)) then
        ! write contributions to norm to file
        write(iunit, '(i6)', advance='no') cas%ind(ia)
        write(iunit, '(99(1x,es15.8))') units_from_atomic(units_out%energy, cas%w(cas%ind(ia))), &
          norm_e/norm, norm_p/norm
      end if
    end do
    if(mpi_grp_is_root(mpi_world)) then
      call io_close(iunit)
    end if
  end if

  ! Calculate and write the transition densities
  if(cas%trandens /= "0") then
    call X(get_transition_densities)(cas, sys)
  end if

  POP_SUB(X(casida_write))
end subroutine X(casida_write)

! ---------------------------------------------------------
subroutine X(write_implied_occupations)(cas, iunit, ind)
  type(casida_t), intent(in) :: cas
  integer,        intent(in) :: iunit
  integer,        intent(in) :: ind
  
  integer :: ik, ast, ist
  FLOAT :: occ
  
  PUSH_SUB(X(write_implied_occupations))
  
  write(iunit, '(a)')
  write(iunit, '(a)') '%Occupations'
  
  do ik = 1, cas%nik
    do ist = 1, cas%n_occ(ik)
      occ = M_ONE * cas%el_per_state
      do ast = cas%n_occ(ik) + 1, cas%nst
        if(cas%index(ist, ast, ik) == 0) cycle  ! we were not using this state
        occ = occ - abs(cas%X(mat)(cas%index(ist, ast, ik), ind))**2
      end do
      write(iunit, '(f8.6,a)', advance='no') occ, ' | '
    end do
    do ast = cas%n_occ(ik) + 1, cas%nst
      occ = M_ZERO
      do ist = 1, cas%n_occ(ik)
        if(cas%index(ist, ast, ik) == 0) cycle  ! we were not using this state
        occ = occ + abs(cas%X(mat)(cas%index(ist, ast, ik), ind))**2
      end do
      write(iunit, '(f8.6)', advance='no') occ
      if(ast < cas%n_occ(ik) + cas%n_unocc(ik)) write(iunit, '(a)', advance='no') ' | '
    end do
    write(iunit, '(a)')
  end do
  
  write(iunit, '(a)') '%'
  
  POP_SUB(X(write_implied_occupations))
end subroutine X(write_implied_occupations)

subroutine X(write_distributed_matrix)(cas, matrix, filename)
  implicit none
  type(casida_t), intent(in) :: cas
  R_TYPE, intent(in) :: matrix(:,:)
  character(len=*), intent(in) :: filename
#ifdef HAVE_MPI
  integer :: outfile, mpistatus, ierr
#endif

  PUSH_SUB(X(write_distributed_matrix))

  if(.not. cas%distributed_matrix) then
    message(1) = "Cannot write distributed matrix if not using ScaLAPACK layout"
    call messages_info(1)
    return
  end if

#ifdef HAVE_MPI
  ! create MPI IO types
  call MPI_Type_create_darray(cas%mpi_grp%size, cas%mpi_grp%rank, 2, (/ cas%n, cas%n /), &
    (/ MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC /), (/ cas%block_size, cas%block_size/), &
    (/ cas%proc_grid%nprow, cas%proc_grid%npcol /), MPI_ORDER_FORTRAN, MPI_FLOAT, &
    cas%darray, ierr)
  call MPI_Type_commit(cas%darray, ierr)

  call MPI_Barrier(cas%mpi_grp%comm, ierr)
  ! write out casida matrix
  call MPI_File_open(cas%mpi_grp%comm, trim(filename), MPI_MODE_CREATE+MPI_MODE_WRONLY, &
    MPI_INFO_NULL, outfile, ierr)
  if(mpi_grp_is_root(cas%mpi_grp)) then
    ! write size of matrix
    call MPI_File_write(outfile, cas%n, 1, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  end if
  ! write matrix with displacement of one integer (size)
  call MPI_File_set_view(outfile, sizeof(cas%n), R_MPITYPE, cas%darray, "native", MPI_INFO_NULL, ierr)
  call MPI_File_write_all(outfile, matrix, cas%nb_rows*cas%nb_cols, R_MPITYPE, mpistatus, ierr)
  call MPI_File_close(outfile, ierr)
#endif
  POP_SUB(X(write_distributed_matrix))
end subroutine X(write_distributed_matrix)

! communication function used for sums over the casida matrix
R_TYPE function X(allreduce_sum)(cas, variable) result(output)
  type(casida_t), intent(in) :: cas
  R_TYPE, intent(in) :: variable
#ifdef HAVE_SCALAPACK
  R_TYPE :: buffer
  integer :: ierr
#endif
  PUSH_SUB(X(allreduce_sum))
  if(.not. cas%distributed_matrix) then
    output = variable
  else
#ifdef HAVE_SCALAPACK
    call MPI_Allreduce(variable, buffer, 1, R_MPITYPE, MPI_SUM, &
      cas%mpi_grp%comm, ierr)
    output = buffer
#endif
  end if
  POP_SUB(X(allreduce_sum))
end function X(allreduce_sum)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
