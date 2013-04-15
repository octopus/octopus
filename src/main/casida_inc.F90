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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

subroutine X(oscillator_strengths)(cas, mesh, st)
  type(casida_t), intent(inout) :: cas
  type(mesh_t), intent(in) :: mesh
  type(states_t), intent(in) :: st

  FLOAT, allocatable :: deltav(:), psi_a(:)
  R_TYPE, allocatable :: xx(:)
  CMPLX, allocatable :: zf(:), zx(:)
  FLOAT, allocatable :: gaus_leg_points(:), gaus_leg_weights(:)
  integer :: ii, jj, ia, ip, idir
  FLOAT :: theta, phi, qlen
  FLOAT :: qvect(MAX_DIM)

  PUSH_SUB(X(oscillator_strengths))

  if(cas%triplet) then
    cas%tm(:,:) = M_ZERO
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
        zf(ip) = exp(M_zI*dot_product(cas%qvector(1:mesh%sb%dim), mesh%x(ip, 1:mesh%sb%dim)))*zf(ip)*psi_a(ip)
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
              zf(ip) = exp(M_zI*dot_product(qvect(1:mesh%sb%dim), mesh%x(ip, 1:mesh%sb%dim)))*zf(ip)*psi_a(ip)
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
      cas%tm(ia, idir) = X(transition_matrix_element)(cas, ia, xx)
    end do
  end do
  SAFE_DEALLOCATE_A(xx)

  ! And the oscillator strengths.
  do ia = 1, cas%n_pairs
    cas%f(ia) = (M_TWO / mesh%sb%dim) * cas%w(ia) * sum( (abs(cas%tm(ia, :)))**2 )
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
  FLOAT :: xx(cas%n_pairs)

  R_TYPE, allocatable :: ff(:)
  R_TYPE, allocatable :: psii(:, :), psia(:, :)
  integer :: ip, ia, idim

  PUSH_SUB(X(ks_matrix_elements))

  SAFE_ALLOCATE(ff(1:mesh%np))
  SAFE_ALLOCATE(psii(1:mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psia(1:mesh%np, 1:st%d%dim))

  do ia = 1, cas%n_pairs
    call states_get_state(st, mesh, cas%pair(ia)%i, cas%pair(ia)%sigma, psii)
    call states_get_state(st, mesh, cas%pair(ia)%a, cas%pair(ia)%sigma, psia)

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
        zz = zz + xx(jb) * (M_ONE/sqrt(cas%s(jb))) * cas%mat(jb, ia)
      end do
      zz = (M_ONE/sqrt(cas%w(ia))) * zz
    else ! TAMM_DANCOFF, VARIATIONAL, PETERSILKA
      do jb = 1, cas%n_pairs
        zz = zz + xx(jb) * cas%mat(jb, ia)
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
  FLOAT,          intent(out) :: rho(:) !< FIXME: should be R_TYPE

  R_TYPE, allocatable :: psi_i(:), psi_a(:)
  integer :: ip

  PUSH_SUB(X(casida_get_rho))

  SAFE_ALLOCATE(psi_i(1:mesh%np))
  SAFE_ALLOCATE(psi_a(1:mesh%np))

  ! FIXME: need to take into account spinor dimension here, not just 1
  call states_get_state(st, mesh, 1, ii, sigma, psi_i)
  call states_get_state(st, mesh, 1, ia, sigma, psi_a)
  
  do ip = 1, mesh%np
    rho(ip) = R_CONJ(psi_i(ip))*psi_a(ip)
  end do
  
  SAFE_DEALLOCATE_A(psi_i)
  SAFE_DEALLOCATE_A(psi_a)

  POP_SUB(X(casida_get_rho))
end subroutine X(casida_get_rho)

! -----------------------------------------------------------------------------

!> one-particle matrix elements of perturbation
subroutine X(casida_lr_hmat1)(sys, hm, pert, hvar, lr_hmat1, st_start, st_end, ik)
  type(system_t),      intent(in)    :: sys
  type(hamiltonian_t), intent(inout) :: hm
  type(pert_t),        intent(in)    :: pert
  FLOAT,               intent(in)    :: hvar(:,:,:)
  R_TYPE,              intent(out)   :: lr_hmat1(:,:,:)
  integer,             intent(in)    :: st_start
  integer,             intent(in)    :: st_end
  integer,             intent(in)    :: ik

  integer :: ist, jst, ispin, idim
  R_TYPE, allocatable :: psi(:,:,:), pert_psi(:,:)

  PUSH_SUB(X(casida_lr_hmat1))

  SAFE_ALLOCATE(psi(1:sys%gr%mesh%np, 1:sys%st%d%dim, st_start:st_end))
  SAFE_ALLOCATE(pert_psi(1:sys%gr%mesh%np, 1:sys%st%d%dim))

  ! could use batches?

  ispin = states_dim_get_spin_index(sys%st%d, ik)
  
  do ist = st_start, st_end
    call states_get_state(sys%st, sys%gr%mesh, ist, ik, psi(:, :, ist))
  enddo
    
  do ist = st_start, st_end
    call X(pert_apply)(pert, sys%gr, sys%geo, hm, ik, psi(:, :, ist), pert_psi(:, :))
    do idim = 1, sys%st%d%dim
      pert_psi(:, idim) = pert_psi(:, idim) + hvar(:, ispin, 1) * psi(:, idim, ist)
    enddo
    
    do jst = ist, st_end
      lr_hmat1(jst, ist, ik) = X(mf_dotp)(sys%gr%mesh, sys%st%d%dim, psi(:, :, jst), pert_psi(:, :))
      if(jst /= ist) lr_hmat1(ist, jst, ik) = R_CONJ(lr_hmat1(jst, ist, ik)) ! Hermiticity
    enddo
  enddo

  POP_SUB(X(casida_lr_hmat1))
end subroutine X(casida_lr_hmat1)

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
!      not sure if this line is appropriate or not. --DAS
!      if(cas%type == CASIDA_PETERSILKA .and. isnt_degenerate(cas, st, ia, jb)) cycle

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

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
