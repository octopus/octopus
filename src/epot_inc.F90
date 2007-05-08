!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!------------------------------------------------------------------------------
! X(project) calculates the action of the sum of the projectors p(1:n_projectors)
! on the psi wavefunction. The result is summed up to ppsi
subroutine X(project)(mesh, p, n_projectors, dim, psi, ppsi, reltype, periodic, ik)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: p(:)
  integer,           intent(in)    :: n_projectors
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: ppsi(:, :)  ! ppsi(1:mesh%np, dim)
  integer,           intent(in)    :: reltype
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik

  integer :: n_s, ip, idim
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)

  call push_sub('epot_inc.project')

  do ip = 1, n_projectors

    n_s = p(ip)%sphere%ns

    if(p(ip)%type == M_LOCAL) then
      do idim = 1, dim
        ppsi(p(ip)%sphere%jxyz(1:n_s), idim) = ppsi(p(ip)%sphere%jxyz(1:n_s), idim) + &
             p(ip)%local_p%v(1:n_s) * psi(p(ip)%sphere%jxyz(1:n_s), idim)
      end do

      cycle

    end if

    if(allocated(lpsi))   deallocate(lpsi)
    if(allocated(plpsi))  deallocate(plpsi)
    ALLOCATE(lpsi(n_s, dim),  n_s*dim)
    ALLOCATE(plpsi(n_s, dim), n_s*dim)

    do idim = 1, dim
      lpsi(1:n_s, idim)  = psi(p(ip)%sphere%jxyz(1:n_s), idim)
      if(periodic) lpsi(1:n_s, idim)  = lpsi(1:n_s, idim) * p(ip)%phases(1:n_s, ik)
    end do

    select case (p(ip)%type)
    case (M_HGH)
      if (periodic) then
        call X(hgh_project)(mesh, p(ip)%sphere, p(ip)%hgh_p, dim, lpsi, plpsi, reltype, p(ip)%phases(:, ik))
      else
        call X(hgh_project)(mesh, p(ip)%sphere, p(ip)%hgh_p, dim, lpsi, plpsi, reltype)
      end if
    case (M_KB)
      if (periodic) then
        call X(kb_project)(mesh, p(ip)%sphere, p(ip)%kb_p, dim, lpsi, plpsi, p(ip)%phases(:, ik))
      else
        call X(kb_project)(mesh, p(ip)%sphere, p(ip)%kb_p, dim, lpsi, plpsi)
      end if
    case (M_RKB)
#ifdef R_TCOMPLEX
      !This can only be aplied to complex spinor wave-functions
      if (periodic) then
        call rkb_project(mesh, p(ip)%sphere, p(ip)%rkb_p, lpsi, plpsi, p(ip)%phases(:, ik))
      else
        call rkb_project(mesh, p(ip)%sphere, p(ip)%rkb_p, lpsi, plpsi)
      end if
#endif
    end select

    do idim = 1, dim
      ppsi(p(ip)%sphere%jxyz(1:n_s), idim) = ppsi(p(ip)%sphere%jxyz(1:n_s), idim) + plpsi(1:n_s, idim)
    end do

  end do

  if(allocated(lpsi))   deallocate(lpsi)
  if(allocated(plpsi))  deallocate(plpsi)

  call pop_sub()
end subroutine X(project)


!------------------------------------------------------------------------------
! X(psidprojectpsi) is used to calculate the contribution of the non-local part
! to the force acting on the ions.
!------------------------------------------------------------------------------
function X(psidprojectpsi)(mesh, p, n_projectors, dim, psi, periodic, ik) result(res)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: p(:)
  integer,           intent(in)    :: n_projectors
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psi(:, :)   ! psi(1:mesh%np, 1:dim)
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik
  R_TYPE :: res(3)


  integer ::  n_s, ip, idim, idir
  R_TYPE, allocatable :: lpsi(:, :)

  call push_sub('epot_inc.psidprojectpsi')

  res = R_TOTYPE(M_ZERO)

  ! index labels the atom
  do ip = 1, n_projectors
    if(p(ip)%type == M_LOCAL) cycle

    n_s = p(ip)%sphere%ns
    if(allocated(lpsi))  deallocate(lpsi)
    ALLOCATE( lpsi(n_s, dim), n_s*dim)
    
    do idim = 1, dim
      lpsi(1:n_s, idim)  = psi(p(ip)%sphere%jxyz(1:n_s), idim)
      if(periodic) lpsi(1:n_s, idim)  = lpsi(1:n_s, idim) * p(ip)%phases(1:n_s, ik)
    end do
    
    select case (p(ip)%type)
    case (M_HGH)
      if (periodic) then
        res = res + X(hgh_dproject)(mesh, p(ip)%sphere, p(ip)%hgh_p, dim, lpsi, p(ip)%phases(:, ik))
      else
        res = res + X(hgh_dproject)(mesh, p(ip)%sphere, p(ip)%hgh_p, dim, lpsi)
      end if
    case (M_KB)
      if (periodic) then
        res = res + X(kb_dproject)(mesh, p(ip)%sphere, p(ip)%kb_p, dim, lpsi, p(ip)%phases(:, ik))
      else
        res = res + X(kb_dproject)(mesh, p(ip)%sphere, p(ip)%kb_p, dim, lpsi)
      end if
    case (M_RKB)
#ifdef R_TCOMPLEX
      !This can only be aplied to complex spinor wave-functions
      if (periodic) then
        res = res + rkb_dproject(mesh, p(ip)%sphere, p(ip)%rkb_p, lpsi, p(ip)%phases(:, ik))
      else
        res = res + rkb_dproject(mesh, p(ip)%sphere, p(ip)%rkb_p, lpsi)
      end if
#endif
    end select

  end do

  if(allocated(lpsi)) deallocate(lpsi)

  call pop_sub()
end function X(psidprojectpsi)


R_TYPE function X(psia_project_psib)(mesh, pj, dim, psia, psib, reltype, periodic, ik) result(apb)
  type(mesh_t),      intent(in)    :: mesh
  type(projector_t), intent(in)    :: pj
  integer,           intent(in)    :: dim
  R_TYPE,            intent(in)    :: psia(:, :)  ! psia(1:mesh%np, dim)
  R_TYPE,            intent(inout) :: psib(:, :)  ! psib(1:mesh%np, dim)
  integer,           intent(in)    :: reltype
  logical,           intent(in)    :: periodic
  integer,           intent(in)    :: ik

  integer ::  n_s, idim
  R_TYPE, allocatable :: lpsi(:, :), plpsi(:,:)

  call push_sub('epot_inc.psia_project_psib')

  n_s = pj%sphere%ns

  ALLOCATE(lpsi(n_s, dim),  n_s*dim)
  ALLOCATE(plpsi(n_s, dim), n_s*dim)
  
  if(pj%type /= M_LOCAL) then
    do idim = 1, dim
      lpsi(1:n_s, idim)  = psib(pj%sphere%jxyz(1:n_s), idim)
      if(periodic) lpsi(1:n_s, idim)  = lpsi(1:n_s, idim) * pj%phases(1:n_s, ik)
    end do
  end if
  
  select case (pj%type)
  case (M_LOCAL)
    do idim = 1, dim
      plpsi(1:n_s, idim) = pj%local_p%v(1:n_s) * psib(pj%sphere%jxyz(1:n_s), idim)
    end do
  case (M_HGH)
    if (periodic) then
      call X(hgh_project)(mesh, pj%sphere, pj%hgh_p, dim, lpsi, plpsi, reltype, pj%phases(:, ik))
    else
      call X(hgh_project)(mesh, pj%sphere, pj%hgh_p, dim, lpsi, plpsi, reltype)
    end if
  case (M_KB)
    if (periodic) then
      call X(kb_project)(mesh, pj%sphere, pj%kb_p, dim, lpsi, plpsi, pj%phases(:, ik))
    else
      call X(kb_project)(mesh, pj%sphere, pj%kb_p, dim, lpsi, plpsi)
    end if
  case (M_RKB)
#ifdef R_TCOMPLEX
    !This can only be aplied to complex spinor wave-functions
    if (periodic) then
      call rkb_project(mesh, pj%sphere, pj%rkb_p, lpsi, plpsi, pj%phases(:, ik))
    else
      call rkb_project(mesh, pj%sphere, pj%rkb_p, lpsi, plpsi)
    end if
#endif
  end select
  
  apb = M_ZERO
  do idim = 1, dim
    plpsi(1:n_s, 1) = R_CONJ(psia(pj%sphere%jxyz(1:n_s), idim)) * plpsi(1:n_s, idim)
  end do
  apb = apb + X(sm_integrate)(mesh, pj%sphere, plpsi(1:n_s, 1))

  call pop_sub()
end function X(psia_project_psib)


subroutine X(calc_forces_from_potential)(gr, geo, ep, st, time)
  type(grid_t), target, intent(inout) :: gr
  type(geometry_t), intent(inout)  :: geo
  type(epot_t),     intent(in)     :: ep
  type(states_t),   intent(inout)     :: st
  FLOAT,            intent(in)     :: time

  integer :: ii, ip, ist, ik, ivnl, idim, idir, ns

  R_TYPE :: psi_proj_gpsi
  R_TYPE :: zz(MAX_DIM)
  R_TYPE, allocatable :: gpsi(:, :, :)
  FLOAT,  allocatable :: grho(:, :), vloc(:), force(:,:)

  type(atom_t), pointer :: atm

  ALLOCATE(force(1:NP, 1:NDIM), NP*NDIM)

  select case(ep%forces)
    
  case(DERIVATE_POTENTIAL)
    
    atm_loop: do ii = 1, geo%natoms
      atm => geo%atom(ii)

      !the local part

      if(.not.simul_box_is_periodic(gr%sb).or.geo%only_user_def) then !we do it in real space

        ns = min(2, st%d%nspin)

        call specie_get_glocal(atm%spec, gr, atm%x, force, time)

        do ip = 1, NP
          force(ip, 1:NDIM) = sum(st%rho(ip, 1:ns))*force(ip, 1:NDIM)
        end do

        do idir = 1, NDIM
          atm%f(idir) = atm%f(idir) - dmf_integrate(gr%m, force(:, idir))
        end do
        
      end if
      
    end do atm_loop

    ASSERT(NDIM == 3)
    
    ! iterate over the projectors
    do ivnl = 1, ep%nvnl

      atm => geo%atom(ep%p(ivnl)%iatom)

      ik_loop: do ik = 1, st%d%nik
        st_loop: do ist = st%st_start, st%st_end
          
          zz = X(psidprojectpsi)(gr%m, ep%p(ivnl:ivnl), 1, st%d%dim, st%X(psi)(:, :, ist, ik), &
               periodic = .false., ik = ik)
          atm%f(1:MAX_DIM) = atm%f(1:MAX_DIM) + M_TWO * R_REAL(st%occ(ist, ik) * zz(1:MAX_DIM))
          
        end do st_loop
      end do ik_loop

    end do

  case(DERIVATE_WAVEFUNCTION)

    ALLOCATE(gpsi(gr%m%np, 1:NDIM, st%d%dim), gr%m%np*NDIM*st%d%dim)
    ALLOCATE(grho(NP, MAX_DIM), NP*MAX_DIM)
    
    grho(1:NP, 1:NDIM) = M_ZERO

    !the non-local part
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        
        ! calculate the gradient of the wave-function
        do idim = 1, st%d%dim
          call X(f_gradient)(gr%sb, gr%f_der, st%X(psi)(:, idim, ist, ik), gpsi(:, :, idim))
          
          !accumulate to calculate the gradient of the density
          do idir = 1, NDIM
            grho(1:NP, idir) = grho(1:NP, idir) + st%d%kweights(ik)*st%occ(ist, ik) * M_TWO * &
                 R_REAL(st%X(psi)(1:NP, idim, ist, ik) * R_CONJ(gpsi(1:NP, idir, idim)))
          end do
        end do

        ! iterate over the projectors
        do ivnl = 1, ep%nvnl
          
          !get the atom corresponding to this projector
          atm => geo%atom(ep%p(ivnl)%iatom)

          do idir = 1, NDIM

            psi_proj_gpsi = X(psia_project_psib)(gr%m, ep%p(ivnl), st%d%dim, &
                 st%X(psi)(:, :, ist, ik), gpsi(:, idir, :), reltype = 0, periodic = .false., ik = ik)
            
            atm%f(idir) = atm%f(idir) - M_TWO * st%occ(ist, ik) * R_REAL(psi_proj_gpsi)

          end do
        
        end do !invl
        
      end do
    end do

    deallocate(gpsi)

    !now add the local part

    if(.not.simul_box_is_periodic(gr%sb).or.geo%only_user_def) then !we do it in real space

      ALLOCATE(vloc(1:NP), NP)
      
      do ii = 1, geo%natoms
        atm => geo%atom(ii)
        
        vloc(1:NP) = M_ZERO
        
        call build_local_part_in_real_space(ep, gr, geo, atm, vloc, time)
        
        do idir = 1, NDIM
          force(1:NP, idir) = grho(1:NP, idir) * vloc(1:NP)
          atm%f(idir) = atm%f(idir) - dmf_integrate(gr%m, force(:, idir))
        end do
      end do
      
      deallocate(vloc)

    end if
    
  end select
  
  deallocate(force)
  
end subroutine X(calc_forces_from_potential)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
