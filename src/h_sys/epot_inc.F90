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

subroutine X(calc_forces_from_potential)(gr, geo, ep, st, time)
  type(grid_t),     intent(inout)  :: gr
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

    !$omp parallel workshare
    grho(1:NP, 1:NDIM) = M_ZERO
    !$omp end parallel workshare    

    !the non-local part
    do ik = 1, st%d%nik
      do ist = st%st_start, st%st_end
        
        ! calculate the gradient of the wave-function
        do idim = 1, st%d%dim
          call X(f_gradient)(gr%sb, gr%f_der, st%X(psi)(:, idim, ist, ik), gpsi(:, :, idim))
          
          !accumulate to calculate the gradient of the density
          do idir = 1, NDIM
            !$omp parallel workshare
            grho(1:NP, idir) = grho(1:NP, idir) + st%d%kweights(ik)*st%occ(ist, ik) * M_TWO * &
                 R_REAL(st%X(psi)(1:NP, idim, ist, ik) * R_CONJ(gpsi(1:NP, idir, idim)))
            !$omp end parallel workshare
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

        !$omp parallel workshare
        vloc(1:NP) = M_ZERO
        !$omp end parallel workshare

        call build_local_part_in_real_space(ep, gr, geo, atm, vloc, time)
        
        do idir = 1, NDIM
          !$omp parallel workshare
          force(1:NP, idir) = grho(1:NP, idir) * vloc(1:NP)
          !$omp end parallel workshare
          atm%f(idir) = atm%f(idir) - dmf_integrate(gr%m, force(:, idir))
        end do
      end do
      
      deallocate(vloc)

    end if
    
  end select
  
  deallocate(force)
  
end subroutine X(calc_forces_from_potential)

!This function calculates |cpsi> = [x,V_nl] |psi>

subroutine X(conmut_vnl_r)(gr, geo, ep, dim, idir, iatom, psi, cpsi, reltype, ik)
  type(grid_t),          intent(in)     :: gr
  type(geometry_t),      intent(in)     :: geo
  type(epot_t), target,  intent(in)     :: ep
  integer,               intent(in)     :: dim
  integer,               intent(in)     :: idir
  integer,               intent(in)     :: iatom
  R_TYPE,                intent(in)     :: psi(:)
  R_TYPE,                intent(out)    :: cpsi(:,:)
  integer,               intent(in)     :: reltype
  integer,               intent(in)     :: ik

  integer ::  n_s, idim, ipj
  R_TYPE, allocatable :: lpsi(:, :), pxlpsi(:,:), xplpsi(:,:)
  logical :: periodic
  integer, pointer :: jxyz(:)

  call push_sub('epot_inc.Xconmut_vnl_r')

  periodic = simul_box_is_periodic(gr%m%sb)

  cpsi(1:gr%m%np, 1:dim) = M_ZERO

  do ipj = ep%atomproj(1, iatom), ep%atomproj(2, iatom)
    
    if (ep%p(ipj)%type == M_LOCAL ) cycle
    
    n_s = ep%p(ipj)%sphere%ns
    jxyz => ep%p(ipj)%sphere%jxyz
    
    ALLOCATE(lpsi(n_s, dim),  n_s*dim)
    ALLOCATE(xplpsi(n_s, dim), n_s*dim)
    ALLOCATE(pxlpsi(n_s, dim), n_s*dim)
    
    do idim = 1, dim
      lpsi(1:n_s, idim) = psi(jxyz(1:n_s))
    end do

    ! x V_nl |psi>
    call X(project_sphere)(gr%m, ep%p(ipj), dim, lpsi, xplpsi, reltype, periodic, ik)
    do idim = 1, dim
      xplpsi(1:n_s, idim) = gr%m%x(jxyz(1:n_s), idir) * xplpsi(1:n_s, idim)
    end do
    
    ! V_nl x |psi>
    do idim = 1, dim
      lpsi(1:n_s, idim) = gr%m%x(jxyz(1:n_s), idir) * lpsi(1:n_s, idim)
    end do
    call X(project_sphere)(gr%m, ep%p(ipj), dim, lpsi, pxlpsi, reltype, periodic, ik)

    ! |cpsi> = x V_nl |psi> - V_nl x |psi> 
    do idim = 1, dim
      cpsi(jxyz(1:n_s), idim) = cpsi(jxyz(1:n_s), idim) + xplpsi(1:n_s, idim) - pxlpsi(1:n_s, idim)
    end do

    deallocate(lpsi, xplpsi, pxlpsi)
  end do
  
  call pop_sub()
  
end subroutine X(conmut_vnl_r)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
