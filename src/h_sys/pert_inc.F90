!! Copyright (C) 2007 the octopus team
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
!! $Id: pert_inc.F90 2548 2006-11-06 21:42:27Z xavier $

! --------------------------------------------------------------------------
subroutine X(pert_apply) (this, gr, geo, h, ik, f_in, f_out)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: h
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:)
  R_TYPE,               intent(out)   :: f_out(:)

  call profiling_in(prof, "PERT_APPLY")

  ASSERT(this%dir.ne.-1)

  select case(this%pert_type)

  case(PERTURBATION_ELECTRIC)
    f_out(1:NP) = f_in(1:NP) * gr%m%x(1:NP, this%dir)

  case(PERTURBATION_MAGNETIC)
    call magnetic()

  case(PERTURBATION_IONIC)
    call ionic()

  case(PERTURBATION_KDOTP)
!    f_out(1:NP) = 0
    call kdotp()

  end select

  call profiling_out(prof)

contains 

  ! --------------------------------------------------------------------------
  subroutine kdotp()
    R_TYPE, allocatable :: f_in_copy(:), grad(:,:)

    ALLOCATE(f_in_copy(1:NP_PART), NP_PART)
    call lalg_copy(NP_PART, f_in, f_in_copy)
    
    ALLOCATE(grad(gr%m%np, gr%sb%dim), gr%m%np*gr%sb%dim)

    call X(derivatives_grad) (gr%f_der%der_discr, f_in_copy, grad) 

    f_out(1:NP) = - M_zI * (grad(1:NP, this%dir)) &
                  + h%d%kpoints(this%dir, ik) * f_in(1:NP)
!    delta_H = (-i*grad + k) . delta_k

!    write(*,*) 'pert_apply: ik = ', ik
!    write(*,*) 'pert_apply: kpoint = ', h%d%kpoints(:,ik)

    deallocate(f_in_copy, grad)
    
  end subroutine kdotp

  ! --------------------------------------------------------------------------
  subroutine magnetic()
    R_TYPE, allocatable :: f_in_copy(:), lf(:,:), vrnl(:,:,:)
    R_TYPE :: cross(1:MAX_DIM), vv(1:MAX_DIM)
    FLOAT :: xx(1:MAX_DIM)
    integer :: iatom, idir, ip

    ALLOCATE(f_in_copy(1:NP_PART), NP_PART)

    call lalg_copy(NP_PART, f_in, f_in_copy)
    
    ALLOCATE(lf(gr%m%np, gr%sb%dim), gr%m%np*gr%sb%dim)

    ! Note that we leave out the term 1/P_c
    call X(f_angular_momentum) (gr%sb, gr%f_der, f_in_copy, lf)
    f_out(1:NP) = M_HALF * lf(1:NP, this%dir)

    deallocate(lf, f_in_copy)

    if(this%gauge == GAUGE_GIPAW .or. this%gauge == GAUGE_ICL) then
      ALLOCATE(vrnl(NP, h%d%dim, gr%sb%dim), NP*h%d%dim*gr%sb%dim)
      vrnl(1:NP, 1:h%d%dim, this%dir) = M_ZERO

      do iatom = 1, geo%natoms

        do idir = 1, gr%sb%dim
          if(this%dir == idir) cycle ! this direction is not used in the cross product
          call X(projector_conmut_r)(h%ep%p(iatom), gr, h%d%dim, idir, ik, f_in, vrnl(:, :, idir))
        end do

        xx(1:MAX_DIM) = geo%atom(iatom)%x(1:MAX_DIM)

        do ip = 1, NP

          if(this%gauge == GAUGE_ICL) xx(1:MAX_DIM) = gr%m%x(ip, 1:MAX_DIM)
         
          vv(1:MAX_DIM) = vrnl(ip, 1, 1:MAX_DIM)

          cross(1) = xx(2)*vv(3) - xx(3)*vv(2)
          cross(2) = xx(3)*vv(1) - xx(1)*vv(3)
          cross(3) = xx(1)*vv(2) - xx(2)*vv(1)

#if !defined(R_TCOMPLEX)
          f_out(ip) = f_out(ip) + M_HALF*cross(this%dir)
#else
          f_out(ip) = f_out(ip) - M_zI*M_HALF*cross(this%dir)
#endif
        end do
      end do

      deallocate(vrnl)
    end if

  end subroutine magnetic

  subroutine ionic
    integer :: iatom, idir
    R_TYPE, allocatable  :: tmp(:)

    ALLOCATE(tmp(1:NP), NP*1)
    
    f_out(1:NP) = M_ZERO
    
    do iatom = 1, geo%natoms
      do idir = 1, NDIM

        if (this%ionic%pure_dir .and. iatom /= this%atom1 .and. idir /= this%dir) cycle

        call X(ionic_perturbation)(this, gr, geo, h, ik, f_in, tmp, iatom, idir)
        
        call lalg_axpy(NP, this%ionic%mix1(iatom, idir), tmp, f_out)

      end do
    end do
      
    deallocate(tmp)

  end subroutine ionic

end subroutine X(pert_apply)

  ! --------------------------------------------------------------------------
subroutine X(ionic_perturbation)(this, gr, geo, h, ik, f_in, f_out, iatom, idir)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: h
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:)
  R_TYPE,               intent(out)   :: f_out(:)
  integer,              intent(in)    :: iatom, idir    

  R_TYPE, allocatable :: grad(:,:), fin(:, :), fout(:, :)
  FLOAT,  allocatable :: vloc(:)
  type(atom_t), pointer :: atm

  atm => geo%atom(iatom)

  ALLOCATE(vloc(1:NP), NP)
  vloc(1:NP) = M_ZERO
  call epot_local_potential(h%ep, gr, gr%m, geo, atm, vloc, CNST(0.0))

  ALLOCATE(fin(1:NP_PART, 1), NP_PART)
  call lalg_copy(NP, f_in, fin(:, 1))    

  !d^T v |f>
  ALLOCATE(fout(1:NP_PART, 1), NP_PART)
  fout(1:NP, 1) = vloc(1:NP) * fin(1:NP, 1)
  call X(project_psi)(gr%m, h%ep%p(iatom:iatom), 1, 1, fin, fout, ik)
  call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, fout(:,1), f_out)

  !v d |f>
  ALLOCATE(grad(1:NP, 1), NP)
  call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, fin(:,1), grad(:,1))
  fout(1:NP, 1) = vloc(1:NP) * grad(1:NP, 1)
  call X(project_psi)(gr%m, h%ep%p(iatom:iatom), 1, 1, grad, fout, ik)
  f_out(1:NP) = -f_out(1:NP) + fout(1:NP, 1)

end subroutine X(ionic_perturbation)

! --------------------------------------------------------------------------
subroutine X(pert_apply_order_2) (this, gr, geo, h, ik, f_in, f_out)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: h
  integer,              intent(in)    :: ik
  R_TYPE,               intent(in)    :: f_in(:)
  R_TYPE,               intent(out)   :: f_out(:)

  select case(this%pert_type)

  case(PERTURBATION_ELECTRIC)
    f_out(1:NP) = R_TOTYPE(M_ZERO)
  case(PERTURBATION_IONIC)
    call ionic()
  case(PERTURBATION_MAGNETIC)
    call magnetic()
  end select

contains

  subroutine magnetic()
    R_TYPE, allocatable :: f_in2(:,:), dnl(:,:), vrnl(:,:), xf(:)
    R_TYPE :: cross1(1:MAX_DIM), bdir(1:MAX_DIM, 2)
    FLOAT  :: rdelta
    R_TYPE :: contr

    integer :: iatom, idir, idir2, ip

    do ip = 1, NP
      rdelta = sum(gr%m%x(ip, 1:MAX_DIM)**2) * ddelta(this%dir, this%dir2)
      f_out(ip) = M_FOURTH*(rdelta - gr%m%x(ip, this%dir)*gr%m%x(ip, this%dir2)) * f_in(ip)
    end do

    ! gauge correction
    apply_gauge: if(this%gauge==GAUGE_GIPAW .or. this%gauge==GAUGE_ICL) then
      bdir(1:MAX_DIM, 1:2) = M_ZERO
      bdir(this%dir,  1)   = M_ONE
      bdir(this%dir2, 2)   = M_ONE

      ALLOCATE(f_in2(NP_PART, NDIM), NP_PART*NDIM)
      ALLOCATE(vrnl(NP_PART, h%d%dim), NP_PART*h%d%dim)
      ALLOCATE(dnl(NP, NDIM), NP*NDIM)
      ALLOCATE(xf(NP), NP)

      f_in2(NP:NP_PART,:) = R_TOTYPE(M_ZERO)
      atoms: do iatom = 1, geo%natoms

        ! This calculates f_in2 = (B x r) f_in
        do ip = 1, NP
          select case(this%gauge)
          case(GAUGE_GIPAW)
            cross1 = X(cross_product)(bdir(:, 2), R_TOTYPE(geo%atom(iatom)%x))
          case(GAUGE_ICL)
            cross1 = X(cross_product)(bdir(:, 2), R_TOTYPE(gr%m%x(ip, :)))
          end select

          f_in2(ip, 1:gr%sb%dim) = cross1(1:gr%sb%dim) * f_in(ip)
        end do

        ! let us now get sum_beta Dnl f_in2
        dnl(1:NP, 1:NDIM) = R_TOTYPE(M_ZERO)
        do idir = 1, gr%sb%dim
          do idir2 = 1, gr%sb%dim
            !calculate dnl |f_in2> = -[x,vnl] |f_in2>
            call X(projector_conmut_r)(h%ep%p(iatom), gr, h%d%dim, idir2, ik, f_in2(:, idir2), vrnl(:, :))

            ! -x vnl |f>
            dnl(1:NP, idir) = dnl(1:NP, idir) - gr%m%x(1:NP, idir) * vrnl(1:NP, 1)

            ! vnl x |f>
            xf(1:NP) = gr%m%x(1:NP, idir) * f_in2(1:NP, idir2)
            call X(projector_conmut_r)(h%ep%p(iatom), gr, h%d%dim, idir2, ik, xf, vrnl(:, :))

            dnl(1:NP, idir) = dnl(1:NP, idir) + vrnl(1:NP, 1)
          end do
        end do

        do ip = 1, NP
          select case(this%gauge)
          case(GAUGE_GIPAW)
            cross1 = X(cross_product)(bdir(:, 1), R_TOTYPE(geo%atom(iatom)%x))
          case(GAUGE_ICL)
            cross1 = X(cross_product)(bdir(:, 1), R_TOTYPE(gr%m%x(ip, :)))
          end select

          contr = M_ZERO
          do idir = 1, gr%sb%dim
            contr = contr + cross1(idir)*dnl(ip, idir)
          end do
          f_out(ip) = f_out(ip) + M_FOURTH*contr
        end do

      end do atoms

      deallocate(f_in2, vrnl, dnl, xf)
    end if apply_gauge

  end subroutine magnetic

  subroutine ionic
    integer :: iatom, idir, jdir
    R_TYPE, allocatable  :: tmp(:)
    
    ALLOCATE(tmp(1:NP), NP*1)
    
    f_out(1:NP) = M_ZERO
    
    do iatom = 1, geo%natoms
      do idir = 1, NDIM
        do jdir = 1, NDIM
          
          if (this%ionic%pure_dir &
               .and. iatom /= this%atom1 .and. idir /= this%dir &
               .and. iatom /= this%atom2 .and. jdir /= this%dir2) cycle
          
          call X(ionic_perturbation_order_2)(this, gr, geo, h, ik, f_in, tmp, iatom, idir, jdir)
          
          call lalg_axpy(NP, this%ionic%mix1(iatom, idir) * this%ionic%mix2(iatom, jdir), tmp, f_out)
          
        end do
      end do
    end do

    deallocate(tmp)

  end subroutine ionic

end subroutine X(pert_apply_order_2)

subroutine X(ionic_perturbation_order_2) (this, gr, geo, h, ik, f_in, f_out, iatom, idir, jdir)
  type(pert_t),        intent(in)    :: this
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: h
  integer,             intent(in)    :: ik
  R_TYPE,              intent(in)    :: f_in(:)
  R_TYPE,              intent(out)   :: f_out(:)
  integer,             intent(in)    :: iatom, idir, jdir

  R_TYPE, allocatable :: fin(:, :)
  R_TYPE, allocatable :: tmp1(:, :), tmp2(:,:)
  FLOAT,  allocatable :: vloc(:)
  type(atom_t), pointer :: atm

  atm => geo%atom(iatom)

  ALLOCATE(fin(1:NP_PART, 1), NP_PART)
  ALLOCATE(tmp1(1:NP_PART, 1), NP_PART)
  ALLOCATE(tmp2(1:NP_PART, 1), NP_PART)
  ALLOCATE(vloc(1:NP), NP)

  vloc(1:NP) = M_ZERO
  call epot_local_potential(h%ep, gr, gr%m, geo, atm, vloc, CNST(0.0))

  call lalg_copy(NP, f_in, fin(:, 1))    

  !di^T dj^T v |f>
  tmp1(1:NP, 1) = vloc(1:NP) * fin(1:NP, 1)
  call X(project_psi)(gr%m, h%ep%p(iatom:iatom), 1, 1, fin, tmp1, ik)
  call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, tmp1(:,1), tmp2(:,1))
  call X(derivatives_oper)(gr%f_der%der_discr%grad(jdir), gr%f_der%der_discr, tmp2(:,1), f_out)

  !di^T v dj |f>
  call X(derivatives_oper)(gr%f_der%der_discr%grad(jdir), gr%f_der%der_discr, fin(:,1), tmp1(:,1))
  tmp2(1:NP, 1) = vloc(1:NP) * tmp1(1:NP, 1)
  call X(project_psi)(gr%m, h%ep%p(iatom:iatom), 1, 1, tmp1, tmp2, ik)
  call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, tmp2(:,1), tmp1(:,1))
  f_out(1:NP) = f_out(1:NP) - tmp1(1:NP, 1)

  !dj^T v di |f>
  call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, fin(:,1), tmp1(:,1))
  tmp2(1:NP, 1) = vloc(1:NP) * tmp1(1:NP, 1)
  call X(project_psi)(gr%m, h%ep%p(iatom:iatom), 1, 1, tmp1, tmp2, ik)
  call X(derivatives_oper)(gr%f_der%der_discr%grad(jdir), gr%f_der%der_discr, tmp2(:,1), tmp1(:,1))
  f_out(1:NP) = f_out(1:NP) - tmp1(1:NP, 1)

  !v di dj |f>
  call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, fin(:,1), tmp1(:,1))
  call X(derivatives_oper)(gr%f_der%der_discr%grad(jdir), gr%f_der%der_discr, tmp1(:,1), tmp2(:,1))
  tmp1(1:NP, 1) = vloc(1:NP) * tmp2(1:NP, 1)
  call X(project_psi)(gr%m, h%ep%p(iatom:iatom), 1, 1, tmp2, tmp1, ik)
  f_out(1:NP) = f_out(1:NP) + tmp1(1:NP, 1)

end subroutine X(ionic_perturbation_order_2)

subroutine X(pert_expectation_density) (this, gr, geo, h, st, psia, psib, density, pert_order)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: h
  type(states_t),       intent(in)    :: st
  R_TYPE,               pointer       :: psia(:, :, :, :)
  R_TYPE,               pointer       :: psib(:, :, :, :)
  R_TYPE,               intent(out)   :: density(:)
  integer, optional,    intent(in)    :: pert_order

  R_TYPE, allocatable :: pertpsib(:)

  integer :: ik, ist, idim, order

  ALLOCATE(pertpsib(1:NP), NP)

  order = 1
  if(present(pert_order)) order = pert_order
  ASSERT(order == 1 .or. order == 2)

  density(1:NP) = R_TOTYPE(M_ZERO)

  do ik = 1, st%d%nik
    do ist  = st%st_start, st%st_end
      do idim = 1, st%d%dim

        if(order == 1) then 
          call X(pert_apply) (this, gr, geo, h, ik, psib(:, idim, ist, ik), pertpsib)
        else
          call X(pert_apply_order_2) (this, gr, geo, h, ik, psib(:, idim, ist, ik), pertpsib)
        end if

        density(1:NP) = density(1:NP) + st%d%kweights(ik)*st%occ(ist, ik)*&
             R_CONJ(psia(1:NP, idim, ist, ik))*pertpsib(1:NP)

      end do
    end do
  end do

  deallocate(pertpsib)

end subroutine X(pert_expectation_density)


R_TYPE function X(pert_expectation_value) (this, gr, geo, h, st, psia, psib, pert_order) result(expval)
  type(pert_t),         intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: h
  type(states_t),       intent(in)    :: st
  R_TYPE,               pointer       :: psia(:, :, :, :)
  R_TYPE,               pointer       :: psib(:, :, :, :)
  integer, optional,    intent(in)    :: pert_order

  R_TYPE, allocatable :: density(:)
#ifdef HAVE_MPI
  R_TYPE :: expval_tmp
#endif
  integer :: order

  order = 1
  if(present(pert_order)) order = pert_order

  ASSERT(order == 1 .or. order == 2)

  ALLOCATE(density(1:NP), NP)

  call X(pert_expectation_density)(this, gr, geo, h, st, psia, psib, density, pert_order = order)

  expval = X(mf_integrate)(gr%m, density)

#ifdef HAVE_MPI
    ! reduce density
    if(st%parallel_in_states) then
        call MPI_Allreduce(expval, expval_tmp, 1, R_MPITYPE, MPI_SUM, st%mpi_grp%comm, mpi_err)
        expval = expval_tmp
    end if
#endif

  deallocate(density)
end function X(pert_expectation_value)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
