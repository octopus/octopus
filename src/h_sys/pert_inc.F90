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
!! $Id: response.F90 2548 2006-11-06 21:42:27Z xavier $

subroutine X(pert_apply) (this, gr, geo, h, f_in, f_out)
  type(pert_t), intent(in)    :: this
  type(grid_t),      intent(inout) :: gr
  type(geometry_t),  intent(in)    :: geo
  type(hamiltonian_t),  intent(in) :: h
  R_TYPE,            intent(in)    :: f_in(:)
  R_TYPE,            intent(out)   :: f_out(:)

  ASSERT(this%dir.ne.-1)

  select case(this%pert_type)

  case(PERTURBATION_ELECTRIC)
    f_out(1:NP) = f_in(1:NP) * gr%m%x(1:NP, this%dir)

  case(PERTURBATION_MAGNETIC)
    call magnetic()

  case(PERTURBATION_IONIC)
    call ionic()

  end select

contains 

  subroutine magnetic()

    R_TYPE, allocatable :: f_in_copy(:), lf(:,:), vrnl(:,:,:)
    R_TYPE :: cross(1:MAX_DIM)
    integer :: iatom, idir, ip

    ALLOCATE(f_in_copy(1:NP_PART), NP_PART)

    call lalg_copy(NP_PART, f_in, f_in_copy)
    
    ALLOCATE(lf(gr%m%np, gr%sb%dim), gr%m%np*gr%sb%dim)

    ! Note that we leave out the term 1/P_c
    call X(f_angular_momentum) (gr%sb, gr%f_der, f_in_copy, lf)
    f_out(1:NP) = M_HALF * lf(1:NP, this%dir)

    deallocate(lf, f_in_copy)

    if(this%gauge==GAUGE_GIPAW .or. this%gauge==GAUGE_ICL) then
      ALLOCATE(vrnl(gr%m%np_part, h%d%dim, gr%sb%dim), gr%m%np_part*h%d%dim*gr%sb%dim)

      do iatom = 1, geo%natoms
        do idir = 1, gr%sb%dim
          call X(conmut_vnl_r)(gr, geo, h%ep, h%d%dim, idir, iatom, f_in, vrnl(:, :, idir), h%reltype, ik=1)
        end do

        do ip = 1, NP
          select case(this%gauge)
          case(GAUGE_GIPAW)
            cross = X(cross_product)(R_TOTYPE(geo%atom(iatom)%x), vrnl(ip, 1, :))
          case(GAUGE_ICL)
            cross = X(cross_product)(R_TOTYPE(gr%m%x(ip, :)), vrnl(ip, 1, :))
          end select

#if !defined(R_TCOMPLEX)
          f_out(ip) = f_out(ip) +        M_HALF * cross(this%dir)
#else
          f_out(ip) = f_out(ip) - M_zI * M_HALF * cross(this%dir)
#endif
        end do
      end do

      deallocate(vrnl)
    end if

  end subroutine magnetic


  subroutine ionic

    R_TYPE, allocatable :: grad(:,:), fin(:, :), fout(:, :)
    FLOAT,  allocatable :: vloc(:)
    type(atom_t), pointer :: atm
    integer :: ipj

    integer :: iatom, idir

    iatom = this%atom1
    idir = this%dir
    atm => geo%atom(iatom)

    ALLOCATE(vloc(1:NP), NP_PART)
    vloc(1:NP) = M_ZERO
    call build_local_part_in_real_space(h%ep, gr, geo, atm, vloc, CNST(0.0))

    ALLOCATE(fin(1:NP_PART, 1), NP_PART)
    call lalg_copy(NP, f_in, fin(:, 1))    

    !d^T v |f>
    ALLOCATE(fout(1:NP_PART, 1), NP_PART)
    fout(1:NP, 1) = vloc(1:NP) * fin(1:NP, 1)
    do ipj = h%ep%atomproj(1, iatom), h%ep%atomproj(2, iatom)
      call X(project_psi)(gr%m, h%ep%p(ipj), 1, fin, fout, 0, .false., ik=1)
    end do
    call X(derivatives_oper)(this%gradt(idir), gr%f_der%der_discr, fout(:,1), f_out)

    !v d |f>
    ALLOCATE(grad(1:NP, 1), NP)
    call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, fin(:,1), grad(:,1))
    fout(1:NP, 1) = vloc(1:NP) * grad(1:NP, 1)
    do ipj = h%ep%atomproj(1, iatom), h%ep%atomproj(2, iatom)
      call X(project_psi)(gr%m, h%ep%p(ipj), 1, grad, fout, 0, .false., ik=1)
    end do
    f_out(1:NP) = f_out(1:NP) + fout(1:NP, 1)

  end subroutine ionic

end subroutine X(pert_apply)


subroutine X(pert_apply_order_2) (this, gr, geo, h, f_in, f_out)
  type(pert_t), intent(in)    :: this
  type(grid_t),      intent(inout) :: gr
  type(geometry_t),  intent(in)    :: geo
  type(hamiltonian_t),  intent(in) :: h
  R_TYPE,            intent(in)    :: f_in(:)
  R_TYPE,            intent(out)   :: f_out(:)

  select case(this%pert_type)

  case(PERTURBATION_ELECTRIC)
    f_out(1:NP) = R_TOTYPE(M_ZERO)

  case(PERTURBATION_IONIC)
    call ionic

  case(PERTURBATION_MAGNETIC)
    call magnetic

  end select

contains

  subroutine magnetic

    R_TYPE, allocatable :: dnl(:,:,:), vrnl(:,:,:), xf(:)
    R_TYPE :: cross1(1:MAX_DIM), cross2(1:MAX_DIM), bdir(1:MAX_DIM, 2)
    FLOAT  :: rdelta
    R_TYPE :: contr

    integer :: iatom, idir, idir2, ip

    do ip = 1, NP
      rdelta = sum(gr%m%x(ip, 1:MAX_DIM)**2) * ddelta(this%dir, this%dir2)
      f_out(ip) = M_FOURTH*(rdelta - gr%m%x(ip, this%dir)*gr%m%x(ip, this%dir2)) * f_in(ip)
    end do

    bdir(1:MAX_DIM, 1:2) = M_ZERO
    bdir(this%dir,  1)   = M_ONE
    bdir(this%dir2, 2)   = M_ONE

    select case(this%gauge)

    case(GAUGE_GIPAW)

      ALLOCATE(vrnl(gr%m%np_part, h%d%dim, gr%sb%dim), gr%m%np_part*h%d%dim*gr%sb%dim)
      ALLOCATE(dnl(gr%m%np_part, gr%sb%dim, gr%sb%dim), gr%m%np_part*gr%sb%dim*gr%sb%dim)
      ALLOCATE(xf(1:NP), NP)

      do iatom = 1, geo%natoms

        !calculate dnl |f> = -[x,vnl] |f>
        do idir = 1, gr%sb%dim
          call X(conmut_vnl_r)(gr, geo, h%ep, h%d%dim, idir, iatom, f_in, vrnl(:, :, idir), h%reltype, ik=1)
        end do

        ! -x vnl |f>
        do idir = 1, gr%sb%dim
          do idir2 = 1, gr%sb%dim
            dnl(1:NP, idir, idir2) = -gr%m%x(1:NP, idir) * vrnl(1:NP, 1, idir2)
          end do
        end do

        ! vnl x |f>
        do idir2 = 1, gr%sb%dim
          do idir = 1, gr%sb%dim
            xf(1:NP) = gr%m%x(1:NP, idir) * f_in(1:NP)
            call X(conmut_vnl_r)(gr, geo, h%ep, h%d%dim, idir2, iatom, xf, vrnl(:, :, idir2), h%reltype, ik=1)
            dnl(1:NP, idir, idir2) = dnl(1:NP, idir, idir2) + vrnl(1:NP, 1, idir2)
          end do
        end do
        
        cross1 = X(cross_product)(bdir(:, 1), R_TOTYPE(geo%atom(iatom)%x))
        cross2 = X(cross_product)(bdir(:, 2), R_TOTYPE(geo%atom(iatom)%x))

        do ip = 1, NP

          contr = M_ZERO
          do idir = 1, gr%sb%dim
            do idir2 = 1, gr%sb%dim
              contr = contr + cross1(idir) * dnl(ip, idir, idir2) * cross2(idir2)
            end do
          end do

          f_out(ip) = f_out(ip) + M_FOURTH*contr

        end do
    
      end do

      deallocate(vrnl)

    end select

  end subroutine magnetic

  subroutine ionic

    R_TYPE, allocatable :: fin(:, :)
    R_TYPE, allocatable :: tmp1(:, :), tmp2(:,:)
    FLOAT,  allocatable :: vloc(:)
    type(atom_t), pointer :: atm
    integer :: ipj

    integer :: iatom, idir, jatom, jdir

    iatom = this%atom1
    idir = this%dir
    jatom = this%atom2
    jdir = this%dir2
    atm => geo%atom(iatom)

    if (iatom /= jatom) then 
      f_out(1:NP) = R_TOTYPE(M_ZERO)
      return
    end if

    ALLOCATE(fin(1:NP_PART, 1), NP_PART)
    ALLOCATE(tmp1(1:NP_PART, 1), NP_PART)
    ALLOCATE(tmp2(1:NP_PART, 1), NP_PART)
    ALLOCATE(vloc(1:NP), NP)

    vloc(1:NP) = M_ZERO
    call build_local_part_in_real_space(h%ep, gr, geo, atm, vloc, CNST(0.0))

    call lalg_copy(NP, f_in, fin(:, 1))    

    !di^T dj^T v |f>
    tmp1(1:NP, 1) = vloc(1:NP) * fin(1:NP, 1)
    do ipj = h%ep%atomproj(1, iatom), h%ep%atomproj(2, iatom)
      call X(project_psi)(gr%m, h%ep%p(ipj), 1, fin, tmp1, 0, .false., ik=1)
    end do    
    call X(derivatives_oper)(this%gradt(idir), gr%f_der%der_discr, tmp1(:,1), tmp2(:,1))
    call X(derivatives_oper)(this%gradt(jdir), gr%f_der%der_discr, tmp2(:,1), f_out)

    !di^T v dj |f>
    call X(derivatives_oper)(gr%f_der%der_discr%grad(jdir), gr%f_der%der_discr, fin(:,1), tmp1(:,1))
    tmp2(1:NP, 1) = vloc(1:NP) * tmp1(1:NP, 1)
    do ipj = h%ep%atomproj(1, iatom), h%ep%atomproj(2, iatom)
      call X(project_psi)(gr%m, h%ep%p(ipj), 1, tmp1, tmp2, 0, .false., ik=1)
    end do    
    call X(derivatives_oper)(this%gradt(idir), gr%f_der%der_discr, tmp2(:,1), tmp1(:,1))
    f_out(1:NP) = f_out(1:NP) + tmp1(1:NP, 1)

    !dj^T v di |f>
    call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, fin(:,1), tmp1(:,1))
    tmp2(1:NP, 1) = vloc(1:NP) * tmp1(1:NP, 1)
    do ipj = h%ep%atomproj(1, iatom), h%ep%atomproj(2, iatom)
      call X(project_psi)(gr%m, h%ep%p(ipj), 1, tmp1, tmp2, 0, .false., ik=1)
    end do    
    call X(derivatives_oper)(this%gradt(jdir), gr%f_der%der_discr, tmp2(:,1), tmp1(:,1))
    f_out(1:NP) = f_out(1:NP) + tmp1(1:NP, 1)

    !v di dj |f>
    call X(derivatives_oper)(gr%f_der%der_discr%grad(idir), gr%f_der%der_discr, fin(:,1), tmp1(:,1))
    call X(derivatives_oper)(gr%f_der%der_discr%grad(jdir), gr%f_der%der_discr, tmp1(:,1), tmp2(:,1))
    tmp1(1:NP, 1) = vloc(1:NP) * tmp2(1:NP, 1)
    do ipj = h%ep%atomproj(1, iatom), h%ep%atomproj(2, iatom)
      call X(project_psi)(gr%m, h%ep%p(ipj), 1, tmp2, tmp1, 0, .false., ik=1)
    end do
    f_out(1:NP) = f_out(1:NP) + tmp1(1:NP, 1)

  end subroutine ionic

end subroutine X(pert_apply_order_2)


subroutine X(pert_expectation_density) (this, gr, geo, h, st, psia, psib, density, pert_order)
  type(pert_t),    intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: h
  type(states_t),       intent(in)    :: st
  R_TYPE,               pointer       :: psia(:, :, :, :)
  R_TYPE,               pointer       :: psib(:, :, :, :)
  R_TYPE,               intent(out)   :: density(:)
  integer, optional,    intent(in)    :: pert_order

  R_TYPE, allocatable :: tmp(:)

  integer :: ik, ist, idim, order

  ALLOCATE(tmp(1:NP), NP)

  order = 1
  if(present(pert_order)) order = pert_order
  ASSERT(order == 1 .or. order == 2)

  density(1:NP) = R_TOTYPE(M_ZERO)

  do ik = 1, st%d%nik
    do ist  = st%st_start, st%st_end
      do idim = 1, st%d%dim

        if(order == 1) then 
          call X(pert_apply) (this, gr, geo, h, psib(:, idim, ist, ik), tmp)
        else
          call X(pert_apply_order_2) (this, gr, geo, h, psib(:, idim, ist, ik), tmp)
        end if

        density(1:NP) = density(1:NP) + st%d%kweights(ik)*st%occ(ist, ik)*&
             R_CONJ(psia(1:NP, idim, ist, ik))*tmp(1:NP)

      end do
    end do
  end do

end subroutine X(pert_expectation_density)


R_TYPE function X(pert_expectation_value) (this, gr, geo, h, st, psia, psib, pert_order) result(expval)
  type(pert_t),    intent(in)    :: this
  type(grid_t),         intent(inout) :: gr
  type(geometry_t),     intent(in)    :: geo
  type(hamiltonian_t),  intent(in)    :: h
  type(states_t),       intent(in)    :: st
  R_TYPE,               pointer       :: psia(:, :, :, :)
  R_TYPE,               pointer       :: psib(:, :, :, :)
  integer, optional,    intent(in)    :: pert_order

  R_TYPE, allocatable :: tmp(:)

  integer :: order

  order = 1
  if(present(pert_order)) order = pert_order

  ASSERT(order == 1 .or. order == 2)

  ALLOCATE(tmp(1:NP), NP)

  call X(pert_expectation_density)(this, gr, geo, h, st, psia, psib, tmp, pert_order = order)

  expval = X(mf_integrate)(gr%m, tmp)

end function X(pert_expectation_value)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
