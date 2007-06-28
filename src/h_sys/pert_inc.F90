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

    call magnetic

  case(PERTURBATION_IONIC)

    call ionic

  end select

contains 

  subroutine magnetic

    R_TYPE, allocatable :: f_in_copy(:), lf(:,:), vrnl(:,:,:)
    R_TYPE :: cross(1:MAX_DIM)
    integer :: iatom, idir, ip

    ALLOCATE(f_in_copy(1:NP_PART), NP_PART)

    call lalg_copy(NP_PART, f_in, f_in_copy)
    
    ALLOCATE(lf(gr%m%np, gr%sb%dim), gr%m%np*gr%sb%dim)

    ! Note that we leave out the term 1/P_c
    call X(f_angular_momentum) (gr%sb, gr%f_der, f_in_copy, lf)
    f_out(1:NP) = lf(1:NP, this%dir)/M_TWO

    deallocate(lf, f_in_copy)

    select case(this%gauge)

    case(GAUGE_GIPAW)

      ALLOCATE(vrnl(gr%m%np_part, h%d%dim, gr%sb%dim), gr%m%np_part*h%d%dim*gr%sb%dim)

      do iatom = 1, geo%natoms

        do idir = 1, gr%sb%dim
          call X(conmut_vnl_r)(gr, geo, h%ep, h%d%dim, idir, iatom, f_in, vrnl(:, :, idir), h%reltype, ik=1)
        end do

        do ip = 1, NP
          cross = X(cross_product)(R_TOTYPE(geo%atom(iatom)%x), vrnl(ip, 1, :))
#if !defined(R_TCOMPLEX)
          f_out(ip) = f_out(ip) + M_HALF * cross(this%dir)
#else
          f_out(ip) = f_out(ip) - M_zI * M_HALF * cross(this%dir)
#endif
        end do
        
      end do

      deallocate(vrnl)

    end select

  end subroutine magnetic

  subroutine ionic

    FLOAT, allocatable :: gv(:,:)
    integer :: iatom, idir

    iatom = this%atom1
    idir = this%dir

    ALLOCATE(gv(NP, NDIM), NP*NDIM)
    call specie_get_glocal( geo%atom(iatom)%spec, gr, geo%atom(iatom)%x, gv)
    f_out(1:NP) = gv(1:NP, idir) * f_in(1:NP)    

    deallocate(gv)

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
      f_out(ip) = CNST(0.25)*(rdelta - gr%m%x(ip, this%dir)*gr%m%x(ip, this%dir2)) * f_in(ip)
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

          f_out(ip) = f_out(ip) + CNST(0.25)*contr

        end do
        
      end do

      deallocate(vrnl)

    end select

  end subroutine magnetic

  subroutine ionic

    FLOAT,  allocatable :: g2v(:,:,:)
    type(atom_t), pointer :: atm
    integer :: iatom, idir, jatom, jdir

    iatom = this%atom1
    idir  = this%dir
    jatom = this%atom2
    jdir = this%dir2
    atm => geo%atom(iatom)

    if( jatom /= iatom ) then 
      f_out(1:NP) = R_TOTYPE(M_ZERO)
    else
      ALLOCATE(g2v(NP, NDIM, NDIM), NP*NDIM**2)
      call specie_get_g2local(atm%spec, gr, atm%x, g2v)
      f_out(1:NP) = g2v(1:NP, idir, jdir) * f_in(1:NP)
    end if

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
