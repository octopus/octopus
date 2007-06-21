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

subroutine X(resp_pert_apply) (this, gr, geo, h, f_in, f_out)
  type(resp_pert_t), intent(in)    :: this
  type(grid_t),      intent(inout) :: gr
  type(geometry_t),  intent(in)    :: geo
  type(hamiltonian_t),  intent(in) :: h
  R_TYPE,            intent(inout) :: f_in(:)
  R_TYPE,            intent(out)   :: f_out(:)

  FLOAT, allocatable :: gv(:,:)
  type(atom_t), pointer :: atm

  R_TYPE :: cross(1:MAX_DIM)

  integer :: ipj, iatom, idir, ip

  ASSERT(this%dir.ne.-1)

  select case(this%resp_type)
  case(RESP_PERTURBATION_ELECTRIC)
    f_out(1:NP) = f_in(1:NP) * gr%m%x(1:NP, this%dir)

  case(RESP_PERTURBATION_MAGNETIC)

    call magnetic

  case(RESP_PERTURBATION_DISPLACE)

    ALLOCATE(gv(NP, NDIM), NP*NDIM)

    atm => geo%atom(this%ion_disp%iatom)

    call specie_get_glocal(atm%spec, gr, atm%x, gv)

    f_out(1:NP) = gv(1:NP, this%ion_disp%idir) * f_in(1:NP)

    deallocate(gv)

  end select

contains 

  subroutine magnetic

    R_TYPE, allocatable :: lf(:,:), vrnl(:,:,:)

    ! I believe this should be NP and not NP_PART
    ALLOCATE(lf(gr%m%np_part, gr%sb%dim), gr%m%np_part*gr%sb%dim)

    ! Note that we leave out the term 1/P_c
    call X(f_angular_momentum) (gr%sb, gr%f_der, f_in, lf)
    f_out(1:NP) = lf(1:NP, this%dir)/M_TWO

    deallocate(lf)

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

end subroutine X(resp_pert_apply)
