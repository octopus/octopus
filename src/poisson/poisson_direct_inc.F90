!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

!-----------------------------------------------------------------
subroutine poisson_solve_direct(this, pot, rho)
  type(poisson_t), intent(in)  :: this
  FLOAT,           intent(out) :: pot(:)
  FLOAT,           intent(in)  :: rho(:)

  FLOAT                :: prefactor, aa1, aa2, aa3, aa4, dressed_factor

  integer              :: ip, jp
  integer              :: dim_ele   !< physical dimensions  (= dimensions of simulation box      if electrons only)
                                    !<                      (= dimensions of simulation box - 1  for dressed electrons)
  integer              :: dim_eff   !< effective dimensions (= dim_ele,     if no soft coulomb)
                                    !<                      (= dim_ele + 1  if soft coulomb)

  FLOAT                :: xx1(1:MAX_DIM+1), xx2(1:MAX_DIM+1), xx3(1:MAX_DIM+1), xx4(1:MAX_DIM+1)
  FLOAT                :: xx(1:MAX_DIM+1), yy(1:MAX_DIM+1) 

  FLOAT                :: lam(1:MAX_DIM) ! for dressed orbitals
  FLOAT                :: qq, pp1, pp2, pp3, pp4 ! for dressed orbitals
  FLOAT                :: lx, lx1, lx2, lx3, lx4, ly

  logical              :: include_diag

  FLOAT                :: xg(MAX_DIM)
  FLOAT                :: pp                        ! for dressed orbitals
  integer, allocatable :: ip_v(:), part_v(:)
  FLOAT, allocatable   :: pvec(:), tmp(:)

  PUSH_SUB(poisson_solve_direct)

  xx  = M_ZERO
  xx1 = M_ZERO
  xx2 = M_ZERO
  xx3 = M_ZERO
  xx4 = M_ZERO
  yy  = M_ZERO
 
  if (.not. this%is_dressed) then
    dim_ele = this%der%mesh%sb%dim
  else
    dim_ele = this%der%mesh%sb%dim - 1
    lam(1:dim_ele) = this%dressed%lambda(1:dim_ele)
    dressed_factor = this%dressed%omega/sqrt(this%dressed%n_electrons)
  end if

  if (this%poisson_soft_coulomb_param**2 > M_ZERO) then
    dim_eff = dim_ele + 1
    xx(dim_eff) = this%poisson_soft_coulomb_param
    xx1(dim_eff) = this%poisson_soft_coulomb_param
    xx2(dim_eff) = this%poisson_soft_coulomb_param
    xx3(dim_eff) = this%poisson_soft_coulomb_param
    xx4(dim_eff) = this%poisson_soft_coulomb_param
    yy(dim_eff) = M_ZERO
    include_diag = .true.
  else
    dim_eff = dim_ele
    include_diag = .false.
  endif

  ASSERT(this%method == POISSON_DIRECT_SUM)

  select case(dim_ele)
  case(3)
    prefactor = M_TWO*M_PI*(M_THREE/(M_PI*M_FOUR))**(M_TWOTHIRD)
  case(2)
    prefactor = M_TWO*sqrt(M_PI)
  case(1)
    prefactor = M_ONE 
  case default
    message(1) = "Internal error: poisson_solve_direct can only be called for 1D, 2D or 3D."
    ! why not? all that is needed is the appropriate prefactors to be defined above, actually. then 1D, 4D etc. can be done
    call messages_fatal(1)
  end select

  if(.not. this%der%mesh%use_curvilinear) then
    prefactor = prefactor / (this%der%mesh%volume_element**(M_ONE/dim_ele))
  end if

  if(this%der%mesh%parallel_in_domains) then
    SAFE_ALLOCATE(pvec(1:this%der%mesh%np))
    SAFE_ALLOCATE(part_v(1:this%der%mesh%np_global))
    SAFE_ALLOCATE(ip_v(1:this%der%mesh%np_global))
    SAFE_ALLOCATE(tmp(1:this%der%mesh%np_global))
    do ip = 1, this%der%mesh%np_global
      ip_v(ip) = ip
    end do
    call partition_get_partition_number(this%der%mesh%inner_partition, this%der%mesh%np_global, ip_v, part_v)
    
    pot = M_ZERO
    do ip = 1, this%der%mesh%np_global
      xg = mesh_x_global(this%der%mesh, ip)
      xx(1:dim_ele) = xg(1:dim_ele)
      if (this%is_dressed) then
        lx = dot_product(lam(1:dim_ele), xx(1:dim_ele))
        pp = xg(dim_ele + 1)
      end if
      if(this%der%mesh%use_curvilinear) then
        do jp = 1, this%der%mesh%np
          if (.not. this%is_dressed) then
            if(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno) == jp .and. .not. include_diag) then
              pvec(jp) = rho(jp)*prefactor**(M_ONE - M_ONE/dim_ele)
            else
              yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
              pvec(jp) = rho(jp)/sqrt(sum((xx(1:dim_eff) - yy(1:dim_eff))**2))
            end if
          else !dressed orbitals
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            ly = dot_product(lam(1:dim_ele), yy(1:dim_ele))
            qq = this%der%mesh%x(jp, dim_ele + 1)
            if (all(yy(1:dim_eff) == xx(1:dim_eff))) then
              pvec(jp) = this%dressed%coulomb*rho(jp)*prefactor**(M_ONE - M_ONE/dim_ele) &
                + rho(jp)*(- dressed_factor*(lx*qq + ly*pp) + lx*ly)
            else
              pvec(jp) = this%dressed%coulomb*rho(jp)/sqrt(sum((xx(1:dim_eff) - yy(1:dim_eff))**2)) &
                + rho(jp)*(- dressed_factor*(lx*qq + ly*pp) + lx*ly)
            end if
          end if
        end do
      else
        do jp = 1, this%der%mesh%np
          if (.not. this%is_dressed) then
            if (vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno) == jp .and. .not. include_diag) then
              pvec(jp) = rho(jp)*prefactor 
            else
              yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
              pvec(jp) = rho(jp)/sqrt(sum((xx(1:dim_eff) - yy(1:dim_eff))**2))
            end if
          else ! dressed orbitals
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            ly = dot_product(lam(1:dim_ele), yy(1:dim_ele))
            qq = this%der%mesh%x(jp, dim_ele + 1)
            if (all(yy(1:dim_eff) == xx(1:dim_eff))) then
              pvec(jp) = this%dressed%coulomb*rho(jp)*prefactor &
                + rho(jp)*(- dressed_factor*(lx*qq + ly*pp) + lx*ly)
            else
              pvec(jp) = this%dressed%coulomb*rho(jp)/sqrt(sum((xx(1:dim_eff) - yy(1:dim_eff))**2)) &
                + rho(jp)*(- dressed_factor*(lx*qq + ly*pp) + lx*ly)
            end if
          end if
        end do
      end if
      tmp(ip) = dmf_integrate(this%der%mesh, pvec, reduce = .false.)
    end do

#ifdef HAVE_MPI
    call comm_allreduce(this%der%mesh%mpi_grp%comm, tmp)
#endif

    do ip = 1, this%der%mesh%np_global
      if (part_v(ip) == this%der%mesh%vp%partno) then
        pot(vec_global2local(this%der%mesh%vp, ip, this%der%mesh%vp%partno)) = tmp(ip)
      end if
    end do

    SAFE_DEALLOCATE_A(pvec)
    SAFE_DEALLOCATE_A(tmp)

  else ! serial mode
    
    do ip = 1, this%der%mesh%np - 4 + 1, 4

      xx1(1:dim_ele) = this%der%mesh%x(ip    , 1:dim_ele)
      xx2(1:dim_ele) = this%der%mesh%x(ip + 1, 1:dim_ele)
      xx3(1:dim_ele) = this%der%mesh%x(ip + 2, 1:dim_ele)
      xx4(1:dim_ele) = this%der%mesh%x(ip + 3, 1:dim_ele)
      if (this%is_dressed) then
        lx1 = dot_product(lam(1:dim_ele), xx1(1:dim_ele))
        lx2 = dot_product(lam(1:dim_ele), xx2(1:dim_ele))
        lx3 = dot_product(lam(1:dim_ele), xx3(1:dim_ele))
        lx4 = dot_product(lam(1:dim_ele), xx4(1:dim_ele))
        pp1 = this%der%mesh%x(ip    , dim_ele + 1)
        pp2 = this%der%mesh%x(ip + 1, dim_ele + 1)
        pp3 = this%der%mesh%x(ip + 2, dim_ele + 1)
        pp4 = this%der%mesh%x(ip + 3, dim_ele + 1)
      end if

      if (this%der%mesh%use_curvilinear) then
      
        if (.not. include_diag) then
          aa1 = prefactor*rho(ip    )*this%der%mesh%vol_pp(ip    )**(M_ONE - M_ONE/dim_ele)
          aa2 = prefactor*rho(ip + 1)*this%der%mesh%vol_pp(ip + 1)**(M_ONE - M_ONE/dim_ele)
          aa3 = prefactor*rho(ip + 2)*this%der%mesh%vol_pp(ip + 2)**(M_ONE - M_ONE/dim_ele)
          aa4 = prefactor*rho(ip + 3)*this%der%mesh%vol_pp(ip + 3)**(M_ONE - M_ONE/dim_ele)
        end if
 
        if (.not. this%is_dressed) then
          !$omp parallel do reduction(+:aa1,aa2,aa3,aa4)
          do jp = 1, this%der%mesh%np
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            if (ip     /= jp .or. include_diag) then
              aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim_eff) - yy(1:dim_eff))**2)) * this%der%mesh%vol_pp(jp)
            end if
            if (ip + 1 /= jp .or. include_diag) then
              aa2 = aa2 + rho(jp)/sqrt(sum((xx2(1:dim_eff) - yy(1:dim_eff))**2)) * this%der%mesh%vol_pp(jp)
            end if
            if (ip + 2 /= jp .or. include_diag) then
              aa3 = aa3 + rho(jp)/sqrt(sum((xx3(1:dim_eff) - yy(1:dim_eff))**2)) * this%der%mesh%vol_pp(jp)
            end if
            if (ip + 3 /= jp .or. include_diag) then
              aa4 = aa4 + rho(jp)/sqrt(sum((xx4(1:dim_eff) - yy(1:dim_eff))**2)) * this%der%mesh%vol_pp(jp)
            end if
          end do
        else !dressed orbitals
          !$omp parallel do reduction(+:aa1,aa2,aa3,aa4)
          do jp = 1, this%der%mesh%np
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            lx1 = dot_product(lam(1:dim_ele), xx1(1:dim_ele))
            qq = this%der%mesh%x(jp, dim_ele + 1)
            ! we do not need the include_diag flag here, because yy and xx1 will differ in the extra dimension 
            ! if include_diag == .true.
            if (.not. all(yy(1:dim_eff) == xx1(1:dim_eff))) then
              aa1 = aa1 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx1(1:dim_eff) - yy(1:dim_eff))**2))*this%der%mesh%vol_pp(jp) &
                + rho(jp)*this%der%mesh%vol_pp(jp)*(- dressed_factor*(lx1*qq + ly*pp1) + lx1*ly)
            end if
            if (.not. all(yy(1:dim_eff) == xx2(1:dim_eff))) then
              aa2 = aa2 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx2(1:dim_eff) - yy(1:dim_eff))**2))*this%der%mesh%vol_pp(jp) &
                + rho(jp)*this%der%mesh%vol_pp(jp)*(- dressed_factor*(lx2*qq + ly*pp2) + lx2*ly)
            end if
            if (.not. all(yy(1:dim_eff) == xx3(1:dim_eff))) then
              aa3 = aa3 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx3(1:dim_eff) - yy(1:dim_eff))**2))*this%der%mesh%vol_pp(jp) &
                + rho(jp)*this%der%mesh%vol_pp(jp)*(- dressed_factor*(lx3*qq + ly*pp3) + lx3*ly)
            end if
            if (.not. all(yy(1:dim_eff) == xx4(1:dim_eff))) then
              aa4 = aa4 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx4(1:dim_eff) - yy(1:dim_eff))**2))*this%der%mesh%vol_pp(jp) &
                + rho(jp)*this%der%mesh%vol_pp(jp)*(- dressed_factor*(lx4*qq + ly*pp4) + lx4*ly)
            end if
          end do
        end if
        
      else

        if(.not. include_diag) then
          aa1 = prefactor*rho(ip    )
          aa2 = prefactor*rho(ip + 1)
          aa3 = prefactor*rho(ip + 2)
          aa4 = prefactor*rho(ip + 3)
        else
          aa1 = M_ZERO
          aa2 = M_ZERO
          aa3 = M_ZERO
          aa4 = M_ZERO
        end if
        
        if (.not. this%is_dressed) then
          !$omp parallel do reduction(+:aa1,aa2,aa3,aa4)
          do jp = 1, this%der%mesh%np
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            if (ip     /= jp .or. include_diag) aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim_eff) - yy(1:dim_eff))**2))
            if (ip + 1 /= jp .or. include_diag) aa2 = aa2 + rho(jp)/sqrt(sum((xx2(1:dim_eff) - yy(1:dim_eff))**2))
            if (ip + 2 /= jp .or. include_diag) aa3 = aa3 + rho(jp)/sqrt(sum((xx3(1:dim_eff) - yy(1:dim_eff))**2))
            if (ip + 3 /= jp .or. include_diag) aa4 = aa4 + rho(jp)/sqrt(sum((xx4(1:dim_eff) - yy(1:dim_eff))**2))
          end do
        else !dressed orbitals
          !$omp parallel do reduction(+:aa1,aa2,aa3,aa4)
          do jp = 1, this%der%mesh%np
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            ly = dot_product(lam(1:dim_ele), yy(1:dim_ele))
            qq = this%der%mesh%x(jp, dim_ele + 1)
            if (.not. all(yy(1:dim_eff) == xx1(1:dim_eff))) then
              aa1 = aa1 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx1(1:dim_eff) - yy(1:dim_eff))**2)) &
                + rho(jp)*(- dressed_factor*(lx1*qq + ly*pp1) + lx1*ly)
            end if
            if (.not. all(yy(1:dim_eff) == xx2(1:dim_eff))) then
              aa2 = aa2 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx2(1:dim_eff) - yy(1:dim_eff))**2)) &
                + rho(jp)*(- dressed_factor*(lx2*qq + ly*pp2) + lx2*ly)
            end if
            if (.not. all(yy(1:dim_eff) == xx3(1:dim_eff))) then
              aa3 = aa3 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx3(1:dim_eff) - yy(1:dim_eff))**2)) &
                + rho(jp)*(- dressed_factor*(lx3*qq + ly*pp3) + lx3*ly)
            end if
            if (.not. all(yy(1:dim_eff) == xx4(1:dim_eff))) then
              aa4 = aa4 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx4(1:dim_eff) - yy(1:dim_eff))**2)) + &
                + rho(jp)*(- dressed_factor*(lx4*qq + ly*pp4) + lx4*ly)
            end if
          end do
        end if
      end if

      pot(ip    ) = this%der%mesh%volume_element*aa1
      pot(ip + 1) = this%der%mesh%volume_element*aa2
      pot(ip + 2) = this%der%mesh%volume_element*aa3
      pot(ip + 3) = this%der%mesh%volume_element*aa4
      
    end do
      
    do ip = ip, this%der%mesh%np

      aa1 = CNST(0.0)
      
      xx1(1:dim_ele) = this%der%mesh%x(ip,1:dim_ele)
      lx1 = dot_product(lam(1:dim_ele), xx1(1:dim_ele))
      if (this%is_dressed) then
        pp1 = this%der%mesh%x(ip, dim_ele + 1)
      end if

      if (this%der%mesh%use_curvilinear) then
        if (.not. this%is_dressed) then
          do jp = 1, this%der%mesh%np
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            if (ip == jp .and. .not. include_diag) then
              aa1 = aa1 + prefactor*rho(ip)*this%der%mesh%vol_pp(jp)**(M_ONE - M_ONE/dim_ele)
            else
              aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim_eff) - yy(1:dim_eff))**2))*this%der%mesh%vol_pp(jp)
            end if
          end do
        else ! dressed orbitals
          do jp = 1, this%der%mesh%np
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            ly = dot_product(lam(1:dim_ele), yy(1:dim_ele))
            qq = this%der%mesh%x(jp, dim_ele + 1)
            if (all( xx1(1:dim_eff) == yy(1:dim_eff))) then
              aa1 = aa1 + this%dressed%coulomb*prefactor*rho(ip)*this%der%mesh%vol_pp(jp)**(M_ONE - M_ONE/dim_ele) &
                + rho(jp)*this%der%mesh%vol_pp(jp)*(- dressed_factor*(lx1*qq + ly*pp1) + lx1*ly)
            else
              aa1 = aa1 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx1(1:dim_eff) - yy(1:dim_eff))**2)) *this%der%mesh%vol_pp(jp) &
                + rho(jp)*this%der%mesh%vol_pp(jp)*(- dressed_factor*(lx1*qq + ly*pp1) + lx1*ly)
            end if
          end do
        end if
      else
        if (.not. this%is_dressed) then
          do jp = 1, this%der%mesh%np
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            if (ip == jp .and. .not. include_diag) then
              aa1 = aa1 + prefactor*rho(ip)
            else
              aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim_eff) - yy(1:dim_eff))**2))
            end if
          end do
        else ! dressed orbitals
          do jp = 1, this%der%mesh%np
            yy(1:dim_ele) = this%der%mesh%x(jp, 1:dim_ele)
            ly = dot_product(lam(1:dim_ele), yy(1:dim_ele))
            qq = this%der%mesh%x(jp, dim_ele + 1)
            if (all( xx1(1:dim_eff) == yy(1:dim_eff))) then
              aa1 = aa1 + this%dressed%coulomb*prefactor*rho(ip) &
                + rho(jp)*(- dressed_factor*(lx1*qq + ly*pp1) + lx1*ly)
            else
              aa1 = aa1 + this%dressed%coulomb*rho(jp)/sqrt(sum((xx1(1:dim_eff) - yy(1:dim_eff))**2)) &
                + rho(jp)*(- dressed_factor*(lx1*qq + ly*pp1) + lx1*ly)
            end if
          end do
        end if
      end if 

      pot(ip) = this%der%mesh%volume_element*aa1
      
    end do
    
  end if

  POP_SUB(poisson_solve_direct) 
end subroutine poisson_solve_direct

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
