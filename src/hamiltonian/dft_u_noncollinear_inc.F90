!! Copyright (C) 2017 N. Tancogne-Dejean 
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

! ---------------------------------------------------------
! TODO: Merge this with the two_body routine in system/output_me_inc.F90
subroutine compute_complex_coulomb_integrals (this, mesh, der, st, psolver, namespace, space)
  type(lda_u_t),       intent(inout) :: this
  type(mesh_t),        intent(in)    :: mesh
  type(derivatives_t), intent(in)    :: der
  type(states_elec_t), intent(in)    :: st
  type(poisson_t),     intent(in)    :: psolver
  type(namespace_t),   intent(in)    :: namespace
  type(space_t),       intent(in)    :: space

  integer :: ist, jst, kst, lst, ijst, klst
  integer :: is1, is2
  integer :: norbs, np_sphere, ios, ip
  integer :: idone, ntodo
  CMPLX, allocatable :: tmp(:), vv(:,:), nn(:,:)
  type(orbitalset_t), pointer :: os
  type(profile_t), save :: prof

  call profiling_in(prof, TOSTRING(X(DFTU_COMPEX_COULOMB_INTEGRALS)))

  PUSH_SUB(compute_complex_coulomb_integrals)

  ASSERT(.not. st%parallel_in_states)
  if(mesh%parallel_in_domains) then
    call messages_not_implemented("Coulomb integrals parallel in domains", namespace=namespace)
  end if

  SAFE_ALLOCATE(nn(1:this%max_np,st%d%dim))
  SAFE_ALLOCATE(vv(1:this%max_np,st%d%dim))
  SAFE_ALLOCATE(tmp(1:this%max_np))
  SAFE_ALLOCATE(this%zcoulomb(1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:this%maxnorbs,1:st%d%dim,1:st%d%dim,1:this%norbsets))
  this%zcoulomb(1:this%maxnorbs, 1:this%maxnorbs, 1:this%maxnorbs, 1:this%maxnorbs, &
                1:st%d%dim, 1:st%d%dim, 1:this%norbsets) = M_ZERO

  !Lets counts the number of orbital to treat, to display a progress bar
  ntodo = 0
  do ios = this%orbs_dist%start, this%orbs_dist%end
    norbs = this%orbsets(ios)%norbs
    ntodo= ntodo + norbs**4*4!((norbs+1)*norbs/2)*((norbs+1)*norbs/2+1)/2
  end do
  idone = 0
  if(mpi_world%rank == 0) call loct_progress_bar(-1, ntodo)


  do ios = this%orbs_dist%start, this%orbs_dist%end
    os => this%orbsets(ios)
    norbs = os%norbs
    np_sphere = os%sphere%np

    call submesh_build_global(os%sphere, space)

    select case (this%sm_poisson)
    case(SM_POISSON_DIRECT)
      call poisson_init_sm(os%poisson, namespace, space, psolver, der, os%sphere, method = POISSON_DIRECT_SUM)
    case(SM_POISSON_ISF)
      call poisson_init_sm(os%poisson, namespace, space, psolver, der, os%sphere, method = POISSON_ISF)
    case(SM_POISSON_PSOLVER)
      call poisson_init_sm(os%poisson, namespace, space, psolver, der, os%sphere, method = POISSON_PSOLVER)
    case(SM_POISSON_FFT)
      call poisson_init_sm(os%poisson, namespace, space, psolver, der, os%sphere, method = POISSON_FFT, force_cmplx=.true.)
    end select

    ijst=0
    do ist = 1, norbs

      do jst = 1, norbs
        ijst=ijst+1

        do is1 = 1, st%d%dim
          !$omp parallel do
          do ip = 1,np_sphere
            nn(ip,is1)  = conjg(os%zorb(ip,is1,ist))*os%zorb(ip,is1,jst)
          end do
          !$omp end parallel do    

          !Here it is important to use a non-periodic poisson solver, e.g. the direct solver
          call zpoisson_solve_sm(os%poisson, os%sphere, vv(1:np_sphere,is1), nn(1:np_sphere,is1))
        end do !is1  

        klst=0
        do kst = 1, norbs
          do lst = 1, norbs
            klst=klst+1

            do is1 = 1, st%d%dim
              do is2 = 1, st%d%dim

                !$omp parallel do
                do ip = 1,np_sphere
                 tmp(ip) = vv(ip,is1)*conjg(os%zorb(ip,is2,kst))*os%zorb(ip,is2,lst)
                end do
                !$omp end parallel do

                this%zcoulomb(ist,jst,kst,lst,is1,is2,ios) = zsm_integrate(mesh, os%sphere, tmp(1:np_sphere))
              end do !is2
            end do !is1

            do is1 = 1, st%d%dim
              do is2 = 1, st%d%dim
                if(abs(this%zcoulomb(ist,jst,kst,lst,is1,is2,ios))<CNST(1.0e-12)) then
                  this%zcoulomb(ist,jst,kst,lst,is1,is2,ios) = M_ZERO
                end if

              end do !is2
            end do !is1

            idone = idone + 1
            if(mpi_world%rank == 0) call loct_progress_bar(idone, ntodo)
          end do !lst
        end do !kst
      end do !jst
    end do !ist
    call poisson_end(os%poisson)
    call submesh_end_cube_map(os%sphere)
    call submesh_end_global(os%sphere)
  end do !iorb

  if(this%orbs_dist%parallel) then
    do ios = 1, this%norbsets
      do is2 = 1, st%d%dim
        do is1 = 1, st%d%dim
          call comm_allreduce(this%orbs_dist%mpi_grp, this%zcoulomb(:,:,:,:,is1,is2,ios))
        end do
      end do
    end do
  end if

  SAFE_DEALLOCATE_A(nn)
  SAFE_DEALLOCATE_A(vv)
  SAFE_DEALLOCATE_A(tmp)

  POP_SUB(compute_complex_coulomb_integrals)
  call profiling_out(prof)
end subroutine compute_complex_coulomb_integrals

! ---------------------------------------------------------
!> This routine computes the effective U in the non-collinear case 
! ---------------------------------------------------------
subroutine compute_ACBNO_U_noncollinear(this, ios, namespace)
  type(lda_u_t),     intent(inout) :: this
  integer,           intent(in)    :: ios
  type(namespace_t), intent(in)    :: namespace

  integer :: im, imp, impp, imppp, ispin1, ispin2, norbs
  CMPLX   :: numU, numJ, tmpU, tmpJ, denomU, denomJ

  PUSH_SUB(compute_ACBNO_U_noncollinear)

  norbs = this%orbsets(ios)%norbs
  numU = M_z0
  numJ = M_z0
  denomU = M_z0
  denomJ = M_z0

  if(norbs > 1) then

    do im = 1, norbs
    do imp = 1,norbs
      do impp = 1, norbs
      do imppp = 1, norbs
        ! We first compute the terms
        ! sum_{alpha,beta} P^alpha_{mmp}P^beta_{mpp,mppp}  
        ! sum_{alpha} P^alpha_{mmp}P^alpha_{mpp,mppp}
        tmpU = M_z0
        tmpJ = M_z0

        do ispin1 = 1, this%spin_channels
          do ispin2 = 1, this%spin_channels
            tmpU = tmpU + this%zn_alt(im,imp,ispin1,ios)*this%zn_alt(impp,imppp,ispin2,ios)&
                         * this%zcoulomb(im,imp,impp,imppp,ispin1,ispin2,ios)
          end do
          tmpJ = tmpJ + this%zn_alt(im,imp,ispin1,ios)*this%zn_alt(impp,imppp,ispin1,ios)&
                        * this%zcoulomb(im,imppp,impp,imp,ispin1,ispin1,ios)
        end do
        tmpJ = tmpJ + this%zn_alt(im,imp,3,ios)*this%zn_alt(impp,imppp,4,ios) &
                           * this%zcoulomb(im,imppp,impp,imp,1,2,ios)                   &
                          +this%zn_alt(im,imp,4,ios)*this%zn_alt(impp,imppp,3,ios) &
                           * this%zcoulomb(im,imppp,impp,imp,2,1,ios)
        ! These are the numerator of the ACBN0 U and J
        numU = numU + tmpU
        numJ = numJ + tmpJ
      end do
      end do

      ! We compute the term
      ! sum_{alpha} sum_{m,mp/=m} N^alpha_{m}N^alpha_{mp}
      tmpJ = M_z0
      tmpU = M_z0
      if(imp/=im) then
        do ispin1 = 1, this%spin_channels
          tmpJ = tmpJ + this%zn(im,im,ispin1,ios)*this%zn(imp,imp,ispin1,ios)
          tmpU = tmpU + this%zn(im,im,ispin1,ios)*this%zn(imp,imp,ispin1,ios)
        end do
        tmpJ = tmpJ + this%zn(im,im,3,ios)*this%zn(imp,imp,4,ios) &
                    + this%zn(im,im,4,ios)*this%zn(imp,imp,3,ios)
      end if
      denomJ = denomJ + tmpJ

      ! We compute the term
      ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
      do ispin1 = 1, this%spin_channels
        do ispin2 = 1, this%spin_channels
          if(ispin1 /= ispin2) then
            tmpU = tmpU + this%zn(im,im,ispin1,ios)*this%zn(imp,imp,ispin2,ios)
          end if
        end do
      end do

      if(im == imp) then
        tmpU = tmpU - (this%zn(im,im,3,ios)*this%zn(im,im,4,ios) &
                          +this%zn(im,im,4,ios)*this%zn(im,im,3,ios))
      end if
      denomU = denomU + tmpU
    end do
    end do
    this%orbsets(ios)%Ueff = TOFLOAT(numU)/TOFLOAT(denomU) - TOFLOAT(numJ)/TOFLOAT(denomJ)
    this%orbsets(ios)%Ubar = TOFLOAT(numU)/TOFLOAT(denomU)
    this%orbsets(ios)%Jbar = TOFLOAT(numJ)/TOFLOAT(denomJ)

  else !In the case of s orbitals, the expression is different
    ! sum_{alpha/=beta} P^alpha_{mmp}P^beta_{mpp,mppp}  
    ! sum_{alpha,beta} sum_{m,mp} N^alpha_{m}N^beta_{mp}
    numU = M_z0
    denomU = M_z0
    do ispin1 = 1, this%spin_channels
      do ispin2 = 1, this%spin_channels
        if(ispin1 /= ispin2) then
          numU = numU + this%zn_alt(1,1,ispin1,ios)*this%zn_alt(1,1,ispin2,ios) &
                          *this%zcoulomb(1,1,1,1,ispin1,ispin2,ios)
          denomU = denomU + this%zn(1,1,ispin1,ios)*this%zn(1,1,ispin2,ios)
        end if
      end do
    end do
    denomU = denomU + this%zn(1,1,3,ios)*this%zn(1,1,4,ios) &
                          +this%zn(1,1,4,ios)*this%zn(1,1,3,ios)

    ! We have to be careful in the case of hydrogen atom for instance 
    if(abs(denomU)> CNST(1.0e-3)) then
      this%orbsets(ios)%Ubar = (TOFLOAT(numU)/TOFLOAT(denomU))
    else
      write(message(1),'(a,a)')' Small denominator value for the s orbital ', this%orbsets(ios)%Ubar
      write(message(2),'(a)')' U is set to zero '
      call messages_warning(2, namespace=namespace)
      this%orbsets(ios)%Ubar = M_ZERO
    end if
    this%orbsets(ios)%Jbar = 0
    this%orbsets(ios)%Ueff = this%orbsets(ios)%Ubar
  end if

  POP_SUB(compute_ACBNO_U_noncollinear)
end subroutine compute_ACBNO_U_noncollinear


