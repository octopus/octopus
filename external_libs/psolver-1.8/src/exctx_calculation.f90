!> @file
!!  Routines which define internal exact exchange operations
!! @author
!! Copyright (C) 2002-2015 BigDFT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~/COPYING file
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the list of contributors, see ~/AUTHORS
subroutine internal_calculation_exctx(istep,factor,pkernel,norb,occup,spinsgn,remote_result,&
     nloc_i,nloc_j,isloc_i,isloc_j,&
     phi_i,phi_j,eexctX,rp_ij)
  use PSbase, wp=>dp
  use overlap_point_to_point
  use Poisson_Solver
  implicit none
  logical, intent(in) :: remote_result
  integer, intent(in) :: istep !<step of the calculation
  integer, intent(in) :: norb
  integer, intent(in) :: nloc_i,nloc_j !<number of local elements to  be treated
  integer, intent(in) :: isloc_i !<starting point of the elements for phi_i
  integer, intent(in) :: isloc_j !<starting point of the elements for phi_j
  real(gp), intent(in) :: factor !<overall factor to treat the data
  real(gp), dimension(norb), intent(in) :: occup,spinsgn !<to treat the data
  type(coulomb_operator), intent(inout) :: pkernel
  type(local_data), intent(inout) :: phi_i,phi_j
  real(gp), intent(inout) :: eexctX
  real(wp), dimension(product(pkernel%ndims)), intent(out) :: rp_ij
  !local variables
  integer :: iorb,jorb,ndim,iorb_glb,jorb_glb,ishift,jshift,ishift_res,jshift_res,i
  real(gp) :: hfaci,hfacj,hfac2,ehart
!loop over all the orbitals
!for the first step do only the upper triangular part
  !do iorb=iorbs,iorbs+norbi-1
!!$  do ind=1,nloc_i
!!$     iorb=

  ndim=product(pkernel%ndims)
  if(pkernel%igpu==1 .and. pkernel%stay_on_gpu /= 1) then
    call synchronize()
 end if
  do iorb=isloc_i,nloc_i+isloc_i-1
  do jorb=isloc_j,nloc_j+isloc_j-1
     !aliasing
     jorb_glb=phi_j%id_glb(jorb)
     iorb_glb=phi_i%id_glb(iorb)
     hfaci=-factor*occup(jorb_glb)
     hfacj=-factor*occup(iorb_glb)
     ishift=phi_i%displ(iorb)
     jshift=phi_j%displ(jorb)
     ishift_res=phi_i%displ_res(iorb)
     jshift_res=phi_j%displ_res(jorb)
     !first cross-check whether the spin indices are the same
     if (spinsgn(iorb_glb) /= spinsgn(jorb_glb)) then
        print *,'temporary',iorb,iorb_glb,jorb,jorb_glb
        stop
     end if
     !do it only for upper triangular results 
     if (istep /= 0 .or. jorb_glb >= iorb_glb) then
        if(pkernel%igpu==1 .and. pkernel%stay_on_gpu /= 1 .and. istep ==0) then
          call reset_gpu_data(ndim,rp_ij,pkernel%w%rho_GPU)
        end if 

        call exctx_pre_computation(iorb, jorb,rp_ij,phi_i,phi_j,pkernel)
!!$        ncalls=ncalls+1
!!$        !Poisson solver in sequential
!!$        if (iproc == iprocref .and. verbose > 1) then
!!$           call yaml_comment('Exact exchange calculation: ' // trim(yaml_toa( &
!!$                nint(real(ncalls,gp)/real(ncalltot,gp)*100.0_gp),fmt='(i3)')) //'%')
!!$        end if

        !call Electrostatic_Solver(
        call H_potential('D',pkernel,rp_ij,rp_ij,ehart,0.0_dp,.false.,&
             quiet='YES')
        !print *,'ext',ehart

        call exctx_accum_eexctX(iorb, jorb, phi_i, phi_j, pkernel, norb,occup,factor,remote_result,istep, ehart, eexctX)

        !accumulate the results for each of the wavefunctions concerned
        call exctx_post_computation(iorb, jorb, rp_ij, phi_i, phi_j, pkernel, norb, occup,factor)
        if ((iorb_glb /= jorb_glb .and. istep==0) .or. remote_result) then
          call exctx_post_computation(jorb, iorb, rp_ij, phi_j, phi_i, pkernel, norb, occup,factor)
        end if
        !write(100+iproc,*)iorb+isorb,jorb+jsorb,igrpr(igroup)
     end if
  end do
  end do

  if(pkernel%igpu==1 .and. pkernel%stay_on_gpu /= 1) then
    call get_gpu_data(ndim,rp_ij,pkernel%w%rho_GPU)
 end if

end subroutine internal_calculation_exctx

subroutine exctx_pre_computation(iorb, jorb, rp_ij, phi1, phi2, pkernel)
  use PSbase, wp=>dp
  use f_precisions, only: f_address
  use overlap_point_to_point
  use dictionaries, only: f_err_throw
  use Poisson_Solver
  implicit none
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(product(pkernel%ndims)), intent(inout) :: rp_ij
  type(local_data), intent(inout) :: phi1,phi2
  real(gp) :: hfac
  integer(f_address) :: myrho_GPU
  integer :: iorb,jorb,ndim,shift1,shift2,i,i_stat
  hfac=1.0_gp/product(pkernel%hgrids)
  shift1=phi1%displ(iorb)
  shift2=phi2%displ(jorb)
  ndim=product(pkernel%ndims)
  if(pkernel%igpu==1 .and. pkernel%stay_on_gpu==1) then
   !   call cudamalloc(ndim,myrho_GPU,i_stat)
    !  if (i_stat /= 0) call f_err_throw('error cudamalloc myrho_GPU (GPU out of memory ?) ')
    !first attempt : send data each time, to remove later.
    !call reset_gpu_data(ndim,rp_ij,myrho_GPU)
    call gpu_pre_computation(pkernel%ndims(1),pkernel%ndims(2),pkernel%ndims(3),&
            pkernel%w%rho_GPU, phi1%data_GPU,shift1,phi2%data_GPU,shift2,hfac)
    !call get_gpu_data(ndim,rp_ij,myrho_GPU)

    !  call cudafree(myrho_GPU)
!myrho_GPU=0

   ! end if
  else
    !$omp parallel do default(shared) private(i)
    do i=1,ndim
      rp_ij(i)=hfac*phi1%data(i+shift1)*phi2%data(i+shift2)
    end do
    !$omp end parallel do
  end if
end subroutine exctx_pre_computation


subroutine exctx_post_computation(orb1, orb2, rp_ij, phi1, phi2, pkernel, norb, occup, factor)
  use PSbase, wp=>dp
  use f_precisions, only: f_address
  use overlap_point_to_point
  use dictionaries, only: f_err_throw
  use Poisson_Solver
  implicit none
  type(coulomb_operator), intent(inout) :: pkernel
  real(wp), dimension(product(pkernel%ndims)), intent(inout) :: rp_ij
  type(local_data), intent(inout) :: phi1,phi2
  integer, intent(in) :: norb
  real(gp), dimension(norb), intent(in) :: occup
  real(gp), intent(in) :: factor !<overall factor to treat the data
  real(gp) :: hfac1
  integer(f_address) :: myrho_GPU
  integer :: orb1,orb2,ndim,orb2_glb,shift2,shift1_res,i,i_stat


  orb2_glb=phi2%id_glb(orb2)

  shift1_res=phi1%displ_res(orb1)
  shift2=phi2%displ(orb2)

  hfac1=-factor*occup(orb2_glb)
  ndim=product(pkernel%ndims)

  if(pkernel%igpu==1 .and. pkernel%stay_on_gpu==1) then
      !call cudamalloc(ndim,myrho_GPU,i_stat)
    !  if (i_stat /= 0) call f_err_throw('error cudamalloc myrho_GPU (GPU out of memory ?) ')

    !first attempt : send data each time, to remove later.
   ! call reset_gpu_data(ndim,rp_ij,myrho_GPU)

    call gpu_post_computation(pkernel%ndims(1),pkernel%ndims(2),pkernel%ndims(3),&
           pkernel%w%rho_GPU, phi1%res_GPU,shift1_res,phi2%data_GPU,shift2,hfac1)

     ! call cudafree(myrho_GPU)
myrho_GPU=0
  else
  !$omp parallel do default(shared) private(i)
  do i=1,ndim
    phi1%res(i+shift1_res)=phi1%res(i+shift1_res)+hfac1*rp_ij(i)*phi2%data(i+shift2)
  end do
  !$omp end parallel do
  end if
end subroutine exctx_post_computation

subroutine exctx_accum_eexctX(orb1, orb2, phi1, phi2, pkernel, norb, occup, factor, remote_result, istep, ehart, eexctX)
  use PSbase, wp=>dp
  use f_precisions, only: f_address
  use overlap_point_to_point
  use dictionaries, only: f_err_throw
  use Poisson_Solver
  use wrapper_MPI
  use iso_c_binding
  implicit none
  type(coulomb_operator), intent(inout) :: pkernel
  type(local_data), intent(inout) :: phi1,phi2
  integer, intent(in) :: norb
  real(gp), dimension(norb), intent(in) :: occup
  real(gp), intent(in) :: factor,ehart
  real(gp) :: hfac2
  real(gp), pointer ::sendbuf
  integer :: orb1,orb2,orb2_glb,orb1_glb, istep
  real(gp), intent(inout) :: eexctX
  logical, intent(in) :: remote_result
  type(c_ptr)::val
  
  orb1_glb=phi1%id_glb(orb1)
  orb2_glb=phi2%id_glb(orb2)
  !this factor is only valid with one k-point
  !can be easily generalised to the k-point case
  hfac2=factor*occup(orb1_glb)*occup(orb2_glb)

  if(pkernel%igpu==1 .and. pkernel%stay_on_gpu==1) then
    !this part is usually computed at the end of h_potential
    hfac2=hfac2*0.5_dp*product(pkernel%hgrids) 
!    val = TRANSFER(pkernel%w%ehart_GPU, C_NULL_PTR)
!    call c_f_pointer(val, sendbuf)
!    call mpiallred(sendbuf,1,MPI_SUM,comm=pkernel%mpi_env%mpi_comm)
    if (orb1_glb == orb2_glb) then
      call gpu_accumulate_eexctX(pkernel%w%ehart_GPU, pkernel%w%eexctX_GPU, hfac2)
    else
      !if the result has to be sent away
       if (remote_result .or. istep==0) then
         call gpu_accumulate_eexctX(pkernel%w%ehart_GPU, pkernel%w%eexctX_GPU, 2.0_gp*hfac2)
       else !otherwise other processors are already calculating it
         call gpu_accumulate_eexctX(pkernel%w%ehart_GPU, pkernel%w%eexctX_GPU, hfac2)
       end if
    end if
  else
        !exact exchange energy
        if (orb1_glb == orb2_glb) then
           eexctX=eexctX+hfac2*real(ehart,gp)
        else
           !if the result has to be sent away
           if (remote_result .or. istep==0) then
              eexctX=eexctX+2.0_gp*hfac2*real(ehart,gp)
           else !otherwise other processors are already calculating it
              eexctX=eexctX+hfac2*real(ehart,gp)
           end if
        end if
  end if
end subroutine exctx_accum_eexctX
