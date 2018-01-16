!> @file
!!    Routines to bind fake argumentf for cuda solver
!! @author
!!    Copyright (C) 2002-2015 BigDFT group  (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
 subroutine cuda_estimate_memory_needs_cu(iproc,n,&
      geo,plansSize,maxPlanSize,freeGPUSize, totalGPUSize )
   use dictionaries
   call f_err_throw('We should not enter into the cuda_estimation_routine')
 end subroutine cuda_estimate_memory_needs_cu

 subroutine pad_data()
   use dictionaries
   call f_err_throw('We should not enter into the pad_data routine')
 end subroutine pad_data

 subroutine unpad_data()
   use dictionaries
   call f_err_throw('We should not enter into the unpad_data routine')
 end subroutine unpad_data

 subroutine finalize_reduction_kernel()
   use dictionaries
   call f_err_throw('We should not enter into the finalize_reduction_kernel routine')
 end subroutine finalize_reduction_kernel

subroutine first_reduction_kernel()
   implicit none
   stop 'reduction1'
 END SUBROUTINE first_reduction_kernel

subroutine second_reduction_kernel()
   implicit none
   stop 'reduction2'
 END SUBROUTINE second_reduction_kernel

 subroutine third_reduction_kernel()
   implicit none
   stop 'reduction3'
 END SUBROUTINE third_reduction_kernel

subroutine cuFFTdestroy()
   implicit none
   stop 'FFTdestroy'
 END SUBROUTINE cuFFTdestroy

subroutine cuda_3d_psolver_general_plan()
   implicit none
   stop 'GPUsolver'
 END SUBROUTINE cuda_3d_psolver_general_plan

subroutine cuda_3d_psolver_general()
   implicit none
   stop 'GPUsolvergen'
 END SUBROUTINE cuda_3d_psolver_general

subroutine cuda_3d_psolver_plangeneral()
   implicit none
   stop 'GPUsolverplan'
 END SUBROUTINE cuda_3d_psolver_plangeneral

subroutine cudacreatestream()
   implicit none
   stop 'cudacreatestream'
 END SUBROUTINE cudacreatestream

subroutine cudadestroystream()
   implicit none
   stop 'cudadestroystream'
 END SUBROUTINE cudadestroystream

subroutine cudacreatecublashandle()
   implicit none
   stop 'cudacreatecublashandle'
 END SUBROUTINE cudacreatecublashandle

subroutine cudadestroycublashandle()
   implicit none
   stop 'cudadestroycublashandle'
 END SUBROUTINE cudadestroycublashandle


subroutine gpu_accumulate_eexctx()
   implicit none
   stop 'gpu_accumulate_eexctx'
 END SUBROUTINE gpu_accumulate_eexctx

subroutine gpu_post_computation()
   implicit none
   stop 'gpu_post_computation'
 END SUBROUTINE gpu_post_computation

subroutine gpu_pre_computation()
   implicit none
   stop 'gpu_pre_computation'
 END SUBROUTINE gpu_pre_computation
