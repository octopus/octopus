!> @file
!! Fake routine in order to compile without CUDA
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

  subroutine set_cpu_gpu_aff()
    implicit none
    stop 'FAKE CPU_GPU_AFF'
  END SUBROUTINE set_cpu_gpu_aff

  subroutine intertamponcGPU()
    implicit none
    stop 'FAKE CUDA Interface'
  END SUBROUTINE intertamponcGPU

  subroutine CUDA_ALLOC_MEM()
    implicit none
    stop 'FAKE CUDA_ALLOC_MEM'
  END SUBROUTINE CUDA_ALLOC_MEM

  subroutine cuda_psi_to_vpsi()
    implicit none
    stop 'fake cuda_psi_to_vpsi'
  END SUBROUTINE cuda_psi_to_vpsi
     
  subroutine cuda_fetch_vpsi()
    implicit none
    stop 'fake cuda_fetch_vpsi'
  END SUBROUTINE cuda_fetch_vpsi

  subroutine CUDA_DEALLOCATE_MEM()
    implicit none
    stop 'fake CUDA_DEALLOCATE_MEM'
  END SUBROUTINE CUDA_DEALLOCATE_MEM

  subroutine GPU_allocate()
    implicit none
    stop 'fake GPU_allocate'
  END SUBROUTINE GPU_allocate

  subroutine GPU_deallocate()
    implicit none
    stop 'fake GPU_deallocate'
  END SUBROUTINE GPU_deallocate

  subroutine GPU_send()
    implicit none
    stop 'fake GPU_send'
  END SUBROUTINE GPU_send

  subroutine GPU_receive()
    implicit none
    stop 'fake GPU_receive'
  END SUBROUTINE GPU_receive

  subroutine localpotential()
    implicit none
    stop 'fake localpotential'
  END SUBROUTINE localpotential

  subroutine localpotentiald()
    implicit none
    stop 'fake localpotentiald'
  END SUBROUTINE localpotentiald

  subroutine kineticterm()
    implicit none
    stop 'fake kineticterm'
  END SUBROUTINE kineticterm

  subroutine kinetictermd()
    implicit none
    stop 'fake kinetictermd'
  END SUBROUTINE kinetictermd

   subroutine prepare_gpu_for_locham()
    implicit none
    stop 'fake prepare_gpu_for_locham'
  END SUBROUTINE prepare_gpu_for_locham

  subroutine gpu_locham()
    implicit none
    stop 'gpu_locham'
  END SUBROUTINE gpu_locham

  subroutine gpu_precond()
    implicit none
    stop 'gpu_locham'
  END SUBROUTINE gpu_precond

  subroutine gpu_locden()
    implicit none
    stop 'gpu_locham'
  END SUBROUTINE gpu_locden

  subroutine free_gpu()
    implicit none
    stop 'free_gpu'
  END SUBROUTINE free_gpu
 
 subroutine preconditionall_gpu()
   implicit none
   stop 'preconditionall_gpu'
 END SUBROUTINE preconditionall_gpu

 subroutine init_gpu_sharing()
   implicit none
   stop 'init_gpu_sharing'
 END SUBROUTINE init_gpu_sharing

 subroutine local_partial_density_gpu()
   implicit none
   stop 'local_partial_density_gpu'
 END SUBROUTINE local_partial_density_gpu

 subroutine local_hamiltonian_gpu()
   implicit none
   stop 'local_hamiltonian_gpu'
 END SUBROUTINE local_hamiltonian_gpu

 subroutine stop_gpu_sharing()
   implicit none
   stop 'stop_gpu_sharing'
 END SUBROUTINE stop_gpu_sharing

 subroutine init_lib()
   implicit none
   stop 'init_lib'
 END SUBROUTINE init_lib

 subroutine sg_init()
   implicit none
   stop 'sg_init'
 END SUBROUTINE sg_init

 subroutine sg_end()
   implicit none
   !stop 'sg_end'
 END SUBROUTINE sg_end

subroutine cudamalloc()
   implicit none
   stop 'allocation'
 END SUBROUTINE cudamalloc

subroutine cudafree()
   implicit none
   stop 'free'
 END SUBROUTINE cudafree

subroutine cuFFTdestroy()
   implicit none
   stop 'FFTdestroy'
 END SUBROUTINE cuFFTdestroy


subroutine reset_gpu_data()
   implicit none
   stop 'reset'
 END SUBROUTINE reset_gpu_data

subroutine get_gpu_data()
   implicit none
   stop 'get'
 END SUBROUTINE get_gpu_data


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
