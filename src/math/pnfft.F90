!! Copyright (C) 2013 Umberto De Giovannini
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


#include "global.h"


!> The includes for the PNFFT
module pnfft_params_m
  use fftw_m
  use pfft_params_m
#ifdef HAVE_PNFFT
  use iso_c_binding
#endif

  implicit none

#ifdef HAVE_PNFFT
  include "pnfft.f03"
#else
!           integer, parameter :: ptrdiff_t_kind = 8
#endif
end module pnfft_params_m
 
!> The low level module to work with the PNFFT library.
!! http://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
module pnfft_m
  use datasets_m
  use global_m
  use io_m
  use loct_math_m
  use math_m
  use messages_m
  use mpi_m
  use parser_m
  use pfft_m
  use pnfft_params_m
  use profiling_m

  implicit none

  private

#ifdef HAVE_PNFFT

  public ::               &
    pnfft_t,              &      
    pnfft_write_info,     &
    pnfft_init_plan,      &
    pnfft_copy_params,    &   
    pnfft_init_procmesh,  &
    pnfft_end,            &
    pnfft_set_sp_nodes,   &
    zpnfft_forward,       &
    zpnfft_backward,      &
    dpnfft_forward,       &
    dpnfft_backward
    
 
  type pnfft_t
  
! Parameters
    integer               :: np(3)           !> Processes
    integer(C_INTPTR_T)   :: N_local(3)     !> Number of Fourier coefficients
    integer(C_INTPTR_T)   :: N(3)            !> Number of Fourier coefficients local
    integer(C_INTPTR_T)   :: Nos(3)          !> FFT grid size
    integer(C_INTPTR_T)   :: M(3)
    integer(C_INTPTR_T)   :: local_M         !> Local number of nodes per process
    integer               :: mm              !> Real space cut-off
    FLOAT                 :: sigma           !> oversampling factor
    integer               :: flags           !> PNFFT initialization options 
    logical               :: set_defaults = .false. !> set default values from the code

    FLOAT                 :: norm       !> Normalization  


! Data 
    type(C_PTR)           :: plan            !> pnfft plan

    complex(C_DOUBLE_COMPLEX), pointer :: f_lin(:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: f(:,:,:) => NULL()
    complex(C_DOUBLE_COMPLEX), pointer :: f_hat(:,:,:) => NULL()
    real(C_DOUBLE),            pointer :: x_lin(:,:) => NULL()
    real(C_DOUBLE),            pointer :: x(:,:,:,:) => NULL()

    real(C_DOUBLE)        :: lower_border(3) !> contains the real-space nodes local borders  
    real(C_DOUBLE)        :: upper_border(3) !> parallelization
    FLOAT                 :: lo_global(3)
    FLOAT                 :: up_global(3)
  
  end type pnfft_t

contains


! ---------------------------------------------------------  
  subroutine pnfft_guru_options(pnfft)
    type(pnfft_t), intent(inout) :: pnfft

    PUSH_SUB(pnfft_guru_options)


    !%Variable PNFFTCutoff
    !%Type integer
    !%Default 6
    !%Section Mesh::FFTs
    !%Description
    !% Cut-off parameter of the window function. 
    !%End
    call parse_integer(datasets_check('PNFFTCutoff'), pnfft%mm, pnfft%mm)


    !%Variable PNFFTOversampling
    !%Type float
    !%Default 2.0
    !%Section Mesh::FFTs
    !%Description
    !% PNFFT oversampling factor (sigma). This will rule the size of the FFT under the hood.
    !%End
    call parse_float(datasets_check('PNFFTOversampling'), pnfft%sigma, pnfft%sigma)



    POP_SUB(pnfft_guru_options)
  end subroutine pnfft_guru_options

  ! ---------------------------------------------------------
  subroutine pnfft_init_params(pnfft, nn, optimize)
    type(pnfft_t),     intent(inout) :: pnfft
    integer,           intent(in)    :: nn(3) !> pnfft bandwidths 
    logical, optional, intent(in)    :: optimize

    integer :: ii, jj, idir, my_nn(3)
    logical :: optimize_
    integer :: pnfft_flags

    PUSH_SUB(pnfft_init_params)

    optimize_ = optional_default(optimize, .true.)

    if(.not. pnfft%set_defaults) then
      !Set defaults
      pnfft%mm = 6 
      pnfft%sigma = M_TWO
    end if
    
    call pnfft_guru_options(pnfft)

    my_nn = 0
    do ii = 1, 3
      my_nn(ii) = nn(ii)*pnfft%sigma
      if(optimize_) call loct_fft_optimize(my_nn(ii), 1) ! ask for an odd number
    end do
    
    pnfft%Nos(1:3) = my_nn(1:3)

    pnfft%flags = PNFFT_MALLOC_X + PNFFT_MALLOC_F_HAT + PNFFT_MALLOC_F + &
                  PNFFT_WINDOW_KAISER_BESSEL




    POP_SUB(pnfft_init_params)
  end subroutine pnfft_init_params


  ! ---------------------------------------------------------
  subroutine pnfft_init_procmesh(pnfft, mpi_grp, comm)
    type(pnfft_t), intent(inout)  :: pnfft
    type(mpi_grp_t),   intent(in) :: mpi_grp
    integer,           intent(out):: comm 
  
    integer :: ierror

    PUSH_SUB(pnfft_init_procmesh)

        call pnfft_init()
      
        pnfft%np(1:3) = 1

        call pfft_decompose(mpi_grp%size, pnfft%np(1), pnfft%np(2))
            
        ierror = pnfft_create_procmesh(2, mpi_grp%comm,  pnfft%np, comm)        

        if (ierror /= 0) then
          message(1) = "The number of rows and columns in PNFFT processor grid is not equal to "
          message(2) = "the number of processor in the MPI communicator."
          message(3) = "Please check it."
          call messages_fatal(3)
        end if
 
    POP_SUB(pnfft_init_procmesh)
  end subroutine pnfft_init_procmesh


  ! ---------------------------------------------------------
  subroutine pnfft_copy_params(in, out)
    type(pnfft_t), intent(in)  :: in
    type(pnfft_t), intent(out) :: out


    PUSH_SUB(pnfft_copy_params)

    out%mm = in%mm
    out%sigma = in%sigma
    out%set_defaults = in%set_defaults       
  
 
    POP_SUB(pnfft_copy_params)
  end subroutine pnfft_copy_params

  ! ---------------------------------------------------------
  subroutine pnfft_write_info(pnfft)
    type(pnfft_t), intent(in) :: pnfft

    integer :: idir

    PUSH_SUB(pnfft_write_info)



    call messages_write("Info: PNFFT parameters")
    call messages_new_line()
    call messages_write("      Fourier coefficients      N = ")
    do idir = 1, 3
      call messages_write(pnfft%N(idir))
      if(idir < 3) call messages_write(" x ")
    end do
    call messages_new_line()
    call messages_write("      Spatial nodes per process   = ")
    call messages_write(pnfft%local_M)
    call messages_new_line()
    call messages_write("      Oversampling factor   sigma = ")
    call messages_write(pnfft%sigma)
    call messages_new_line()
    call messages_write("      FFT grid size             n = ")
    do idir = 1, 3
      call messages_write(pnfft%Nos(idir))
      if(idir < 3) call messages_write(" x ")
    end do
    call messages_new_line()
    call messages_write("      Real Space cutoff           = ")
    call messages_write(pnfft%mm)  
    call messages_new_line()
    call messages_write("      Process mesh             np = ")
    do idir = 1, 3
      call messages_write(pnfft%np(idir))
      if(idir < 3) call messages_write(" x ")
    end do
    call messages_info()



 
    POP_SUB(pnfft_write_info)
  end subroutine pnfft_write_info



  ! ---------------------------------------------------------
  subroutine pnfft_init_plan(pnfft, mpi_comm, fs_n_global, fs_n, fs_istart, rs_n, rs_istart)
    type(pnfft_t), intent(inout) :: pnfft
    integer,         intent(in)  :: mpi_comm         !< MPI comunicator
    integer,         intent(in)  :: fs_n_global(1:3) !< The general number of elements in each dimension in Fourier space
    integer,         intent(out) :: fs_n(1:3)        !< Local number of elements in each direction in Fourier space
    integer,         intent(out) :: fs_istart(1:3)   !< Where does the local portion of the function start in Fourier space
    integer,         intent(out) :: rs_n(1:3)        !< Local number of elements in each direction in real space
    integer,         intent(out) :: rs_istart(1:3)   !< Where does the local portion of the function start in real space

    real(C_DOUBLE) :: lower_border(3), upper_border(3)
    real(C_DOUBLE) :: x_max(3)
    integer(C_INTPTR_T) :: local_N(3), local_N_start(3), d=3, local_M
    type(C_PTR) :: cf_hat, cf, cx


    PUSH_SUB(pnfft_init_plan)

    pnfft%N(1:3) = fs_n_global(1:3)
    
    call pnfft_init_params(pnfft, fs_n_global(1:3), optimize = .true.)
    
    x_max(:) = CNST(0.5)
         
    call pnfft_local_size_guru(3, pnfft%N, pnfft%Nos, x_max, pnfft%mm, mpi_comm, &
         PNFFT_TRANSPOSED_F_HAT, local_N, local_N_start, lower_border, upper_border)
!          PNFFT_TRANSPOSED_NONE, local_N, local_N_start, lower_border, upper_border)

    pnfft%lower_border = lower_border
    pnfft%upper_border = upper_border
    pnfft%N_local(1:3)   = local_N(1:3) 
    pnfft%M(1:3)   = pnfft%N_local(1:3) 
! 
!     pnfft%M(1)   = pnfft%N_local(1) 
!     pnfft%M(2)   = pnfft%N_local(3) 
!     pnfft%M(3)   = pnfft%N_local(2) 
    local_M = pnfft%M(1) * pnfft%M(2) * pnfft%M(3)


! In order to stick with octopus convention that the forward FFT maps rs to fs 
! we have to use the PNFFT library with swapped rs and fs definitions.
! This might seem wrong because, while for ordinary FFTs rs and fs are symmetric,
! PNFFT (and NFFT) breaks this symmetry by introducing different grids for rs
! and fs. 
! For some reason I don`t understand however this procedure of exchanging the definitions
! is working fine for NFFT. Hopefully it will also do for PNFFT 
    
    fs_n(1:3)      = local_N(1:3) 
    fs_istart(1:3) = pnfft%N(1:3)/2 + local_N_start(1:3) + 1 

!     if(iand(pnfft%flags, PNFFT_TRANSPOSED_F_HAT) /= 0 ) then
      ! we use the same decomposition 
      ! but since we use PNFFT_TRANSPOSED_F_HAT the indices are transposed
!     rs_n(1)      = fs_n(2) 
!     rs_n(2)      = fs_n(1) 
!     rs_n(3)      = fs_n(3) 
    rs_n(1:3) = pnfft%M(1:3) 
    rs_istart(1:3) = fs_istart(1:3) 

!     rs_istart(1) = fs_istart(1) 
!     rs_istart(2) = fs_istart(3) 
!     rs_istart(3) = fs_istart(2) 

!     else
!       fs_n(:)      = rs_n(:)
!       fs_istart(:) = rs_istart(:)
!     end if
    

    pnfft%plan = pnfft_init_guru(3, pnfft%N, pnfft%Nos, x_max, local_M, pnfft%mm, &
                 pnfft%flags, PFFT_ESTIMATE, mpi_comm)

    pnfft%local_M=local_M
     
    ! Get data pointers in C format
    cf_hat = pnfft_get_f_hat(pnfft%plan)
    cf     = pnfft_get_f(pnfft%plan)
    cx     = pnfft_get_x(pnfft%plan)

    ! Convert data pointers to Fortran format
!     call c_f_pointer(cf_hat, pnfft%f_hat, [local_N(3),local_N(2),local_N(1)])
    call c_f_pointer(cf_hat, pnfft%f_hat, [local_N(1),local_N(3),local_N(2)])
    call c_f_pointer(cf,     pnfft%f_lin, [pnfft%local_M])
    call c_f_pointer(cf,     pnfft%f,     [pnfft%M(1),pnfft%M(2),pnfft%M(3)])
    call c_f_pointer(cx,     pnfft%x_lin, [d, pnfft%local_M])
    call c_f_pointer(cx,     pnfft%x,     [d, pnfft%M(1),pnfft%M(2),pnfft%M(3)])




    print *, mpi_world%rank, "local_N(3)      ", local_N
    print *, mpi_world%rank, "local_N_start(3)", local_N_start
    print *, mpi_world%rank, "fs_istart(1:3)  ", fs_istart
    print *, mpi_world%rank, "fs_n(1:3)       ", fs_n
    print *, mpi_world%rank, "rs_istart(1:3)  ", rs_istart
    print *, mpi_world%rank, "rs_n(1:3)       ", rs_n
    print *, mpi_world%rank, "lower_border    ", lower_border
    print *, mpi_world%rank, "upper_border    ", upper_border
    print *, mpi_world%rank, "rs_range        ", upper_border(:) - lower_border(:)
    print *, mpi_world%rank, "local_M         ", local_M
    print *, mpi_world%rank, "pnfft%N_local   ", pnfft%N_local 
    print *, mpi_world%rank, "pnfft%M         ", pnfft%M 


    POP_SUB(pnfft_init_plan)
  end subroutine pnfft_init_plan

  ! ---------------------------------------------------------
  subroutine pnfft_end(pnfft)
    type(pnfft_t), intent(inout) :: pnfft

    PUSH_SUB(pnfft_end)

    call pnfft_finalize(pnfft%plan, PNFFT_FREE_X + PNFFT_FREE_F_HAT + PNFFT_FREE_F)
    call pnfft_cleanup()
   
    pnfft%f_lin => NULL()
    pnfft%f => NULL()
    pnfft%f_hat => NULL()
    pnfft%x => NULL()
    pnfft%x_lin => NULL()
   
    POP_SUB(pnfft_end)
  end subroutine pnfft_end
  
  
  ! ---------------------------------------------------------  
  ! Set the coordinates for the spatial nodes rescaled to 
  ! the 3D torus [-0.5,0.5)
  ! ---------------------------------------------------------  
  subroutine pnfft_set_sp_nodes(pnfft, X, rs_istart)
    type(pnfft_t),    intent(inout) :: pnfft
    FLOAT,            intent(in)    :: X(:,:) !X(i, dim)
    integer,          intent(in)    :: rs_istart(1:3)

    FLOAT   :: len(3), cc(3), eps
    integer :: ii, idir, i1, i2, i3, size(3)
    FLOAT, allocatable ::  dX(:,:)

    PUSH_SUB(pnfft_set_sp_nodes)
 
    eps = 1.000001 ! the sample nodes must be in [0.5,0.5)
    
    do idir=1,3  
      len(:) = (maxval(X(:,idir))-minval(X(:,idir)))*eps
      cc(:) = (minval(X(:,idir))+maxval(X(:,idir)))/M_TWO
    end do
    
    
    
    do i1 = 1, pnfft%M(1)
      do i2 = 1, pnfft%M(2)
        do i3 = 1, pnfft%M(3)
          pnfft%x(1, i1,i2,i3) = (X(rs_istart(1)+i1-1, 1)  - cc(1))/len(1)
          pnfft%x(2, i1,i2,i3) = (X(rs_istart(2)+i2-1, 2)  - cc(2))/len(2)
          pnfft%x(3, i1,i2,i3) = (X(rs_istart(3)+i3-1, 3)  - cc(3))/len(3)
!           pnfft%x_lin(1, pnfft_idx_3to1(pnfft,ix,iy,iz)) = (X(rs_istart(1)+ix-1, 1)  - cc(1))/len(1)
!           pnfft%x_lin(2, pnfft_idx_3to1(pnfft,ix,iy,iz)) = (X(rs_istart(2)+iy-1, 2)  - cc(2))/len(2)
!           pnfft%x_lin(3, pnfft_idx_3to1(pnfft,ix,iy,iz)) = (X(rs_istart(3)+iz-1, 3)  - cc(3))/len(3)

        end do
      end do
    end do


    call pnfft_messages_debug(pnfft)



!     x1_ =(X(:,1)-cc)/length
!     length = (maxval(X(:,2))-minval(X(:,2)))*eps
!     cc = (minval(X(:,2))+maxval(X(:,2)))/M_TWO
!     x2_ =(X(:,2)-cc)/length
!     length = (maxval(X(:,3))-minval(X(:,3)))*eps
!     cc = (minval(X(:,3))+maxval(X(:,3)))/M_TWO
!     x3_ =(X(:,3)-cc)/length
    
     
!     ! Set the normalization factor  
!     do ii = 1, pnfft%M-1
!       dX(ii,1)= abs(x1_(ii+1)-x1_(ii))
!       dX(ii,2)= abs(x2_(ii+1)-x2_(ii))
!       dX(ii,3)= abs(x3_(ii+1)-x3_(ii))
!     end do
! 
!     pnfft%norm = M_ONE/(minval(dX(:,1)) * minval(dX(:,2)) * minval(dX(:,3)))


    POP_SUB(pnfft_set_sp_nodes)
  end subroutine pnfft_set_sp_nodes
  
  ! ---------------------------------------------------------
  subroutine pnfft_messages_debug(pnfft)
    type(pnfft_t), intent(in) :: pnfft 

    integer          :: nn, i1, i2, i3 
    integer          :: npart, ierr
    integer          :: iunit          !< For debug output to files.
    character(len=3) :: filenum

    PUSH_SUB(pnfft_messages_debug)

    if(in_debug_mode) then
  
      if(mpi_grp_is_root(mpi_world)) then
        call io_mkdir('debug/PNFFT')
      end if      
      call MPI_Barrier(mpi_world%comm, ierr)

      nn = mpi_world%rank
      write(filenum, '(i3.3)') nn

      iunit = io_open('debug/PNFFT/rs_partition.'//filenum, &
           action='write')
           
      do i1 = 1, pnfft%M(1)
       do i2 = 1, pnfft%M(2)
         do i3 = 1, pnfft%M(3)
           write(iunit, '(3f18.8)') pnfft%x(1, i1,i2,i3), pnfft%x(2, i1,i2,i3), pnfft%x(3, i1,i2,i3) 
         end do
       end do
      end do     
      call io_close(iunit)

    end if

    POP_SUB(pnfft_messages_debug)

  end subroutine pnfft_messages_debug

  !---------------------------------------------------------------------------------
  integer function pnfft_idx_3to1(pnfft, ix , iy, iz) result(idx)
    type(pnfft_t),  intent(in) :: pnfft
    integer,        intent(in) :: ix 
    integer,        intent(in) :: iy 
    integer,        intent(in) :: iz 

    idx =  (ix-1)*pnfft%M(2)*pnfft%M(3) + (iy-1)*pnfft%M(3) + (iz-1) + 1

  end function pnfft_idx_3to1




  #include "undef.F90"
  #include "real.F90"
  #include "pnfft_inc.F90"

  #include "undef.F90"
  #include "complex.F90"
  #include "pnfft_inc.F90"


#endif

end module pnfft_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
