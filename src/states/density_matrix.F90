!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: states.F90 5022 2009-03-03 17:47:58Z nitsche $

#include "global.h"

module density_matrix_m

  use datasets_m
  use embedded_particles_m
  use global_m
  use grid_m
  use io_function_m
  use io_m
  use lalg_adv_m
  use loct_m
  use loct_parser_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use mpi_lib_m
  use profiling_m
  use states_m

!  use derivatives_m
!  use calc_mode_m
!  use crystal_m
!  use datasets_m
!  use distributed_m
!  use blas_m
!  use geometry_m
!  use hardware_m
!  use loct_parser_m
!  use math_m
!  use mesh_m
!  use multicomm_m
!  use simul_box_m
!  use smear_m
!  use states_dim_m
!  use units_m
!  use varinfo_m
!  use lalg_adv_m

  implicit none

  private

  public :: density_matrix_write

contains

  ! ---------------------------------------------------------
  subroutine density_matrix_write(dir, gr, st, embeddedparticles)
    character(len=*), intent(in) :: dir
    type(grid_t),     intent(in) :: gr
    type(states_t),   intent(in) :: st
    type(embedded_particle_t),   intent(in) :: embeddedparticles

    integer :: mm, jj, ll, j, err_code, iunit, ndims
    integer :: ipart,ndensmat_to_calculate,ncols
    integer, allocatable :: n1(:), n2(:), npoints(:)
    logical :: bof
    FLOAT :: origin
    character(len=200) :: dirname, filename
    CMPLX, allocatable :: densmatr(:, :), evectors(:, :), wavef(:,:)
    FLOAT, allocatable :: evalues(:), sqrdensity(:), density(:), graddens(:)
    FLOAT, allocatable :: hartreep(:), potential(:)

    type(block_t) :: blk
    character(80), allocatable :: labels_densmat(:)
    integer, allocatable :: particle_kept_densmat(:)
    integer, allocatable :: nnatorb_prt_densmat(:)

    call push_sub('states.density_matrix_write')

    !%Variable DensityMatricestoCalc
    !%Type block
    !%Section States
    !%Description
    !% choice of which particle density matrices will be calculated and output, in the
    !%  embedded particles scheme
    !%
    !% <tt>%DensityMatricestoCalc
    !% <br>&nbsp;&nbsp; proton   | 1 | 10
    !% <br>&nbsp;&nbsp; electron | 2 | 15
    !% <br>%</tt>
    !%
    !% would ask octopus to calculate the density matrix corresponding to the 1st
    !% particle (whose coordinates correspond to dimensions 1 to ndim_embedded),
    !% which is an proton, then that corresponding to the 2nd particle
    !% (electron with dimensions ndim_embedded+1 to 2*ndim_embedded), printing
    !% 10 natural orbitals for the first and 15 for the second.
    !%
    !%End
   
    if(loct_parse_block(datasets_check('DensityMatricestoCalc'), blk)/=0) then
     message(1) = 'To print out density matrices, you must specify the DensityMatricestoCalc block in input'
     call write_fatal(1)
    end if
   
    ncols = loct_parse_block_cols(blk, 0)
    if(ncols /= 3 ) then
      call input_error("DensityMatricestoCalc")
    end if
    ndensmat_to_calculate=loct_parse_block_n(blk)
    if (ndensmat_to_calculate < 0 .or. &
        ndensmat_to_calculate > embeddedparticles%nparticle_embedded) then
      call input_error("DensityMatricestoCalc")
    end if

    allocate (labels_densmat(ndensmat_to_calculate))
    allocate (particle_kept_densmat(ndensmat_to_calculate))
    allocate (nnatorb_prt_densmat(ndensmat_to_calculate))
   
    do ipart=1,ndensmat_to_calculate
      call loct_parse_block_string(blk, ipart-1, 0, labels_densmat(ipart))
      call loct_parse_block_int(blk, ipart-1, 1, particle_kept_densmat(ipart))
      call loct_parse_block_int(blk, ipart-1, 2, nnatorb_prt_densmat(ipart))
   
      write (message(1),'(a,a)') 'labels_densmat = ', labels_densmat(ipart)
      write (message(2),'(a,i6)') 'particle_kept_densmat = ', particle_kept_densmat(ipart)
      write (message(3),'(a,i6)') 'nnatorb_prt_densmat = ', nnatorb_prt_densmat(ipart)
      call write_info(3)
    end do
    call loct_parse_block_end(blk)


    ! This is the root directory where everything will be written.
    dirname = trim(dir)//'/density-matrix'
    call loct_mkdir(trim(dirname))

    ! The algorithm should consider how many dimensions the wavefunction has (ndims),
    ! and how many (and which) dimensions should be integrated away. For the moment
    ! being, we assume there are only two dimensions and we trace away the second one.
    ndims = gr%sb%dim


    ! Allocatation of the arrays that store the limiting indexes for each direction,
    ! and the total number of points.
    ALLOCATE(n1(ndims), ndims)
    ALLOCATE(n2(ndims), ndims)
    ALLOCATE(npoints(ndims), ndims)
    do j = 1, ndims
      n1(j) = gr%mesh%idx%nr(1, j)
      n2(j) = gr%mesh%idx%nr(2, j)
      npoints(j) = gr%mesh%idx%ll(j)
    end do


    ! not always the real origin if the box is shifted, no?
    !  which happens to be my case...
    !  only important for printout, so it is ok
    origin=(npoints(1)/2+1)*gr%mesh%h(1)

    ALLOCATE(densmatr(npoints(1), npoints(1)), npoints(1)*npoints(1) )
    ALLOCATE(evectors(npoints(1), npoints(1)), npoints(1)*npoints(1) )
    ALLOCATE(evalues(npoints(1)), npoints(1))
!    ALLOCATE(wavef(npoints(1), npoints(2)), npoints(1)*npoints(2) )
    ALLOCATE(sqrdensity(npoints(1)), npoints(1))   
    ALLOCATE(graddens(npoints(1)), npoints(1))
    ALLOCATE(potential(npoints(1)), npoints(1))
    ALLOCATE(density(npoints(1)), npoints(1))
    ALLOCATE(hartreep(npoints(1)), npoints(1))

    write (*,*) 'shape(densmatr) ', shape(densmatr)
    write (*,*) ' st%nst = ', st%nst

    states_loop: do mm = 1, st%nst
      densmatr  = M_z0
      graddens  = M_ZERO
      potential = M_ZERO
      hartreep  = M_ZERO

      if(states_are_real(st)) then
        call dmf_calculate_gamma(gr%mesh, st%dpsi(:, 1, mm, 1), densmatr)
      else
        call zmf_calculate_gamma(gr%mesh, st%zpsi(:, 1, mm, 1), densmatr)
      end if

      ! Only node zero writes.
      if(mpi_grp_is_root(mpi_world)) then

        do ll = 1, npoints(1)
          density(ll)= real(densmatr(ll,ll))
          sqrdensity(ll)=sqrt(density(ll))
        end do
        
        do ll = 1, npoints(1)
          do jj = 1, npoints(1)
            hartreep(ll) = hartreep(ll)+density(jj)/(sqrt(((ll-jj)*gr%mesh%h(1))**2+1))
          end do
          hartreep(ll) = 0.5*hartreep(ll)*gr%mesh%h(1)
        end do

        ! WARNING: This gradient probably could be improved.
        do ll = 3, npoints(1) - 2
          if(sqrdensity(ll) > CNST(1.0e-12) ) then
            graddens(ll) = (-sqrdensity(ll+2)-sqrdensity(ll-2) & 
              & +16d0*(sqrdensity(ll+1)+sqrdensity(ll-1))&
              & -30d0*sqrdensity(ll))/(12d0*gr%mesh%h(1)**2)
            potential(ll)=graddens(ll)/(2d0*sqrdensity(ll))
          endif
        end do

        !Diagonalize the density matrix
        call lalg_eigensolve(npoints(1), densmatr, evectors, evalues, bof, err_code)
      
        !Write everything into files
        !NOTE: The highest eigenvalues are the last ones not the first!!!
        !      Writing is therefore in reverse order
        evectors = evectors/sqrt(gr%mesh%h(1))
        evalues  = evalues*gr%mesh%h(1)

        write(filename,'(a,i2.2)') trim(dirname)//'/occnumb_',mm
        iunit = io_open(trim(filename), action='write')

        do jj = npoints(1), 1, -1
          write(iunit,'(i4.4,es11.3)') npoints(1)-jj+1, evalues(jj)
        end do
            
        call io_close(iunit)

        do jj = npoints(1)-10, npoints(1)
          write(filename,'(a,i2.2,a,i4.4)') trim(dirname)//'/natorb_', mm, '_', npoints(1)-jj+1
          iunit = io_open(filename, action='write')
          do ll = 1, npoints(1)
            write(iunit,'(es11.3,es11.3,es11.3)') ll*gr%mesh%h(1)-origin, real(evectors(ll,jj)), & 
              & aimag(evectors(ll,jj))
          end do
        call io_close(iunit)
        end do

        write(filename,'(a,i2.2)') trim(dirname)//'/densmatr_', mm
        iunit = io_open(filename,action='write')
        do jj = 1, npoints(1)
          do ll = 1, npoints(1)
            write(iunit,'(es11.3,es11.3)') jj*gr%mesh%h(1)-origin, &
              & ll*gr%mesh%h(1)-origin, real(densmatr(jj,ll)), aimag(densmatr(jj,ll))
          end do
        end do
        call io_close(iunit)

        write(filename,'(a,i2.2)') trim(dirname)//'/potential_', mm
        iunit = io_open(filename,action='write')
        do ll=1, npoints(1)
          write(iunit,'(es11.3,es11.3,es11.3,es11.3)') ll*gr%mesh%h(1)-origin, potential(ll), &
            & hartreep(ll), potential(ll)-hartreep(ll)
        end do
        call io_close(iunit)

      end if

      !Diagonalize 2 particle wave function as well
        
      !  call lalg_eigensolve(npointsx,wavef,evectors,evalues,bof,err_code)
         
      !  evectors=evectors/sqrt(gr%mesh%h(1))
      !  evalues=evalues*gr%mesh%h(1)
        
      !  write(filename,'(a,i6.6,a,i2.2)') 'Tdstep_',i,'/wavefeva_',mm
      !  iunit=io_open(filename,action='write')
      
      !  do jj=npointsx, 1, -1
      !   write(iunit,'(i4.4,es11.3)') npointsx-jj+1, evalues(jj)
      !  enddo
            
      !  call io_close(iunit)
            
      !  do jj=1, npointsx
      !   write(filename,'(a,i6.6,a,i2.2,a,i4.4)') 'Tdstep_',i,'/wavefeve_', mm, '_', npointsx-jj+1
      !   iunit=io_open(filename,action='write')
      
      !   do ll=1, npointsx
      !    write(iunit,'(es11.3,es11.3,es11.3)') ll*gr%mesh%h(1)-origin, real(evectors(ll,jj)), & 
      !                                    & aimag(evectors(ll,jj))
             
      !   enddo
      !   call io_close(iunit)
      !  enddo
         
      !   write(filename,'(a,i6.6,a,i2.2)') 'Tdstep_',i,'/wavef_',mm
      !   iunit=io_open(filename,action='write')
         
      !   do jj=1, npointsx
      !    do ll=1, npointsx
      !write(iunit,'(i6.4,i6.4,es11.3, es11.3)') jj,ll,real(wavef(jj,ll)),aimag(wavef(jj,ll)) 
      !    enddo
      !   enddo

    end do states_loop

    deallocate(evectors)
    deallocate(evalues)
!    deallocate(wavef)
    deallocate(sqrdensity)
    deallocate(graddens)
    deallocate(potential)
    deallocate(density)
    deallocate(hartreep)
    deallocate(densmatr)

    deallocate (n1, n2, npoints)
    deallocate (labels_densmat)
    deallocate (particle_kept_densmat)
    deallocate (nnatorb_prt_densmat)
    call pop_sub()
  end subroutine density_matrix_write
  ! ---------------------------------------------------------

end module density_matrix_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
