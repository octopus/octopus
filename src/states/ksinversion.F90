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
  use modelmb_particles_m
  use global_m
  use grid_m
  use hypercube_m
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
  subroutine kspotential_inversion_write(dir, gr, st, modelmbparticles)
    character(len=*), intent(in) :: dir
    type(grid_t),     intent(in) :: gr
    type(states_t),   intent(in) :: st
    type(modelmb_particle_t),   intent(in) :: modelmbparticles

    integer :: mm, jj, ll, j, err_code, iunit, ndims, ndim1part
    integer :: ipart,nkspot_to_calculate,ncols
    integer :: ikeeppart,npt_1part, idir, irealdir
    integer :: ikspot
    integer, allocatable :: npoints(:)
    integer, allocatable :: nr_1part(:,:)
    integer, allocatable :: ix_1part(:), ix_1part_p(:)
    integer, allocatable :: enlarge_1part(:)
    logical :: bof
    character(len=200) :: dirname, filename
    FLOAT :: vol_elem_1part
    FLOAT, allocatable :: origin(:), h_1part(:)
    CMPLX, allocatable :: densmatr(:, :), evectors(:, :), wavef(:,:)
    FLOAT, allocatable :: evalues(:), sqrdensity(:), density(:), graddens(:)
    FLOAT, allocatable :: hartreep(:), potential(:)

    type(hypercube_t) :: hypercube_1part

    type(block_t) :: blk
    character(80), allocatable :: labels_kspot(:)
    integer, allocatable :: particle_kept_dens(:)
    
    call push_sub('states.density_matrix_write')

    !%Variable KSPotentialtoCalc
    !%Type block
    !%Section States
    !%Description
    !% choice of which KS potential will be calculated and output, in the
    !%  modelmb particles scheme (the corresponding density is also output)
    !%
    !% <tt>%KSPotentialtoCalc
    !% <br>&nbsp;&nbsp; proton   | 1 
    !% <br>&nbsp;&nbsp; electron | 2 
    !% <br>%</tt>
    !%
    !% would ask octopus to calculate the KS potential corresponding to the 1st
    !% particle (whose coordinates correspond to dimensions 1 to ndim_modelmb),
    !% which is a proton, then that corresponding to the 2nd particle
    !% (electron with dimensions ndim_modelmb+1 to 2*ndim_modelmb)
    !%
    !%End
   
    if(loct_parse_block(datasets_check('KSPotentialtoCalc'), blk)/=0) then
     message(1) = 'To print out KS potentials, you must specify the KSPotentialtoCalc block in input'
     call write_fatal(1)
    end if
   
    ncols = loct_parse_block_cols(blk, 0)
    if(ncols /= 2 ) then
      call input_error("KSPotentialtoCalc")
    end if
    nkspot_to_calculate=loct_parse_block_n(blk)
    if (nkspot_to_calculate < 0 .or. &
        nkspot_to_calculate > modelmbparticles%nparticle_modelmb) then
      call input_error("KSPotentialtoCalc")
    end if

    ALLOCATE (labels_kspot(nkspot_to_calculate),nkspot_to_calculate)
    ALLOCATE (particle_kept_dens(nkspot_to_calculate),nkspot_to_calculate)
   
    do ipart=1,nkspot_to_calculate
      call loct_parse_block_string(blk, ipart-1, 0, labels_kspot(ipart))
      call loct_parse_block_int(blk, ipart-1, 1, particle_kept_dens(ipart))
   
      write (message(1),'(a,a)') 'labels_kspot = ', labels_kspot(ipart)
      write (message(2),'(a,i6)') 'particle_kept_dens = ', particle_kept_dens(ipart)
      call write_info(2)
    end do
    call loct_parse_block_end(blk)
! END reading in of input var block

    ! This is the root directory where everything will be written.
    dirname = trim(dir)//'/KS-potential'
    call loct_mkdir(trim(dirname))

    ! The algorithm should consider how many dimensions the wavefunction has (ndims),
    ! and how many (and which) dimensions should be integrated away.
    ndims = gr%sb%dim

    ! Allocatation of the arrays that store the limiting indexes for each direction
    ALLOCATE(npoints(ndims), ndims)
    do j = 1, ndims
      npoints(j) = gr%mesh%idx%ll(j)
      write (*,*) 'npoints, nr = ', npoints(j), gr%mesh%idx%nr(:,j)
    end do

    ndim1part=modelmbparticles%ndim_modelmb

    ALLOCATE(origin(ndim1part), ndim1part)
    ALLOCATE(ix_1part(ndim1part), ndim1part)
    ALLOCATE(ix_1part_p(ndim1part), ndim1part)
    ALLOCATE(nr_1part(2,ndim1part), 2*ndim1part)
    ALLOCATE(h_1part(ndim1part), ndim1part)
    ALLOCATE(enlarge_1part(ndim1part), ndim1part)

! for the moment force keeping the 1st particle coordinate(s)
!  this will be made into a loop over desired density matrices
    kspot_loop: do ikspot = 1, nkspot_to_calculate
      ikeeppart = particle_kept_dens(idensmat)

!   get full size of arrays for 1 particle only in ndim_modelmb dimensions
      npt_1part=1
      do idir = 1, ndim1part
        npt_1part = npt_1part*npoints((ikeeppart - 1)*ndim1part + idir)
      end do

      ALLOCATE(density(npt_1part), npt_1part)
      ALLOCATE(hartreep(npt_1part), npt_1part)
      ALLOCATE(potential(npt_1part), npt_1part)

      write (*,*) 'shape(density) ', shape(density)
      write (*,*) ' st%nst = ', st%nst

!   volume element for the chosen particle
      vol_elem_1part=1.0d0
      do idir=1,ndim1part
        irealdir=(ikeeppart-1)*ndim1part + idir
        vol_elem_1part=vol_elem_1part*gr%mesh%h(irealdir)
        h_1part(idir) = gr%mesh%h(irealdir)
      end do

!   store start and end positions for the relevant dimensions for this particle
      nr_1part(:,:)=gr%mesh%idx%nr(:,(ikeeppart-1)*ndim1part+1:ikeeppart*ndim1part)

!   initialize a hypercube for just this particle
!   NB: hypercube_* presume that enlarge is the same for all dimensions!
      enlarge_1part=gr%mesh%idx%enlarge((ikeeppart-1)*ndim1part+1:ikeeppart*ndim1part)
      call hypercube_init(hypercube_1part, ndim1part, nr_1part, enlarge_1part(1))

      ! not always the real origin if the box is shifted, no?
      !  which happens to be Matthieu's case...
      !  only important for printout, so it is ok
      do idir=1,ndim1part
        irealdir=(ikeeppart-1)*ndim1part + idir
        origin(idir)=(npoints(irealdir)/2+1)*gr%mesh%h(irealdir)
      end do

        density   = M_ZERO
        potential = M_ZERO
        hartreep  = M_ZERO

!   calculate the 1 particle density matrix for this Many Body state, and for the chosen
!   particle being the free coordinate
        if(states_are_real(st)) then
          call dmf_calculate_gamma(ikeeppart, ndim1part,&
                hypercube_1part, npt_1part, nr_1part, enlarge_1part(1),&
                gr%mesh, st%dpsi(:, 1, mm, 1), densmatr)
        else
          call zmf_calculate_gamma(ikeeppart, ndim1part,&
                hypercube_1part, npt_1part, nr_1part, enlarge_1part(1),&
                gr%mesh, st%zpsi(:, 1, mm, 1), densmatr)
        end if

        ! Only node zero writes.
        if(.not. mpi_grp_is_root(mpi_world)) cycle

!        do ll = 1, npt_1part
!          do jj = 1, npt_1part
!            hartreep(ll) = hartreep(ll)+density(jj)/(sqrt(((ll-jj)*gr%mesh%h(1))**2+1))
!          end do
!          hartreep(ll) = 0.5*hartreep(ll)*gr%mesh%h(1)
!        end do

        write(filename,'(a,i3.3,a,i2.2)') trim(dirname)//'/density_', ikeeppart,'_', mm
        iunit = io_open(filename,action='write')
        do jj = 1, npt_1part
          call hypercube_i_to_x(hypercube_1part, ndim1part, nr_1part, 0, jj, ix_1part)
          do idir=1,ndim1part
            write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir)*h_1part(idir)-origin(idir)
          end do
          write(iunit,'(es11.3,es11.3)') real(densmatr(jj,jj)), aimag(densmatr(jj,jj))
        end do
        call io_close(iunit)

!        write(filename,'(a,i3.3,a,i2.2)') trim(dirname)//'/potential_', ikeeppart,'_', mm
!        iunit = io_open(filename,action='write')
!        do ll=1, npt_1part
!          call hypercube_i_to_x(hypercube_1part, ndim1part, nr_1part, 0, ll, ix_1part)
!          do idir=1,ndim1part
!            write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir)*h_1part(idir)-origin(idir)
!          end do
!          write(iunit,'(es11.3,es11.3,es11.3)') potential(ll), hartreep(ll), potential(ll)-hartreep(ll)
!        end do
!        call io_close(iunit)

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

      deallocate(potential)
      deallocate(density)
      deallocate(hartreep)
      deallocate(densmatr)
      
    end do kspot_loop ! loop over densmats to output

    deallocate(origin)
    deallocate(ix_1part)
    deallocate(ix_1part_p)
    deallocate(nr_1part)
    deallocate(enlarge_1part)

    deallocate (npoints)
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
