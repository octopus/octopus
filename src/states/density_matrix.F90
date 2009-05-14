!! Copyright (C) 2009 N. Helbig and M. Verstraete
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
  use modelmb_particles_m
  use modelmb_1part_m
  use mpi_m
  use mpi_lib_m
  use profiling_m
  use states_m

  implicit none

  private

  public :: density_matrix_write, &
            density_matrix_init, &
            density_matrix_end, &
            density_matrix_nullify, &
            modelmb_denmat_t

  type modelmb_denmat_t
    integer :: ndensmat_to_calculate
    character(len=200) :: dirname

    character(80), pointer :: labels_densmat(:)
    integer, pointer :: particle_kept_densmat(:)
    integer, pointer :: nnatorb_prt_densmat(:)
  end type modelmb_denmat_t

contains

  subroutine density_matrix_init(dir, st, denmat)

    implicit none

    character(len=*),       intent(in) :: dir
    type(states_t),         intent(in) :: st
    type(modelmb_denmat_t), intent(out) :: denmat

    ! local vars
    integer :: ncols, ipart
    type(block_t) :: blk

    call push_sub('states.density_matrix_init')

    !%Variable DensityMatricestoCalc
    !%Type block
    !%Section States
    !%Description
    !% choice of which particle density matrices will be calculated and output, in the
    !%  modelmb particles scheme (the corresponding density is also output)
    !%
    !% <tt>%DensityMatricestoCalc
    !% <br>&nbsp;&nbsp; proton   | 1 | 10
    !% <br>&nbsp;&nbsp; electron | 2 | 15
    !% <br>%</tt>
    !%
    !% would ask octopus to calculate the density matrix corresponding to the 1st
    !% particle (whose coordinates correspond to dimensions 1 to ndim_modelmb),
    !% which is an proton, then that corresponding to the 2nd particle
    !% (electron with dimensions ndim_modelmb+1 to 2*ndim_modelmb), printing
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
    denmat%ndensmat_to_calculate=loct_parse_block_n(blk)
    if (denmat%ndensmat_to_calculate < 0 .or. &
        denmat%ndensmat_to_calculate > st%modelMBparticles%nparticle_modelmb) then
      call input_error("DensityMatricestoCalc")
    end if

    SAFE_ALLOCATE(denmat%labels_densmat(1:denmat%ndensmat_to_calculate))
    SAFE_ALLOCATE(denmat%particle_kept_densmat(1:denmat%ndensmat_to_calculate))
    SAFE_ALLOCATE(denmat%nnatorb_prt_densmat(1:denmat%ndensmat_to_calculate))
   
    do ipart=1,denmat%ndensmat_to_calculate
      call loct_parse_block_string(blk, ipart-1, 0, denmat%labels_densmat(ipart))
      call loct_parse_block_int(blk, ipart-1, 1, denmat%particle_kept_densmat(ipart))
      call loct_parse_block_int(blk, ipart-1, 2, denmat%nnatorb_prt_densmat(ipart))
   
      write (message(1),'(a,a)') 'labels_densmat = ', denmat%labels_densmat(ipart)
      write (message(2),'(a,i6)') 'particle_kept_densmat = ', denmat%particle_kept_densmat(ipart)
      write (message(3),'(a,i6)') 'nnatorb_prt_densmat = ', denmat%nnatorb_prt_densmat(ipart)
      call write_info(3)
    end do
    call loct_parse_block_end(blk)
    ! END reading in of input var block DensityMatricestoCalc

    ! ensure parent directory exists
    call io_mkdir(dir)

    ! This is the root directory where everything will be written.
    denmat%dirname = trim(dir)//'/density-matrix'
    call loct_mkdir(trim(denmat%dirname))
    
    call pop_sub()

  end subroutine density_matrix_init


  ! ---------------------------------------------------------
  subroutine density_matrix_write(gr, st, wf, mm, denmat)
    type(grid_t),     intent(in) :: gr
    type(states_t),   intent(in) :: st
    CMPLX,            intent(in) :: wf(1:gr%mesh%np_part_global)
    integer,          intent(in) :: mm
    type(modelmb_denmat_t), intent(in) :: denmat

    integer :: jj, ll, j, err_code, iunit, ndims, ndim1part
    integer :: ikeeppart, idir
    integer :: idensmat, nparticles_densmat
    integer, allocatable :: npoints(:)
    integer, allocatable :: ix_1part(:), ix_1part_p(:)
    logical :: bof
    character(len=200) :: filename
    CMPLX, allocatable :: densmatr(:, :), evectors(:, :)
    CMPLX, allocatable :: densmatr_tmp(:, :)
    FLOAT, allocatable :: evalues(:), density(:)

    type(modelmb_1part_t) :: mb_1part
    FLOAT, allocatable :: dipole_moment(:)

    call push_sub('states.density_matrix_write')


    ! The algorithm should consider how many dimensions the wavefunction has (ndims),
    ! and how many (and which) dimensions should be integrated away.
    ndims = gr%sb%dim

    ndim1part=st%modelMBparticles%ndim_modelmb

    call modelmb_1part_nullify(mb_1part)
    SAFE_ALLOCATE(  ix_1part(1:ndim1part))
    SAFE_ALLOCATE(ix_1part_p(1:ndim1part))
    SAFE_ALLOCATE(dipole_moment(1:ndim1part))

    ! Allocatation of the arrays that store the limiting indexes for each direction
    SAFE_ALLOCATE(npoints(1:ndims))
    do j = 1, ndims
      npoints(j) = gr%mesh%idx%ll(j)
    end do


    ! loop over desired density matrices
    densmat_loop: do idensmat=1,denmat%ndensmat_to_calculate
      ikeeppart=denmat%particle_kept_densmat(idensmat)
      nparticles_densmat = st%modelMBparticles%nparticles_per_type(st%modelMBparticles%particletype_modelMB(ikeeppart))

      call modelmb_1part_init(mb_1part, gr%mesh, ikeeppart, ndim1part, gr%sb%box_offset)

      SAFE_ALLOCATE(densmatr(1:mb_1part%npt_1part, 1:mb_1part%npt_1part))
      SAFE_ALLOCATE(evectors(1:mb_1part%npt_1part, 1:mb_1part%npt_1part))
      SAFE_ALLOCATE(evalues(1:mb_1part%npt_1part))
      SAFE_ALLOCATE(density(1:mb_1part%npt_1part))

      
      densmatr  = M_z0

      !   calculate the 1 particle density matrix for this Many Body state, and for the chosen
      !   particle being the free coordinate
      call zmf_calculate_gamma(ikeeppart, mb_1part, nparticles_densmat, &
            gr%mesh, wf, densmatr)

      ! Only node zero writes.
      ! mjv 14/3/2009: is this still at the right place in the file? None of
      ! this works in parallel yet...
      if(.not. mpi_grp_is_root(mpi_world)) cycle

      !Diagonalize the density matrix
      bof=.true.
      SAFE_ALLOCATE(densmatr_tmp(1:mb_1part%npt_1part, 1:mb_1part%npt_1part))
      densmatr_tmp=densmatr
      ! CHECK: should we only be diagonalizing the main grid points, as
      ! opposed to the full mb_1part%npt_1part?
      call lalg_eigensolve(mb_1part%npt_1part, densmatr_tmp, evectors, evalues, bof, err_code)
      SAFE_DEALLOCATE_A(densmatr_tmp)
    
      !NOTE: The highest eigenvalues are the last ones not the first!!!
      !      Writing is therefore in reverse order
      evectors = evectors/sqrt(mb_1part%vol_elem_1part)
      evalues  = evalues*mb_1part%vol_elem_1part

      !Write everything into files
      write(filename,'(a,i3.3,a,i2.2)') trim(denmat%dirname)//'/occnumb_ip',ikeeppart,'_iMB',mm
      iunit = io_open(trim(filename), action='write')

      do jj = mb_1part%npt_1part, 1, -1
        write(iunit,'(i4.4,es11.3)') mb_1part%npt_1part-jj+1, evalues(jj)
      end do
          
      call io_close(iunit)

      do jj = mb_1part%npt_1part-denmat%nnatorb_prt_densmat(ikeeppart)+1, mb_1part%npt_1part
        write(filename,'(a,i3.3,a,i2.2,a,i4.4)') trim(denmat%dirname)//'/natorb_ip', ikeeppart,'_iMB', mm, '_', mb_1part%npt_1part-jj+1
        iunit = io_open(filename, action='write')
        do ll = 1, mb_1part%npt_1part
          call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, mb_1part%enlarge_1part(1), ll, ix_1part)
          do idir=1,ndim1part
            write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
          end do
          write(iunit,'(es11.3,es11.3)') real(evectors(ll,jj)), aimag(evectors(ll,jj))
        end do
      call io_close(iunit)
      end do

      write(filename,'(a,i3.3,a,i2.2)') trim(denmat%dirname)//'/densmatr_ip', ikeeppart,'_iMB', mm
      iunit = io_open(filename,action='write')
      do jj = 1, mb_1part%npt_1part
        call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, mb_1part%enlarge_1part(1), jj, ix_1part)
        do ll = 1, mb_1part%npt_1part
          call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, mb_1part%enlarge_1part(1), ll, ix_1part_p)
          do idir=1,ndim1part
            write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
          end do
          do idir=1,ndim1part
            write(iunit,'(es11.3)', ADVANCE='no') ix_1part_p(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
          end do
          write(iunit,'(es11.3,es11.3)') real(densmatr(jj,ll)), aimag(densmatr(jj,ll))
        end do
        write(iunit,*)
      end do
      call io_close(iunit)

      write(filename,'(a,i3.3,a,i2.2)') trim(denmat%dirname)//'/density_ip', ikeeppart,'_iMB', mm
      iunit = io_open(filename,action='write')
      do jj = 1, mb_1part%npt_1part
        call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, mb_1part%enlarge_1part(1), jj, ix_1part)
        do idir=1,ndim1part
          write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir)*mb_1part%h_1part(idir)+mb_1part%origin(idir)
        end do
        write(iunit,'(es18.10)') real(densmatr(jj,jj))
      end do
      call io_close(iunit)


      ! calculate dipole moment from density for this particle
      dipole_moment(:) = 0.0d0
      do jj = 1,mb_1part%npt_1part
        call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, mb_1part%enlarge_1part(1), jj, ix_1part)
        dipole_moment = dipole_moment+(ix_1part(:)*mb_1part%h_1part(:)+mb_1part%origin(:))*real(densmatr(jj,jj))*&
                      st%modelMBparticles%charge_particle_modelMB(ikeeppart)
      end do
      ! note: for eventual multiple particles in 4D (eg 8D total) this would fail to give the last values of dipole_moment
      write (message(1),'(a,I6,a,I6,a,I6)') 'For particle ', ikeeppart, ' of MB state ', mm
      write (message(2),'(a,3E20.10)') 'The dipole moment is (in a.u. = e bohr):     ', dipole_moment(1:min(3,ndim1part))
      write (message(3),'(a,E15.3)') '     with intrinsic numerical error usually <= ', 1.e-6*mb_1part%npt_1part
      call write_info(3)

      SAFE_DEALLOCATE_A(evectors)
      SAFE_DEALLOCATE_A(evalues)
      SAFE_DEALLOCATE_A(density)
      SAFE_DEALLOCATE_A(densmatr)
    
      call modelmb_1part_end(mb_1part)
    
    end do densmat_loop ! loop over densmats to output


    SAFE_DEALLOCATE_A(ix_1part)
    SAFE_DEALLOCATE_A(ix_1part_p)
    SAFE_DEALLOCATE_A(npoints)
    SAFE_DEALLOCATE_A(dipole_moment)

    call pop_sub()
  end subroutine density_matrix_write
  ! ---------------------------------------------------------


  subroutine density_matrix_nullify(this)
    implicit none
    type(modelmb_denmat_t), intent(in) :: this
    call push_sub('states.density_matrix_nullify')
    nullify(this%labels_densmat)
    nullify(this%particle_kept_densmat)
    nullify(this%nnatorb_prt_densmat)
    call pop_sub()
  end subroutine density_matrix_nullify

  ! ---------------------------------------------------------
  subroutine density_matrix_end(this)
    implicit none
    type(modelmb_denmat_t), intent(in) :: this
    call push_sub('states.density_matrix_end')
    SAFE_DEALLOCATE_P(this%labels_densmat)
    SAFE_DEALLOCATE_P(this%particle_kept_densmat)
    SAFE_DEALLOCATE_P(this%nnatorb_prt_densmat)
    call pop_sub()
  end subroutine density_matrix_end

end module density_matrix_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
