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

module modelmb_density_m

  use comm_m
  use datasets_m
  use global_m
  use grid_m
  use hypercube_m
  use index_m
  use io_m
  use lalg_adv_m
  use loct_m
  use mesh_m
  use mesh_function_m
  use messages_m
  use modelmb_particles_m
  use modelmb_1part_m
  use mpi_m
  use mpi_lib_m
  use par_vec_m
  use parser_m
  use profiling_m
  use states_m

  implicit none

  private

  public :: modelmb_density_write,      &
            modelmb_density_init,       &
            modelmb_density_end,        &
            modelmb_density_nullify,    &
            dmodelmb_density_calculate, &
            zmodelmb_density_calculate, &
            modelmb_density_t

  type modelmb_density_t
    integer :: ndensities_to_calculate
    character(200) :: dirname
    character(80), pointer :: labels(:)
    integer, pointer :: particle_kept(:)
  end type modelmb_density_t

contains

  subroutine modelmb_density_init(dir, st, den)
    character(len=*),        intent(in)  :: dir
    type(states_t),          intent(in)  :: st
    type(modelmb_density_t), intent(out) :: den

    integer :: ncols, ipart
    type(block_t) :: blk

    PUSH_SUB(modelmb_density_init)

    ! the description for this variable is in modelmb_particles.F90
    if(parse_block(datasets_check('DensitiestoCalc'), blk) /= 0) then
      message(1) = 'To print out density, you must specify the DensitiestoCalc block in input'
      call messages_fatal(1)
    end if
   
    ncols = parse_block_cols(blk, 0)
    if(ncols /= 2 ) then
      call input_error("DensitiestoCalc")
    end if
    den%ndensities_to_calculate = parse_block_n(blk)
    if (den%ndensities_to_calculate < 0 .or. &
        den%ndensities_to_calculate > st%modelmbparticles%nparticle) then
      call input_error("DensitiestoCalc")
    end if

    SAFE_ALLOCATE(den%labels(1:den%ndensities_to_calculate))
    SAFE_ALLOCATE(den%particle_kept(1:den%ndensities_to_calculate))
   
    do ipart = 1, den%ndensities_to_calculate
      call parse_block_string(blk, ipart-1, 0, den%labels(ipart))
      call parse_block_integer(blk, ipart-1, 1, den%particle_kept(ipart))
   
      write (message(1),'(a,a)') 'labels_densities = ', den%labels(ipart)
      write (message(2),'(a,i6)') 'particle_kept_densities = ', den%particle_kept(ipart)
      call messages_info(2)
    end do
    call parse_block_end(blk)
    ! END reading in of input var block DensitiestoCalc

    den%dirname = trim(dir)

    POP_SUB(modelmb_density_init)

  end subroutine modelmb_density_init


  ! ---------------------------------------------------------
  subroutine modelmb_density_write(gr, st, wf, mm, den)
    type(grid_t),            intent(in) :: gr
    type(states_t),          intent(in) :: st
    CMPLX,                   intent(in) :: wf(1:gr%mesh%np_part_global)
    integer,                 intent(in) :: mm
    type(modelmb_density_t), intent(in) :: den

    integer :: jj, idir, iunit, ndims, ndim1part
    integer :: ikeeppart, idensities, nparticles_density
    integer, allocatable :: npoints(:)
    integer, allocatable :: ix_1part(:), ix_1part_p(:)
    character(len=200) :: filename
    FLOAT, allocatable :: density(:)

    type(modelmb_1part_t) :: mb_1part
    FLOAT, allocatable :: dipole_moment(:)

! FIXME: this should be separated into calculating and writing the density
    PUSH_SUB(modelmb_density_write)


    ! The algorithm should consider how many dimensions the wavefunction has (ndims),
    ! and how many (and which) dimensions should be integrated away.
    ndims = gr%sb%dim

    ndim1part = st%modelmbparticles%ndim

    call modelmb_1part_nullify(mb_1part)
    SAFE_ALLOCATE(  ix_1part(1:ndim1part))
    SAFE_ALLOCATE(ix_1part_p(1:ndim1part))
    SAFE_ALLOCATE(dipole_moment(1:ndim1part))

    ! Allocation of the arrays that store the limiting indices for each direction
    SAFE_ALLOCATE(npoints(1:ndims))
    do idir = 1, ndims
      npoints(idir) = gr%mesh%idx%ll(idir)
    end do


    ! loop over desired density matrices
    dens_loop: do idensities = 1, den%ndensities_to_calculate
      ikeeppart = den%particle_kept(idensities)
      nparticles_density = st%modelmbparticles%nparticles_per_type(st%modelmbparticles%particletype(ikeeppart))

      call modelmb_1part_init(mb_1part, gr%mesh, ikeeppart, ndim1part, gr%sb%box_offset)

      SAFE_ALLOCATE(density(1:mb_1part%npt_part))

      density = M_ZERO

      !   calculate the 1-particle density for this many-body state, and for the chosen
      !   particle being the free coordinate
      call zmodelmb_density_calculate(ikeeppart, mb_1part, nparticles_density, &
            gr%mesh, wf, density)

      ! Only node zero writes.
      ! mjv 14/3/2009: is this still at the right place in the file? None of
      ! this works in parallel yet...
      if(.not. mpi_grp_is_root(mpi_world)) cycle

      write(filename,'(a,i3.3,a,i2.2)') trim(den%dirname)//'/density_ip', ikeeppart,'_imb', mm
      iunit = io_open(filename,action='write')
      do jj = 1, mb_1part%npt
        call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, &
             mb_1part%enlarge_1part(1), jj, ix_1part)
        do idir = 1, ndim1part
          write(iunit,'(es11.3)', ADVANCE='no') ix_1part(idir) * mb_1part%h_1part(idir) + mb_1part%origin(idir)
        end do
        write(iunit,'(es18.10)') real(density(jj))
      end do
      call io_close(iunit)


      ! calculate dipole moment from density for this particle
      dipole_moment(:) = M_ZERO
      do jj = 1, mb_1part%npt_part
        call hypercube_i_to_x(mb_1part%hypercube_1part, ndim1part, mb_1part%nr_1part, &
             mb_1part%enlarge_1part(1), jj, ix_1part)
        dipole_moment = dipole_moment + (ix_1part(:) * mb_1part%h_1part(:) + mb_1part%origin(:)) &
                      * real(density(jj)) &
                      * st%modelmbparticles%charge_particle(ikeeppart)
      end do
      ! note: for eventual multiple particles in 4D (eg 8D total) this would fail to give the last values of dipole_moment
      write (message(1),'(a,I6,a,I6,a,I6)') 'For particle ', ikeeppart, ' of mb state ', mm
      write (message(2),'(a,3E20.10)') 'The dipole moment is (in a.u. = e bohr):     ', dipole_moment(1:min(3, ndim1part))
      write (message(3),'(a,E15.3)') '     with intrinsic numerical error usually <= ', 1.e-6 * mb_1part%npt
      call messages_info(3)
      ! dipole should be expressed in units_out%length

      SAFE_DEALLOCATE_A(density)
    
      call modelmb_1part_end(mb_1part)
    
    end do dens_loop ! loop over densities to output


    SAFE_DEALLOCATE_A(ix_1part)
    SAFE_DEALLOCATE_A(ix_1part_p)
    SAFE_DEALLOCATE_A(npoints)
    SAFE_DEALLOCATE_A(dipole_moment)

    POP_SUB(modelmb_density_write)
  end subroutine modelmb_density_write
  ! ---------------------------------------------------------


  subroutine modelmb_density_nullify(this)
    type(modelmb_density_t), intent(out) :: this

    PUSH_SUB(modelmb_density_nullify)
    nullify(this%labels)
    nullify(this%particle_kept)

    POP_SUB(modelmb_density_nullify)
  end subroutine modelmb_density_nullify

  ! ---------------------------------------------------------
  subroutine modelmb_density_end(this)
    type(modelmb_density_t), intent(inout) :: this

    PUSH_SUB(modelmb_density_end)
    SAFE_DEALLOCATE_P(this%labels)
    SAFE_DEALLOCATE_P(this%particle_kept)

    POP_SUB(modelmb_density_end)
  end subroutine modelmb_density_end

#include "undef.F90"
#include "real.F90"
#include "modelmb_density_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "modelmb_density_inc.F90"

end module modelmb_density_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
