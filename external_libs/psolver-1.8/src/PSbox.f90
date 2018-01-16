!> @file
!!    Modulefile for handling of the Parallelization of the Simulation box of the Poisson Solver
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2002-2015 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module PSbox
  use wrapper_MPI
  use PStypes, only: coulomb_operator, PSolver_energies
  use PSbase
  implicit none
  private

  interface PS_reduce
     module procedure reduce_scalar,reduce_array,reduce_energies
  end interface PS_reduce
  

  public :: PS_reduce,PS_gather

contains

  !>gather a distributed array to have a full array
  !!if only src is present this is assumed to be a full array
  !!otherwise it is assumed to be a distributed array
  subroutine PS_gather(src,kernel,dest,nsrc)
    use dynamic_memory, only: f_memcpy
    implicit none
    !> input array. If dest is present, the values are assumed to be distributed
    !!otherwise the values are not modified and are gathered in the dest
    !!array
    real(dp), dimension(*), intent(in) :: src
    type(coulomb_operator), intent(in) :: kernel
    real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3),*), intent(out), optional :: dest
    integer, intent(in), optional :: nsrc !< number of copies of the array src (useful for spin-polarized)
    !local variables
    integer :: ispin,nspin,isrc

    nspin=1
    if (present(nsrc)) nspin=nsrc
    if (present(dest)) then
       if (kernel%mpi_env%nproc > 1) then
          isrc=1
          do ispin=1,nspin
             call mpiallgather(src(isrc),recvbuf=dest(1,1,1,ispin),&
               recvcounts=kernel%counts,&
               displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
             isrc=isrc+kernel%grid%m1*kernel%grid%m3*kernel%grid%n3p
          end do
       else
          call f_memcpy(n=product(kernel%ndims)*nspin,src=src(1),dest=dest(1,1,1,1))
       end if
    else
       if (kernel%mpi_env%nproc > 1) then
          isrc=1
          do ispin=1,nspin
             call mpiallgather(src(isrc),recvcounts=kernel%counts,&
                  displs=kernel%displs,comm=kernel%mpi_env%mpi_comm)
             isrc=isrc+kernel%grid%m1*kernel%grid%m3*kernel%grid%n3p
          end do
       end if
    end if
  end subroutine PS_gather

  !> reduce all the given information 
  !! MPI_SUM is applied in the case of unspecified op
  subroutine reduce_scalar(val,kernel,op)
    implicit none
    real(dp), intent(inout) :: val
    type(coulomb_operator), intent(in) :: kernel
    integer, intent(in), optional :: op !< operation to be done
    !local variables
    integer :: mpi_op

    mpi_op=MPI_SUM
    if (present(op)) mpi_op=op

    if (kernel%mpi_env%nproc > 1) &
         call mpiallred(val,1,op=mpi_op,comm=kernel%mpi_env%mpi_comm)
    
  end subroutine reduce_scalar

  subroutine reduce_array(val,kernel,op)
    implicit none
    real(dp), dimension(:), intent(inout) :: val
    type(coulomb_operator), intent(in) :: kernel
    integer, intent(in), optional :: op !< operation to be done
    !local variables
    integer :: mpi_op

    mpi_op=MPI_SUM
    if (present(op)) mpi_op=op

    if (kernel%mpi_env%nproc > 1) &
         call mpiallred(val(1),size(val),op=mpi_op,comm=kernel%mpi_env%mpi_comm)

  end subroutine reduce_array

  !>this is of course to do the sum
  subroutine reduce_energies(e,kernel)
    type(PSolver_energies), intent(inout) :: e
    type(coulomb_operator), intent(in) :: kernel
    !local variables
    integer, parameter :: energsize=10
    real(gp), dimension(energsize) :: vals

    if (kernel%mpi_env%nproc > 1) then
       vals(1)   =e%hartree   
       vals(2)   =e%elec      
       vals(3)   =e%eVextra   
       vals(4)   =e%cavitation
       vals(5:10)=e%strten     
       call PS_reduce(vals,kernel)
       e%hartree   =vals(1)   
       e%elec      =vals(2)   
       e%eVextra   =vals(3)   
       e%cavitation=vals(4)   
       e%strten    =vals(5:10)
    end if

  end subroutine reduce_energies

  !>the invers of gathering. Transforms a full array in a distributed one
  subroutine PS_scatter(src,dest,kernel)
    implicit none
    type(coulomb_operator), intent(in) :: kernel
    real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(in) :: src
    real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), intent(out), optional :: dest

    call mpiscatterv(sendbuf=src,sendcounts=kernel%counts,&
         displs=kernel%displs,recvbuf=dest,comm=kernel%mpi_env%mpi_comm)
    
  end subroutine PS_scatter

  !>read a field in a distributed way
  !! only good to read the nspin==1 fields
  subroutine PS_read_field(filename,kernel,pot,npot)
    use yaml_output
    use dynamic_memory
    use IObox
    use dictionaries, only: f_err_throw
    implicit none
    character(len=*), intent(in) :: filename
    type(coulomb_operator), intent(in) :: kernel
    real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p,*), intent(out) :: pot
    integer, intent(in), optional :: npot
    !local variables
    logical :: wrong_spin
    integer :: nspin_dim, nspin,ispin
    character(len=1) :: geocode
    integer, dimension(3) :: ndims
    real(gp), dimension(3) :: hgrids
    real(dp), dimension(:,:,:,:), allocatable :: pot_from_disk

    if (kernel%mpi_env%iproc==0) then
       call yaml_map('Reading local potential from file:',filename)
       call read_field_dimensions(filename,geocode,ndims,nspin_dim)
       !allocate the potential in full
       pot_from_disk=f_malloc([ndims(1),ndims(2),ndims(3),nspin_dim],id='pot_from_disk')
       !> Read a density file using file format depending on the extension.
       call read_field(filename,&
            geocode,ndims,hgrids,nspin,product(ndims),nspin_dim,pot_from_disk)
       if (nspin/=nspin_dim) call f_err_throw('nspin/=nspin_dim')
    else
       pot_from_disk=f_malloc([1,1,1,1],id='pot_from_disk')
    end if

    wrong_spin= nspin/=1
    if (present(npot)) wrong_spin= npot /= nspin
    if (wrong_spin) call f_err_throw('Error in read_file: the npot is not correct')


    !then scatter the result
    do ispin=1,nspin
       call PS_scatter(pot_from_disk(1,1,1,ispin),pot(1,1,ispin),kernel)
    end do

    call f_free(pot_from_disk)
  end subroutine PS_read_field

  !> write a distributed array in the disk
  subroutine PS_dump_field(kernel,filename,src_dist,src_full,ixyz0)
    use dynamic_memory
    use IObox
    implicit none
    character(len=*), intent(in) :: filename
    type(coulomb_operator), intent(in) :: kernel
    real(dp), dimension(kernel%grid%m1,kernel%grid%m3*kernel%grid%n3p), optional :: src_dist
    real(dp), dimension(kernel%ndims(1),kernel%ndims(2),kernel%ndims(3)), intent(in), optional :: src_full
    integer, dimension(3), intent(in), optional :: ixyz0
    !local variables
    real(dp), dimension(:,:,:), allocatable :: work_full

    if (present(src_dist)) then
       work_full=f_malloc(kernel%ndims,id='work_full')
       call PS_gather(src=src_dist,dest=work_full,kernel=kernel)
       if (present(ixyz0)) then
          call dump_field(filename,kernel%geocode,kernel%ndims,kernel%hgrids,1,work_full,ixyz0=ixyz0)
       else
          call dump_field(filename,kernel%geocode,kernel%ndims,kernel%hgrids,1,work_full)
       end if
          call f_free(work_full)
    else if (present(src_full)) then
       if (present(ixyz0)) then
          call dump_field(filename,kernel%geocode,kernel%ndims,kernel%hgrids,1,src_full,ixyz0=ixyz0)
       else
          call dump_field(filename,kernel%geocode,kernel%ndims,kernel%hgrids,1,src_full)
       end if
    end if

  end subroutine PS_dump_field

!!$  !> Read the densit and put the values in the rhopot arrays according to the parallelization indicated by
!!$  !! nscatterarr array
!!$  subroutine read_potential_from_disk(iproc,nproc,filename,geocode,ngatherarr,n1i,n2i,n3i,n3p,nspin,hxh,hyh,hzh,pot)
!!$    use dynamic_memory
!!$    implicit none
!!$    integer, intent(in) :: iproc,nproc,n1i,n2i,n3i,n3p,nspin
!!$    real(gp), intent(in) :: hxh,hyh,hzh
!!$    character(len=*), intent(in) :: filename
!!$    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
!!$    integer, dimension(0:nproc-1,2), intent(in) :: ngatherarr
!!$    real(dp), dimension(n1i,n2i,max(n3p,1),nspin), intent(out) :: pot
!!$    !local variables
!!$    character(len=*), parameter :: subname='read_potential_from_disk'
!!$    integer :: n1t,n2t,n3t,nspint,ierror,ierr,ispin
!!$    real(gp) :: hxt,hyt,hzt
!!$    real(dp), dimension(:,:,:,:), pointer :: pot_from_disk
!!$
!!$    !only the first processor should read this
!!$    if (iproc == 0) then
!!$       write(*,'(1x,a)')'Reading local potential from file:'//trim(filename)
!!$       call read_density(trim(filename),geocode,&
!!$            n1t,n2t,n3t,nspint,hxt,hyt,hzt,pot_from_disk)
!!$       if (abs(hxt-hxh) <= 1.e-5_gp .and. abs(hyt-hyh) <= 1.e-5_gp .and. abs(hzt-hzh) <= 1.e-5_gp .and. &
!!$            nspint == nspin .and. &
!!$            n1i  == n1t  .and. n2i == n2t .and. n3i == n3t) then
!!$       else
!!$          write(*,*)'ERROR (to be documented): some of the parameters do not coincide'
!!$          write(*,*)hxh,hyh,hzh,hxt,hyt,hzt,nspin,nspint,n1i,n2i,n3i,n1t,n2t,n3t
!!$       end if
!!$    else
!!$       pot_from_disk = f_malloc_ptr((/ 1, 1, 1, nspin /),id='pot_from_disk')
!!$    end if
!!$
!!$    if (nproc > 1) then
!!$       do ispin=1,nspin
!!$          call MPI_SCATTERV(pot_from_disk(1,1,1,ispin),&
!!$               ngatherarr(0,1),ngatherarr(0,2),mpidtypd, &
!!$               pot(1,1,1,ispin),&
!!$               n1i*n2i*n3p,mpidtypd,0,%mpi_comm,ierr)
!!$       end do
!!$    else
!!$       call vcopy(n1i*n2i*n3i*nspin,pot_from_disk(1,1,1,1),1,pot(1,1,1,1),1)
!!$    end if
!!$
!!$    call f_free_ptr(pot_from_disk)
!!$
!!$  end subroutine read_potential_from_disk


end module PSbox
