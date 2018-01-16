module IOBoxETSF

  implicit none

  private

  public :: write_etsf_density
  public :: read_etsf

contains

  !> Write a field in the ISF basis in the ETSF format
  subroutine write_etsf_density(filename,message,geocode,ndims,hgrids,&
       x,nspin,nat,rxyz,iatype,ntypes,nzatom)
    use PSbase
    implicit none
    character(len=*), intent(in) :: filename,message
    character(len = 1), intent(in) :: geocode
    integer, intent(in) :: nat, ntypes, nspin
    integer, dimension(3), intent(in) :: ndims
    real(gp), dimension(3), intent(in) :: hgrids
    real(dp), dimension(ndims(1),ndims(2),ndims(3),nspin), target, intent(in) :: x
    real(gp), dimension(3,nat), intent(in) :: rxyz
    integer, dimension(nat), target, intent(in) :: iatype
    integer, dimension(ntypes), intent(in) :: nzatom

    !local variables

    write(0, "(A)") "Illegal call to write_etsf_density(), not compiled with ETSF_IO support."
    stop

    !To avoid warnings from the compiler
    write(*,*) filename,message,ndims,hgrids,x(1,1,1,1),rxyz(1,1)
  END SUBROUTINE write_etsf_density
  
  !> @file
  !!    File for handling the read-write of a given simulation in ETSF format, fake version
  !!
  !! @author
  !!    G. Fisicaro, L. Genovese (September 2015)
  !!    Copyright (C) 2002-2015 BigDFT group 
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS 
  !> Read a field in the ISF basis in the ETSF format
  subroutine read_etsf(filename,geocode,n1i,n2i,n3i,nspin,hxh,hyh,hzh,ldrho,nrho,rho,&
       nat,rxyz,iatypes,znucl)
    use PSbase
    implicit none
    character(len=*), intent(in) :: filename
    character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
    integer, intent(in) :: ldrho,nrho !<dimensions of the rho array
    integer, intent(out) :: nspin
    integer, intent(out) ::  n1i,n2i,n3i
    real(gp), intent(out) :: hxh,hyh,hzh
    real(dp), dimension(ldrho,nrho) :: rho
    real(gp), dimension(:,:), pointer :: rxyz
    integer, intent(out) ::  nat
    integer, dimension(:), pointer :: iatypes, znucl
    write(0, "(A)") "Illegal call to read_etsf(), not compiled with ETSF_IO support."
    stop

    !To avoid warnings from the compiler
    write(*,*) filename,geocode,nspin,n1i,n2i,n3i
    n1i=0
    n2i=0
    n3i=0
    hxh=0.0_gp
    hyh=0.0_gp
    hzh=0.0_gp
    rho(1,1)=0.0_gp
    rxyz(1,1)=0.0_gp
    nat=0
    nspin=0
  END SUBROUTINE read_etsf

end module IOBoxETSF
