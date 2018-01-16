!> @file
!!    Modulefile for handling the read-write of a given simulation
!!    box
!!
!! @author
!!    G. Fisicaro, L. Genovese (September 2015)
!!    Copyright (C) 2002-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module IObox
  use PSbase
  implicit none

  integer, parameter :: UNKNOWN=0
  integer, parameter :: CUBE=1
  integer, parameter :: ETSF=2
  integer, parameter :: POT=3

  private


  public :: read_field,read_field_dimensions,dump_field,field_multipoles

  contains

    !>Calculate a set of multipoles from a given box iterator.
    !! Defined up to quadrupoles
    subroutine field_multipoles(bit,field,nfield,mp)
      use box
      use f_harmonics
      implicit none
      integer, intent(in) :: nfield !<number of components of the field to reduce on (useful for spin)
      type(box_iterator), intent(inout) :: bit
      real(dp), dimension(bit%mesh%ndims(1)*bit%mesh%ndims(2)*(bit%i3e-bit%i3s+1),nfield), intent(in) :: field
      type(f_multipoles), intent(inout) :: mp !<multipole structure

      !local variables
      integer :: ispin
      real(dp) :: q

      !the multipoles have to be created and initialized before.
      do ispin=1,nfield
         do while(box_next_point(bit))
            q= field(bit%ind,ispin)*bit%mesh%volume_element
            bit%tmp=closest_r(bit%mesh,bit%rxyz,mp%rxyz)
            call f_multipoles_accumulate(mp%Q,mp%lmax,bit%tmp,q)
         end do
      end do

    end subroutine field_multipoles
    
    pure subroutine cube_dimensions(geocode,ndims,nc1,nc2,nc3)
      implicit none
      character(len=1), intent(in) :: geocode
      integer, dimension(3), intent(in) :: ndims
      integer, intent(out) :: nc1,nc2,nc3
      !value of the buffer in the x and z direction
      if (geocode /= 'F') then
         nc1=ndims(1)
         nc3=ndims(3)
      else
         nc1=ndims(1)-31
         nc3=ndims(3)-31
      end if
      !value of the buffer in the y direction
      if (geocode == 'P') then
         nc2=ndims(2)
      else
         nc2=ndims(2)-31
      end if
    end subroutine cube_dimensions

    pure subroutine startend_buffers(geocode,nl1,nl2,nl3,nbx,nby,nbz)
      implicit none
      character(len=1), intent(in) :: geocode
      integer, intent(out) :: nl1,nl2,nl3,nbx,nby,nbz

      if (geocode /= 'F') then
         nl1=1
         nl3=1
         nbx = 1
         nbz = 1
      else
         nl1=15
         nl3=15
         nbx = 0
         nbz = 0
      end if
      !value of the buffer in the y direction
      if (geocode == 'P') then
         nl2=1
         nby = 1
      else
         nl2=15
         nby = 0
      end if

    end subroutine startend_buffers

    function get_file_format(filename,isuffix) result(fformat)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(out) :: isuffix
      integer :: fformat
      !local variables
      integer :: ipos,iext
      character(len=8), dimension(4), parameter :: &
           exts=['.cube   ','.etsf   ','.etsf.nc','.pot    ']

      ! Format = 1 -> cube (default)
      ! Format = 2 -> ETSF
      ! ...
      fformat = UNKNOWN
      identify: do iext=1,size(exts)
         isuffix = index(filename, trim(exts(iext)), back = .true.)
         if (isuffix >0) exit identify
      end do identify
      if (isuffix > 0) then
         isuffix = isuffix - 1
         select case(iext)
         case(1)
            fformat=CUBE
         case(2,3)
            fformat=ETSF
         case(4)
            fformat=POT
         end select
      else
         isuffix = len(filename)
      end if

!!$      isuffix = index(filename, ".cube", back = .true.)
!!$      if (isuffix > 0) then
!!$         fformat = CUBE
!!$      else
!!$         isuffix = index(filename, ".etsf", back = .true.)
!!$         if (isuffix <= 0) isuffix = index(filename, ".etsf.nc", back = .true.)
!!$         if (isuffix > 0) fformat = ETSF
!!$      end if
!!$      if (isuffix > 0) then
!!$         isuffix = isuffix - 1
!!$      else
!!$         isuffix = len(filename)
!!$      end if

      !fallback to CUBE format if there is no other extension for the file
      if (fformat == UNKNOWN) then
         !eliminate slashes
         ipos=index(filename, "/",back=.true.)
         ipos=max(ipos,1)
         if (index(filename(ipos:), ".") == 0) fformat=CUBE
      end if

    end function get_file_format

    subroutine read_cube_header(filename,geocode,ndims,hgrids,&
         nat,rxyz, iatypes, znucl)
      use f_utils
      use dynamic_memory
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, dimension(3), intent(out) :: ndims
      real(gp), dimension(3), intent(out) :: hgrids
      integer, intent(out) :: nat
      real(gp), dimension(:,:), pointer :: rxyz
      integer, dimension(:), pointer :: iatypes, znucl
      !local variables
      character(len=*), parameter :: subname='read_cube_header'
      integer :: n1t,n2t,n3t,n1,n2,n3,idum,iat,j,unt
      integer :: nl1,nl2,nl3,nbx,nby,nbz,n1i,n2i,n3i
      real(gp) :: dum1,dum2,dum3,hxh,hyh,hzh
      integer, dimension(:), allocatable :: znucl_

      call startend_buffers(geocode,nl1,nl2,nl3,nbx,nby,nbz)
      unt=22
      call f_open_file(unit=unt,file=trim(filename)//".cube",status='old')
      read(unt,*)! 'CUBE file for charge density'
      read(unt,*)! 'Case for '//trim(message)

      read(unt,'(i5,3(f12.6))') nat, dum1, dum2, dum3
      read(unt,'(i5,3(f12.6))') n1t , hxh  , dum1 , dum2
      read(unt,'(i5,3(f12.6))') n2t , dum1 , hyh  , dum2
      read(unt,'(i5,3(f12.6))') n3t , dum1 , dum2 , hzh

      !grid positions
      n1=n1t/2-nbx
      n1i=2*n1+(1-nbx)+2*nl1
      n2=n2t/2-nby
      n2i=2*n2+(1-nby)+2*nl2
      n3=n3t/2-nbz
      n3i=2*n3+(1-nbz)+2*nl3

      !atomic positions
      rxyz = f_malloc_ptr((/ 3, nat /),id='rxyz')
      iatypes = f_malloc_ptr(nat,id='iatypes')
      znucl_ = f_malloc(nat,id='znucl_')
      znucl_(:) = -1

      do iat=1,nat
         read(unt,'(i5,4(f12.6))') idum , dum1 , (rxyz(j,iat),j=1,3)
         do j = 1, nat, 1
            if (znucl_(j) == idum .or. znucl_(j) == -1) then
               znucl_(j) = idum
               exit
            end if
         end do
         iatypes(iat) = j
      end do

      do j = 1, nat, 1
         if (znucl_(j) == -1) then
            exit
         end if
      end do
      znucl = f_malloc_ptr(j-1,id='znucl')
      znucl(1:j-1) = znucl_(1:j-1)

      call f_free(znucl_)

      call f_close(unt)

      ndims(1)=n1i
      ndims(2)=n2i
      ndims(3)=n3i
      hgrids(1)=hxh
      hgrids(2)=hyh
      hgrids(3)=hzh

    END SUBROUTINE read_cube_header

    !>   Read a cube field which have been plotted previously by write_cube_fields
    subroutine read_cube_field(filename,geocode,ndims,rho)
      use f_utils
      use dynamic_memory
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, dimension(3), intent(in) :: ndims
      real(dp), dimension(ndims(1),ndims(2),ndims(3)) :: rho
      !local variables
      !n(c) character(len=*), parameter :: subname='read_cube_field'
      character(len=3) :: advancestring
      integer :: n1t,n2t,n3t,n1,n2,n3,i1,i2,i3,nat,iat,unt
      integer :: nl1,nl2,nl3,nbx,nby,nbz,icount,ind,n1i,n2i,n3i
      real(gp) :: dum1,dum2,dum3,tt

      call startend_buffers(geocode,nl1,nl2,nl3,nbx,nby,nbz)

      !aliasing
      n1i=ndims(1)
      n2i=ndims(2)
      n3i=ndims(3)

      unt=f_get_free_unit(22)
      call f_open_file(unit=unt,file=trim(filename)//'.cube',status='old')
      read(unt,*)! 'CUBE file for charge density'
      read(unt,*)! 'Case for '//trim(message)

      read(unt,'(i5,3(f12.6),a)')  nat , dum1, dum2, dum3
      read(unt,'(i5,3(f12.6))') n1t , dum3,   dum1 ,dum2
      read(unt,'(i5,3(f12.6))') n2t ,dum1 , dum3  ,  dum2
      read(unt,'(i5,3(f12.6))') n3t ,dum1 , dum2 , dum3

      !grid positions
      n1=n1t/2-nbx
      if (n1i /= 2*n1+(1-nbx)+2*nl1) stop 'n1i not valid'
      n2=n2t/2-nby
      if (n2i /= 2*n2+(1-nby)+2*nl2) stop 'n2i not valid'
      n3=n3t/2-nbz
      if (n3i /= 2*n3+(1-nbz)+2*nl3) stop 'n3i not valid'

      !zero the buffer
      call f_zero(rho)

      do iat=1,nat
         !read(unt,'(i5,4(f12.6))')! idum , dum1 , (rxyz(j,iat),j=1,3)
         read(unt,*)! idum , dum1 , (rxyz(j,iat),j=1,3)
      end do

      !the loop is reverted for a cube file
      !charge normalised to the total charge
      do i1=0,2*(n1+nbx) - 1
         do i2=0,2*(n2+nby) - 1
            icount=0
            do i3=0,2*(n3+nbz) - 1
               icount=icount+1
               if (icount == 6 .or. i3==2*(n3+nbz) - 1) then
                  advancestring='yes'
                  icount=0
               else
                  advancestring='no'
               end if
               ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
               read(unt,'(1x,1pe13.6)',advance=advancestring) tt !rho(ind)
               !rho(ind)=tt
               rho(i1+nl1,i2+nl2,i3+nl3)=tt
               !           write(16,*)i1,i2,i3,ind,rho(ind)
               !read(unt,*)',advance=advancestring) rho(ind)
            end do
         end do
      end do
      call f_close(unt)

      ! write(14,*)rho

    END SUBROUTINE read_cube_field

    !> routine to be used to estimate the dimension of the array to be allocated
    subroutine read_field_dimensions(filename,geocode,ndims,nspin)
      use dictionaries, only: f_err_throw
      use dynamic_memory
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(out) :: nspin
      integer, dimension(3), intent(out) ::  ndims
      !local variables
      integer :: nat,fformat,isuffix
      real(gp), dimension(1) :: rho_fake
      real(gp), dimension(3) :: hgrids
      real(gp), dimension(:,:), pointer   :: rxyz2
      integer, dimension(:), pointer   :: iatypes2, znucl2

      fformat=get_file_format(filename,isuffix)

      select case(fformat)
      case(CUBE)
         call read_cube(filename(1:isuffix),geocode,ndims,hgrids,nspin,1,1,rho_fake,&
              nat,rxyz2, iatypes2, znucl2,dry_run=.true.)
         call f_free_ptr(rxyz2)
         call f_free_ptr(iatypes2)
         call f_free_ptr(znucl2)
      case(ETSF)
         call f_err_throw('Size estimator not (yet) possible for ETSF format')
      end select

    end subroutine read_field_dimensions

    !> Read a density file using file format depending on the extension.
    subroutine read_field(filename,geocode,ndims,hgrids,nspin,ldrho,nrho,rho,&
         nat,rxyz,iatypes,znucl)
      use dynamic_memory
      use dictionaries, only: f_err_throw
      use IOboxETSF, only: read_etsf
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: ldrho,nrho
      integer, intent(out) :: nspin
      integer, dimension(3), intent(out) ::  ndims
      real(gp), dimension(3), intent(out) :: hgrids
      real(dp), dimension(ldrho,nrho), intent(inout) :: rho
      real(gp), dimension(:,:), pointer, optional :: rxyz
      integer, intent(out), optional ::  nat
      integer, dimension(:), pointer, optional :: iatypes, znucl
      !local variables
      character(len = *), parameter :: subname = "read_field"
      logical, dimension(4) :: optargs
      integer :: isuffix,fformat,nat_read
      real(gp), dimension(:,:), pointer :: rxyz_read
      integer, dimension(:), pointer :: iatypes_read, znucl_read

      optargs=[present(rxyz),present(nat),present(iatypes),present(znucl)]
      !check the arguments
      if ((.not. all(optargs)) .and. (any(optargs))) &
           call f_err_throw('Wrong usage of read_field, rxyz, znucl and iatypes should be _all_ present')

      call f_routine(id=subname)

      fformat=get_file_format(filename,isuffix)

      select case(fformat)
      case(CUBE)
         call read_cube(filename(1:isuffix),geocode,ndims,hgrids,nspin,ldrho,nrho,rho,&
              nat_read,rxyz_read, iatypes_read, znucl_read)
      case(ETSF)
         call read_etsf(filename(1:isuffix),geocode,&
              ndims(1),ndims(2),ndims(3),nspin,hgrids(1),hgrids(2),hgrids(3),&
              ldrho,nrho,rho,&
              nat_read,rxyz_read, iatypes_read, znucl_read)
         if (ldrho < product(ndims) .or. nrho < nspin) &
              call f_err_throw('Severe error, the sizes of the rho array have revealed not to be sufficient. '//&
              'The etsf reading  might have caused a boundary error')
      end select

      if (all(optargs)) then
         rxyz => rxyz_read
         iatypes => iatypes_read
         znucl => znucl_read
         nat=nat_read
      else
         call f_free_ptr(rxyz_read)
         call f_free_ptr(iatypes_read)
         call f_free_ptr(znucl_read)
      end if
      call f_release_routine()
    END SUBROUTINE read_field

    !> Read density or potential in cube format
    subroutine read_cube(filename,geocode,ndims,hgrids,nspin,ldrho,nrho,rho,&
         nat,rxyz, iatypes, znucl,dry_run)
      use PSbase
      use f_utils
      use yaml_strings, only: operator(+),f_strcpy
      use dictionaries, only: f_err_throw
      use dynamic_memory
      implicit none
      character(len=*), intent(in) :: filename
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: ldrho,nrho !<dimensions of the rho array
      integer, intent(out) :: nspin
      integer, dimension(3), intent(out) ::  ndims
      real(gp), dimension(3), intent(out) :: hgrids
      real(dp), dimension(ldrho,nrho), intent(inout) :: rho
      real(gp), dimension(:,:), pointer   :: rxyz
      integer, intent(out)   ::  nat
      integer, dimension(:), pointer   :: iatypes, znucl
      logical, intent(in), optional :: dry_run !<only retrieve dimensions, do not fill the array
      !local variables
      !n(c) character(len=*), parameter :: subname='read_cube'
      character(len=5) :: suffix
      integer, dimension(3) :: na,nb
      real(gp), dimension(3) :: ha,hb
      integer :: nat2
      logical :: exists,drr
      real(gp), dimension(:,:), pointer   :: rxyz2
      integer, dimension(:), pointer   :: iatypes2, znucl2

      ! Test if we have up and down densities.
      call f_file_exists(file=trim(filename)//"-up.cube",exists=exists)
      if (exists) then
         call f_file_exists(file=trim(filename)//"-down.cube",exists=exists)
         if (.not.exists) then
            call f_err_throw("found a "+filename+"-up.cube file, but no -down.cube...")
            nspin = 1
         else
            nspin = 2
         end if
      else
         call f_file_exists(file=trim(filename)//".cube",exists=exists)
         if (.not. exists) call f_err_throw('The file '+filename+' does not exists')
         nspin = 1
      end if

      drr=.false.
      if (present(dry_run)) drr=dry_run

      !read the header of the files and verify it is coherent for nspin==2
      if (nspin /= 2) then
         call f_strcpy(src='',dest=suffix)
         call read_cube_header(filename//trim(suffix),geocode,ndims,hgrids,&
              nat,rxyz, iatypes, znucl)
      else
         call f_strcpy(src='-up',dest=suffix)
         call read_cube_header(filename//trim(suffix),geocode,na,ha,&
              nat,rxyz, iatypes, znucl)
         call f_strcpy(src='-down',dest=suffix)
         call read_cube_header(filename//trim(suffix),geocode,nb,hb,&
              nat2,rxyz2, iatypes2, znucl2)
         if (any(na /= nb .or. ha /= hb)) then
            call f_err_throw('Error in reading .cube file, the dimensions between spin up and down are not coherent')
         else
            ndims=na
            hgrids=ha
         end if
         call f_free_ptr(rxyz2)
         call f_free_ptr(iatypes2)
         call f_free_ptr(znucl2)
      end if

      !now the information to allocate the array has been set up
      if (drr) then
         return
      else
         !check if the given dimension for rho is compatible with the array
         if (ldrho < product(ndims) .or. nrho /= nspin) &
              call f_err_throw('The dimension of the rho array is not coherent with the file')
      end if

      !read the header of the files and verify it is coherent for nspin==2
      if (nspin /= 2) then
         call f_strcpy(src='',dest=suffix)
         call read_cube_field(filename//trim(suffix),geocode,ndims,rho)
      else
         call f_strcpy(src='-up',dest=suffix)
         call read_cube_field(filename//trim(suffix),geocode,ndims,rho(1,1))
         call f_strcpy(src='-down',dest=suffix)
         call read_cube_field(filename//trim(suffix),geocode,ndims,rho(1,2))
      end if

    END SUBROUTINE read_cube

    !> Write a (sum of two) field in the ISF basis in the cube format
    subroutine write_cube_fields(prefix,message,geocode,ndims,ns,hgrids,&
         factor,a,x,nexpo,b,y,nat,rxyz,iatype,nzatom,nelpsp,ixyz0)
      use f_utils
      use dynamic_memory
      implicit none
      !integer,intent(in) :: fileunit0,fileunitx,fileunity,fileunitz
      character(len=1), intent(in) :: geocode
      character(len=*), intent(in) :: message,prefix
      integer, intent(in) :: nexpo,nat
      real(gp), intent(in) :: a,b,factor
      integer, dimension(3), intent(in) :: ndims,ns
      real(gp), dimension(3), intent(in) :: hgrids
      !type(atoms_data), intent(in) :: at
      real(dp), dimension(ndims(1),ndims(2),ndims(3)), intent(in) :: x,y
      real(gp), dimension(3,nat), intent(in) :: rxyz
      integer, dimension(nat), intent(in) :: iatype
      integer, dimension(*), intent(in) :: nzatom !< of dimension ntypes
      integer, dimension(*), intent(in) :: nelpsp !< of dimension ntypes
      integer, dimension(3), intent(in) :: ixyz0
      !local variables
      character(len=3) :: advancestring
      integer :: fileunit0,fileunitx,fileunity,fileunitz,n1i,n2i,n3i,n1s,n2s,n3s
      integer :: nl1,nl2,nl3,nbx,nby,nbz,i1,i2,i3,icount,j,iat,nc1,nc2,nc3
      real(dp) :: later_avg,xx,yy,zz
      real(gp), dimension(3) :: cell_dim

      call f_routine(id='write_cube_fields')

      call startend_buffers(geocode,nl1,nl2,nl3,nbx,nby,nbz)
      call cube_dimensions(geocode,ndims,nc1,nc2,nc3)

      n1i=ndims(1)
      n2i=ndims(2)
      n3i=ndims(3)
      n1s=ns(1)
      n2s=ns(2)
      n3s=ns(3)
      fileunit0=f_get_free_unit(22)
      call f_open_file(unit=fileunit0,file=trim(prefix)//'.cube',status='unknown')

      ! A nonstandard .CUBE file where the field is written with the maximum number of
      ! decimal places can be obtained by uncommenting the writes to unit 23
      !open(unit=22,file=trim(filename)//'.cube',status='unknown')
      !  open(unit=23,file=trim(filename)//'.CUBE',status='unknown')
      write(fileunit0,*)'CUBE file for ISF field'
      write(fileunit0,*)'Case for '//trim(message)
      write(fileunit0,'(i5,3(f12.6))') nat,0.0_gp,0.0_gp,0.0_gp
      !  write(23,*)'CUBE file for ISF field'
      !  write(23,*)'Case for '//trim(message)
      !  write(23,'(i5,3(f12.6))') at%astruct%nat,0.0_gp,0.0_gp,0.0_gp
      !grid and grid spacings
      write(fileunit0,'(i5,3(f12.6))') nc1,hgrids(1),0.0_gp,0.0_gp
      write(fileunit0,'(i5,3(f12.6))') nc2,0.0_gp,hgrids(2),0.0_gp
      write(fileunit0,'(i5,3(f12.6))') nc3,0.0_gp,0.0_gp,hgrids(3)
      !  write(23,'(i5,3(f12.6))') nc1,hxh,0.0_gp,0.0_gp
      !  write(23,'(i5,3(f12.6))') nc2,0.0_gp,hyh,0.0_gp
      !  write(23,'(i5,3(f12.6))') nc3,0.0_gp,0.0_gp,hzh
      !atomic number and positions
      do iat=1,nat
         write(fileunit0,'(i5,4(f12.6))') nzatom(iatype(iat)), real(nelpsp(iatype(iat)),gp) &
              ,(rxyz(j,iat),j=1,3)
         !     write(23,'(i5,f12.6,3(f19.12))') at%nzatom(at%astruct%iatype(iat)), at%nelpsp(at%astruct%iatype(iat))*1. &
         !          ,(rxyz(j,iat),j=1,3)
      end do

      !the loop is reverted for a cube file
      !charge normalised to the total charge
      do i1=0,nc1 - 1
         do i2=0,nc2 - 1
            icount=0
            do i3=0,nc3 - 1
               icount=icount+1
               if (icount == 6 .or. i3==nc3 - 1) then
                  advancestring='yes'
                  icount=0
               else
                  advancestring='no'
               end if
               !ind=i1+nl1+(i2+nl2-1)*n1i+(i3+nl3-1)*n1i*n2i
               write(fileunit0,'(1x,1pe13.6)',advance=advancestring)&
                    a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
               !           write(23,'(1x,e24.17)',advance=advancestring)&
               !                a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
            end do
         end do
      end do
      !close(22)
      call f_close(fileunit0)
      !  close(23)

      cell_dim(1)=nc1*hgrids(1)
      cell_dim(2)=nc2*hgrids(2)
      cell_dim(3)=nc3*hgrids(3)


      fileunitx=f_get_free_unit(23)
      fileunity=f_get_free_unit(24)
      fileunitz=f_get_free_unit(25)

      call f_open_file(unit=fileunitx,file=trim(prefix)//'_avg_x.dat',status='unknown')
      call f_open_file(unit=fileunity,file=trim(prefix)//'_avg_y.dat',status='unknown')
      call f_open_file(unit=fileunitz,file=trim(prefix)//'_avg_z.dat',status='unknown')

      !average in x direction
      !open(unit=23,file=trim(filename)//'_avg_x',status='unknown')
      !open(unit=24,file=trim(filename)//'_centre_x',status='unknown')
      !  do i1=0,2*n1+1
      do i1=0,nc1 - 1
         later_avg=0.0_dp
         do i3=0,nc3 -1
            do i2=0,nc2 - 1
               later_avg=later_avg+&
                    a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
            end do
         end do
         later_avg=later_avg/real(nc2*nc3,dp) !2D integration/2D Volume
         !to be checked with periodic/isolated BC
         write(fileunitx,*)i1+n1s,cell_dim(1)/real(factor*nc1,dp)*(i1+n1s),later_avg
         !write(24,*)i1+n1s,at%astruct%cell_dim(1)/real(factor*nc1,dp)*(i1+n1s),&
         !     a*x(i1+nl1,nc2/2+nl2,nc3/2+nl3)**nexpo+b*y(i1+nl1,nc2/2+nl2,nc3/2+nl3)
      end do
      !close(23)
      !close(24)
      !average in y direction
      !open(unit=23,file=trim(filename)//'_avg_y',status='unknown')
      !open(unit=24,file=trim(filename)//'_centre_y',status='unknown')
      do i2=0,nc2 - 1
         later_avg=0.0_dp
         do i3=0,nc3 - 1
            do i1=0,nc1 -1
               later_avg=later_avg+&
                    a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
            end do
         end do
         later_avg=later_avg/real(nc1*nc3,dp) !2D integration/2D Volume
         !to be checked with periodic/isolated BC
         write(fileunity,*)i2+n2s,cell_dim(2)/real(factor*nc2,dp)*(i2+n2s),later_avg
         !write(24,*)i2+n2s,at%astruct%cell_dim(2)/real(factor*nc2,dp)*(i2+n2s),&
         !    a*x(nc1/2+nl1,i2+nl2,nc3/2+nl3)**nexpo+b*y(nc1/2+nl1,i2+nl2,nc3/2+nl3)
      end do
      !close(23)
      !close(24)
      !average in z direction
      !open(unit=23,file=trim(filename)//'_avg_z',status='unknown')
      !open(unit=24,file=trim(filename)//'_centre_z',status='unknown')
      do i3=0,nc3 - 1
         later_avg=0.0_dp
         do i2=0,nc2 - 1
            do i1=0,nc1 -1
               later_avg=later_avg+&
                    a*x(i1+nl1,i2+nl2,i3+nl3)**nexpo+b*y(i1+nl1,i2+nl2,i3+nl3)
            end do
         end do
         later_avg=later_avg/real(nc1*nc2,dp) !2D integration/2D Volume
         !to be checked with periodic/isolated BC
         write(fileunitz,*)i3+n3s,cell_dim(3)/real(factor*nc3,dp)*(i3+n3s),later_avg
         !write(24,*)i3+n3s,at%astruct%cell_dim(3)/real(factor*nc3,dp)*(i3+n3s),&
         !     a*x(nc1/2+nl1,nc2/2+nl2,i3+nl3)**nexpo+b*y(nc1/2+nl1,nc2/2+nl2,i3+nl3)
      end do

      call f_close(fileunitx)
      call f_close(fileunity)
      call f_close(fileunitz)
      !close(23)
      !close(24)

      !plot one single line of the file
      if ( .not. all(ixyz0 == -1)) then
         fileunitx=f_get_free_unit(26)
         fileunity=f_get_free_unit(27)
         fileunitz=f_get_free_unit(28)

         call f_open_file(unit=fileunitx,file=trim(prefix)//'_line_x.dat',status='unknown')
         call f_open_file(unit=fileunity,file=trim(prefix)//'_line_y.dat',status='unknown')
         call f_open_file(unit=fileunitz,file=trim(prefix)//'_line_z.dat',status='unknown')

         do i1=1,n1i
            xx=i1*hgrids(1)
            write(fileunitx,'(1x,I8,4(1x,e22.15))')i1,xx,a*x(i1,ixyz0(2),ixyz0(3))**nexpo+b*y(i1,ixyz0(2),ixyz0(3))
         end do

         do i2=1,n2i
            yy=i2*hgrids(2)
            write(fileunity,'(1x,I8,3(1x,e22.15))')i2,yy,a*x(ixyz0(1),i2,ixyz0(3))**nexpo+b*y(ixyz0(1),i2,ixyz0(3))
         end do

         do i3=1,n3i
            zz=i3*hgrids(3)
            write(fileunitz,'(1x,I8,3(1x,e22.15))')i3,zz,a*x(ixyz0(1),ixyz0(2),i3)**nexpo+b*y(ixyz0(1),ixyz0(2),i3)
         end do

         call f_close(fileunitx)
         call f_close(fileunity)
         call f_close(fileunitz)
      end if

      call f_release_routine()

    END SUBROUTINE write_cube_fields

    subroutine dump_field(filename,geocode,ndims,hgrids,nspin,rho,&
         rxyz,iatype,nzatom,nelpsp,ixyz0)
      use dynamic_memory
      use dictionaries, only: f_err_throw
      use f_utils
      use IOboxETSF, only: write_etsf_density
      use yaml_strings
      implicit none
      !integer,intent(in) :: fileunit0,fileunitx,fileunity,fileunitz
      integer, intent(in) :: nspin
      character(len=1), intent(in) :: geocode
      character(len=*), intent(in) :: filename
      integer, dimension(3), intent(in) :: ndims
      real(gp), dimension(3), intent(in) :: hgrids
      !type(atoms_data), intent(in) :: at
      real(dp), dimension(ndims(1),ndims(2),ndims(3),nspin), intent(in) :: rho
      real(gp), dimension(:,:), intent(in), optional, target :: rxyz
      integer, dimension(:), intent(in), optional, target  :: iatype
      integer, dimension(:), intent(in), optional, target  :: nzatom !< of dimension ntypes
      integer, dimension(:), intent(in), optional, target  :: nelpsp !< of dimension ntypes
      integer, dimension(3), intent(in), optional ::  ixyz0 !< points that have to be plot as lines
      !local variables
      character(len=5) :: suffix
      character(len=65) :: message
      integer :: ia,ib,isuffix,fformat,nat!ierr,n1i,n2i,n3i
      real(gp), dimension(:,:), pointer :: rxyz_
      integer, dimension(:), pointer :: iatype_
      integer, dimension(:), pointer :: nzatom_
      integer, dimension(:), pointer :: nelpsp_
      !local variables
      real(dp) :: a,b
      !real(gp) :: hxh,hyh,hzh
      integer, dimension(3) :: ns,ixyz0_ !< starting points (zero in this case)
!      real(dp), dimension(:,:), pointer :: pot_ion
      integer,parameter :: unit0 = 22
      integer,parameter :: unitx = 23
      integer,parameter :: unity = 24
      integer,parameter :: unitz = 25

      call f_routine(id='dump_field')

      nat=0
      if (present(iatype)) then
         nat=size(iatype) !should also check the arguments
         !associate the pointers
         iatype_=>iatype
         rxyz_=>rxyz
         nzatom_=>nzatom
         nelpsp_=>nelpsp
      else
         iatype_=f_malloc_ptr(1,id='iatype_')
         nzatom_=f_malloc_ptr(1,id='nzatom_')
         nelpsp_=f_malloc_ptr(1,id='nelpsp_')
         rxyz_=f_malloc_ptr([1,1],id='rxyz_')
      end if

      ixyz0_=-1
      if (present(ixyz0)) then
         if (any(ixyz0 < 1) .or. any(ixyz0 > ndims)) then
            call f_err_throw('The values of ixyz0='+yaml_toa(ixyz0)+&
                 ' should be within the size of the box (1 to'+&
                 yaml_toa(ndims)+')') !,&
                   !err_name='BIGDFT_RUNTIME_ERROR')
         end if
         ixyz0_=ixyz0
      end if

      call f_zero(ns)

      fformat=get_file_format(filename,isuffix)



!!$      if (iproc == 0) then


!!$         open(unit=unit0,file=trim(filename(:isuffix))//'.cube',status='unknown')
!!$         open(unit=unitx,file=trim(filename(:isuffix))//'_avg_x',status='unknown')
!!$         open(unit=unity,file=trim(filename(:isuffix))//'_avg_y',status='unknown')
!!$         open(unit=unitz,file=trim(filename(:isuffix))//'_avg_z',status='unknown')

         if (nspin /=2) then
            message='total spin'
            if (fformat == CUBE) then
               suffix=''
               a=1.0_dp
               ia=1
               b=0.0_dp
               ib=1
               call write_cube_fields(filename(:isuffix),message,geocode,&
                    ndims,ns,hgrids,&
                    1.0_dp,a,rho(1,1,1,ia),1,b,rho(1,1,1,ib),nat,rxyz_,iatype_,nzatom_,nelpsp_,ixyz0_)
            else
               call write_etsf_density(filename(:isuffix),message,geocode,&
                    ndims,hgrids,&
                    rho, 1,nat,rxyz_,iatype_,size(nzatom_),nzatom_)
            end if
         else
            if (fformat == CUBE) then
               suffix=''
               message='total spin'
               a=1.0_dp
               ia=1
               b=0.0_dp
               ib=2
               call write_cube_fields(trim(filename(:isuffix))//trim(suffix),message,geocode,&
                    ndims,ns,hgrids,&
                    1.0_dp,a,rho(1,1,1,ia),1,b,rho(1,1,1,ib),nat,rxyz_,iatype_,nzatom_,nelpsp_,ixyz0_)
               suffix='-down'
               message='spin down'
               a=0.0_dp
               ia=1
               b=1.0_dp
               ib=2
               call write_cube_fields(trim(filename(:isuffix))//trim(suffix),message,geocode,&
                    ndims,ns,hgrids,&
                    1.0_dp,a,rho(1,1,1,ia),1,b,rho(1,1,1,ib),nat,rxyz_,iatype_,nzatom_,nelpsp_,ixyz0_)
               suffix='-u-d'
               message='spin difference'
               a=1.0_dp
               ia=1
               b=-2.0_dp
               ib=2
               call write_cube_fields(trim(filename(:isuffix))//trim(suffix),message,geocode,&
                    ndims,ns,hgrids,&
                    1.0_dp,a,rho(1,1,1,ia),1,b,rho(1,1,1,ib),nat,rxyz_,iatype_,nzatom_,nelpsp_,ixyz0_)
               suffix='-up'
               message='spin up'
               a=1.0_dp
               ia=1
               b=-1.0_dp
               ib=2
               call write_cube_fields(trim(filename(:isuffix))//trim(suffix),message,geocode,&
                    ndims,ns,hgrids,&
                    1.0_dp,a,rho(1,1,1,ia),1,b,rho(1,1,1,ib),nat,rxyz_,iatype_,nzatom_,nelpsp_,ixyz0_)
            else
               message = 'spin up, down, total, difference'
               call write_etsf_density(filename(:isuffix),message,geocode,&
                    ndims,hgrids,&
                    rho, 2,nat,rxyz_,iatype_,size(nzatom_),nzatom_)
            end if

         end if

         close(unit=unit0)
         close(unit=unitx)
         close(unit=unity)
         close(unit=unitz)

!!$      end if

         if (present(iatype)) then
            nullify(iatype_)
            nullify(rxyz_)
            nullify(nzatom_)
            nullify(nelpsp_)
         else
            call f_free_ptr(iatype_)
            call f_free_ptr(rxyz_)
            call f_free_ptr(nzatom_)
            call f_free_ptr(nelpsp_)
         end if

         call f_release_routine()

       END SUBROUTINE dump_field

end module IObox
