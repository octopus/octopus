!> @file
!!  Provide the test functions for the different examples of Poisson Solver
!! @author
!!    Copyright (C) 2002-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
subroutine PS_Check_options(parser)
  use yaml_parse
  use dictionaries, only: dict_new,operator(.is.)
  implicit none
  type(yaml_cl_parse), intent(inout) :: parser

  call yaml_cl_parse_option(parser,'ndim','None',&
       'Domain Sizes','n',&
       dict_new('Usage' .is. &
       'Sizes of the simulation domain of the check',&
       'Allowed values' .is. &
       'Yaml list of integers. If a scalar integer is given, all the dimensions will have this size.'),first_option=.true.)

  call yaml_cl_parse_option(parser,'geocode','F',&
       'Boundary conditions','g',&
       dict_new('Usage' .is. &
       'Set the boundary conditions of the run',&
       'Allowed values' .is. &
       'String scalar. "F","S","W","P" boundary conditions are allowed'))

  call yaml_cl_parse_option(parser,'method','None',&
       'Embedding method','m',&
       dict_new('Usage' .is. &
       'Set the embedding method used. A non present value implies vacuum treatment.',&
       'Allowed values' .is. &
       dict_new("PI" .is. 'Polarization iteration Method',&
       "PCG" .is. 'Preconditioned Conjugate Gradient')))

  call yaml_cl_parse_option(parser,'seteps','4',&
       'Epsilon determination method','e',&
       dict_new('Usage' .is. &
       'Set the dielectric constant determination method.',&
       'Allowed values' .is. &
       dict_new('1' .is. 'Analytical epsilon' ,&
       '2' .is. 'analytical electron dependence',&
       '3' .is. 'real electron density from cube file (need electroninc_density.cube)',&
       '4' .is. 'calculate the cavity and dump it on disk',&
       '5' .is. 'Solves GPe with PCG customized (should be identical to 4 + PCG)',&
       '6' .is. 'Modified Poisson Botzmann Equation solver',&
       '7' .is. 'Solves GPe with PCG customized and a restart is implemented',&
       '8' .is. 'Solves GPe with PI customized (should be identical to 4 + PI)',&
       '9' .is. 'Solves GPe with PSD',&
       '10' .is. 'Solves GPe with PCG customized and is coupled with PSD')))

  call yaml_cl_parse_option(parser,'accel','No',&
       'GPU Acceleration','a',&
       dict_new('Usage' .is. &
       'Boolean, set the GPU acceleration'))

  call yaml_cl_parse_option(parser,'logfile','Yes',&
       'Write logfile','l',&
       dict_new('Usage' .is. &
       'Boolean, set the logfile as log.yaml'))

  call yaml_cl_parse_option(parser,'deltacav','None',&
       'Delta rigid cavity','d',&
       dict_new('Usage' .is. &
       'Sizes of the delta for error function',&
       'Allowed values' .is. &
       'Real value'))

  call yaml_cl_parse_option(parser,'input','None',&
       'Inputfile of Poisson solver','i',&
       dict_new('Usage' .is. &
       'Put the dictionary of PS inputfile to preload some parameters',&
       'Allowed values' .is. &
       'yaml Dictionary'))


end subroutine PS_Check_options


!!!!> This subroutine builds some analytic functions that can be used for 
!!!!! testing the poisson solver.
!!!!! The default choice is already well-tuned for comparison.
!!!!! WARNING: not all the test functions can be used for all the boundary conditions of
!!!!! the poisson solver, in order to have a reliable analytic comparison.
!!!!! The parameters of the functions must be adjusted in order to have a sufficiently localized
!!!!! function in the isolated direction and an explicitly periodic function in the periodic ones.
!!!!! Beware of the high-frequency components that may false the results when hgrid is too high.
!!!recursive subroutine test_functions_box(mesh,nspden,a_gauss,&
!!!     density,potential,rhopot,pot_ion,offset)
!!!  use box
!!!  use f_utils
!!!  use f_precisions
!!!  use PSbase
!!!  use f_jmp
!!!  use time_profiling
!!!  use numerics
!!!  use dynamic_memory
!!!  implicit none
!!!  type(cell), intent(in) :: mesh !<definition of the cell
!!!  integer, intent(in) :: nspden
!!!  real(kind=8), intent(in) :: a_gauss
!!!  real(kind=8), intent(out) :: offset
!!!  real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(out) :: pot_ion,potential
!!!  real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3),nspden), intent(out) :: density,rhopot
!!!  !local variables
!!!  integer, parameter :: nrep=10
!!!  logical ::  separable=.true.,old=.false. !save attribute for the profiling
!!!  integer :: i1,i2,i3,ifx,ify,ifz,i,n01,n02,n03
!!!  
!!!  real(kind=8) :: x1,x2,x3,length,denval,a2,derf_tt,factor,r,r2,hx,hy,hz
!!!  real(kind=8) :: fx,fx2,fy,fy2,fz,fz2,a,ax,ay,az,bx,by,bz,tt,acell
!!!  integer(f_long) :: t1,t0
!!!  type(box_iterator) :: bit
!!!  type(f_jmpbuf), save :: jmpbuf
!!!
!!!  !backward compatibility
!!!  acell=mesh%ndims(1)*mesh%hgrids(1)
!!!  hx=mesh%hgrids(1)
!!!  hy=mesh%hgrids(2)
!!!  hz=mesh%hgrids(3)
!!!  n01=mesh%ndims(1)
!!!  n02=mesh%ndims(2)
!!!  n03=mesh%ndims(3)
!!!
!!!  !select the specifications for the loop
!!!  select case (cell_geocode(mesh))
!!!  case('P')
!!!     !parameters for the test functions
!!!     length=acell
!!!     a=0.5d0/a_gauss**2
!!!     !test functions in the three directions
!!!     ifx=5
!!!     ify=5
!!!     ifz=5
!!!     !parameters of the test functions
!!!     ax=length
!!!     ay=length
!!!     az=length
!!!     bx=2.d0!real(nu,kind=8)
!!!     by=2.d0!real(nu,kind=8)
!!!     bz=2.d0
!!!     !other factors
!!!     factor = 1.0_dp
!!!  case('S')
!!!     !parameters for the test functions
!!!     length=acell
!!!     a=0.5d0/a_gauss**2
!!!     !test functions in the three directions
!!!     ifx=5
!!!     ifz=5
!!!     !non-periodic dimension
!!!     ify=6
!!!     !parameters of the test functions
!!!     ax=length
!!!     az=length
!!!     bx=2.d0!real(nu,kind=8)
!!!     bz=2.d0!real(nu,kind=8)
!!!     !non-periodic dimension
!!!     ay=length
!!!     by=a
!!!     factor = oneofourpi
!!!  case('F')
!!!     a2 = a_gauss**2
!!!     !Normalization
!!!     factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
!!!     separable=.false.
!!!  case default
!!!     print *,'geometry code not admitted',cell_geocode(mesh)
!!!     stop
!!!  end select
!!!
!!!
!!!  t0=f_time()
!!!  !     do irep=1,nrep 
!!!  entry test_new()
!!!  if (separable .and. .not. old) then
!!!     denval=0.d0 !value for keeping the density positive
!!!     bit=box_iter(mesh,centered=.true.)
!!!     !here separability is restored
!!!     do while(box_next_z(bit))
!!!        call functions(bit%rxyz(3),az,bz,fz,fz2,ifz)
!!!        do while(box_next_y(bit))
!!!           call functions(bit%rxyz(2),ay,by,fy,fy2,ify)
!!!           do while(box_next_x(bit))
!!!              call functions(bit%rxyz(1),ax,bx,fx,fx2,ifx)
!!!              do i=1,nspden
!!!                 density(bit%i,bit%j,bit%k,i) = factor/real(nspden,kind=8)*(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)
!!!              end do
!!!              potential(bit%i,bit%j,bit%k) = -factor*fourpi*fx*fy*fz
!!!              denval=max(denval,-density(bit%i,bit%j,bit%k,1))
!!!           end do
!!!        end do
!!!     end do
!!!  end if
!!!  call f_profile_end(test_new,jmpbuf)
!!!
!!!  !     end do
!!!  t1=f_time()
!!!!!$  print *,f_humantime(t1-t0)
!!!!!$  print *,'there2',denval
!!!  t0=f_time()
!!!  !     do irep=1,nrep 
!!!  entry test_old()
!!!
!!!!print *,'here',separable,old,nrep,associated(jmpbuf%jmp_buf)
!!!  if (separable .and. old) then
!!!     denval=0.d0 !value for keeping the density positive
!!!     do i3=1,n03
!!!        x3 = hz*real(i3-n03/2-1,kind=8)
!!!        call functions(x3,az,bz,fz,fz2,ifz)
!!!        do i2=1,n02
!!!           x2 = hy*real(i2-n02/2-1,kind=8)
!!!           call functions(x2,ay,by,fy,fy2,ify)
!!!           do i1=1,n01
!!!              x1 = hx*real(i1-n01/2-1,kind=8)
!!!              call functions(x1,ax,bx,fx,fx2,ifx)
!!!              print *,'i1,i2,i3',x1,x2,x3
!!!              do i=1,nspden
!!!                 density(i1,i2,i3,i) = factor/real(nspden,kind=8)*(fx2*fy*fz+fx*fy2*fz+fx*fy*fz2)
!!!              end do
!!!              potential(i1,i2,i3) = -fourpi*factor*fx*fy*fz
!!!              denval=max(denval,-density(i1,i2,i3,1))
!!!           end do
!!!        end do
!!!     end do
!!!  end if
!!!  call f_profile_end(test_old,jmpbuf)
!!!  !     end do
!!!  t1=f_time()
!!!
!!!!!$  print *,f_humantime(t1-t0)
!!!!!$  print *,'there3',denval
!!!
!!!  if (.not. separable .and. .not. old) then
!!!     bit=box_iter(mesh,centered=.true.)
!!!     do while(box_next_point(bit))
!!!        r2=square(mesh,bit%rxyz)
!!!        do i=1,nspden
!!!           density(bit%i,bit%j,bit%k,i) = &
!!!                1.d0/real(nspden,kind=8)*max(factor*safe_exp(-r2/a2),1d-24)
!!!        end do
!!!        r = sqrt(r2)
!!!        !Potential from a gaussian
!!!        if (r == 0.0_dp) then
!!!           potential(bit%i,bit%j,bit%k) = 2.d0/(sqrt(pi)*a_gauss)
!!!        else
!!!           call derf_local(derf_tt,r/a_gauss)
!!!           potential(bit%i,bit%j,bit%k) = derf_tt/r
!!!        end if
!!!     end do
!!!  end if
!!!  if (.not. separable .and. old) then
!!!     !gaussian function
!!!     do i3=1,n03
!!!        x3 = hz*real(i3-n03/2,kind=8)
!!!        do i2=1,n02
!!!           x2 = hy*real(i2-n02/2,kind=8)
!!!           do i1=1,n01
!!!              x1 = hx*real(i1-n01/2,kind=8)
!!!              r2 = x1*x1+x2*x2+x3*x3
!!!              do i=1,nspden
!!!                 density(i1,i2,i3,i) = 1.d0/real(nspden,kind=8)*max(factor*safe_exp(-r2/a2),1d-24)
!!!              end do
!!!              r = sqrt(r2)
!!!              !Potential from a gaussian
!!!              if (r == 0.d0) then
!!!                 potential(i1,i2,i3) = 2.d0/(sqrt(pi)*a_gauss)
!!!              else
!!!                 call derf_local(derf_tt,r/a_gauss)
!!!                 potential(i1,i2,i3) = derf_tt/r
!!!              end if
!!!           end do
!!!        end do
!!!     end do
!!!  end if
!!!
!!!  !here we can profile both cases
!!!  !call f_profile(test_old,'Test separable iterator',&
!!!  !repeat=nrep,jmpbuf=jmpbuf,dump_results=.true.)
!!!  !call f_profile(test_new,'Original loop',&
!!!  !          repeat=nrep,jmpbuf=jmpbuf,dump_results=.true.)
!!!
!!!
!!!  denval=0.d0
!!!     
!!!  ! For ixc/=0 the XC potential is added to the solution, and an analytic comparison is no more
!!!  ! possible. In that case the only possible comparison is between the serial and the parallel case
!!!  ! To ease the comparison between the serial and the parallel case we add a random pot_ion
!!!  ! to the potential.
!!!  call f_memcpy(src=density,dest=rhopot)
!!!
!!!  offset=0.d0
!!!  do i3=1,n03
!!!     do i2=1,n02
!!!        do i1=1,n01
!!!           tt=abs(dsin(real(i1+i2+i3,kind=8)+.7d0))
!!!           pot_ion(i1,i2,i3)=tt
!!!           offset=offset+potential(i1,i2,i3)
!!!           !add the case for offset in the surfaces case 
!!!           !(for periodic case it is absorbed in offset)
!!!           if (cell_geocode(mesh) == 'S' .and. denval /= 0.d0) then
!!!              x2 = hy*real(i2-1,kind=8)-0.5d0*acell+0.5d0*hy
!!!              potential(i1,i2,i3)=potential(i1,i2,i3)&
!!!                   -8.d0*datan(1.d0)*denval*real(nspden,kind=8)*(x2**2+0.25d0*acell**2)
!!!              !this stands for
!!!              !denval*2pi*Lx*Lz/Ly^2(y^2-Ly^2/4), less accurate in hgrid
!!!           end if
!!!        end do
!!!     end do
!!!  end do
!!!  if (denval /= 0.d0) density=rhopot
!!!  offset=offset*hx*hy*hz
!!!
!!!  !print *,'offset',offset
!!!
!!!END SUBROUTINE test_functions_box
!!!
!!$subroutine functions(x,a,b,f,f2,whichone)
!!$  implicit none
!!$  integer, intent(in) :: whichone
!!$  real(kind=8), intent(in) :: x,a,b
!!$  real(kind=8), intent(out) :: f,f2
!!$  !local variables
!!$  real(kind=8) :: r,r2,y,yp,ys,factor,pi,g,h,g1,g2,h1,h2
!!$  real(kind=8) :: length,frequency,nu,sigma,agauss
!!$
!!$  factor=0.d0
!!$  pi = 4.d0*datan(1.d0)
!!$  select case(whichone)
!!$  case(1)
!!$     !constant
!!$     f=1.d0
!!$     f2=0.d0
!!$  case(2)
!!$     !gaussian of sigma s.t. a=1/(2*sigma^2)
!!$     r2=a*x**2
!!$     f=dexp(-r2)
!!$     f2=(-2.d0*a+4.d0*a*r2)*dexp(-r2)
!!$  case(3)
!!$     !gaussian "shrinked" with a=length of the system
!!$     length=a
!!$     r=pi*x/length
!!$     y=dtan(r)
!!$     yp=pi/length*1.d0/(dcos(r))**2
!!$     ys=2.d0*pi/length*y*yp
!!$     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
!!$     f2=factor*dexp(-y**2)
!!$     f=dexp(-y**2)
!!$  case(4)
!!$     !cosine with a=length, b=frequency
!!$     length=a
!!$     frequency=b
!!$     r=frequency*pi*x/length
!!$     f=dcos(r)
!!$     f2=-(frequency*pi/length)**2*dcos(r)
!!$  case(5)
!!$     !exp of a cosine, a=length
!!$     nu=2.d0
!!$     r=pi*nu/a*x
!!$     y=dcos(r)
!!$     yp=dsin(r)
!!$     f=dexp(y)
!!$     factor=(pi*nu/a)**2*(-y+yp**2)
!!$     f2=factor*f
!!$  case(6)
!!$     !gaussian times "shrinked" gaussian, sigma=length/10
!!$     length=a
!!$     r=pi*x/length
!!$     y=dtan(r)
!!$     yp=pi/length*1.d0/(dcos(r))**2
!!$     ys=2.d0*pi/length*y*yp
!!$     factor=-2.d0*ys*y-2.d0*yp**2+4.d0*yp**2*y**2
!!$     g=dexp(-y**2)
!!$     g1=-2.d0*y*yp*g
!!$     g2=factor*dexp(-y**2)
!!$
!!$     sigma=length/10
!!$     agauss=0.5d0/sigma**2
!!$     r2=agauss*x**2
!!$     h=dexp(-r2)
!!$     h1=-2.d0*agauss*x*h
!!$     h2=(-2.d0*agauss+4.d0*agauss*r2)*dexp(-r2)
!!$     f=g*h
!!$     f2=g2*h+g*h2+2.d0*g1*h1
!!$  case(7)
!!$     !sine with a=length, b=frequency
!!$     length=a
!!$     frequency=b
!!$     r=frequency*pi*x/length
!!$     f=dsin(r)
!!$     f2=-(frequency*pi/length)**2*dsin(r)
!!$  end select
!!$
!!$END SUBROUTINE functions

subroutine fill_functions_arrays(separable,mesh,funcs,factor,density,potential)
  use PSbase, only: dp
  use numerics, only: fourpi
  use box
  use f_functions
  implicit none
  logical, intent(in) :: separable
  real(dp), intent(in) :: factor
  type(cell), intent(in) :: mesh !<definition of the cell
  type(f_function), dimension(3), intent(in) :: funcs
  real(dp), dimension(mesh%ndim), intent(out) :: density,potential
  !local variables
  integer, parameter :: DENSITY_=1,POTENTIAL_=2 !to improve readability
  type(box_iterator) :: bit

  bit=box_iter(mesh,centered=.true.)
  if (separable) then
     call separable_3d_function(bit,funcs,-factor*fourpi,potential)
     call separable_3d_laplacian(bit,funcs,factor,density)
  else
     call radial_3d_function(bit,funcs(DENSITY_),factor,density)
     call radial_3d_function(bit,funcs(POTENTIAL_),1.0_dp,potential)
  end if
end subroutine fill_functions_arrays
  
subroutine test_functions_new(mesh,nspden,a_gauss,&
     density,potential,rhopot,pot_ion,offset)
  use box
  use f_utils
  use f_precisions
  use PSbase
  use numerics
  use dynamic_memory
  use f_functions
  implicit none
  type(cell), intent(in) :: mesh !<definition of the cell
  integer, intent(in) :: nspden
  real(kind=8), intent(in) :: a_gauss
  real(kind=8), intent(out) :: offset
  real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(out) :: pot_ion,potential
  real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3),nspden), intent(out) :: density,rhopot
  !local variables
  integer, parameter :: DENSITY_=1,POTENTIAL_=2
  logical, dimension(3) :: pers
  integer :: i
  real(dp) :: factor,a2
  type(f_function), dimension(3) :: funcs

  pot_ion=0.d0
  pers=cell_periodic_dims(mesh)
  do i=1,3
     if (pers(i)) then
        funcs(i)=f_function_new(f_exp_cosine,&
             length=mesh%ndims(1)*mesh%hgrids(1),frequency=2.0_dp)
     else
        funcs(i)=f_function_new(f_shrink_gaussian,&
             length=mesh%ndims(1)*mesh%hgrids(1))
     end if
  end do

  !select the specifications for the loop
  select case (cell_geocode(mesh))
  case('P')
     !parameters for the test functions
     factor =1.0_dp
  case('S')
     factor =oneofourpi
  case('F')
     a2 = a_gauss**2
     !Normalization
     factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
     funcs(DENSITY_)=f_function_new(f_gaussian,exponent=1.0_dp/a2)
     funcs(POTENTIAL_)=f_function_new(f_erf,scale=a_gauss/sqrt(2.0_dp))
     
     !test call
     call radial_3d_function_mp(mesh,0.5_dp*a2,density)
  end select

  !laplacian potential = -4pi density
  call fill_functions_arrays(cell_geocode(mesh) /= 'F',mesh,funcs,factor,&
       density,potential)
  
  !treatment for rhopot and pot_ion
  call f_memcpy(src=density,dest=rhopot)

  offset=sum(potential)*mesh%volume_element

END SUBROUTINE test_functions_new

!> This subroutine builds some analytic functions that can be used for 
!! testing the poisson solver.
!! The default choice is already well-tuned for comparison.
!! WARNING: not all the test functions can be used for all the boundary conditions of
!! the poisson solver, in order to have a reliable analytic comparison.
!! The parameters of the functions must be adjusted in order to have a sufficiently localized
!! function in the isolated direction and an explicitly periodic function in the periodic ones.
!! Beware of the high-frequency components that may falsify the results when hgrid is too high.
subroutine test_functions_new2(mesh,acell,a_gauss,mu0,density,potential)
  use yaml_output
  use f_utils
  use f_precisions, only: dp=> f_double
  use box
  use f_functions
  use numerics, only: pi,fourpi
  use f_blas, only: f_axpy
  implicit none
  real(dp), intent(in) :: acell,a_gauss,mu0
  type(cell), intent(in) :: mesh
  real(dp), dimension(mesh%ndim), intent(out) :: density,potential
  !local variables
  integer, parameter :: DENSITY_=1,POTENTIAL_=2
  logical :: separable
  !integer :: i
  real(dp) :: a2,factor
  type(f_function), dimension(3) :: funcs 

  !non-orthorhombic lattice compatibility, it will replace the test_functioncs_new routine

  separable=.false.
  select case(trim(cell_geocode(mesh)))
  case('P')
     !parameters for the test functions
!!$     do i=1,3
!!$        funcs(i)=f_function_new(f_shrink_gaussian,&
!!$             length=acell)
!!$        !funcs(i)=f_function_new(f_exp_cosine,&
!!$        !     length=acell,frequency=2.0_dp)
!!$     end do
     funcs(1)=f_function_new(f_exp_cosine,&
          length=acell,frequency=2.0_dp)
     funcs(2)=f_function_new(f_shrink_gaussian,&
          length=acell)
     funcs(3)=funcs(1)
     factor=1.0_dp
     separable=.true.
  case('S')
     funcs(1)=f_function_new(f_exp_cosine,&
          length=acell,frequency=2.0_dp)
     funcs(2)=f_function_new(f_shrink_gaussian,&
          length=acell)
     funcs(3)=funcs(1)
     factor=1.0_dp
     separable=.true.
  case('F')
     a2 = a_gauss**2
     !Normalization
     factor = 1.d0/(a_gauss*a2*pi*sqrt(pi))
     funcs(DENSITY_)=f_function_new(f_gaussian,exponent=1.0_dp/a2)
     funcs(POTENTIAL_)=f_function_new(f_erf,scale=a_gauss/sqrt(2.0_dp))
     separable=.false.
  case('W')
     funcs(1)=f_function_new(f_shrink_gaussian,&
          length=acell)
     funcs(2)=funcs(1)
     funcs(3)=f_function_new(f_exp_cosine,&
          length=acell,frequency=2.0_dp)
     separable=.true.
     factor=1.0_dp/fourpi
     !the different modes for the wires-like bc are not used here
  end select

  call fill_functions_arrays(separable,mesh,funcs,factor,&
       density,potential)

  !add the screened term to the density if present
  call f_axpy(a=-mu0**2/fourpi*factor,x=potential,y=density)
  
end subroutine test_functions_new2

subroutine radial_3d_function_mp(mesh,rloc,f)
  use f_functions
  use PSbase, only: dp
  use box
  use multipole_preserving
  use dynamic_memory
  implicit none
  real(dp), intent(in) :: rloc
  type(cell), intent(in) :: mesh
  real(dp), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(out) :: f
  !local variables
  integer, parameter :: mp_isf=16
  type(box_iterator) :: bit
  real(dp), dimension(:), allocatable :: mpx,mpy,mpz

  mpx=f_malloc(mesh%ndims(1),id='mpx')
  mpy=f_malloc(mesh%ndims(2),id='mpy')
  mpz=f_malloc(mesh%ndims(3),id='mpz')

  !test of the multipole preserving routine
  !initialize the work arrays needed to integrate with isf
  !names of the routines to be redefined
  call initialize_real_space_conversion(isf_m=mp_isf)

  bit=box_iter(mesh,centered=.true.)
  call gaussian_density(bit%oxyz,rloc,1,.true.,.false.,bit,&
       mp_isf,mesh%ndims(1)-1,mesh%ndims(2)-1,mesh%ndims(3)-1,&
       mpx, mpy, mpz,int(mesh%ndim),f)

!!$  do while(box_next_point(bit))
!!$     r2=square(mesh,bit%rxyz)
!!$     r = sqrt(r2)
!!$     f(bit%i,bit%j,bit%k) =factor*eval(func,r)
!!$  end do

  call f_free(mpx,mpy,mpz)

  call finalize_real_space_conversion()

end subroutine radial_3d_function_mp
