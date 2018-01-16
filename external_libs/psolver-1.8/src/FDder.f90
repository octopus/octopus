  !> @file
  !!   Module containing finite difference derivative operations
  !!
  !! @author
  !!    G. Fisicaro, L. Genovese (September 2015)
  !!    Copyright (C) 2002-2016 BigDFT group
  !!    This file is distributed under the terms of the
  !!    GNU General Public License, see ~/COPYING file
  !!    or http://www.gnu.org/copyleft/gpl.txt .
  !!    For the list of contributors, see ~/AUTHORS
module FDder
  use PSbase
  implicit none
  private

  public :: nabla_u_square,update_rhopol,Delta_GPe_operator
  public :: nabla_u,div_u_i,nabla_u_and_square,nonvacuum_projection

  contains
    !>This routine computes 'nord' order accurate first derivatives
    !! on a equally spaced grid with coefficients from 'Mathematica' program.
    !!
    !! input:
    !! ngrid       = number of points in the grid,
    !! u(ngrid)    = function values at the grid points
    !!
    !! output:
    !! du(ngrid)   = first derivative values at the grid points
    subroutine nabla_u_and_square(geocode,n01,n02,n03,u,du,du2,nord,hgrids)
      use dictionaries, only: f_err_throw
      implicit none


      !c..declare the pass
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: n01,n02,n03,nord
      real(kind=8), dimension(3), intent(in) :: hgrids
      real(kind=8), dimension(n01,n02,n03), intent(in) :: u
      real(kind=8), dimension(n01,n02,n03,3), intent(out) :: du !<nabla of u
      real(kind=8), dimension(n01,n02,n03) :: du2 !<square module of nabla u

      !c..local variables
      integer :: n,m,n_cell
      integer :: i,j,i1,i2,i3,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
      real(kind=8) :: hx,hy,hz,d
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      hx = hgrids(1)!acell/real(n01,kind=8)
      hy = hgrids(2)!acell/real(n02,kind=8)
      hz = hgrids(3)!acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
         call f_err_throw('Ngrid in has to be setted > than n=nord + 1')
         !stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
         !O.K.
      case default
         write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
         stop
      end select

      do i=-m,m
         do j=-m,m
            c1D(i,j)=0.d0
         end do
      end do

      include 'FiniteDiffCorff.inc'
      c1D_1 = c1D/hx
      c1D_2 = c1D/hy
      c1D_3 = c1D/hz
!!!default(shared) private(i1,i2,i3,j,ii, d)
      !$omp parallel do default(none) &
      !$omp private(i3,i2,i1,d,ii,j) &
      !$omp shared(du2,perx,m,n01,n02,n03,du,c1D_1,u)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d = 0.0d0
               du2(i1,i2,i3) = 0.0d0

               if (i1.le.m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                     end do
                  end if
               else if (i1.gt.n01-m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                  end do
               end if
               du(i1,i2,i3,1)= d
               du2(i1,i2,i3) = d*d

            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d = 0.0d0

               if (i2.le.m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                     end do
                  end if
               else if (i2.gt.n02-m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                  end do
               end if
               du(i1,i2,i3,2)= d
               du2(i1,i2,i3) = du2(i1,i2,i3) + d*d
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0

               if (i3.le.m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                     end do
                  end if
               else if (i3.gt.n03-m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                  end do
               end if
               du(i1,i2,i3,3)=d
               du2(i1,i2,i3) = du2(i1,i2,i3) + d*d
               !du2(i1,i2,i3) = square(mesh,du(i1,i2,i3,1:3))

               !der(1:3)=0.d0
               !do i=1,3
               ! do j=1,3
               !  der(i) =  der(i) + mesh%gu(i,j)*du(i1,i2,i3,j)
               ! end do
               !end do
               !du(i1,i2,i3,1:3) = der(1:3)


            end do
         end do
      end do
      !$omp end parallel do
    end subroutine nabla_u_and_square

    !>only calculate the square of the gradient (in orthorhombic cells)
    subroutine nabla_u_square(geocode,n01,n02,n03,u,du2,nord,hgrids)
      use dictionaries, only: f_err_throw
      implicit none

      !c..declare the pass
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: n01,n02,n03,nord
      real(kind=8), dimension(3), intent(in) :: hgrids
      real(kind=8), dimension(n01,n02,n03), intent(in) :: u
      real(kind=8), dimension(n01,n02,n03), intent(out) :: du2 !<square module of nabla u

      !c..local variables
      integer :: n,m,n_cell
      integer :: i,j,i1,i2,i3,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
      real(kind=8) :: hx,hy,hz,d
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      hx = hgrids(1)!acell/real(n01,kind=8)
      hy = hgrids(2)!acell/real(n02,kind=8)
      hz = hgrids(3)!acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
         call f_err_throw('Ngrid in has to be set >  n=nord + 1')
         !stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
         !O.K.
      case default
         write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
         stop
      end select

      do i=-m,m
         do j=-m,m
            c1D(i,j)=0.d0
         end do
      end do

      include 'FiniteDiffCorff.inc'
      c1D_1 = c1D/hx
      c1D_2 = c1D/hy
      c1D_3 = c1D/hz
!!!default(shared) private(i1,i2,i3,j,ii, d)
      !$omp parallel do default(none) &
      !$omp private(i3,i2,i1,d,ii,j) &
      !$omp shared(du2,perx,m,n01,n02,n03,c1D_1,u)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d = 0.0d0
               !du2(i1,i2,i3) = 0.0d0

               if (i1.le.m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                     end do
                  end if
               else if (i1.gt.n01-m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                  end do
               end if
               !du(i1,i2,i3,1)= d
               du2(i1,i2,i3) = d*d

            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d = 0.0d0

               if (i2.le.m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                     end do
                  end if
               else if (i2.gt.n02-m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                  end do
               end if
               !du(i1,i2,i3,2)= d
               du2(i1,i2,i3) = du2(i1,i2,i3) + d*d
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii,d)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0

               if (i3.le.m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                     end do
                  end if
               else if (i3.gt.n03-m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                  end do
               end if
               !du(i1,i2,i3,3)=d
               du2(i1,i2,i3) = du2(i1,i2,i3) + d*d

            end do
         end do
      end do
      !$omp end parallel do
    end subroutine nabla_u_square

    !>this routine computes 'nord' order accurate first derivatives
    !!on a equally spaced grid with coefficients from 'Matematica' program.
    !!
    !!input:
    !!ngrid       = number of points in the grid,
    !!u(ngrid)    = function values at the grid points
    !!
    !!output:
    !!du(ngrid)   = first derivative values at the grid points
    !!
    !!declare the pass
    subroutine div_u_i(geocode,n01,n02,n03,u,du,nord,hgrids,cc)
      implicit none

      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: n01,n02,n03,nord
      real(kind=8), dimension(3), intent(in) :: hgrids
      real(kind=8), dimension(n01,n02,n03,3), intent(in) :: u !<vector field u_i
      real(kind=8), dimension(n01,n02,n03), intent(out) :: du !<divergence d_i u_i
      real(kind=8), dimension(n01,n02,n03), intent(out), optional :: cc !< u_i u_j d_i u_j , needed for the surface term where u_i=d_i rho

      !c..local variables
      integer :: n,m,n_cell
      integer :: i,j,i1,i2,i3,ii
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
      real(kind=8) :: hx,hy,hz,d1,uxy,uyz,uxz
      real(kind=8), parameter :: zero = 0.d0! 1.0d-11
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      hx = hgrids(1)!acell/real(n01,kind=8)
      hy = hgrids(2)!acell/real(n02,kind=8)
      hz = hgrids(3)!acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')

      ! Beware that n_cell has to be > than n.
      if (n_cell < n) then
         write(*,*)'ngrid in has to be set > than n=nord + 1'
         stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
         !O.K.
      case default
         write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
         stop
      end select

      do i=-m,m
         do j=-m,m
            c1D(i,j)=0.d0
         end do
      end do

      include 'FiniteDiffCorff.inc'

      c1D_1 = c1D/hx
      c1D_2 = c1D/hy
      c1D_3 = c1D/hz

!!$  if (present(cc)) then
!!$     cc(i1,i2,i3) = (u(i1,i2,i3,1)**2)*d1+(u(i1,i2,i3,2)**2)*d2+(u(i1,i2,i3,3)**2)*d3+&
!!$          2.d0*u(i1,i2,i3,1)*u(i1,i2,i3,2)*uxy+2.d0*u(i1,i2,i3,2)*u(i1,i2,i3,3)*uyz+&
!!$          2.d0*u(i1,i2,i3,1)*u(i1,i2,i3,3)*uxz
!!$  end if


      if (present(cc)) then
         !$omp parallel do default(none) &
         !$omp private(i1,i2,i3,j,ii, d1,uxz,uxy) &
         !$omp shared(du,c1D_1,u,perx,m,hx,n01,n02,n03,cc)
         do i3=1,n03
            do i2=1,n02
               do i1=1,n01

                  du(i1,i2,i3) = 0.0d0

                  d1 = 0.d0
                  uxy = 0.d0
                  uxz = 0.d0
                  if (i1.le.m) then
                     if (perx) then
                        do j=-m,m
                           ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                           d1 = d1 + c1D_1(j,0)*u(ii,i2,i3,1)
                           uxy = uxy + c1D_1(j,0)*u(ii,i2,i3,2)!/hx
                           uxz = uxz + c1D_1(j,0)*u(ii,i2,i3,3)!/hx
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3,1)
                           uxy = uxy + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3,2)!/hx
                           uxz = uxz + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3,3)!/hx
                        end do
                     end if
                  else if (i1.gt.n01-m) then
                     if (perx) then
                        do j=-m,m
                           ii=modulo(i1 + j - 1, n01 ) + 1
                           d1 = d1 + c1D_1(j,0)*u(ii,i2,i3,1)
                           uxy = uxy + c1D_1(j,0)*u(ii,i2,i3,2)!/hx
                           uxz = uxz + c1D_1(j,0)*u(ii,i2,i3,3)!/hx
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3,1)
                           uxy = uxy + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3,2)!/hx
                           uxz = uxz + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3,3)!/hx
                        end do
                     end if
                  else
                     do j=-m,m
                        d1 = d1 + c1D_1(j,0)*u(i1 + j,i2,i3,1)
                        uxy = uxy + c1D_1(j,0)*u(i1 + j,i2,i3,2)!/hx
                        uxz = uxz + c1D_1(j,0)*u(i1 + j,i2,i3,3)!/hx
                     end do
                  end if
                  !uxy=uxy/hx
                  !uxz=uxz/hx

                  du(i1,i2,i3) =d1
                  cc(i1,i2,i3) = (u(i1,i2,i3,1)**2)*d1 + &
                       2.d0*u(i1,i2,i3,1)*u(i1,i2,i3,2)*uxy + &
                       2.d0*u(i1,i2,i3,1)*u(i1,i2,i3,3)*uxz

               end do
            end do
         end do
         !$omp end parallel do

         !shared)
         !$omp parallel do default(none) &
         !$omp private(i1,i2,i3,j,ii,d1,uyz) &
         !$omp shared(n01,n02,n03,pery,m,c1D_2,u,du,cc)
         do i3=1,n03
            do i2=1,n02
               do i1=1,n01
                  d1=0.d0
                  uyz = 0.d0

                  if (i2.le.m) then
                     if (pery) then
                        do j=-m,m
                           ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                           d1 = d1 + c1D_2(j,0)*u(i1,ii,i3,2)
                           uyz = uyz + c1D_2(j,0)*u(i1,ii,i3,3)!/hy
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3,2)
                           uyz = uyz + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3,3)!/hy
                        end do
                     end if
                  else if (i2.gt.n02-m) then
                     if (pery) then
                        do j=-m,m
                           ii=modulo(i2 + j - 1, n02 ) + 1
                           d1 = d1 + c1D_2(j,0)*u(i1,ii,i3,2)
                           uyz = uyz + c1D_2(j,0)*u(i1,ii,i3,3)!/hy
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3,2)
                           uyz = uyz + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3,3)!/hy
                        end do
                     end if
                  else
                     do j=-m,m
                        d1 = d1 + c1D_2(j,0)*u(i1,i2 + j,i3,2)
                        uyz = uyz + c1D_2(j,0)*u(i1,i2 + j,i3,3)!/hy
                     end do
                  end if
                  du(i1,i2,i3) = du(i1,i2,i3) + d1
                  !uyz=uyz/hy
                  cc(i1,i2,i3) = cc(i1,i2,i3) + (u(i1,i2,i3,2)**2)*d1+ &
                       2.d0*u(i1,i2,i3,2)*u(i1,i2,i3,3)*uyz
               end do
            end do
         end do
         !$omp end parallel do


         !(shared) private(i1,i2,i3,j,ii, d1)
         !$omp parallel do default(none) &
         !$omp private(i1,i2,i3,j,ii,d1) &
         !$omp shared(n01,n02,n03,perz,m,c1D_3,u,du,cc)
         do i3=1,n03
            do i2=1,n02
               do i1=1,n01
                  d1=0.d0
                  if (i3.le.m) then
                     if (perz) then
                        do j=-m,m
                           ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                           d1 = d1 + c1D_3(j,0)*u(i1,i2,ii,3)
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1,3)
                        end do
                     end if
                  else if (i3.gt.n03-m) then
                     if (perz) then
                        do j=-m,m
                           ii=modulo(i3 + j - 1, n03 ) + 1
                           d1 = d1 + c1D_3(j,0)*u(i1,i2,ii,3)
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m,3)
                        end do
                     end if
                  else
                     do j=-m,m
                        d1 = d1 + c1D_3(j,0)*u(i1,i2,i3 + j,3)
                     end do
                  end if
                  du(i1,i2,i3) = du(i1,i2,i3)+d1
                  cc(i1,i2,i3) = cc(i1,i2,i3) + (u(i1,i2,i3,3)**2)*d1
               end do
            end do
         end do
         !$omp end parallel do
      else

         !$omp parallel do default(none) &
         !$omp private(i1,i2,i3,j,ii, d1) &
         !$omp shared(du,c1D_1,u,perx,m,hx,n01,n02,n03)
         do i3=1,n03
            do i2=1,n02
               do i1=1,n01

                  du(i1,i2,i3) = 0.0d0

                  d1 = 0.d0
                  if (i1.le.m) then
                     if (perx) then
                        do j=-m,m
                           ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                           d1 = d1 + c1D_1(j,0)*u(ii,i2,i3,1)
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3,1)
                        end do
                     end if
                  else if (i1.gt.n01-m) then
                     if (perx) then
                        do j=-m,m
                           ii=modulo(i1 + j - 1, n01 ) + 1
                           d1 = d1 + c1D_1(j,0)*u(ii,i2,i3,1)
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3,1)
                        end do
                     end if
                  else
                     do j=-m,m
                        d1 = d1 + c1D_1(j,0)*u(i1 + j,i2,i3,1)
                     end do
                  end if
                  du(i1,i2,i3) =d1
               end do
            end do
         end do
         !$omp end parallel do

         !default(shared) private(i1,i2,i3,j,ii,d1,uyz)
         !$omp parallel do default(none) &
         !$omp private(i1,i2,i3,j,ii,d1) &
         !$omp shared(c1D_2,u,n01,n02,n03,du,m,pery)
         do i3=1,n03
            do i2=1,n02
               do i1=1,n01
                  d1=0.d0
                  if (i2.le.m) then
                     if (pery) then
                        do j=-m,m
                           ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                           d1 = d1 + c1D_2(j,0)*u(i1,ii,i3,2)
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3,2)
                        end do
                     end if
                  else if (i2.gt.n02-m) then
                     if (pery) then
                        do j=-m,m
                           ii=modulo(i2 + j - 1, n02 ) + 1
                           d1 = d1 + c1D_2(j,0)*u(i1,ii,i3,2)
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3,2)
                        end do
                     end if
                  else
                     do j=-m,m
                        d1 = d1 + c1D_2(j,0)*u(i1,i2 + j,i3,2)
                     end do
                  end if
                  du(i1,i2,i3) = du(i1,i2,i3) + d1
               end do
            end do
         end do
         !$omp end parallel do

         !$omp parallel do default(none) &
         !$omp private(i1,i2,i3,j,ii,d1) &
         !$omp shared(c1D_3,u,n01,n02,n03,du,m,perz)
         do i3=1,n03
            do i2=1,n02
               do i1=1,n01
                  d1=0.d0
                  if (i3.le.m) then
                     if (perz) then
                        do j=-m,m
                           ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                           d1 = d1 + c1D_3(j,0)*u(i1,i2,ii,3)
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1,3)
                        end do
                     end if
                  else if (i3.gt.n03-m) then
                     if (perz) then
                        do j=-m,m
                           ii=modulo(i3 + j - 1, n03 ) + 1
                           d1 = d1 + c1D_3(j,0)*u(i1,i2,ii,3)
                        end do
                     else
                        do j=-m,m
                           d1 = d1 + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m,3)
                        end do
                     end if
                  else
                     do j=-m,m
                        d1 = d1 + c1D_3(j,0)*u(i1,i2,i3 + j,3)
                     end do
                  end if
                  du(i1,i2,i3) = du(i1,i2,i3)+d1
               end do
            end do
         end do
         !$omp end parallel do
      end if

    end subroutine div_u_i

    subroutine nabla_u_epsilon(geocode,n01,n02,n03,u,du,nord,hgrids,eps)
      implicit none

      !c..this routine computes 'nord' order accurate first derivatives
      !c..on a equally spaced grid with coefficients from 'Matematica' program.

      !c..input:
      !c..ngrid       = number of points in the grid,
      !c..u(ngrid)    = function values at the grid points

      !c..output:
      !c..du(ngrid)   = first derivative values at the grid points

      !c..declare the pass
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: n01,n02,n03,nord
      real(kind=8), dimension(3), intent(in) :: hgrids
      real(kind=8), dimension(n01,n02,n03), intent(in) :: u !<scalar field
      real(kind=8), dimension(n01,n02,n03,3), intent(out) :: du !< nabla d_i u
      real(gp), dimension(n01,n02,n03), intent(in) :: eps !< array of epsilon, multiply the result with that

      !c..local variables
      integer :: n,m,n_cell,ii
      integer :: i,j,i1,i2,i3
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
      real(kind=8) :: hx,hy,hz, d,e
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      hx = hgrids(1)!acell/real(n01,kind=8)
      hy = hgrids(2)!acell/real(n02,kind=8)
      hz = hgrids(3)!acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')


      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
         write(*,*)'ngrid in has to be setted > than n=nord + 1'
         stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
         !O.K.
      case default
         write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
         stop
      end select

      do i=-m,m
         do j=-m,m
            c1D(i,j)=0.d0
         end do
      end do

      include 'FiniteDiffCorff.inc'

      c1D_1 = c1D/hx
      c1D_2 = c1D/hy
      c1D_3 = c1D/hz
      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d,e)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d= 0.0d0
               e=eps(i1,i2,i3)


               if (i1.le.m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                     end do
                  end if
               else if (i1.gt.n01-m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                  end do
               end if
               du(i1,i2,i3,1) = d*e
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d,e)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0
               e=eps(i1,i2,i3)
               if (i2.le.m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                     end do
                  end if
               else if (i2.gt.n02-m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                  end do
               end if
               du(i1,i2,i3,2)=d*e
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d,e)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0
               e=eps(i1,i2,i3)
               if (i3.le.m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                     end do
                  end if
               else if (i3.gt.n03-m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                  end do
               end if
               du(i1,i2,i3,3)=d*e
            end do
         end do
      end do
      !$omp end parallel do

    end subroutine nabla_u_epsilon

    subroutine nabla_u(geocode,n01,n02,n03,u,du,nord,hgrids)
      implicit none

      !c..this routine computes 'nord' order accurate first derivatives
      !c..on a equally spaced grid with coefficients from 'Matematica' program.

      !c..input:
      !c..ngrid       = number of points in the grid,
      !c..u(ngrid)    = function values at the grid points

      !c..output:
      !c..du(ngrid)   = first derivative values at the grid points

      !c..declare the pass
      character(len=1), intent(in) :: geocode !< @copydoc poisson_solver::doc::geocode
      integer, intent(in) :: n01,n02,n03,nord
      real(kind=8), dimension(3), intent(in) :: hgrids
      real(kind=8), dimension(n01,n02,n03), intent(in) :: u !<scalar field
      real(kind=8), dimension(n01,n02,n03,3), intent(out) :: du !< nabla d_i u

      !c..local variables
      integer :: n,m,n_cell,ii
      integer :: i,j,i1,i2,i3
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D, c1D_1, c1D_2, c1D_3
      real(kind=8) :: hx,hy,hz, d
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      hx = hgrids(1)!acell/real(n01,kind=8)
      hy = hgrids(2)!acell/real(n02,kind=8)
      hz = hgrids(3)!acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      perx=(geocode /= 'F')
      pery=(geocode == 'P')
      perz=(geocode /= 'F')


      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
         write(*,*)'ngrid in has to be setted > than n=nord + 1'
         stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16

      select case(nord)
      case(2,4,6,8,16)
         !O.K.
      case default
         write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
         stop
      end select

      do i=-m,m
         do j=-m,m
            c1D(i,j)=0.d0
         end do
      end do

      include 'FiniteDiffCorff.inc'

      c1D_1 = c1D/hx
      c1D_2 = c1D/hy
      c1D_3 = c1D/hz
      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01

               d= 0.0d0

               if (i1.le.m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-m-1)*u(j+m+1,i2,i3)
                     end do
                  end if
               else if (i1.gt.n01-m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j - 1, n01 ) + 1
                        d = d + c1D_1(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_1(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_1(j,0)*u(i1 + j,i2,i3)
                  end do
               end if
               du(i1,i2,i3,1) = d
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0

               if (i2.le.m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-m-1)*u(i1,j+m+1,i3)
                     end do
                  end if
               else if (i2.gt.n02-m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j - 1, n02 ) + 1
                        d = d + c1D_2(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_2(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_2(j,0)*u(i1,i2 + j,i3)
                  end do
               end if
               du(i1,i2,i3,2)=d
            end do
         end do
      end do
      !$omp end parallel do

      !$omp parallel do default(shared) private(i1,i2,i3,j,ii, d)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               d = 0.0d0
               if (i3.le.m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-m-1)*u(i1,i2,j+m+1)
                     end do
                  end if
               else if (i3.gt.n03-m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j - 1, n03 ) + 1
                        d = d + c1D_3(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        d = d + c1D_3(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                     end do
                  end if
               else
                  do j=-m,m
                     d = d + c1D_3(j,0)*u(i1,i2,i3 + j)
                  end do
               end if
               du(i1,i2,i3,3)=d
            end do
         end do
      end do
      !$omp end parallel do
    end subroutine nabla_u

    !> Like fssnord3DmatNabla but corrected such that the index goes at the beginning
    !! Multiplies also times (nabla epsilon)/(4pi*epsilon)= nabla (log(epsilon))/(4*pi)
    subroutine update_rhopol(mesh,u,nord,eta,eps,dlogeps,rhopol,rhores2)
      use numerics, only: oneofourpi
      use box
      implicit none

      !c..this routine computes 'nord' order accurate first derivatives
      !c..on a equally spaced grid with coefficients from 'Matematica' program.

      !c..input:
      !c..ngrid       = number of points in the grid,
      !c..u(ngrid)    = function values at the grid points

      !c..output:
      !c..du(ngrid)   = first derivative values at the grid points

      !c..declare the pass
      integer, intent(in) :: nord
      real(kind=8), intent(in) :: eta
      type(cell), intent(in) :: mesh
      real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(in) :: u
      real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(in) :: eps
      real(kind=8), dimension(3,mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(in) :: dlogeps
      real(kind=8), dimension(mesh%ndims(1),mesh%ndims(2),mesh%ndims(3)), intent(inout) :: rhopol
      real(kind=8), intent(out) :: rhores2

      !c..local variables
      integer :: n,m,n_cell,n01,n02,n03
      integer :: i,j,i1,i2,i3,ii
      !real(kind=8), parameter :: oneo4pi=0.25d0/pi_param
      real(kind=8), dimension(-nord/2:nord/2,-nord/2:nord/2) :: c1D,c1DF
      real(kind=8) :: hx,hy,hz,dx,dy,dz,res,rho
      real(kind=8) :: rpoints
      real(dp), dimension(3) :: tmp
      logical, dimension(3) :: peri
      logical :: perx,pery,perz

      n = nord+1
      m = nord/2
      n01=mesh%ndims(1)
      n02=mesh%ndims(2)
      n03=mesh%ndims(3)
      hx = mesh%hgrids(1)!acell/real(n01,kind=8)
      hy = mesh%hgrids(2)!acell/real(n02,kind=8)
      hz = mesh%hgrids(3)!acell/real(n03,kind=8)
      n_cell = max(n01,n02,n03)
      rpoints=product(real([n01,n02,n03],dp))

      !buffers associated to the geocode
      !conditions for periodicity in the three directions
      peri=cell_periodic_dims(mesh)
      perx=peri(1)
      pery=peri(2)
      perz=peri(3)

      ! Beware that n_cell has to be > than n.
      if (n_cell.lt.n) then
         write(*,*)'ngrid in has to be set > n=nord + 1'
         stop
      end if

      ! Setting of 'nord' order accurate first derivative coefficient from 'Matematica'.
      !Only nord=2,4,6,8,16
      if (all(nord /=[2,4,6,8,16])) then
         write(*,*)'Only nord-order 2,4,6,8,16 accurate first derivative'
         stop
      end if

      do i=-m,m
         do j=-m,m
            c1D(i,j)=0.d0
            c1DF(i,j)=0.d0
         end do
      end do

      include 'FiniteDiffCorff.inc'

      rhores2=0.d0
      !$omp parallel do default(shared) &
      !$omp private(i1,i2,i3,j,ii, dx,dy,dz,res,rho)&
      !$omp reduction(+:rhores2)
      do i3=1,n03
         do i2=1,n02
            do i1=1,n01
               dx=0.d0
               if (i1.le.m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j + n01 - 1, n01 ) + 1
                        dx = dx + c1D(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        dx = dx + c1D(j,i1-m-1)*u(j+m+1,i2,i3)
                     end do
                  end if
               else if (i1.gt.n01-m) then
                  if (perx) then
                     do j=-m,m
                        ii=modulo(i1 + j - 1, n01 ) + 1
                        dx = dx + c1D(j,0)*u(ii,i2,i3)
                     end do
                  else
                     do j=-m,m
                        dx = dx + c1D(j,i1-n01+m)*u(n01 + j - m,i2,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     dx = dx + c1D(j,0)*u(i1 + j,i2,i3)
                  end do
               end if
               tmp(1)=dx/hx
               dy = 0.0d0
               if (i2.le.m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j + n02 - 1, n02 ) + 1
                        dy = dy + c1D(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        dy = dy + c1D(j,i2-m-1)*u(i1,j+m+1,i3)
                     end do
                  end if
               else if (i2.gt.n02-m) then
                  if (pery) then
                     do j=-m,m
                        ii=modulo(i2 + j - 1, n02 ) + 1
                        dy = dy + c1D(j,0)*u(i1,ii,i3)
                     end do
                  else
                     do j=-m,m
                        dy = dy + c1D(j,i2-n02+m)*u(i1,n02 + j - m,i3)
                     end do
                  end if
               else
                  do j=-m,m
                     dy = dy + c1D(j,0)*u(i1,i2 + j,i3)
                  end do
               end if
               tmp(2)=dy/hy
               dz = 0.0d0
               if (i3.le.m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j + n03 - 1, n03 ) + 1
                        dz = dz + c1D(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        dz = dz + c1D(j,i3-m-1)*u(i1,i2,j+m+1)
                     end do
                  end if
               else if (i3.gt.n03-m) then
                  if (perz) then
                     do j=-m,m
                        ii=modulo(i3 + j - 1, n03 ) + 1
                        dz = dz + c1D(j,0)*u(i1,i2,ii)
                     end do
                  else
                     do j=-m,m
                        dz = dz + c1D(j,i3-n03+m)*u(i1,i2,n03 + j - m)
                     end do
                  end if
               else
                  do j=-m,m
                     dz = dz + c1D(j,0)*u(i1,i2,i3 + j)
                  end do
               end if
               tmp(3)=dz/hz

               !retrieve the previous treatment
               res = dotp(mesh,dlogeps(1,i1,i2,i3),tmp)
               res = res*oneofourpi
               rho=rhopol(i1,i2,i3)
               res=res-rho
               res=eta*res
               !rhores2=rhores2+res*res*eps(i1,i2,i3)*eps(i1,i2,i3)
               rhores2=rhores2+res*res
               rhopol(i1,i2,i3)=res+rho

            end do
         end do
      end do
      !$omp end parallel do

    end subroutine update_rhopol

    !subroutine Apply_GPe_operator(nord,geocode,ndims,hgrids,eps,pot,a_pot,work)
    subroutine Delta_GPe_operator(nord,geocode,ndims,hgrids,eps,pot,a_pot)
      use wrapper_linalg, only: axpy
      use numerics, only: oneofourpi
      implicit none
      integer, intent(in) :: nord
      character(len=1), intent(in) :: geocode
      integer, dimension(3), intent(in) :: ndims
      real(gp), dimension(3), intent(in) :: hgrids
      !>dielectric function
      real(dp), dimension(ndims(1),ndims(2),ndims(3)), intent(in) :: eps
      real(dp), dimension(ndims(1),ndims(2),ndims(3)), intent(in) :: pot
      real(dp), dimension(ndims(1),ndims(2),ndims(3)), intent(inout) :: a_pot
      !>work array containing in output the gradient of the potential
      !!multiplied by the dielectric function
      !real(dp), dimension(ndims(1),ndims(2),ndims(3),3), intent(out) :: work
      !LG @GF: work arrays should not be used that way, there creates stack overflows
      real(dp), dimension(ndims(1),ndims(2),ndims(3),3) :: work
      real(dp), dimension(ndims(1),ndims(2),ndims(3)) :: work1
!      real(dp) :: pi
      integer :: i1,i2,i3
     !local variables
      !first calculate the derivative of the potential

      call nabla_u_epsilon(geocode,ndims(1),ndims(2),ndims(3),pot,work,nord,hgrids,eps)
      call div_u_i(geocode,ndims(1),ndims(2),ndims(3),work,work1,nord,hgrids)

      do i3=1,ndims(3)
       do i2=1,ndims(2)
        do i1=1,ndims(1)
         a_pot(i1,i2,i3)=a_pot(i1,i2,i3) + oneofourpi*work1(i1,i2,i3)
        end do
       end do
      end do

    end subroutine Delta_GPe_operator

    !> verify that the density is considerably zero in the region where epsilon is different from one
    subroutine nonvacuum_projection(n1,n23,rho,oneoeps,norm)
      implicit none
      integer, intent(in) :: n1,n23 !< parallelized box dimensions
      real(dp), dimension(n1,n23), intent(in) :: rho !<charge density
      !>inverse of epsilon (might also be inverse of sqrt(eps))
      real(dp), dimension(n1,n23), intent(in) :: oneoeps
      real(dp), intent(out) :: norm !< \int of rho where epsilon /=1
      !local variables
      integer :: i1,i23
      real(dp), parameter :: tol= 5.d-1

      norm=0.0_dp
      !$omp parallel do default(shared) private(i1,i23)&
      !$omp reduction(+:norm)
      do i23=1,n23
         do i1=1,n1
            if (abs(oneoeps(i1,i23) - 1.0_dp) > tol) norm=norm+rho(i1,i23)
         end do
      end do
      !$omp end parallel do

    end subroutine nonvacuum_projection

    !>calculate polarization charge and epsilon for plotting purposes.
    !! no need of gathering the results as the arrays are already given in full form
    subroutine polarization_charge(geocode,ndims,hgrids,nord,rho,pot,nabla_pot,rho_pol)
      use numerics, only: oneofourpi
      implicit none
      integer, intent(in) :: nord
      character(len=1), intent(in) :: geocode
      integer, dimension(3), intent(in) :: ndims
      real(dp), dimension(3), intent(in) :: hgrids
      real(dp), dimension(ndims(1),ndims(2),ndims(3)), intent(in) :: rho,pot
      real(dp), dimension(ndims(1),ndims(2),ndims(3),3), intent(out) :: nabla_pot
      real(dp), dimension(ndims(1),ndims(2),ndims(3)), intent(out) :: rho_pol
      !local variables
      integer :: i1,i2,i3

      call nabla_u(geocode,ndims(1),ndims(2),ndims(3),pot,nabla_pot,nord,hgrids)
      call div_u_i(geocode,ndims(1),ndims(2),ndims(3),nabla_pot,rho_pol,nord,hgrids)

      do i3=1,ndims(3)
         do i2=1,ndims(2)
            do i1=1,ndims(1)
               !this section has to be inserted into a optimized calculation of the
               !derivative
               rho_pol(i1,i2,i3)=-oneofourpi*rho_pol(i1,i2,i3)-rho(i1,i2,i3)
            end do
         end do
      end do

    end subroutine polarization_charge

end module FDder
