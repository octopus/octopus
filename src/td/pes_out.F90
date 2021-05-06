!! Copyright (C) 2016 U. De Giovannini
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module pes_out_oct_m
  use cube_function_oct_m
  use cube_oct_m
  use global_oct_m
  use io_oct_m
  use io_function_oct_m
  use loct_oct_m
  use math_oct_m
  use messages_oct_m
  use namespace_oct_m
#if defined(HAVE_NETCDF)
  use netcdf
#endif    
  use profiling_oct_m
  use qshep_oct_m
  use space_oct_m
  use sort_oct_m
  use string_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use vtk_oct_m
    
  implicit none

  private

  public ::                             &
    pes_out_velocity_map,               &
    pes_out_velocity_map_cut,           &
    pes_out_arpes_cut

  integer, parameter ::   &
    INTEGRATE_NONE    = -1,       &
    INTEGRATE_PHI     =  1,       &
    INTEGRATE_THETA   =  2,       &
    INTEGRATE_R       =  3,       &
    INTEGRATE_KX      =  4,       &
    INTEGRATE_KY      =  5,       &
    INTEGRATE_KZ      =  6



contains

  ! ---------------------------------------------------------
  subroutine pes_out_velocity_map(pesK, file, namespace, space, Lk, ll, how, pmesh)
    FLOAT,             intent(in) :: pesK(:,:,:)
    character(len=*),  intent(in) :: file
    type(namespace_t), intent(in) :: namespace
    type(space_t),     intent(in) :: space
    FLOAT,             intent(in) :: Lk(:,:)
    integer,           intent(in) :: ll(:)  
    integer(8),        intent(in) :: how
    FLOAT, optional,   intent(in) :: pmesh(:,:,:,:)  
  
    integer :: ierr
    character(len=512) :: filename
    type(cube_t) :: cube
    type(cube_function_t) :: cf
    FLOAT :: dk(3)  

    PUSH_SUB(pes_out_velocity_map)

    call cube_init(cube, ll, namespace, space)
    call dcube_function_alloc_RS(cube, cf, force_alloc = .true.)
    cf%dRS = pesK
  
  
    if (.not. present(pmesh) ) then
      ! Ignore Lk and use pmesh
      dk(:) = M_ZERO
      dk(1:space%dim) = abs(Lk(2, 1:space%dim) - Lk(1,1:space%dim))
    end if
  
#if defined(HAVE_NETCDF)  
  
    if(bitand(how, OPTION__OUTPUTFORMAT__NETCDF) /= 0) then
      filename = trim(file)//".ncdf"
      write(message(1), '(a)') 'Writing netcdf format file: '
      call messages_info(1)
  
      call dout_cf_netcdf(filename, ierr, cf, cube, space, dk(:), .false., unit_one/units_out%length)

    end if

#endif
  
    if(bitand(how, OPTION__OUTPUTFORMAT__VTK) /= 0)  then
      filename = trim(file)//".vtk"
      write(message(1), '(a)') 'Writing vtk format file: '
      call messages_info(1)
    
      if (present(pmesh)) then          
        call dvtk_out_cf_structured(filename, namespace, 'PES_vel_map', ierr, cf, cube,& 
           unit_one/units_out%length, pmesh, ascii = .false.)
      else 
        call dvtk_out_cf(filename, namespace, 'PES_vel_map', ierr, cf, cube, dk(:),& 
          unit_one/units_out%length)
      end if        
    end if
  
    call cube_end(cube)
    call dcube_function_free_RS(cube, cf)

    POP_SUB(pes_out_velocity_map)
    
  end subroutine pes_out_velocity_map


  ! ---------------------------------------------------------
  subroutine pes_out_arpes_cut(namespace, arpes, file, dim, ll, pmesh, Ekin)
    type(namespace_t), intent(in) :: namespace
    FLOAT,             intent(in) :: arpes(:,:,:)
    character(len=*),  intent(in) :: file
    integer,           intent(in) :: dim
    integer,           intent(in) :: ll(:)
    FLOAT,             intent(in) :: pmesh(:,:,:,:)
    FLOAT,             intent(in) :: Ekin(:,:,:)
    
    integer :: iunit, ip, ie
    FLOAT   :: dp, length, pp(1:2), pdiff(1:2)
    
    PUSH_SUB(pes_out_arpes_cut)
    
    iunit = io_open(file, namespace, action='write')
    write(iunit, '(a)') '##################################################'
    write(iunit, '(a1,a18,2x,a18,2x,a18,2x,a18,2x, a18,2x,a18)') '#', &
                                      str_center("Ppath", 18), &
                                      str_center("px", 18), str_center("py", 18), &
                                      str_center("E", 18),  str_center("P[Ppath,E]", 18)

    write(iunit, '(a1,a18,2x,a18,2x,a18,2x,a18,2x, a18,2x,a18)') '#', &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &
                                      str_center('[hbar/'//trim(units_abbrev(units_out%length)) // ']', 18), &    
                                      str_center('['//trim(units_abbrev(units_out%energy)) // ']', 18), &
                                      str_center('[1/' //trim(units_abbrev(units_out%energy))//']', 18)
    write(iunit, '(a)') '##################################################'
    
    
    length = M_ZERO
    pdiff(:) = M_ZERO
    
    do ip = 1, ll(1)
      pp(1:2) = pmesh(ip,1,1,1:2)
      if (ip > 1) pdiff = pp(1:2) - pmesh(ip-1,1,1,1:2)
      dp = sqrt(sum(pdiff(1:2)**2))
      length = length + dp 

      select case (dim)
        case (2)
          do ie = 1, ll(2) 
            write(iunit, '(es19.12,2x,es19.12,2x,es19.12,2x,es19.12,2x,es19.12,2x,es19.12)')   &
                                            units_from_atomic(unit_one/units_out%length, length),&
                                            units_from_atomic(unit_one/units_out%length, pp(1)), &
                                            units_from_atomic(unit_one/units_out%length, pp(2)), &
                                            units_from_atomic(units_out%energy, Ekin(ip,ie,1)), &
                                            arpes(ip,ie,1)
          end do      

        case(3)
          do ie = 1, ll(3)
            write(iunit, '(es19.12,2x,es19.12,2x,es19.12,2x,es19.12,2x,es19.12,2x,es19.12)')   &
                                            units_from_atomic(unit_one/units_out%length, length),&
                                            units_from_atomic(unit_one/units_out%length, pp(1)), &
                                            units_from_atomic(unit_one/units_out%length, pp(2)), &
                                            units_from_atomic(units_out%energy, Ekin(ip,1,ie)), &
                                            arpes(ip,1,ie)
          end do
          
      end select
      

      write(iunit, *)
      
       
    end do
    
    
    call io_close(iunit)
    
    POP_SUB(pes_out_arpes_cut)    
  end subroutine pes_out_arpes_cut


  ! ---------------------------------------------------------
  subroutine pes_out_velocity_map_cut(namespace, pesK, file, ll, dim, pol, dir, integrate, pos, Lk, pmesh)
    type(namespace_t), intent(in) :: namespace
    FLOAT,             intent(in) :: pesK(:,:,:)
    character(len=*),  intent(in) :: file
    integer,           intent(in) :: ll(:)
    integer,           intent(in) :: dim
    FLOAT,             intent(in) :: pol(3)
    integer,           intent(in) :: dir
    integer,           intent(in) :: integrate
    integer, optional, intent(in) :: pos(3)
    FLOAT, optional,   intent(in) :: Lk(:,:)
    FLOAT, optional,   intent(in) :: pmesh(:,:,:,:)

    integer              :: ii, ix, iy, iunit,ldir(2), icut(3)
    FLOAT                :: KK(3),temp
    integer, allocatable :: idx(:,:)
    FLOAT, allocatable   :: Lk_(:,:)
    FLOAT                :: rotation(1:dim,1:dim)
    logical              :: aligned_axis
  ! integration
    FLOAT                :: K, KKK(3), theta, phi, Dphi, Dk(3)
    integer              :: iph, Nphi
  ! progress
    integer              :: idone, ntodo

    FLOAT, allocatable :: cube_f(:)
    type(qshep_t) :: interp

    PUSH_SUB(pes_out_velocity_map_cut)

    iunit = io_open(file, namespace, action='write')



    ASSERT(size(pesK, 1) == ll(1))

    if (.not. present(pmesh)) then
      ASSERT(present(Lk))

      SAFE_ALLOCATE(idx(1:maxval(ll(:)), 1:3))
      SAFE_ALLOCATE(Lk_(1:maxval(ll(:)), 1:3))

      Dk(:) = M_ZERO
      Dk(1:dim) = abs(Lk(2,1:dim)-Lk(1,1:dim))

      do ii = 1, 3
        Lk_(:,ii) = Lk(:,ii)
        call sort(Lk_(1:ll(ii), ii), idx(1:ll(ii), ii)) !We need to sort the k-vectors in order to dump in gnuplot format
      end do
    end if

    aligned_axis = sum((pol-(/0 ,0 ,1/))**2)  <= M_EPSILON  .or. &
                   sum((pol-(/0 ,1 ,0/))**2)  <= M_EPSILON  .or. &
                   sum((pol-(/1 ,0 ,0/))**2)  <= M_EPSILON

    if (present(pos)) then
      icut(1:3) = pos(1:3)
    else
      icut(1:3) = ll(1:3)/2 + 1
    end if

    if (aligned_axis .and. integrate == INTEGRATE_NONE) then !no need to rotate and interpolate

      select case (dir)
        case (1)
          ldir(:) =(/2,3/)
        case (2)
          ldir(:) =(/1,3/)
        case (3)
          ldir(:) =(/1,2/)

      end select


      if (present(pmesh)) then
        do ix = 1, ll(ldir(1))
          do iy = 1, ll(ldir(2))

            select case (dir)
              case (1)
                temp = pesK(icut(dir), ix, iy)
                KK(1) = pmesh(icut(dir), ix, iy, ldir(1))
                KK(2) = pmesh(icut(dir), ix, iy, ldir(2))
              case (2)
                temp = pesK(ix, icut(dir), iy)
                KK(1) = pmesh(ix, icut(dir), iy, ldir(1))
                KK(2) = pmesh(ix, icut(dir), iy, ldir(2))
              case (3)
                temp = pesK(ix, iy, icut(dir))
                KK(1) = pmesh(ix, iy, icut(dir), ldir(1))
                KK(2) = pmesh(ix, iy, icut(dir), ldir(2))
            end select

            write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
                    units_from_atomic(sqrt(units_out%energy), KK(1)),&
                    units_from_atomic(sqrt(units_out%energy), KK(2)),&
                    temp
          end do
          write(iunit, *)
        end do

      else
        do ix = 1, ll(ldir(1))
          KK(1) = Lk_(ix, ldir(1))
          do iy = 1, ll(ldir(2))
            KK(2) = Lk_(iy, ldir(2))

            select case (dir)
              case (1)
                temp = pesK(idx(icut(dir), 1), idx(ix, 2), idx(iy, 3))
              case (2)
                temp = pesK(idx(ix, 1), idx(icut(dir), 2), idx(iy, 3))
              case (3)
                temp = pesK(idx(ix, 1), idx(iy, 2), idx(icut(dir), 3))
            end select

            write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
                    units_from_atomic(sqrt(units_out%energy), KK(1)),&
                    units_from_atomic(sqrt(units_out%energy), KK(2)),&
                    temp
          end do
          write(iunit, *)
        end do
      end if

    else
      ! We set the z-axis along the pol vector
      call generate_rotation_matrix(rotation, (/M_ZERO, M_ZERO, M_ONE/), pol )

      if(debug%info) then
        print *,"Rotate z-axis over the zenith axis"
        print *,rotation(1,:)
        print *,rotation(2,:)
        print *,rotation(3,:)
      end if

      call pes_out_interpolator_init(namespace, pesK, Lk, ll, dim, cube_f, interp, pmesh)

      ntodo = product(ll(1:2))
      idone = 0
      call loct_progress_bar(-1, ntodo)

      do ix = 1, ll(1)
       do iy = 1, ll(2)

         !cut
         select case (dir)
           case (1)
             KK(1) = M_ZERO
             if (present(pmesh)) then
               KK(2) = pmesh(ix, 1, 1, 1)
               KK(3) = pmesh(1, iy, 1, 2)
             else
               KK(2) = Lk_(ix, 1)
               KK(3) = Lk_(iy, 2)
             end if
           case (2)
             KK(2) = M_ZERO
             if(present(pmesh)) then
               KK(1) = pmesh(ix,  1, 1, 1)
               KK(3) = pmesh(1,  iy, 1, 2)
             else
               KK(1) = Lk_(ix, 1)
               KK(3) = Lk_(iy, 2)
             end if

           case (3)
             KK(3) = M_ZERO
             if(present(pmesh)) then
               KK(1) = pmesh(ix,  1, 1, 1)
               KK(2) = pmesh(1,  iy, 1, 2)
             else
               KK(1) = Lk_(ix, 1)
               KK(2) = Lk_(iy, 2)
             end if

         end select

         temp = qshep_interpolate(interp, cube_f, matmul(rotation,KK(1:3)) )

         select case (integrate)
           case (INTEGRATE_PHI)
             temp = M_ZERO
             K = sqrt(KK(1)**2 + KK(2)**2 + KK(3)**2)

             Nphi = 360
             Dphi = M_TWO * M_PI/Nphi

             do iph = 0, Nphi
               phi = iph * Dphi
               theta = atan2(sqrt(KK(1)**2+KK(2)**2),KK(3))

               KKK(1) = K *sin(theta)*cos(phi)
               KKK(2) = K *sin(theta)*sin(phi)
               KKK(3) = K *cos(theta)

               temp = temp + &
                             abs(qshep_interpolate(interp, cube_f, matmul(rotation,KKK(1:3)) ))
             end do
             temp = temp * Dphi

           case (INTEGRATE_KX)
             temp = M_ZERO
             do ii =1, ll(1)
                KKK(:) = KK(:) + (/Lk(ii,1), M_ZERO, M_ZERO/)
                temp = temp + &
                              abs(qshep_interpolate(interp, cube_f, matmul(rotation,KKK(1:3)) ))
             end do
             temp = temp * Dk(1)

           case (INTEGRATE_KY)
             temp = M_ZERO
             do ii =1, ll(2)
                KKK(:) = KK(:) + (/M_ZERO, Lk(ii,2), M_ZERO/)
                temp = temp + &
                              abs(qshep_interpolate(interp, cube_f, matmul(rotation,KKK(1:3)) ))
             end do
             temp = temp * Dk(2)

           case (INTEGRATE_KZ)
             temp = M_ZERO
             do ii =1, ll(3)
                KKK(:) = KK(:) + (/M_ZERO, M_ZERO, Lk(ii,3)/)
                temp = temp + &
                              abs(qshep_interpolate(interp, cube_f, matmul(rotation,KKK(1:3)) ))
             end do
             temp = temp * Dk(3)

         end select

         write(iunit, '(es19.12,2x,es19.12,2x,es19.12)') &
                 units_from_atomic(sqrt(units_out%energy), KK(1)),&
                 units_from_atomic(sqrt(units_out%energy), KK( 2)),&
                 temp


         idone = idone +1
         call loct_progress_bar(idone, ntodo)

       end do
       write(iunit, *)
      end do

      write(stdout, '(1x)')


    end if

    call io_close(iunit)

    call pes_out_interpolator_end(cube_f, interp)

    SAFE_DEALLOCATE_A(idx)
    SAFE_DEALLOCATE_A(Lk_)

    POP_SUB(pes_out_velocity_map_cut)
  end subroutine pes_out_velocity_map_cut

  ! --------------------------------------------------------
  !
  !>  Qshep interpolation helper function initialization.
  !!  Generates the linearized version of pesK (cube_f) and the associated
  !!  qshep interpolator opbject (interp).
  !
  ! ---------------------------------------------------------
  subroutine pes_out_interpolator_init(namespace, pesK, Lk, ll, dim, cube_f, interp, pmesh)
    type(namespace_t),  intent(in)    :: namespace
    FLOAT,              intent(in)    :: pesK(:,:,:)
    FLOAT,              intent(in)    :: Lk(:,:)
    integer,            intent(in)    :: ll(:)
    integer,            intent(in)    :: dim
    FLOAT, allocatable, intent(out)   :: cube_f(:)
    type(qshep_t),      intent(out)   :: interp
    FLOAT, optional,    intent(in)    :: pmesh(:,:,:,:)

    integer :: np, ii, ix, iy, iz
    FLOAT   :: KK(3)
    FLOAT, allocatable ::  kx(:),ky(:),kz(:)

    PUSH_SUB(pes_out_interpolator_init)


    call messages_write("Initializing Qshep interpolator. Be patient it may take a while... ")
    call messages_info()

    np = ll(1)*ll(2)*ll(3)

    !check dim
    if (dim  <  2 .or. dim > 3) then
      message(1) = "This interpolator works only for 2 <= dim <= 3."
      call messages_fatal(1, namespace=namespace)
    end if

    SAFE_ALLOCATE(cube_f(1:np))

    SAFE_ALLOCATE(kx(1:np))
    SAFE_ALLOCATE(ky(1:np))
    SAFE_ALLOCATE(kz(1:np))

    cube_f = M_ZERO
    kx = M_ZERO
    ky = M_ZERO
    kz = M_ZERO


    ii=1
    if (present(pmesh)) then
      do ix = 1, ll(1)
        do iy = 1, ll(2)
          do iz = 1, ll(3)

            cube_f(ii) =  pesK(ix,iy,iz)

            kx(ii) = pmesh(ix,iy,iz,1)
            ky(ii) = pmesh(ix,iy,iz,2)
            kz(ii) = pmesh(ix,iy,iz,3)

            ii = ii +1
          end do
        end do
      end do

    else

      do ix = 1, ll(1)
        KK(1) = Lk(ix, 1)
        do iy = 1, ll(2)
          KK(2) = Lk(iy, 2)
          do iz = 1, ll(3)
            KK(3) = Lk(iz, 3)

            cube_f(ii) =  pesK(ix,iy,iz)

            kx(ii) = KK(1)
            ky(ii) = KK(2)
            kz(ii) = KK(3)

            ii = ii +1
          end do
        end do
      end do
    end if


    select case(dim)
      case (2)
        call qshep_init(interp, np, cube_f, kx, ky)
      case (3)
        call qshep_init(interp, np, cube_f, kx, ky, kz)
    end select

    SAFE_DEALLOCATE_A(kx)
    SAFE_DEALLOCATE_A(ky)
    SAFE_DEALLOCATE_A(kz)

    call messages_write("done")
    call messages_new_line()
    call messages_info()


    POP_SUB(pes_out_interpolator_init)
  end subroutine pes_out_interpolator_init

  ! ---------------------------------------------------------
  !>  Destroy the interpolation objects
  ! ---------------------------------------------------------
  subroutine pes_out_interpolator_end(cube_f, interp)
    FLOAT, allocatable, intent(inout) :: cube_f(:)
    type(qshep_t),      intent(inout) :: interp

    PUSH_SUB(pes_out_interpolator_end)

    call qshep_end(interp)

    SAFE_DEALLOCATE_A(cube_f)

    POP_SUB(pes_out_interpolator_end)
  end subroutine pes_out_interpolator_end

end module pes_out_oct_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
