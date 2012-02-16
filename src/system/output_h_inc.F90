!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch.
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
!! $Id: output_linear_h.F90 5765 2009-07-24 20:47:06Z dstrubbe $

  ! ---------------------------------------------------------
  subroutine output_hamiltonian(hm, der, dir, outp, geo)
    type(hamiltonian_t),   intent(in) :: hm
    type(derivatives_t),   intent(in) :: der
    character(len=*),      intent(in) :: dir
    type(output_t),        intent(in) :: outp
    type(geometry_t),      intent(in) :: geo

    integer :: is, err, idir
    character(len=80) :: fname
    FLOAT, allocatable :: v0(:,:), nxc(:)
    
    PUSH_SUB(output_hamiltonian)

    if(iand(outp%what, C_OUTPUT_POTENTIAL).ne.0) then
      SAFE_ALLOCATE(v0(1:der%mesh%np, 1:hm%d%dim))
      v0(1:der%mesh%np, 1) = hm%ep%vpsl(1:der%mesh%np)
      call dio_function_output(outp%how, dir, "v0", der%mesh, v0(:, 1), units_out%energy, err, geo = geo)
      SAFE_DEALLOCATE_A(v0)

      if(hm%ep%classical_pot > 0) then
        call dio_function_output(outp%how, dir, "vc", der%mesh, hm%ep%Vclassical, units_out%energy, err, geo = geo)
      end if

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call dio_function_output(outp%how, dir, 'vh', der%mesh, hm%vhartree, units_out%energy, err, geo = geo)
        do is = 1, min(hm%d%ispin, 2)
          if(hm%d%ispin == 1) then
            write(fname, '(a)') 'vxc'
          else
            write(fname, '(a,i1)') 'vxc-sp', is
          endif
          call dio_function_output(outp%how, dir, fname, der%mesh, hm%vxc(:, is), units_out%energy, err, geo = geo)

          ! finally the full KS potential (without non-local PP contributions)
          if(hm%d%ispin == 1) then
            write(fname, '(a)') 'vks'
          else
            write(fname, '(a,i1)') 'vks-sp', is
          endif
          if (hm%ep%classical_pot > 0) then
            call dio_function_output(outp%how, dir, fname, der%mesh, &
              hm%ep%vpsl + hm%ep%Vclassical + hm%vhxc(:, is), units_out%energy, err, geo = geo)
          else
            call dio_function_output(outp%how, dir, fname, der%mesh, &
              hm%ep%vpsl + hm%vhxc(:, is), units_out%energy, err, geo = geo)
          end if
        end do
      end if

      if(hm%self_induced_magnetic) then
        ! unit of magnetic field is same as of electric field, and same as force (since e = 1)
        select case(der%mesh%sb%dim)
        case(3)
          do idir = 1, der%mesh%sb%dim
            call dio_function_output(outp%how, dir, 'Bind_'//index2axis(idir), der%mesh, hm%b_ind(:, idir), &
              units_out%force, err, geo = geo)
          enddo
        case(2)
          call dio_function_output(outp%how, dir, 'Bind_z', der%mesh, hm%b_ind(:, 1), units_out%force, err, geo = geo)
        end select
      end if
    end if

    if(iand(outp%what, C_OUTPUT_XC_DENSITY) /= 0 .and. hm%theory_level /= INDEPENDENT_PARTICLES) then
      SAFE_ALLOCATE(v0(1:der%mesh%np_part, 1))
      SAFE_ALLOCATE(nxc(1:der%mesh%np))

      do is = 1, min(hm%d%ispin, 2)
        if(hm%d%ispin == 1) then
          write(fname, '(a)') 'nxc'
        else
          write(fname, '(a,i1)') 'nxc-sp', is
        endif
                
        v0(1:der%mesh%np, 1) = hm%vxc(1:der%mesh%np, is)

        call dderivatives_lapl(der, v0(:, 1), nxc)

        call dio_function_output(outp%how, dir, fname, der%mesh, nxc, units_out%energy, err, geo = geo)
        
      end do

      SAFE_DEALLOCATE_A(v0)
      SAFE_DEALLOCATE_A(nxc)
    end if


    POP_SUB(output_hamiltonian)
  end subroutine output_hamiltonian


  ! ---------------------------------------------------------
  subroutine output_scalar_pot(outp, gr, geo, hm, dir)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(hamiltonian_t),  intent(inout) :: hm
    type(output_t),       intent(in)    :: outp
    character(len=*),     intent(in)    :: dir

    integer :: is, err
    character(len=80) :: fname
    FLOAT, allocatable :: scalar_pot(:)

    PUSH_SUB(output_scalar_pot)

    if(iand(outp%what, C_OUTPUT_TD_POTENTIAL).ne.0) then
      SAFE_ALLOCATE(scalar_pot(1:gr%mesh%np))
      do is = 1, hm%ep%no_lasers
        write(fname, '(a,i1)') 'scalar_pot-', is
        scalar_pot = M_ZERO
        call laser_potential(hm%ep%lasers(is), gr%mesh, scalar_pot)
        write(0, *) trim(fname)
        call dio_function_output(outp%how, dir, fname, gr%mesh, scalar_pot, units_out%energy, err, geo = geo)
      end do
      SAFE_DEALLOCATE_A(scalar_pot)
    end if

    POP_SUB(output_scalar_pot)
  end subroutine output_scalar_pot


  ! ---------------------------------------------------------
  subroutine output_kick(outp, gr, geo, hm, dir)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(in)    :: geo
    type(hamiltonian_t),  intent(inout) :: hm
    type(output_t),       intent(in)    :: outp
    character(len=*),     intent(in)    :: dir

    integer :: err
    CMPLX, allocatable :: kick_function(:)
    
    PUSH_SUB(output_kick)

    if(iand(outp%what, C_OUTPUT_KICK_FUNCTION).ne.0) then
      SAFE_ALLOCATE(kick_function(1:gr%mesh%np))
      call kick_function_get(gr, hm%ep%kick, kick_function)
      call zio_function_output(outp%how, dir, "kick_function", gr%mesh, kick_function(:), &
        units_out%energy, err, geo = geo)
      SAFE_DEALLOCATE_A(kick_function)
    end if
    
    POP_SUB(output_kick)
  end subroutine output_kick


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
