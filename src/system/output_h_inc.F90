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
  subroutine output_hamiltonian(hm, mesh, dir, outp, geo)
    type(hamiltonian_t),   intent(in) :: hm
    type(mesh_t),          intent(in) :: mesh
    character(len=*),      intent(in) :: dir
    type(output_t),        intent(in) :: outp
    type(geometry_t),      intent(in) :: geo

    integer :: is, err, idir
    character(len=80) :: fname

    FLOAT, allocatable :: v0(:,:)
    
    PUSH_SUB(output_hamiltonian)

    if(iand(outp%what, C_OUTPUT_POTENTIAL).ne.0) then
      SAFE_ALLOCATE(v0(1:mesh%np, 1:hm%d%dim))
      v0(1:mesh%np, 1) = hm%ep%vpsl(1:mesh%np)
      call dio_function_output(outp%how, dir, "v0", mesh, v0(:, 1), units_out%energy, err, geo = geo)
      SAFE_DEALLOCATE_A(v0)

      if(hm%ep%classical_pot > 0) then
        call dio_function_output(outp%how, dir, "vc", mesh, hm%ep%Vclassical, units_out%energy, err, geo = geo)
      end if

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call dio_function_output(outp%how, dir, 'vh', mesh, hm%vhartree, units_out%energy, err, geo = geo)
        do is = 1, min(hm%d%ispin, 2)
          if(hm%d%ispin == 1) then
            write(fname, '(a)') 'vxc'
          else
            write(fname, '(a,i1)') 'vxc-sp', is
          endif
          call dio_function_output(outp%how, dir, fname, mesh, hm%vxc(:, is), units_out%energy, err, geo = geo)

          ! finally the full KS potential (without non-local PP contributions)
          if(hm%d%ispin == 1) then
            write(fname, '(a)') 'vks'
          else
            write(fname, '(a,i1)') 'vks-sp', is
          endif
          if (hm%ep%classical_pot > 0) then
            call dio_function_output(outp%how, dir, fname, mesh, &
              hm%ep%vpsl + hm%ep%Vclassical + hm%vhxc(:, is), units_out%energy, err, geo = geo)
          else
            call dio_function_output(outp%how, dir, fname, mesh, &
              hm%ep%vpsl + hm%vhxc(:, is), units_out%energy, err, geo = geo)
          end if
        end do
      end if

      if(hm%self_induced_magnetic) then
        ! unit of magnetic field is same as of electric field, and same as force (since e = 1)
        select case(mesh%sb%dim)
        case(3)
          do idir = 1, mesh%sb%dim
            call dio_function_output(outp%how, dir, 'Bind_'//index2axis(idir), mesh, hm%b_ind(:, idir), &
              units_out%force, err, geo = geo)
          enddo
        case(2)
          call dio_function_output(outp%how, dir, 'Bind_z', mesh, hm%b_ind(:, 1), units_out%force, err, geo = geo)
        end select
      end if
    end if

    POP_SUB(output_hamiltonian)
  end subroutine output_hamiltonian

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
