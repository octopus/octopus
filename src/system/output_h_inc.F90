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
  subroutine h_sys_output_hamiltonian(hm, m, sb, dir, outp, geo)
    type(hamiltonian_t),   intent(in) :: hm
    type(mesh_t),          intent(in) :: m
    type(simul_box_t),     intent(in) :: sb
    character(len=*),      intent(in) :: dir
    type(h_sys_output_t),  intent(in) :: outp
    type(geometry_t),      intent(in) :: geo

    integer :: is, err
    character(len=80) :: fname

    FLOAT, allocatable :: v0(:,:)
    
    call push_sub('output_h.h_sys_output_hamiltonian')

    if(iand(outp%what, output_potential).ne.0) then
      SAFE_ALLOCATE(v0(1:m%np, 1:hm%d%dim))
      v0(1:m%np, 1) = hm%ep%vpsl(1:m%np)
      call doutput_function(outp%how, dir, "v0", m, sb, v0(:, 1), units_out%energy, err, geo = geo)
      SAFE_DEALLOCATE_A(v0)

      if(hm%ep%classical_pot > 0) then
        call doutput_function(outp%how, dir, "vc", m, sb, hm%ep%Vclassical, units_out%energy, err, geo = geo)
      end if

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call doutput_function(outp%how, dir, 'vh', m, sb, hm%vhartree, units_out%energy, err, geo = geo)
        do is = 1, min(hm%d%ispin, 2)
          write(fname, '(a,i1)') 'vxc-', is
          call doutput_function(outp%how, dir, fname, m, sb, hm%vxc(:, is), units_out%energy, err, geo = geo)

          ! finally the full KS potential (without non-local PP contributions)
          write(fname, '(a,i1)') 'vks-', is
          if (hm%ep%classical_pot > 0) then
            call doutput_function(outp%how, dir, fname, m, sb, &
              hm%ep%vpsl + hm%ep%Vclassical + hm%vhxc(:, is), units_out%energy, err, geo = geo)
          else
            call doutput_function(outp%how, dir, fname, m, sb, &
              hm%ep%vpsl + hm%vhxc(:, is), units_out%energy, err, geo = geo)
          end if
        end do
      end if

      if(hm%self_induced_magnetic) then
        select case(sb%dim)
        case(3)
          call doutput_function(outp%how, dir, 'Bind_x', m, sb, hm%b_ind(:, 1), unit_one, err, geo = geo)
          call doutput_function(outp%how, dir, 'Bind_y', m, sb, hm%b_ind(:, 2), unit_one, err, geo = geo)
          call doutput_function(outp%how, dir, 'Bind_z', m, sb, hm%b_ind(:, 3), unit_one, err, geo = geo)
        case(2)
          call doutput_function(outp%how, dir, 'Bind_z', m, sb, hm%b_ind(:, 1), unit_one, err, geo = geo)
        end select
      end if
    end if

    call pop_sub()
  end subroutine h_sys_output_hamiltonian

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
