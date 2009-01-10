  ! ---------------------------------------------------------
  subroutine h_sys_output_hamiltonian(hm, m, sb, dir, outp)
    type(hamiltonian_t),   intent(in) :: hm
    type(mesh_t),          intent(in) :: m
    type(simul_box_t),     intent(in) :: sb
    character(len=*),      intent(in) :: dir
    type(h_sys_output_t),  intent(in) :: outp

    integer :: is, err
    character(len=80) :: fname
    FLOAT :: u
    FLOAT, allocatable :: v0(:,:)
    
    call push_sub('h.hamiltonian_output')

    u = units_out%energy%factor
    if(iand(outp%what, output_potential).ne.0) then
      ALLOCATE(v0(1:m%np, 1:hm%d%dim), m%np*hm%d%dim)
      v0(1:m%np, 1) = hm%ep%vpsl(1:m%np)
      call doutput_function(outp%how, dir, "v0", m, sb, v0(:, 1), u, err)
      deallocate(v0)

      if(hm%ep%classical_pot > 0) then
        call doutput_function(outp%how, dir, "vc", m, sb, hm%ep%Vclassical, u, err)
      end if

      if(hm%theory_level.ne.INDEPENDENT_PARTICLES) then
        call doutput_function(outp%how, dir, 'vh', m, sb, hm%vhartree, u, err)
        do is = 1, min(hm%d%ispin, 2)
          write(fname, '(a,i1)') 'vxc-', is
          call doutput_function(outp%how, dir, fname, m, sb, hm%vxc(:, is), u, err)

          ! finally the full KS potential (without non-local PP contributions)
          write(fname, '(a,i1)') 'vks-', is
          if (hm%ep%classical_pot > 0) then
            call doutput_function(outp%how, dir, fname, m, sb, &
              hm%ep%vpsl + hm%ep%Vclassical + hm%vhartree + hm%vxc(:, is), u, err)
          else
            call doutput_function(outp%how, dir, fname, m, sb, &
              hm%ep%vpsl + hm%vhartree + hm%vxc(:, is), u, err)
          end if
        end do
      end if

      if(hm%self_induced_magnetic) then
        select case(sb%dim)
        case(3)
          call doutput_function(outp%how, dir, 'Bind_x', m, sb, hm%b_ind(:, 1), M_ONE, err)
          call doutput_function(outp%how, dir, 'Bind_y', m, sb, hm%b_ind(:, 2), M_ONE, err)
          call doutput_function(outp%how, dir, 'Bind_z', m, sb, hm%b_ind(:, 3), M_ONE, err)
        case(2)
          call doutput_function(outp%how, dir, 'Bind_z', m, sb, hm%b_ind(:, 1), M_ONE, err)
        end select
      end if
    end if

    call pop_sub()
  end subroutine h_sys_output_hamiltonian

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
