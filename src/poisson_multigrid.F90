!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id$

subroutine poisson_multigrid(gr, pot, rho)
  type(grid_type), intent(inout) :: gr
  FLOAT, intent(out) :: pot(:)
  FLOAT, intent(in)  :: rho(:)

  type(mesh_type), pointer :: m_fine, m_coarse
  FLOAT, allocatable :: test(:), test2(:), ltest(:), ltest2(:)
  integer :: ierr

  m_fine   => gr%mgrid%level(0)%m
  m_coarse => gr%mgrid%level(1)%m

  ! test transfer between coarse and fine grids
  allocate(test(m_coarse%np), ltest(m_coarse%np), test2(m_fine%np), ltest2(m_fine%np))
  call multigrid_fine2coarse(gr%mgrid, 1, rho, test)
  call multigrid_coarse2fine(gr%mgrid, 1, test, test2)

  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
     "lixo", "fine", m_fine, m_fine%sb, rho, M_ONE, ierr)
  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
     "lixo", "coarse", m_coarse, m_coarse%sb, test, M_ONE, ierr)
  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
     "lixo", "fine2", m_fine, m_fine%sb, test2, M_ONE, ierr)

  ! test laplacians
  call df_laplacian(gr%sb, gr%mgrid%level(0)%f_der, rho, ltest2)
  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
     "lixo", "l-fine", m_fine, m_fine%sb, ltest2, M_ONE, ierr)

  call df_laplacian(gr%sb, gr%mgrid%level(1)%f_der, test, ltest)
  call doutput_function(output_fill_how("AxisX_and_PlaneX_Gnuplot"), &
     "lixo", "l-coarse", m_coarse, m_coarse%sb, ltest, M_ONE, ierr)
  
  stop "aqui"
end subroutine poisson_multigrid
