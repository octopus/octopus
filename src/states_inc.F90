! Orthonormalizes nst orbital in mesh m
subroutine R_FUNC(states_gram_schmidt)(nst, m, dim, psi)
  integer, intent(in) :: nst, dim
  type(mesh_type), intent(IN) :: m
  R_TYPE, intent(inout) :: psi(0:m%np, dim, nst)

  integer :: p, q, id
  real(r8) :: nrm2
  R_TYPE :: ss

  do p = 1, nst
    do q = 1, p - 1
      ss = R_FUNC(states_ddot)(m, dim, psi(1:m%np, :, q), psi(1:m%np, :, p))
      do id = 1, dim
        call R_FUNC(axpy) (m%np, -ss, psi(1:m%np, id, q), 1, psi(1:m%np, id, p), 1)
      end do
    enddo
    nrm2 = R_FUNC(states_ddot)(m, dim, psi(1:m%np, :, p), psi(1:m%np, :, p))
    ss = REALORCOMPLEX(1.0_r8/sqrt(nrm2))
    do id = 1, dim
      call R_FUNC(scal) (m%np, ss, psi(1:m%np, id, p), 1)
    end do
  end do

  return
end subroutine R_FUNC(states_gram_schmidt)

function R_FUNC(states_ddot)(m, dim, f1, f2)
  type(mesh_type), intent(IN) :: m
  integer, intent(in) :: dim
  R_TYPE, intent(IN) :: f1(1:m%np, dim), f2(1:m%np, dim)
  real(r8) :: R_FUNC(states_ddot)

  integer i
  real(r8) :: d

  d = 0._r8
  do i = 1, dim
    d = d + R_FUNC(mesh_dp)(m, f1, f2)
  end do

  R_FUNC(states_ddot) = d
end function R_FUNC(states_ddot)
