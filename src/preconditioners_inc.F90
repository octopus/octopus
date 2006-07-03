subroutine X(apply_preconditioner_smoothing)(this, a, b)
  type(preconditioner_smoothing_t), intent(inout) :: this
  R_TYPE, intent(inout) :: a(:)
  R_TYPE, intent(inout) :: b(:)
  
  call X(nl_operator_operate) (this%op, a(:), b(:))
  
end subroutine X(apply_preconditioner_smoothing)


subroutine X(apply_preconditioner_smoothing_wfs)(this, a, b, n)
  type(preconditioner_smoothing_t), intent(inout) :: this
  R_TYPE, intent(inout) :: a(:,:)
  R_TYPE, intent(inout) :: b(:,:)
  integer, intent(in)   :: n

  integer :: i

  do i=1, n
    call X(nl_operator_operate) (this%op, a(:,i), b(:,i))
  end  do
  
end subroutine X(apply_preconditioner_smoothing_wfs)
