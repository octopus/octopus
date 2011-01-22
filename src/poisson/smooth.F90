!!***h* BigDFT/restrict
!! NAME
!!   restrict
!!
!! FUNCTION
!!   An wavelet analysis where only the scaling function coefficients are calculated
!!   The input array x is not overwritten
!!
!! SOURCE
!!

subroutine restrict(n01,n02,n03,n1,n2,n3,x,y)
  implicit none
  !integer, parameter :: order=8,omin=order/2,omax=order-1-omin
  integer, parameter :: order=8,omin=4,omax=3
  !Arguments
  integer, intent(in) :: n01,n02,n03,n1,n2,n3
  real(kind=8), intent(in) :: x(-1:n01-2,-1:n02-2,-1:n03-2)
  !real(kind=8), intent(out) :: y(-4:n1+3,-4:n2+3,-4:n3+3)
  real(kind=8), intent(out) :: y(-omin:n1+omax,-omin:n2+omax,-omin:n3+omax)
  !Local variables
  real(kind=8), allocatable :: w1(:,:,:) , w2(:,:,:)
  real(kind=8) :: t1,t2
  integer :: i1,i2,i3,nt
  !Dimensions n01 should be >= 2*n1+1
  if (n01 < 2*n1+1) then
     print "(1x,a,i0,a,i0,a)","restrict: n01 (",n01,") < 2*n1+1 (",2*n1+1,")"
     stop
  end if
  if (n02 < 2*n2+1) then
     print "(1x,a,i0,a,i0,a)","restrict: n02 (",n02,") < 2*n2+1 (",2*n2+1,")"
     stop
  end if
  if (n03 < 2*n3+1) then
     print "(1x,a,i0,a,i0,a)","restrict: n03 (",n03,") < 2*n3+1 (",2*n3+1,")"
     stop
  end if
  !Allocations
  allocate(w1(-1:2*n2-1,-1:2*n3-1,-omin:n1+omax))
  allocate(w2(-1:2*n3-1,-omin:n1+omax,-omin:n2+omax))
  t1=0.d0
  do i3=-1,2*n3-1
    do i2=-1,2*n2-1
      do i1=-1,2*n1-1
        t1=t1+x(i1,i2,i3)
      enddo
    enddo
  enddo
        
! i1,i2,i3 -> i2,i3,I1
  nt=(2*n2+1)*(2*n3+1)
  call  restrict_rot_grow(n1,nt,x,w1)
! i2,i3,I1 -> i3,I1,I2
  nt=(2*n3+1)*(n1+order)
  call  restrict_rot_grow(n2,nt,w1,w2)
! i3,I1,I2  -> I1,I2,I3
  nt=(n1+order)*(n2+order)
  call  restrict_rot_grow(n3,nt,w2,y)
  
  t2=0.d0
  do i3=-omin,n3+omax
    do i2=-omin,n2+omax
      do i1=-omin,n1+omax
        t2=t2+y(i1,i2,i3)
      enddo
    enddo
  enddo

  if (abs((t1-8*t2)/t1).gt.1.d-10) then
     write(*,*) 'WARNING Restrict Monopoles',t1,t2*8
  end if

  deallocate (w1,w2)

end subroutine restrict
!!***


!!***h* BigDFT/restrict_rot_grow
!! NAME
!!   restrict_rot_grow
!!
!! RESTRICTIONS
!!   8-th order lifted interpolating wavelets
!!
!! SOURCE
!!
subroutine restrict_rot_grow(n,nt,x,y)
  implicit none
  !integer, parameter :: order=8,omin=order/2,omax=order-1-omin
  integer, parameter :: order=8,omin=4,omax=3
  !Arguments
  integer, intent(in) :: n,nt
  real(kind=8), intent(in) :: x(-1:2*n-1,nt)
  real(kind=8), intent(out) :: y(nt,-omin:n+omax)
  !Local variables  
  real(kind=8) ::  cht(-order:order)
  real(kind=8) :: tt
  integer :: i,j,l
!  8-th order lifted interpolating wavelets
  cht(0)=2871.d0/4096.d0     
  cht(1)=1.d0/4.d0            
  cht(2)=-245.d0/2048.d0     
  cht(3)=0.d0             
  cht(4)=49.d0/2048.d0       
  cht(5)=0.d0              
  cht(6)=-11.d0/2048.d0       
  cht(7)=0.d0              
  cht(8)=5.d0/8192.d0        
  do i=1,order
     cht(-i)=cht(i)
  enddo
  
  do j=1,nt
    do i=-omin,n+omax
      tt=0.d0
      do l=max(-order,-2*i-1),min(order,2*n-2*i-1)
        tt=tt+cht(l)*x(2*i+l,j)
        !write(*,*) 'i,l,2*i+l',i,l,2*i+l
      enddo
      y(j,i)=tt
    end do
  end do
end subroutine restrict_rot_grow
!!***


!!***h* BigDFT/extrapolate
!! NAME
!!   extrapolate
!!
!! FUNCTIONS
!!   An wavelet analysis where only the scaling function coefficients are calculated
!!   The input array x is not overwritten
!!
!! RESTRICTIONS
!!   8-th order lifted interpolating wavelets
!!
!! SOURCE
!!
subroutine extrapolate(n01,n02,n03,n1,n2,n3,x,y)
  implicit none
  !integer, parameter :: order=8,omin=order/2,omax=order-1-o2
  integer, parameter :: order=8,omin=4,omax=3
  !Arguments
  integer, intent(in) :: n01,n02,n03,n1,n2,n3
!  real(kind=8), intent(in) ::  x(-4:n1+3,-4:n2+3,-4:n3+3)
  real(kind=8), intent(in) ::  x(-omin:n1+omax,-omin:n2+omax,-omin:n3+omax)
  real(kind=8), intent(out) :: y(-1:n01-2,-1:n02-2,-1:n03-2)
  !Local variables
  real(kind=8), allocatable :: w1(:) , w2(:)
  integer :: nt
  !Allocations
  allocate(w1((n2+order)*(n3+order)*(2*n1+1)))
  allocate(w2((n3+order)*(2*n1+1)*(2*n2+1)))
  
! i1,i2,i3 -> i2,i3,I1
  nt=(n2+order)*(n3+order)
  call  extrapolate_rot_shrink(n1,nt,x,w1)
! i2,i3,I1 -> i3,I1,I2
  nt=(n3+order)*(2*n1+1)
  call  extrapolate_rot_shrink(n2,nt,w1,w2)
! i3,I1,I2  -> I1,I2,I3
  nt=(2*n1+1)*(2*n2+1)
  call  extrapolate_rot_shrink(n3,nt,w2,y)
  
  deallocate (w1,w2)
  
end subroutine extrapolate
!!***


!!***h* BigDFT/extrapolate_rot_shrink
!! NAME
!!   extrapolate_rot_shrink
!!
!! RESTRICTIONS
!!   8-th order lifted interpolating wavelets
!!
!! SOURCE
!!
subroutine extrapolate_rot_shrink(n,nt,x,y)
  implicit none
  !integer, parameter :: order=8,omin=order/2,omax=order-1-omin
  integer, parameter :: order=8,omin=4,omax=3
  !Arguments
  integer, intent(in) :: n,nt
  real(kind=8), intent(in) :: x(-omin:n+omax,nt)
  real(kind=8), intent(out) :: y(nt,-1:2*n-1)
  !Local variables
  real(kind=8) :: ch1,ch2,ch3,ch4
  integer :: i,j
!  8-th order lifted interpolating wavelets
  ch1=1225.d0/2048.d0
  ch2=-245.d0/2048.d0
  ch3=49.d0/2048.d0
  ch4=-5.d0/2048.d0

  do j=1,nt
     i=-1
     y(j,2*i+1)=ch4*(x(i-3,j)+x(i+4,j))+ch3*(x(i-2,j)+x(i+3,j))+ &
                ch2*(x(i-1,j)+x(i+2,j))+ch1*(x(i  ,j)+x(i+1,j))
     do i=0,n-1
        y(j,2*i+0)=x(i,j)
        y(j,2*i+1)=ch4*(x(i-3,j)+x(i+4,j))+ch3*(x(i-2,j)+x(i+3,j))+ &
                   ch2*(x(i-1,j)+x(i+2,j))+ch1*(x(i  ,j)+x(i+1,j))
     enddo
  end do
end subroutine extrapolate_rot_shrink
!!***
