! Copyright 2009
! Wenjie Tang, Andri Arnaldsson, Samuel T. Chill, and Graeme Henkelman
!
! Bader is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! A copy of the GNU General Public License is available at
! http://www.gnu.org/licenses/

!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module containing matrix functions
!-----------------------------------------------------------------------------------!
MODULE matrix_mod
  USE kind_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: matrix_mult,matrix_vector,vector_matrix
  PUBLIC :: matrix_transpose,matrix_volume
  PUBLIC :: matrix_3x3_inverse

  CONTAINS

!-----------------------------------------------------------------------------------!
! matrix_mult:  Multiply the matricies A and B and return the product C
!-----------------------------------------------------------------------------------!

  SUBROUTINE matrix_mult(A,B,C)

    REAL(q2),INTENT(IN),DIMENSION(:,:) :: A,B
    REAL(q2),INTENT(OUT),DIMENSION(:,:) :: C
    INTEGER :: m,n,j,k

    m=SIZE(A,2)
    n=SIZE(B,2)

! comment this check to optimize
!    IF(m .ne. SIZE(B,1)) THEN
!      WRITE(*,*) 'ERROR: matrix multiplication dimensions do not match'
!      STOP
!    END IF

    DO j=1,n
      C(:,j)=0._q2
      DO k=1,m
        C(:,j)=C(:,j)+B(k,j)*A(:,k)
      END DO
    END DO

  RETURN
  END SUBROUTINE matrix_mult

!-----------------------------------------------------------------------------------!
! matrix_vector:  Multiply the matrix M with the vector V and return the product Vp
!-----------------------------------------------------------------------------------!

  SUBROUTINE matrix_vector(M,V,Vp)

    REAL(q2),INTENT(IN),DIMENSION(:,:) :: M
    REAL(q2),INTENT(IN),DIMENSION(:) :: V
    REAL(q2),INTENT(OUT),DIMENSION(:) :: Vp

    INTEGER :: i,n

    n=SIZE(V)

!    IF(n .ne. SIZE(M,2)) THEN
!      WRITE(*,*) 'ERROR: matrix-vector multiplication dimensions do not match'
!      STOP
!    END IF

    Vp=0._q2
    DO i=1,n
      Vp=Vp+V(i)*M(:,i)
    END DO

  RETURN
  END SUBROUTINE matrix_vector

!-----------------------------------------------------------------------------------!
! vector_matrix:  Multiply the vector V with the matrix M and return the product Vp
!-----------------------------------------------------------------------------------!

  SUBROUTINE vector_matrix(V,M,Vp)

    REAL(q2),INTENT(IN),DIMENSION(:) :: V
    REAL(q2),INTENT(IN),DIMENSION(:,:) :: M
    REAL(q2),INTENT(OUT),DIMENSION(:) :: Vp

    INTEGER :: i,n

    n=SIZE(V)

!    IF(n .ne. SIZE(M,1)) THEN
!      WRITE(*,*) 'ERROR: matrix-vector multiplication dimensions do not match'
!      STOP
!    END IF

    Vp=0._q2
    DO i=1,n
      Vp=Vp+V(i)*M(i,:)
    END DO

  RETURN
  END SUBROUTINE vector_matrix

!-----------------------------------------------------------------------------------!
! matrix_transpose:  Set matrix B to be the transpose of A
!-----------------------------------------------------------------------------------!

  SUBROUTINE matrix_transpose(A,B)

    REAL(q2),INTENT(IN),DIMENSION(:,:) :: A
    REAL(q2),INTENT(INOUT),DIMENSION(:,:) :: B

    INTEGER :: i,j,n,m

    n=SIZE(A,1)
    m=SIZE(A,2)

    DO i=1,n
      DO j=1,m
        B(j,i)=A(i,j)
      END DO
    END DO

  RETURN
  END SUBROUTINE matrix_transpose

!-----------------------------------------------------------------------------------!
! matrix_3x3_inverse:  Set matrix B to be the inverse of A
!-----------------------------------------------------------------------------------!

  SUBROUTINE matrix_3x3_inverse(A,B)

    REAL(q2),INTENT(IN),DIMENSION(3,3) :: A
    REAL(q2),INTENT(OUT),DIMENSION(3,3) :: B
    REAL(q2) :: det
    INTEGER :: i,j,it,jt

    det=0
    DO i=1,3
      it=i-1
      DO j=1,3
        jt=j-1
        B(j,i) = & 
        &      A(mod(it+1,3)+1,mod(jt+1,3)+1)*A(mod(it+2,3)+1,mod(jt+2,3)+1)  &
        &     -A(mod(it+1,3)+1,mod(jt+2,3)+1)*A(mod(it+2,3)+1,mod(jt+1,3)+1)
      END DO
      det=det+A(i,1)*B(1,i)
    END DO

    DO i=1,3
      DO j=1,3
        B(i,j)=B(i,j)/det
      END DO
    END DO

  RETURN
  END SUBROUTINE matrix_3x3_inverse

!-----------------------------------------------------------------------------------!
! matrix_volume: Function returning the triple product of the lattice vectors.
!-----------------------------------------------------------------------------------!
    
  FUNCTION matrix_volume(h)

    REAL(q2),INTENT(IN),DIMENSION(3,3) :: h
    REAL(q2) :: matrix_volume
    
    matrix_volume = h(1,1)*(h(2,2)*h(3,3)-h(2,3)*h(3,2))  &
    &              -h(1,2)*(h(2,1)*h(3,3)-h(3,1)*h(2,3))  &
    &              +h(1,3)*(h(2,1)*h(3,2)-h(3,1)*h(2,2))

    matrix_volume = abs(matrix_volume)

  RETURN
  END FUNCTION matrix_volume

!-----------------------------------------------------------------------------------!

END MODULE matrix_mod

