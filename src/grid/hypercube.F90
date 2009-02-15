!! Copyright (C) 2009 N. Helbig
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
!! $Id: hypercube.F90 4907 2009-01-31 11:21:32Z xavier $

subroutine x_to_i(ndim,nr,enlarge,boxdim,coord,icoord)

implicit none

integer, intent(in):: ndim, enlarge
integer, intent(in):: coord(1:ndim), boxdim(1:ndim+1)
integer, intent(in):: nr(1:2,1:ndim)

integer, intent(out):: icoord

integer:: boxnumb
integer, allocatable:: border(:),npoints(:), lowerb(:), tempcoord(:)
integer:: ii,jj

allocate(border(1:ndim), npoints(1:ndim), lowerb(1:ndim), tempcoord(1:ndim))

do ii=1, ndim
  npoints(ii)=nr(2,ii)-nr(1,ii)+1
  border(ii)=nr(1,ii)+2*enlarge
  lowerb(ii)=nr(1,ii)
enddo

!move coordinates such that inner box is in the upper right corner
do ii=1, ndim
  tempcoord(ii)=coord(ii)
  tempcoord(ii)=tempcoord(ii)+enlarge-nr(1,ii)
  tempcoord(ii)=mod(tempcoord(ii),npoints(ii))
  tempcoord(ii)=tempcoord(ii)+nr(1,ii)  
enddo

!determine which box we are in
boxnumb=1
do ii=1, ndim
  if(tempcoord(ii)<border(ii)) then
    boxnumb=ii+1
    exit
  endif
enddo 

!transform coordinates
icoord=0

if(boxnumb==1) then
  npoints=npoints-2*enlarge
  do ii=ndim, 1, -1
    icoord=icoord*npoints(ii)
    icoord=icoord+tempcoord(ii)-border(ii)
  enddo
  icoord=icoord+1
  if(icoord>boxdim(1).or.icoord<1) then
    write(*,*) "Error, box point outside box"
    stop
  endif
else
  do jj=1, boxnumb-2
    npoints(jj)=npoints(jj)-2*enlarge
    lowerb(jj)=nr(1,jj)+2*enlarge
  enddo
  npoints(boxnumb-1)=2*enlarge
  do jj=ndim, 1, -1
    icoord=icoord*npoints(jj)
    icoord=icoord+(tempcoord(jj)-lowerb(jj))
  enddo
  icoord=icoord+1    
  if(icoord>boxdim(boxnumb).or.icoord<1) then
    write(*,*) "Error, box point outside box"
    stop
  else
    do jj=1, boxnumb-1
      icoord=icoord+boxdim(jj)
    enddo
  endif
endif  

deallocate(border, npoints)

end subroutine



subroutine i_to_x(ndim,nr,enlarge,boxdim,coord,icoord)

implicit none

integer, intent(in):: ndim, enlarge
integer, intent(in):: icoord, boxdim(1:ndim+1)
integer, intent(in):: nr(1:2,1:ndim)

integer, intent(out):: coord(1:ndim)

integer:: boxnumb, tempcoord
integer, allocatable:: border(:),npoints(:), lowerb(:)
integer:: ii,jj

allocate(border(1:ndim), npoints(1:ndim), lowerb(1:ndim))

boxnumb=1

jj=0

do ii=1, ndim
  jj=jj+boxdim(ii)  
  if(icoord>jj) then
    boxnumb=ii+1
  endif
enddo
  
do ii=1, ndim
  npoints(ii)=nr(2,ii)-nr(1,ii)+1
  border(ii)=nr(1,ii)+2*enlarge
  lowerb(ii)=nr(1,ii)
enddo

tempcoord=icoord

if(boxnumb==1) then
  npoints=npoints-2*enlarge
  tempcoord=tempcoord-1
  do ii=1, ndim
    coord(ii)=mod(tempcoord,npoints(ii))
    tempcoord=tempcoord-coord(ii)
    tempcoord=tempcoord/npoints(ii)
    coord(ii)=coord(ii)+border(ii)
  enddo   
else
  do ii=1, boxnumb-2
    npoints(ii)=npoints(ii)-2*enlarge
    lowerb(ii)=nr(1,ii)+2*enlarge
    tempcoord=tempcoord-boxdim(ii)
  enddo
  npoints(boxnumb-1)=2*enlarge
  tempcoord=tempcoord-boxdim(boxnumb-1)
  tempcoord=tempcoord-1
  do ii=1, ndim
    coord(ii)=mod(tempcoord,npoints(ii))
    tempcoord=tempcoord-coord(ii)
    tempcoord=tempcoord/npoints(ii)
    coord(ii)=coord(ii)+lowerb(ii)
  enddo
endif

do ii=1, ndim
  npoints(ii)=nr(2,ii)-nr(1,ii)+1
enddo
    
!move inner box back to the middle

do ii=1, ndim
  coord(ii)=coord(ii)-nr(1,ii)-enlarge
  if(coord(ii)<0) then
    coord(ii)=coord(ii)+npoints(ii)
  endif
  coord(ii)=coord(ii)+nr(1,ii)
enddo

deallocate(border, npoints)

end subroutine
