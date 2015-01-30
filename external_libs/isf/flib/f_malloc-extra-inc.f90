!> @file
!! Other fortran file for f_malloc routines
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  !guess the rank
  m%rank=0

  if (present(lbounds)) then
     m%rank=size(lbounds)
     m%lbounds(1:m%rank)=lbounds
  end if

  if (present(sizes)) then
     if (m%rank == 0) then
        m%rank=size(sizes)
     else if (f_err_raise(m%rank/=size(sizes),&
          'sizes not conformal with lbounds',ERR_INVALID_MALLOC)) then
        return
     end if
     m%shape(1:m%rank)=sizes
     do i=1,m%rank
        m%ubounds(i)=m%lbounds(i)+m%shape(i)-1
     end do
     if (present(ubounds)) then
        if (f_err_raise(m%rank/=size(ubounds),'sizes not conformal with ubounds')) return
        do i=1,m%rank
           if (f_err_raise(m%ubounds(i) /=ubounds(i),&
                'ubounds not conformal with sizes and lbounds',ERR_INVALID_MALLOC)) return
        end do
     end if
  else
     if (present(ubounds)) then
        if (m%rank == 0) then
           m%rank=size(ubounds)
        else if (f_err_raise(m%rank/=size(ubounds),&
             'ubounds not conformal with lbounds',ERR_INVALID_MALLOC)) then
           return
        end if
        m%ubounds(1:m%rank)=ubounds
        do i=1,m%rank
           m%shape(i)=m%ubounds(i)-m%lbounds(i)+1
        end do
     else
        call f_err_throw('at least sizes or ubounds should be defined',ERR_INVALID_MALLOC)
        return
     end if
  end if

