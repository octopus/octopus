/*
 Copyright (C) 2013-2020 X. Andrade, M. Lueders

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <cl_global.h>
#include <cl_complex.h>

__kernel void boundaries_periodic(const int nper, 
               __global const int * __restrict per_points,
               __global double * __restrict ff,
               const int ldff){
  const int ist  = get_global_id(0);
  const int iper = get_global_id(1);
  
  if(iper >= nper) return;
  
  const int ip_bnd = per_points[iper*2    ] - 1;
  const int ip_inn = per_points[iper*2 + 1] - 1;

  ff[(ip_bnd<<ldff) + ist] = ff[(ip_inn<<ldff) + ist];

}

__kernel void boundaries_periodic_corr(const int nper, 
               __global const int * __restrict per_points,
               __global double2 * __restrict ff,
               const int ldff,
               __global double2 * __restrict phase_correction,
               const int np){
  const int ist  = get_global_id(0);
  const int iper = get_global_id(1);
  
  if(iper >= nper) return;
  
  const int ip_bnd = per_points[iper*2    ] - 1;
  const int ip_inn = per_points[iper*2 + 1] - 1;

  ff[(ip_bnd<<ldff) + ist] = complex_mul( ff[(ip_inn<<ldff) + ist], phase_correction[ip_bnd - np]);

}

__kernel void boundaries_periodic_send(const int maxsend,
               __global const int * __restrict nsend,
               __global const int * __restrict per_send,
               __global const double * __restrict ff,
               const int ldff, 
               __global double * __restrict sendbuffer){

  const int ist   = get_global_id(0);
  const int ip    = get_global_id(1);
  const int ipart = get_global_id(2);
  
  const int np = nsend[ipart];

  if(ip >= np) return;
  
  const int ip_send = per_send[maxsend*ipart + ip] - 1;

  sendbuffer[((maxsend*ipart + ip)<<ldff) + ist] = ff[(ip_send<<ldff) + ist];
  
}

__kernel void boundaries_periodic_recv(const int maxrecv,
               __global const int * __restrict nrecv,
               __global const int * __restrict per_recv,
               const int ldper_recv,
               __global const double * __restrict recvbuffer,
               __global double * __restrict ff,
               const int ldff){

  const int ist   = get_global_id(0);
  const int ip    = get_global_id(1);
  const int ipart = get_global_id(2);
  
  const int np = nrecv[ipart];

  if(ip >= np) return;
  
  const int ip_recv = per_recv[ldper_recv*ipart + ip] - 1;

  ff[(ip_recv<<ldff) + ist] = recvbuffer[((maxrecv*ipart + ip)<<ldff) + ist];
  
}

__kernel void boundaries_periodic_recv_corr(const int maxrecv,
               __global const int * __restrict nrecv,
               __global const int * __restrict per_recv,
               const int ldper_recv,
               __global const double2 * __restrict recvbuffer,
               __global double2 * __restrict ff,
               const int ldff,
               __global double2 * __restrict phase_correction,
               const int np){

  const int ist   = get_global_id(0);
  const int ip    = get_global_id(1);
  const int ipart = get_global_id(2);
  
  const int np_local = nrecv[ipart];

  if(ip >= np_local) return;
  
  const int ip_recv = per_recv[ldper_recv*ipart + ip] - 1;

  ff[(ip_recv<<ldff) + ist] = complex_mul(recvbuffer[((maxrecv*ipart + ip)<<ldff) + ist], phase_correction[ip_recv-np]);
  
}


/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
