/*
 Copyright (C) 2010 X. Andrade

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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id: projector.cl 2146 2006-05-23 17:36:00Z xavier $
*/

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void dprojector_gather(const int npoints,
				__global const int * map,
				__global const double * psi, const int ldpsi,
				__global double * lpsi, const int ldlpsi){

  int ip = get_global_id(0);
  int ist = get_global_id(1);
  
  if(ip < npoints) lpsi[ldlpsi*ist + ip] = psi[ldpsi*(map[ip] - 1) + ist];
}

__kernel void zprojector_gather(const int npoints,
				__global const int * map,
				__global const double2 * psi, const int ldpsi,
				__global double2 * lpsi, const int ldlpsi){

  int ip = get_global_id(0);
  int ist = get_global_id(1);
  
  if(ip < npoints) lpsi[ldlpsi*ist + ip] = psi[ldpsi*(map[ip] - 1) + ist];
}

__kernel void dprojector_scatter(const int npoints,
				 __global const int * map,
				 __global const double * lpsi, const int ldlpsi,
				 __global double * psi, const int ldpsi){
				

  int ip = get_global_id(0);
  int ist = get_global_id(1);
  
  if(ip < npoints) psi[ldpsi*(map[ip] - 1) + ist] += lpsi[ldlpsi*ist + ip];
}

__kernel void zprojector_scatter(const int npoints,
				 __global const int * map,
				 __global const double2 * lpsi, const int ldlpsi,
				 __global double2 * psi, const int ldpsi){
				

  int ip = get_global_id(0);
  int ist = get_global_id(1);
  
  if(ip < npoints) psi[ldpsi*(map[ip] - 1) + ist] += lpsi[ldlpsi*ist + ip];
}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/

