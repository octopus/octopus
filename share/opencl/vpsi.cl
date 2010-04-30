#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void dvpsi(const int nn, const int npsi, const int offset, const __global double * vv, const __global double * psi, __global double * vpsi){
  int i = get_global_id(0);

  double vi = vv[offset + i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*nn + i] = vi*psi[j*nn + i];
  }

}

__kernel void zvpsi(const int nn, const int npsi, const int offset, const __global double * vv, const __global double2 * psi, __global double2 * vpsi){
  int i = get_global_id(0);

  double vi = vv[offset + i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*nn + i] = vi*psi[j*nn + i];
  }

}

__kernel void zvpsi_spinors(const int nn, const int npsi, const __global double * vv, const __global double2 * psi,  __global double2 * vpsi){
  int i = get_global_id(0);

  double vi1 = vv[       i];
  double vi2 = vv[nn   + i];
  double vi3 = vv[2*nn + i];
  double vi4 = vv[3*nn + i];
  for(int j = 0; j < 2*npsi; j+=2){
    double2 psi1 = psi[      j*nn + i];
    double2 psi2 = psi[(j + 1)*nn + i];

    vpsi[      j*nn + i] = vi1*psi1 + (double2)(vi3*psi2.x - vi4*psi2.y, vi3*psi2.y + vi4*psi2.x);
    vpsi[(j + 1)*nn + i] = vi2*psi2 + (double2)(vi3*psi1.x + vi4*psi1.y, vi3*psi1.y - vi4*psi1.x);
  }

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
