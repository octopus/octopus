#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void svpsi(const int nn, const int npsi, const __global float * vv, const __global float * psi, __global float * vpsi){
  int i = get_global_id(0);

  float vi = vv[i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*nn + i] = vi*psi[j*nn + i];
  }
}

__kernel void dvpsi(const int nn, const int npsi, const __global double * vv, const __global double * psi, __global double * vpsi){
  int i = get_global_id(0);
  
  double vi = vv[i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*nn + i] = vi*psi[j*nn + i];
  }

}

__kernel void cvpsi(const int nn, const int npsi, const __global float * vv, const __global float2 * psi, __global float2 * vpsi){
  int i = get_global_id(0);

  float vi = vv[i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*nn + i] = vi*psi[j*nn + i];
  }

}

__kernel void zvpsi(const int nn, const int npsi, const __global double * vv, const __global double2 * psi, __global double2 * vpsi){
  int i = get_global_id(0);

  double vi = vv[i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*nn + i] = vi*psi[j*nn + i];
  }

}


/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
