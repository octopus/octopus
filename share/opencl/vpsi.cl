#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void dvpsi(const int npsi, const int offset, const __global double * vv, 
		    const __global double * psi, const int ldpsi,
		    __global double * vpsi, const int ldvpsi){
  int i = get_global_id(0);

  double vi = vv[offset + i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*ldvpsi + i] = vi*psi[j*ldpsi + i];
  }

}

__kernel void zvpsi(const int npsi, const int offset, const __global double * vv, 
		    const __global double2 * psi, const int ldpsi,
		    __global double2 * vpsi, const int ldvpsi){
  int i = get_global_id(0);

  double vi = vv[offset + i];
  for(int j = 0; j < npsi; j++){
    vpsi[j*ldvpsi + i] = vi*psi[j*ldpsi + i];
  }

}

__kernel void zvpsi_spinors(const int npsi, 
			    const __global double * vv, const int ldvv,
			    const __global double2 * psi, const int ldpsi,
			    __global double2 * vpsi, const int ldvpsi){
  int i = get_global_id(0);

  double vi1 = vv[         i];
  double vi2 = vv[ldvv   + i];
  double vi3 = vv[2*ldvv + i];
  double vi4 = vv[3*ldvv + i];
  for(int j = 0; j < 2*npsi; j+=2){
    double2 psi1 = psi[      j*ldpsi + i];
    double2 psi2 = psi[(j + 1)*ldpsi + i];

    vpsi[      j*ldvpsi + i] = vi1*psi1 + (double2)(vi3*psi2.x - vi4*psi2.y, vi3*psi2.y + vi4*psi2.x);
    vpsi[(j + 1)*ldvpsi + i] = vi2*psi2 + (double2)(vi3*psi1.x + vi4*psi1.y, vi3*psi1.y - vi4*psi1.x);
  }

}

/*
 Local Variables:
 mode: c
 coding: utf-8
 End:
*/
