#include <config.h>

#ifdef HAVE_CLAMDFFT

#include <clAmdFft.h>

void FC_FUNC_(clamdfftgetversion_low, CLAMDFFTGETVERSION_LOW)(int * major, int * minor, int * patch, int * status){
  cl_uint cl_major, cl_minor, cl_patch;

  *status = clAmdFftGetVersion(&cl_major, &cl_minor, &cl_patch);
  *major = cl_major;
  *minor = cl_minor;
  *patch = cl_patch;
}

void FC_FUNC_(clamdfftsetup_low, CLAMDFFTSETUP_LOW)(int * status){
  clAmdFftSetupData setup_data;

  *status = clAmdFftSetup(&setup_data);
}

void FC_FUNC_(clamdfftteardown_low, CLAMDFFTTEARDOWN_LOW)(){
  clAmdFftTeardown();
}

#endif
