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

/**************************************************/

void FC_FUNC_(clamdfftsetup_low, CLAMDFFTSETUP_LOW)(int * status){
  clAmdFftSetupData setup_data;

  *status = clAmdFftInitSetupData(&setup_data);
  *status = clAmdFftSetup(&setup_data);
}

/**************************************************/

void FC_FUNC_(clamdfftteardown_low, CLAMDFFTTEARDOWN_LOW)(){
  clAmdFftTeardown();
}

/**************************************************/

void FC_FUNC_(clamdfftcreatedefaultplan_low, CLAMDFFTCREATEDEFAULTPLAN_LOW)
     (clAmdFftPlanHandle ** plHandle, cl_context * context, const int * dim, const cl_long * clLengths, int * status){
  size_t * lengths_size_t;
  int i;

  lengths_size_t = (size_t *) malloc(sizeof(size_t)*(*dim));
  for(i = 0; i < *dim; i++){
    lengths_size_t[i] = (size_t) clLengths[i];
    printf("%d %d\n", clLengths[i], lengths_size_t[i]);
  }

  *plHandle = (clAmdFftPlanHandle *) malloc(sizeof(clAmdFftPlanHandle));

  *status = clAmdFftCreateDefaultPlan(*plHandle, *context, *dim, lengths_size_t);

  free(lengths_size_t);
  
}

/**************************************************/

void FC_FUNC_(clamdfftdestroyplan_low, CLAMDFFTDESTROYPLAN_LOW)
     (clAmdFftPlanHandle ** plHandle, int * status){

  *status = clAmdFftDestroyPlan(*plHandle);
  free(*plHandle);
}

/**************************************************/

#endif
