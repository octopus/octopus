#include <config.h>

#ifdef HAVE_CLFFT

#include <clFFT.h>

void FC_FUNC_(clfftgetversion_low, CLFFTGETVERSION_LOW)(int * major, int * minor, int * patch, int * status){
  cl_uint cl_major, cl_minor, cl_patch;

  *status = clfftGetVersion(&cl_major, &cl_minor, &cl_patch);
  *major = cl_major;
  *minor = cl_minor;
  *patch = cl_patch;
}

/**************************************************/

void FC_FUNC_(clfftsetup_low, CLFFTSETUP_LOW)(int * status){
  clfftSetupData setup_data;

  *status = clfftInitSetupData(&setup_data);
  *status = clfftSetup(&setup_data);
}

/**************************************************/

void FC_FUNC_(clfftteardown_low, CLFFTTEARDOWN_LOW)(){
  clfftTeardown();
}

/**************************************************/

void FC_FUNC_(clfftcreatedefaultplan_low, CLFFTCREATEDEFAULTPLAN_LOW)
     (clfftPlanHandle ** plHandle, cl_context * context, const int * dim, const cl_long * clLengths, int * status){
  size_t * lengths_size_t;
  int i;

  lengths_size_t = (size_t *) malloc(sizeof(size_t)*(*dim));
  for(i = 0; i < *dim; i++){
    lengths_size_t[i] = (size_t) clLengths[i];
    /* printf("%d %d\n", clLengths[i], lengths_size_t[i]);*/
  }

  *plHandle = (clfftPlanHandle *) malloc(sizeof(clfftPlanHandle));

  *status = clfftCreateDefaultPlan(*plHandle, *context, *dim, lengths_size_t);

  free(lengths_size_t);
  
}

/**************************************************/

void FC_FUNC_(clfftdestroyplan_low, CLFFTDESTROYPLAN_LOW)
     (clfftPlanHandle ** plHandle, int * status){

  *status = clfftDestroyPlan(*plHandle);
  free(*plHandle);
}

/**************************************************/

void FC_FUNC_(clfftenqueuetransform_low, CLFFTENQUEUETRANSFORM_LOW)(clfftPlanHandle ** plHandle, 
									  const int * dir, 
									  cl_command_queue * commQueues,
									  cl_mem * inputBuffers, 
									  cl_mem * outputBuffers, 
									  int * status){
  
  *status = clfftEnqueueTransform(**plHandle, *dir, 1, commQueues, 0, NULL, NULL, inputBuffers, outputBuffers, NULL);

}

/**************************************************/

void FC_FUNC_(clfftenqueuetransform_tmpbuf_low, CLFFTENQUEUETRANSFORM_TMPBUF_LOW)(clfftPlanHandle ** plHandle, 
											const int * dir, 
											cl_command_queue * commQueues,
											cl_mem * inputBuffers, 
											cl_mem * outputBuffers, 
											cl_mem * tmpBuffer, 
											int * status){
  
  *status = clfftEnqueueTransform(**plHandle, *dir, 1, commQueues, 0, NULL, NULL, inputBuffers, outputBuffers, *tmpBuffer);

}

/**************************************************/

void FC_FUNC_(clfftsetplanprecision_low, CLFFTSETPLANPRECISION_LOW)(clfftPlanHandle ** plHandle, 
									  const int * precision, 
									  int * status){

  *status = clfftSetPlanPrecision (**plHandle, *precision);

}


/**************************************************/


void FC_FUNC_(clfftsetplanbatchsize_low, CLFFTSETPLANBATCHSIZE_LOW)(clfftPlanHandle ** plHandle, 
									  const cl_long * batchSize, 
									  int * status){

  *status = clfftSetPlanBatchSize(**plHandle, (size_t) *batchSize);

}

/**************************************************/

void FC_FUNC_(clfftsetlayout_low, CLFFTSETLAYOUT_LOW)(clfftPlanHandle ** plHandle, 
							    const int * iLayout,
							    const int * oLayout,
							    int * status){

  
  *status = clfftSetLayout(**plHandle, *iLayout, *oLayout);

}

/**************************************************/

void FC_FUNC_(clfftsetresultlocation_low, CLFFTSETRESULTLOCATION_LOW)(clfftPlanHandle ** plHandle, 
									    const int * placeness, 
									    int * status){

  *status = clfftSetResultLocation (**plHandle, *placeness);
  
}


/**************************************************/

void FC_FUNC_(clfftsetplaninstride_low, CLFFTSETPLANINSTRIDE_LOW)(clfftPlanHandle ** plHandle, 
									const int * dim, 
									const cl_long * clStrides, 
									int * status){
  size_t * strides;
  int i;

  strides = (size_t *) malloc(sizeof(size_t)*(*dim));

  for(i = 0; i < *dim; i++){
    strides[i] = (size_t) clStrides[i];
    /* printf("%d %d\n", clStrides[i], strides[i]);*/
  }

  *status = clfftSetPlanInStride(**plHandle, *dim, strides);
  
  free(strides);

}

/**************************************************/

void FC_FUNC_(clfftsetplanoutstride_low, CLFFTSETPLANOUTSTRIDE_LOW)(clfftPlanHandle ** plHandle, 
									const int * dim, 
									const cl_long * clStrides, 
									int * status){
  size_t * strides;
  int i;

  strides = (size_t *) malloc(sizeof(size_t)*(*dim));

  for(i = 0; i < *dim; i++){
    strides[i] = (size_t) clStrides[i];
    /* printf("%d %d\n", clStrides[i], strides[i]);*/
  }

  *status = clfftSetPlanOutStride(**plHandle, *dim, strides);
  
  free(strides);

}

/**************************************************/

void FC_FUNC_(clfftgetplaninstride_low, CLFFTGETPLANINSTRIDE_LOW)(clfftPlanHandle ** plHandle, 
									const int * dim, 
									cl_long * clStrides, 
									int * status){
  size_t * strides;
  int i;

  strides = (size_t *) malloc(sizeof(size_t)*(*dim));

  *status = clfftGetPlanInStride(**plHandle, *dim, strides);

  for(i = 0; i < *dim; i++){
    clStrides[i] = (cl_long) strides[i];
    /*    printf("%d %d\n", clStrides[i], strides[i]);*/
  }

  free(strides);

}

/**************************************************/

void FC_FUNC_(clfftgetplanoutstride_low, CLFFTGETPLANOUTSTRIDE_LOW)(clfftPlanHandle ** plHandle, 
									  const int * dim, 
									  cl_long * clStrides, 
									  int * status){
  size_t * strides;
  int i;

  strides = (size_t *) malloc(sizeof(size_t)*(*dim));

  *status = clfftGetPlanOutStride(**plHandle, *dim, strides);

  for(i = 0; i < *dim; i++){
    clStrides[i] = (cl_long) strides[i];
    /*    printf("%d %d\n", clStrides[i], strides[i]);*/
  }

  free(strides);

}

/**************************************************/

void FC_FUNC_(clfftsetplanscale_low, CLFFTSETPLANSCALE_LOW)(clfftPlanHandle ** plHandle, 
								  const int * dir, 
								  const double * scale,
								  int * status){

  *status = clfftSetPlanScale(**plHandle, *dir, (cl_float) *scale);

}

/**************************************************/

void FC_FUNC_(clfftgetplanscale_low, CLFFTGETPLANSCALE_LOW)(clfftPlanHandle ** plHandle, 
								  const int * dir, 
								  double * scale,
								  int * status){

  cl_float fscale;

  *status = clfftGetPlanScale(**plHandle, *dir, &fscale);
  *scale = fscale;

}

/**************************************************/

void FC_FUNC_(clfftbakeplan_low, CLFFTBAKEPLAN_LOW)(clfftPlanHandle ** plHandle,
							  cl_command_queue * commQueue,
							  int * status){

  *status = clfftBakePlan(**plHandle, 1, commQueue, NULL, NULL);

}


/**************************************************/

void FC_FUNC_(clfftgettmpbufsize_low, CLFFTGETTMPBUFSIZE_LOW)(const clfftPlanHandle ** plHandle, 
								    cl_long * buffersize,
								    int * status){

  size_t buffersize_size_t;
  
  *status = clfftGetTmpBufSize(**plHandle, &buffersize_size_t);
  *buffersize = (cl_long) buffersize_size_t;
  
}

#endif
