#include <config.h>

#ifdef HAVE_CLAMDFFT

#include <clFFT.h>

void FC_FUNC_(clamdfftgetversion_low, CLAMDFFTGETVERSION_LOW)(int * major, int * minor, int * patch, int * status){
  cl_uint cl_major, cl_minor, cl_patch;

  *status = clfftGetVersion(&cl_major, &cl_minor, &cl_patch);
  *major = cl_major;
  *minor = cl_minor;
  *patch = cl_patch;
}

/**************************************************/

void FC_FUNC_(clamdfftsetup_low, CLAMDFFTSETUP_LOW)(int * status){
  clfftSetupData setup_data;

  *status = clfftInitSetupData(&setup_data);
  *status = clfftSetup(&setup_data);
}

/**************************************************/

void FC_FUNC_(clamdfftteardown_low, CLAMDFFTTEARDOWN_LOW)(){
  clfftTeardown();
}

/**************************************************/

void FC_FUNC_(clamdfftcreatedefaultplan_low, CLAMDFFTCREATEDEFAULTPLAN_LOW)
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

void FC_FUNC_(clamdfftdestroyplan_low, CLAMDFFTDESTROYPLAN_LOW)
     (clfftPlanHandle ** plHandle, int * status){

  *status = clfftDestroyPlan(*plHandle);
  free(*plHandle);
}

/**************************************************/

void FC_FUNC_(clamdfftenqueuetransform_low, CLAMDFFTENQUEUETRANSFORM_LOW)(clfftPlanHandle ** plHandle, 
									  const int * dir, 
									  cl_command_queue * commQueues,
									  cl_mem * inputBuffers, 
									  cl_mem * outputBuffers, 
									  int * status){
  
  *status = clfftEnqueueTransform(**plHandle, *dir, 1, commQueues, 0, NULL, NULL, inputBuffers, outputBuffers, NULL);

}

/**************************************************/

void FC_FUNC_(clamdfftenqueuetransform_tmpbuf_low, CLAMDFFTENQUEUETRANSFORM_TMPBUF_LOW)(clfftPlanHandle ** plHandle, 
											const int * dir, 
											cl_command_queue * commQueues,
											cl_mem * inputBuffers, 
											cl_mem * outputBuffers, 
											cl_mem * tmpBuffer, 
											int * status){
  
  *status = clfftEnqueueTransform(**plHandle, *dir, 1, commQueues, 0, NULL, NULL, inputBuffers, outputBuffers, *tmpBuffer);

}

/**************************************************/

void FC_FUNC_(clamdfftsetplanprecision_low, CLAMDFFTSETPLANPRECISION_LOW)(clfftPlanHandle ** plHandle, 
									  const int * precision, 
									  int * status){

  *status = clfftSetPlanPrecision (**plHandle, *precision);

}


/**************************************************/


void FC_FUNC_(clamdfftsetplanbatchsize_low, CLAMDFFTSETPLANBATCHSIZE_LOW)(clfftPlanHandle ** plHandle, 
									  const cl_long * batchSize, 
									  int * status){

  *status = clfftSetPlanBatchSize(**plHandle, (size_t) *batchSize);

}

/**************************************************/

void FC_FUNC_(clamdfftsetlayout_low, CLAMDFFTSETLAYOUT_LOW)(clfftPlanHandle ** plHandle, 
							    const int * iLayout,
							    const int * oLayout,
							    int * status){

  
  *status = clfftSetLayout(**plHandle, *iLayout, *oLayout);

}

/**************************************************/

void FC_FUNC_(clamdfftsetresultlocation_low, CLAMDFFTSETRESULTLOCATION_LOW)(clfftPlanHandle ** plHandle, 
									    const int * placeness, 
									    int * status){

  *status = clfftSetResultLocation (**plHandle, *placeness);
  
}


/**************************************************/

void FC_FUNC_(clamdfftsetplaninstride_low, CLAMDFFTSETPLANINSTRIDE_LOW)(clfftPlanHandle ** plHandle, 
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

void FC_FUNC_(clamdfftsetplanoutstride_low, CLAMDFFTSETPLANOUTSTRIDE_LOW)(clfftPlanHandle ** plHandle, 
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

void FC_FUNC_(clamdfftgetplaninstride_low, CLAMDFFTGETPLANINSTRIDE_LOW)(clfftPlanHandle ** plHandle, 
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

void FC_FUNC_(clamdfftgetplanoutstride_low, CLAMDFFTGETPLANOUTSTRIDE_LOW)(clfftPlanHandle ** plHandle, 
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

void FC_FUNC_(clamdfftsetplanscale_low, CLAMDFFTSETPLANSCALE_LOW)(clfftPlanHandle ** plHandle, 
								  const int * dir, 
								  const double * scale,
								  int * status){

  *status = clfftSetPlanScale(**plHandle, *dir, (cl_float) *scale);

}

/**************************************************/

void FC_FUNC_(clamdfftgetplanscale_low, CLAMDFFTGETPLANSCALE_LOW)(clfftPlanHandle ** plHandle, 
								  const int * dir, 
								  double * scale,
								  int * status){

  cl_float fscale;

  *status = clfftGetPlanScale(**plHandle, *dir, &fscale);
  *scale = fscale;

}

/**************************************************/

void FC_FUNC_(clamdfftbakeplan_low, CLAMDFFTBAKEPLAN_LOW)(clfftPlanHandle ** plHandle,
							  cl_command_queue * commQueue,
							  int * status){

  *status = clfftBakePlan(**plHandle, 1, commQueue, NULL, NULL);

}


/**************************************************/

void FC_FUNC_(clamdfftgettmpbufsize_low, CLAMDFFTGETTMPBUFSIZE_LOW)(const clfftPlanHandle ** plHandle, 
								    cl_long * buffersize,
								    int * status){

  size_t buffersize_size_t;
  
  *status = clfftGetTmpBufSize(**plHandle, &buffersize_size_t);
  *buffersize = (cl_long) buffersize_size_t;
  
}

#endif
