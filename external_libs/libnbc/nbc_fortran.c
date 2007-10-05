/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#include "nbc.h"

/* fortran bindings */
#include "nbc_fortran.h"

#define NBC_HANDLE_SIZE_START 1024
#define NBC_REQUEST_NULL -1

static int handle_size;
static NBC_Handle **handles;
static int F_initialized;

int NBC_Create_fortran_handle(int *fhandle, NBC_Handle **handle) {
  NBC_Handle *newhandle;
  int i;
  
  if(!F_initialized) {
    handle_size=NBC_HANDLE_SIZE_START;
    handles = malloc(handle_size*sizeof(NBC_Handle*));
    /* initialize new handles */
    for(i=0; i<handle_size; i++) {
      handles[i] = NULL;
    }
    F_initialized = 1;
  }

  for(i=0; i<handle_size; i++) {
    if(handles[i] == NULL) break;
  }

  /* we reached the last one and did not find a free entry */
  if(handles[i] != NULL) {
    /* double number of handles */
    handle_size = handle_size*2;
    handles = realloc(handles, handle_size*sizeof(NBC_Handle*));
    /* initialize new handles */
    for(i=handle_size/2; i<handle_size; i++) {
      handles[i] = NULL;
    }
  }
  i=handle_size/2;

  newhandle = malloc(sizeof(NBC_Handle));
  handles[i] = newhandle;
  *handle = newhandle;
  *fhandle = i;

  return 0;
}

void NBC_F77_FUNC_(nbc_wait,NBC_WAIT)(int *fhandle, int *ierr) {

  if((*fhandle == -1) || (handles[*fhandle] == NULL)) return;

  *ierr = NBC_Wait(handles[*fhandle]);
  free(handles[*fhandle]);
  handles[*fhandle] = NULL;
  *fhandle = NBC_REQUEST_NULL;
}

/* TODO: this could be a problem, when a handle is actually finished and
 * we call nbc_test or wait later on it and the same index has already
 * been given to another request ... but this is the same when I use a
 * memory address as handle (like ompi does???) -> ask Jeff */
void NBC_F77_FUNC_(nbc_test,NBC_TEST)(int *fhandle, int *ierr) {

  if((*fhandle == -1) || (handles[*fhandle] == NULL)) return;

  *ierr = NBC_Progress(handles[*fhandle]);
  if(*ierr == NBC_OK) {
    free(handles[*fhandle]);
    handles[*fhandle] = NULL;
    *fhandle = NBC_REQUEST_NULL;
  }
}
