/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#include "nbc.h"

/* only used in this file */
static __inline__ int NBC_Start_round(NBC_Handle *handle);

/* the keyval (global) */
static int gkeyval=MPI_KEYVAL_INVALID; 

static int NBC_Key_copy(MPI_Comm oldcomm, int keyval, void *extra_state, void *attribute_val_in, void *attribute_val_out, int *flag) {
  /* delete the attribute in the new comm  - it will be created at the
   * first usage */
  *flag = 0;

  return MPI_SUCCESS;
}

static int NBC_Key_delete(MPI_Comm comm, int keyval, void *attribute_val, void *extra_state) {
  if(keyval == gkeyval) {
    free(attribute_val);
  } else {
    printf("Got wrong keyval!(%i)\n", keyval); 
  }

  return MPI_SUCCESS;
}

/* allocates a new schedule array */
int NBC_Sched_create(NBC_Schedule* schedule) {
  
  *schedule=malloc(2*sizeof(int));
  if(*schedule == NULL) { return NBC_OOR; }
  *(int*)*schedule=2*sizeof(int);
  *(((int*)*schedule)+1)=0;

  return NBC_OK;
}

/* frees a schedule array */
int NBC_Free(NBC_Handle* handle) {
  
  free(*handle->schedule);
  /* the schedule pointer itself is also malloc'd */
  free(handle->schedule);
  /* if the nbc_I<collective> attached some data */
  if(NULL != handle->tmpbuf) free(handle->tmpbuf);

  return NBC_OK;
}

/* this function puts a send into the schedule */
int NBC_Sched_send(void* buf, int count, MPI_Datatype datatype, int dest, NBC_Schedule *schedule) {
  int size;
  NBC_Args_send* send_args;
  
  /* get size of actual schedule */
  NBC_GET_SIZE(*schedule, size);
  /*printf("schedule is %i bytes\n", size);*/
  *schedule = realloc(*schedule, size+sizeof(NBC_Args_send)+sizeof(NBC_Fn_type));
  if(*schedule == NULL) { printf("Error in realloc()\n"); return NBC_OOR; }
  
  /* adjust the function type */
  *(NBC_Fn_type*)((char*)*schedule+size)=SEND;
  
  /* store the passed arguments */
  send_args = (NBC_Args_send*)((char*)*schedule+size+sizeof(NBC_Fn_type));
  send_args->buf=buf;
  send_args->count=count;
  send_args->datatype=datatype;
  send_args->dest=dest;

  /* increase number of elements in schedule */
  NBC_INC_NUM_ROUND(*schedule);
  DEBUG(10, "adding send - ends at byte %i\n", (int)(size+sizeof(NBC_Args_send)+sizeof(NBC_Fn_type)));

  /* increase size of schedule */
  NBC_INC_SIZE(*schedule, sizeof(NBC_Args_send)+sizeof(NBC_Fn_type));

  return NBC_OK;
}

/* this function puts a receive into the schedule */
int NBC_Sched_recv(void* buf, int count, MPI_Datatype datatype, int source, NBC_Schedule *schedule) {
  int size;
  NBC_Args_recv* recv_args;
  
  /* get size of actual schedule */
  NBC_GET_SIZE(*schedule, size);
  /*printf("schedule is %i bytes\n", size);*/
  *schedule = realloc(*schedule, size+sizeof(NBC_Args_recv)+sizeof(NBC_Fn_type));
  if(*schedule == NULL) { printf("Error in realloc()\n"); return NBC_OOR; }
  
  /* adjust the function type */
  *(NBC_Fn_type*)((char*)*schedule+size)=RECV;

  /* store the passed arguments */
  recv_args=(NBC_Args_recv*)((char*)*schedule+size+sizeof(NBC_Fn_type));
  recv_args->buf=buf;
  recv_args->count=count;
  recv_args->datatype=datatype;
  recv_args->source=source;

  /* increase number of elements in schedule */
  NBC_INC_NUM_ROUND(*schedule);
  DEBUG(10, "adding receive - ends at byte %i\n", (int)(size+sizeof(NBC_Args_recv)+sizeof(NBC_Fn_type)));

  /* increase size of schedule */
  NBC_INC_SIZE(*schedule, sizeof(NBC_Args_recv)+sizeof(NBC_Fn_type));

  return NBC_OK;
}

/* this function puts an operation into the schedule */
int NBC_Sched_op(void *buf3, void* buf1, void* buf2, int count, MPI_Datatype datatype, MPI_Op op, NBC_Schedule *schedule) {
  int size;
  NBC_Args_op* op_args;
  
  /* get size of actual schedule */
  NBC_GET_SIZE(*schedule, size);
  /*printf("schedule is %i bytes\n", size);*/
  *schedule = realloc(*schedule, size+sizeof(NBC_Args_op)+sizeof(NBC_Fn_type));
  if(*schedule == NULL) { printf("Error in realloc()\n"); return NBC_OOR; }
  
  /* adjust the function type */
  *(NBC_Fn_type*)((char*)*schedule+size)=OP;

  /* store the passed arguments */
  op_args=(NBC_Args_op*)((char*)*schedule+size+sizeof(NBC_Fn_type));
  op_args->buf1=buf1;
  op_args->buf2=buf2;
  op_args->buf3=buf3;
  op_args->count=count;
  op_args->op=op;
  op_args->datatype=datatype;

  /* increase number of elements in schedule */
  NBC_INC_NUM_ROUND(*schedule);
  DEBUG(10, "adding op - ends at byte %i\n", (int)(size+sizeof(NBC_Args_op)+sizeof(NBC_Fn_type)));

  /* increase size of schedule */
  NBC_INC_SIZE(*schedule, sizeof(NBC_Args_op)+sizeof(NBC_Fn_type));
  
  return NBC_OK;
}

/* this function puts a copy into the schedule */
int NBC_Sched_copy(void *src, int srccount, MPI_Datatype srctype, void *tgt, int tgtcount, MPI_Datatype tgttype, NBC_Schedule *schedule) {
  int size;
  NBC_Args_copy* copy_args;
  
  /* get size of actual schedule */
  NBC_GET_SIZE(*schedule, size);
  /*printf("schedule is %i bytes\n", size);*/
  *schedule = realloc(*schedule, size+sizeof(NBC_Args_copy)+sizeof(NBC_Fn_type));
  if(*schedule == NULL) { printf("Error in realloc()\n"); return NBC_OOR; }
  
  /* adjust the function type */
  *(NBC_Fn_type*)((char*)*schedule+size)=COPY;
  
  /* store the passed arguments */
  copy_args = (NBC_Args_copy*)((char*)*schedule+size+sizeof(NBC_Fn_type));
  copy_args->src=src;
  copy_args->srccount=srccount;
  copy_args->srctype=srctype;
  copy_args->tgt=tgt;
  copy_args->tgtcount=tgtcount;
  copy_args->tgttype=tgttype;

  /* increase number of elements in schedule */
  NBC_INC_NUM_ROUND(*schedule);
  DEBUG(10, "adding copy - ends at byte %i\n", (int)(size+sizeof(NBC_Args_copy)+sizeof(NBC_Fn_type)));

  /* increase size of schedule */
  NBC_INC_SIZE(*schedule, sizeof(NBC_Args_copy)+sizeof(NBC_Fn_type));

  return NBC_OK;
}

/* this function ends a round of a schedule */
int NBC_Sched_barrier(NBC_Schedule *schedule) {
  int size;
  
  /* get size of actual schedule */
  NBC_GET_SIZE(*schedule, size);
  /*printf("round terminated at %i bytes\n", size);*/
  *schedule = realloc(*schedule, size+sizeof(char)+sizeof(int));
  if(*schedule == NULL) { printf("Error in realloc()\n"); return NBC_OOR; }
  
  /* add the barrier char (1) because another round follows */
  *(char*)((char*)*schedule+size)=1;
  
  /* set round count elements = 0 for new round */
  *(int*)((char*)*schedule+size+sizeof(char))=0;
  DEBUG(10, "ending round at byte %i\n", (int)(size+sizeof(char)+sizeof(int)));
  
  /* increase size of schedule */
  NBC_INC_SIZE(*schedule, sizeof(char)+sizeof(int));

  return NBC_OK;
}

/* this function ends a schedule */
int NBC_Sched_commit(NBC_Schedule *schedule) {
  int size;
 
  /* get size of actual schedule */
  NBC_GET_SIZE(*schedule, size);
  /*printf("schedule terminated at %i bytes\n", size);*/
  *schedule = realloc(*schedule, size+sizeof(char));
  if(*schedule == NULL) { printf("Error in realloc()\n"); return NBC_OOR; }
 
  /* add the barrier char (0) because this is the last round */
  *(char*)((char*)*schedule+size)=0;
  DEBUG(10, "closing schedule %p at byte %i\n", *schedule, (int)(size+sizeof(char)));

  /* increase size of schedule */
  NBC_INC_SIZE(*schedule, sizeof(char));
 
  return NBC_OK;
}

int NBC_Progress(NBC_Handle *handle) {
  int flag, res;
  long size;
  char *delim;

  if((handle->req_count > 0) && (handle->req_array != NULL)) {
    DEBUG(50, "NBC_Progress: testing for %i requests\n", handle->req_count);
    res = MPI_Testall(handle->req_count, handle->req_array, &flag, MPI_STATUS_IGNORE);
    if(res != MPI_SUCCESS) { printf("MPI Error in MPI_Testall() (%i)\n", res); return res; }
  } else { 
    flag = 1; /* we had no open requests -> proceed to next round */
  }

  /* a round is finished */
  if(flag) {
    /* adjust delim to start of current round */
    DEBUG(10, "NBC_Progress: going in schedule %p to row-offset: %li\n", *handle->schedule, handle->row_offset);
    delim = (char*)*handle->schedule + handle->row_offset;
    DEBUG(10, "delim: %p\n", delim);
    NBC_GET_ROUND_SIZE(delim, size);
    DEBUG(10, "size: %li\n", size);
    /* adjust delim to end of current round -> delimiter */
    delim = delim + size;

    if(handle->req_array != NULL) {
      /* free request array */
      free(handle->req_array);
      handle->req_array = NULL;
    }
    handle->req_count = 0;

    if(*delim == 0) {
      /* this was the last round - we're done */
      DEBUG(5, "NBC_Progress last round finished - we're done\n");
      return NBC_OK;
    } else {
      DEBUG(5, "NBC_Progress round finished - goto next round\n");
      /* move delim to start of next round */
      delim = delim+1;
      /* initializing handle for new virgin round */
      handle->row_offset = (long)delim - (long)*handle->schedule;
      /* kick it off */
      res = NBC_Start_round(handle);
      if(NBC_OK != res) { printf("Error in NBC_Start_round() (%i)\n", res); return res; }
    }
  }
 
  return NBC_CONTINUE;
}

int NBC_Progress_block(NBC_Handle *handle) {
  int res;
  long size;
  char *delim;

  do {
    if((handle->req_count > 0) && (handle->req_array != NULL)) {
      DEBUG(50, "NBC_Progress_block: waiting for %i requests\n", handle->req_count);
      res = MPI_Waitall(handle->req_count, handle->req_array, MPI_STATUS_IGNORE);
      if(res != MPI_SUCCESS) { printf("MPI Error in MPI_Waitall() (%i)\n", res); return res; }
    }

    /* a round is finished - adjust delim to start of current round */
    DEBUG(10, "NBC_Progress_block: going in schedule %p to row-offset: %li\n", *handle->schedule, handle->row_offset);
    delim = (char*)*handle->schedule + handle->row_offset;
    DEBUG(10, "delim: %p\n", delim);
    NBC_GET_ROUND_SIZE(delim, size);
    DEBUG(10, "size: %li\n", size);
    /* adjust delim to end of current round -> delimiter */
    delim = delim + size;
    
    if(handle->req_array != NULL) {
      /* free request array */
      free(handle->req_array);
      handle->req_array = NULL;
    }
    handle->req_count = 0;

    if(*delim == 0) {
      /* this was the last round - we're done */
      DEBUG(5, "NBC_Progress_block: last round finished - we're done\n");
      break;
    } else {
      DEBUG(5, "NBC_Progress_block: round finished - goto next round\n");
      /* move delim to start of next round */
      delim = delim+1;
      /* initializing handle for new virgin round */
      handle->row_offset = (long)delim - (long)*handle->schedule;
      /* kick it off */
      res = NBC_Start_round(handle);
      if(NBC_OK != res) { printf("Error in NBC_Start_round() (%i)\n", res); return res; }
    }
  } while(1);
 
  return NBC_OK;
}

static __inline__ int NBC_Start_round(NBC_Handle *handle) {
  int *numptr; /* number of operations */
  int i, res;
  NBC_Fn_type *typeptr;
  NBC_Args_send *sendargs; 
  NBC_Args_recv *recvargs; 
  NBC_Args_op *opargs; 
  NBC_Args_copy *copyargs; 
  NBC_Schedule myschedule;

  /* get schedule address */
  myschedule = (NBC_Schedule*)((char*)*handle->schedule + handle->row_offset);

  numptr = (int*)myschedule;
  DEBUG(10, "start_round round at address %p : posting %i operations\n", myschedule, *numptr);

  /* typeptr is increased by sizeof(int) bytes to point to type */
  typeptr = (NBC_Fn_type*)(numptr+1);
  for (i=0; i<*numptr; i++) {
    /* go sizeof op-data forward */
    switch(*typeptr) {
      case SEND:
        DEBUG(5,"  SEND (offset %li) ", (long)typeptr-(long)myschedule);
        sendargs = (NBC_Args_send*)(typeptr+1);
        DEBUG(5,"*buf: %p, count: %i, type: %lu, dest: %i, tag: %i)\n", sendargs->buf, sendargs->count, (unsigned long)sendargs->datatype, sendargs->dest, handle->tag);
        typeptr = (NBC_Fn_type*)(((NBC_Args_send*)typeptr)+1);
        /* get an additional request - TODO: req_count NOT thread safe */
        handle->req_count++;
        handle->req_array = realloc(handle->req_array, (handle->req_count)*sizeof(MPI_Request));
        CHECK_NULL(handle->req_array);
        res = MPI_Isend(sendargs->buf, sendargs->count, sendargs->datatype, sendargs->dest, handle->tag, handle->mycomm, handle->req_array+handle->req_count-1);
        if(MPI_SUCCESS != res) { printf("Error in MPI_Isend() (%i)\n", res); return res; }
        break;
      case RECV:
        DEBUG(5, "  RECV (offset %li) ", (long)typeptr-(long)myschedule);
        recvargs = (NBC_Args_recv*)(typeptr+1);
        DEBUG(5, "*buf: %p, count: %i, type: %lu, source: %i, tag: %i)\n", recvargs->buf, recvargs->count, (unsigned long)recvargs->datatype, recvargs->source, handle->tag);
        typeptr = (NBC_Fn_type*)(((NBC_Args_recv*)typeptr)+1);
        /* get an additional request - TODO: req_count NOT thread safe */
        handle->req_count++;
        handle->req_array = realloc(handle->req_array, (handle->req_count)*sizeof(MPI_Request));
        CHECK_NULL(handle->req_array);
        res = MPI_Irecv(recvargs->buf, recvargs->count, recvargs->datatype, recvargs->source, handle->tag, handle->mycomm, handle->req_array+handle->req_count-1);
        if(MPI_SUCCESS != res) { printf("Error in MPI_Irecv() (%i)\n", res); return res; }
        break;
      case OP:
        DEBUG(5, "  OP   (offset %li) ", (long)typeptr-(long)myschedule);
        opargs = (NBC_Args_op*)(typeptr+1);
        DEBUG(5, "*buf1: %lu, buf2: %lu, count: %i, type: %lu)\n", (unsigned long)opargs->buf1, (unsigned long)opargs->buf2, opargs->count, (unsigned long)opargs->datatype);
        typeptr = (NBC_Fn_type*)((NBC_Args_op*)typeptr+1);
        res = NBC_Operation(opargs->buf3, opargs->buf1, opargs->buf2, opargs->op, opargs->datatype, opargs->count);
        if(res != NBC_OK) { printf("NBC_Operation() failed (code: %i)\n", res); return res; }
        break;
      case COPY:
        DEBUG(5, "  COPY   (offset %li) ", (long)typeptr-(long)myschedule);
        copyargs = (NBC_Args_copy*)(typeptr+1);
        DEBUG(5, "*src: %lu, srccount: %i, srctype: %lu, *tgt: %lu, tgtcount: %i, tgttype: %lu)\n", (unsigned long)copyargs->src, copyargs->srccount, (unsigned long)copyargs->srctype, (unsigned long)copyargs->tgt, copyargs->tgtcount, (unsigned long)copyargs->tgttype);
        typeptr = (NBC_Fn_type*)((NBC_Args_copy*)typeptr+1);
        res = NBC_Copy(copyargs->src, copyargs->srccount, copyargs->srctype, copyargs->tgt, copyargs->tgtcount, copyargs->tgttype, handle->mycomm);
        if(res != NBC_OK) { printf("NBC_Copy() failed (code: %i)\n", res); return res; }
        break;
      default:
        printf("NBC_Start_round: bad type %li at offset %li\n", (long)*typeptr, (long)typeptr-(long)myschedule);
        return NBC_BAD_SCHED;
    }
    /* increase ptr by size of fn_type enum */
    typeptr = (NBC_Fn_type*)((NBC_Fn_type*)typeptr+1);
  }

  /* check if we can make progress */
  res = NBC_Progress(handle);
  if((NBC_OK != res) && (NBC_CONTINUE != res)) { printf("Error in NBC_Progress() (%i)\n", res); return res; }
  
  return NBC_OK;
}

int NBC_Start(NBC_Handle *handle, MPI_Comm comm, NBC_Schedule *schedule) {
  int res, flag;
  NBC_Comminfo *comminfo;
  
  /* create a new state and return handle to it */
  handle->req_array = NULL;
  handle->req_count = 0;
  handle->comm = comm;
  /* first int is the schedule size */
  handle->row_offset = sizeof(int);
  handle->schedule = schedule;

  /******************** Do the tag and shadow comm administration ...  ***************/
  
  /* keyval is not initialized yet, we have to init it */
  if(MPI_KEYVAL_INVALID == gkeyval) {
    res = MPI_Keyval_create(NBC_Key_copy, NBC_Key_delete, &(gkeyval), NULL); 
    if((MPI_SUCCESS != res)) { printf("Error in MPI_Keyval_create() (%i)\n", res); return res; }
  } 

  res = MPI_Attr_get(comm, gkeyval, &comminfo, &flag);
  if((MPI_SUCCESS != res)) { printf("Error in MPI_Attr_get() (%i)\n", res); return res; }
  if (flag) {
    /* we found it */
    comminfo->tag++;
    handle->tag=comminfo->tag;
    handle->mycomm=comminfo->mycomm;
  } else {
    /* we have to create a new one */
    comminfo = malloc(sizeof(NBC_Comminfo));
    if(comminfo == NULL) { printf("Error in malloc()\n"); return NBC_OOR; }
    
    /* set tag to 1 */
    handle->tag=1;
    comminfo->tag=1;

    /* dup and save shadow communicator */
    res = MPI_Comm_dup(comm, &(handle->mycomm));
    if((MPI_SUCCESS != res)) { printf("Error in MPI_Comm_dup() (%i)\n", res); return res; }
    comminfo->mycomm=handle->mycomm;
    DEBUG(1, "created a shadow communicator for %lu ... %lu\n", (unsigned long)handle->comm, (unsigned long)handle->mycomm); \
  
    /* put the new attribute to the comm */
    res = MPI_Attr_put(comm, gkeyval, comminfo); 
    if((MPI_SUCCESS != res)) { printf("Error in MPI_Attr_put() (%i)\n", res); return res; }
  }
  
  /* reset counter ... */ 
  if(handle->tag == 32767) {  
    handle->tag=1;
    comminfo->tag=1;
    DEBUG(2,"resetting tags ...\n"); 
  } 
  
  /******************** end of tag and shadow comm administration ...  ***************/
  
  DEBUG(3, "got tag %i\n", handle->tag);

  /* kick off first round */
  res = NBC_Start_round(handle);
  if((NBC_OK != res)) { printf("Error in NBC_Start_round() (%i)\n", res); return res; }
  
  return NBC_OK;
}

int NBC_Wait_poll(NBC_Handle *handle) {
  int res;

  while(NBC_OK != NBC_Progress(handle));

  /* this should not be done here - schedules should be reused */
  res = NBC_Free(handle);
  if((NBC_OK != res)) { printf("Error in NBC_Free() (%i)\n", res); return res; }

  return NBC_OK;
}

int NBC_Wait(NBC_Handle *handle) {
  int res;
  
  NBC_Progress_block(handle);
  
  /* this should not be done here - schedules should be reused */
  res = NBC_Free(handle);
  if((NBC_OK != res)) { printf("Error in NBC_Free() (%i)\n", res); return res; }
  
  return NBC_OK;
}

