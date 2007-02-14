/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#include "nbc.h"

/* Dissemination implementation of MPI_Ibarrier */
int NBC_Ibarrier(MPI_Comm comm, NBC_Handle* handle) {
  int round, rank, p, maxround, res, recvpeer, sendpeer;
  NBC_Schedule *schedule;
  
  res = MPI_Comm_rank(comm, &rank);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_rank() (%i)\n", res); return res; }
  res = MPI_Comm_size(comm, &p);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_size() (%i)\n", res); return res; }
  
  schedule = malloc(sizeof(NBC_Schedule));
  if (NULL == schedule) { printf("Error in malloc()\n"); return res; }
  
  round = -1;
  handle->tmpbuf=NULL;

  res = NBC_Sched_create(schedule);
  if(res != NBC_OK) { printf("Error in NBC_Sched_create (%i)\n", res); return res; }

  maxround = (int)ceil((log(p)/LOG2)-1);

  do {
    round++;
    sendpeer = (rank + (1<<round)) % p;
    /* add p because modulo does not work with negative values */
    recvpeer = ((rank - (1<<round))+p) % p;

    /* send msg to sendpeer */
    res = NBC_Sched_send(NULL, 0, MPI_BYTE, sendpeer, schedule);
    if (NBC_OK != res) { printf("Error in NBC_Sched_send() (%i)\n", res); return res; }

    /* recv msg from recvpeer */
    res = NBC_Sched_recv(NULL, 0, MPI_BYTE, recvpeer, schedule);
    if (NBC_OK != res) { printf("Error in NBC_Sched_recv() (%i)\n", res); return res; }
    /* end communication round */
    if(round < maxround){
      res = NBC_Sched_barrier(schedule);
      if (NBC_OK != res) { printf("Error in NBC_Sched_barrier() (%i)\n", res); return res; }
    }
  } while (round < maxround);

  res = NBC_Sched_commit(schedule);
  if (NBC_OK != res) { printf("Error in NBC_Sched_commit() (%i)\n", res); return res; }
  
  res = NBC_Start(handle, comm, schedule);
  if (NBC_OK != res) { printf("Error in NBC_Start() (%i)\n", res); return res; }
  
  return NBC_OK;
}

/*void NBC_IBARRIER(MPI_Fint *comm, MPI_Fint *ierr) {
 *ierr = NBC_Ibarrier(MPI_Comm comm, NBC_Handle* handle);
}*/
