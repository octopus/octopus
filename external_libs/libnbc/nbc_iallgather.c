/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#include "nbc.h"

/* simple linear MPI_Iallgather 
 * the algorithm uses p-1 rounds
 * each node sends the packet it received last round (or has in round 0) to it's right neighbor (modulo p)
 * each node receives from it's left (modulo p) neighbor */
int NBC_Iallgather(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, MPI_Comm comm, NBC_Handle *handle) {
  int rank, p, res, r;
  MPI_Aint rcvext, sndext;
  NBC_Schedule *schedule;
  char *rbuf;
  
  res = MPI_Comm_rank(comm, &rank);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_rank() (%i)\n", res); return res; }
  res = MPI_Comm_size(comm, &p);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_size() (%i)\n", res); return res; }
  res = MPI_Type_extent(sendtype, &sndext);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Type_extent() (%i)\n", res); return res; }
  res = MPI_Type_extent(recvtype, &rcvext);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Type_extent() (%i)\n", res); return res; }

  schedule = malloc(sizeof(NBC_Schedule));
  if (NULL == schedule) { printf("Error in malloc()\n"); return res; }

  handle->tmpbuf = NULL;
 
  res = NBC_Sched_create(schedule);
  if(NBC_OK != res) { printf("Error in NBC_Sched_create, (%i)\n", res); return res; }
  
  /* do p-1 rounds */
  for(r=0;r<p;r++) {
    if(r == rank) {
      /* it's me :) - copy my data to receive buffer */
      rbuf = ((char *)recvbuf) + (rank*recvcount*rcvext);
      res = NBC_Copy(sendbuf, sendcount, sendtype, rbuf, recvcount, recvtype, comm);
      if (NBC_OK != res) { printf("Error in NBC_Copy() (%i)\n", res); return res; }
    } else {
      /* recv from rank r */
      rbuf = ((char *)recvbuf) + r*(recvcount*rcvext);
      res = NBC_Sched_recv(rbuf, recvcount, recvtype, r, schedule);
      if (NBC_OK != res) { printf("Error in NBC_Sched_recv() (%i)\n", res); return res; }
      /* send to rank r */
      res = NBC_Sched_send(sendbuf, sendcount, sendtype, r, schedule);
      if (NBC_OK != res) { printf("Error in NBC_Sched_send() (%i)\n", res); return res; }
    }
  }

  res = NBC_Sched_commit(schedule);
  if (NBC_OK != res) { printf("Error in NBC_Sched_commit() (%i)\n", res); return res; }

  /*NBC_PRINT_SCHED(*schedule);*/
 
  res = NBC_Start(handle, comm, schedule);
  if (NBC_OK != res) { printf("Error in NBC_Start() (%i)\n", res); return res; }
 
  return NBC_OK;
}
