/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#include "nbc.h"

/* simple linear MPI_Iallgatherv
 * the algorithm uses p-1 rounds
 * first round:
 *   each node sends to it's left node (rank+1)%p sendcount elements 
 *   each node begins with it's right node (rank-11)%p and receives from it recvcounts[(rank+1)%p] elements
 * second round: 
 *   each node sends to node (rank+2)%p sendcount elements 
 *   each node receives from node (rank-2)%p recvcounts[(rank+2)%p] elements */
int NBC_Iallgatherv(void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int *recvcounts, int *displs, MPI_Datatype recvtype, MPI_Comm comm, NBC_Handle *handle) {
  int rank, p, res, r, speer, rpeer;
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
  if (NULL == schedule) { printf("Error in malloc() (%i)\n", res); return res; }

  handle->tmpbuf=NULL;
 
  res = NBC_Sched_create(schedule);
  if(res != NBC_OK) { printf("Error in NBC_Sched_create, (%i)\n", res); return res; }
  
  /* copy my data to receive buffer */
  rbuf = ((char *)recvbuf) + (displs[rank]*rcvext);
  NBC_Copy(sendbuf, sendcount, sendtype, rbuf, recvcounts[rank], recvtype, comm);
  if (NBC_OK != res) { printf("Error in NBC_Copy() (%i)\n", res); return res; }

  /* do p-1 rounds */
  for(r=1;r<p;r++) {
    speer = (rank+r)%p;
    rpeer = (rank-r+p)%p;
    rbuf = ((char *)recvbuf) + (displs[rpeer]*rcvext);
    
    res = NBC_Sched_recv(rbuf, recvcounts[rpeer], recvtype, rpeer, schedule);
    if (NBC_OK != res) { printf("Error in NBC_Sched_recv() (%i)\n", res); return res; }
    res = NBC_Sched_send(sendbuf, sendcount, sendtype, speer, schedule);
    if (NBC_OK != res) { printf("Error in NBC_Sched_send() (%i)\n", res); return res; }
  }

  res = NBC_Sched_commit(schedule);
  if (NBC_OK != res) { printf("Error in NBC_Sched_commit() (%i)\n", res); return res; }
 
  res = NBC_Start(handle, comm, schedule);
  if (NBC_OK != res) { printf("Error in NBC_Start() (%i)\n", res); return res; }
 
  return NBC_OK;
}
