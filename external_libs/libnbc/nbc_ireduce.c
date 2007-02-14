/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
#include "nbc.h"

/* binomial reduce
 * working principle:
 * - each node gets a virtual rank vrank
 * - the 'root' node get vrank 0 
 * - node 0 gets the vrank of the 'root'
 * - all other ranks stay identical (they do not matter)
 *
 * Algorithm:
 * pairwise exchange
 * round r: 
 *  grp = rank % 2^r
 *  if grp == 0: receive from rank + 2^(r-1) if it exists and reduce value
 *  if grp == 1: send to rank - 2^(r-1) and exit function
 *  
 * do this for R=log_2(p) rounds
 *    
 */
#define RANK2VRANK(rank, vrank, root) \
{ \
  vrank = rank; \
  if (rank == 0) vrank = root; \
  if (rank == root) vrank = 0; \
}
#define VRANK2RANK(rank, vrank, root) \
{ \
  rank = vrank; \
  if (vrank == 0) rank = root; \
  if (vrank == root) rank = 0; \
}
int NBC_Ireduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm, NBC_Handle* handle) {
  int vrank, peer, vpeer, rank, maxr, p, r, res, firstred;
  MPI_Aint ext;
  NBC_Schedule *schedule;
  char *redbuf=NULL;
  
  res = MPI_Comm_rank(comm, &rank);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_rank() (%i)\n", res); return res; }
  res = MPI_Comm_size(comm, &p);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Comm_size() (%i)\n", res); return res; }
  res = MPI_Type_extent(datatype, &ext);
  if (MPI_SUCCESS != res) { printf("MPI Error in MPI_Type_extent() (%i)\n", res); return res; }
  
  schedule = malloc(sizeof(NBC_Schedule));
  if (NULL == schedule) { printf("Error in malloc() (%i)\n", res); return res; }

  if(rank == root) {
    /* root reduces in receivebuffer */
    handle->tmpbuf = malloc(ext*count);
  } else {
    /* recvbuf may not be valid on non-root nodes */
    handle->tmpbuf = malloc(ext*count*2);
    redbuf = ((char*)handle->tmpbuf)+(ext*count);
  }
  if (NULL == handle->tmpbuf) { printf("Error in malloc() (%i)\n", res); return res; }

  res = NBC_Sched_create(schedule);
  if(res != NBC_OK) { printf("Error in NBC_Sched_create (%i)\n", res); return res; }

  RANK2VRANK(rank, vrank, root);
  maxr = (int)ceil((log(p)/LOG2));

  /* only one node -> copy data */
  if(p == 1) {
    res = NBC_Copy(sendbuf, count, datatype, recvbuf, count, datatype, comm);
    if (NBC_OK != res) { printf("Error in NBC_Copy() (%i)\n", res); return res; }
  }

  firstred = 1;
  for(r=1; r<=maxr; r++) {
    if((vrank % (1<<r)) == 0) {
      /* we have to receive this round */
      vpeer = vrank + (1<<(r-1));
      VRANK2RANK(peer, vpeer, root)
      if(peer<p) {
        res = NBC_Sched_recv(handle->tmpbuf, count, datatype, peer, schedule);
        if (NBC_OK != res) { free(handle->tmpbuf); printf("Error in NBC_Sched_recv() (%i)\n", res); return res; }
        /* we have to wait until we have the data */
        res = NBC_Sched_barrier(schedule);
        if (NBC_OK != res) { free(handle->tmpbuf); printf("Error in NBC_Sched_barrier() (%i)\n", res); return res; }
        /* perform the reduce in my local buffer */
        if(firstred) {
          if(rank == root) {
            /* root is the only one who reduces in the receivebuffer 
             * take data from sendbuf in first round - save copy */
            res = NBC_Sched_op(recvbuf, sendbuf, handle->tmpbuf, count, datatype, op, schedule);
          } else {
            /* all others may not have a receive buffer 
             * take data from sendbuf in first round - save copy */
            res = NBC_Sched_op(redbuf, sendbuf, handle->tmpbuf, count, datatype, op, schedule);
          }
          firstred = 0;
        } else {
          if(rank == root) {
            /* root is the only one who reduces in the receivebuffer */
            res = NBC_Sched_op(recvbuf, recvbuf, handle->tmpbuf, count, datatype, op, schedule);
          } else {
            /* all others may not have a receive buffer */
            res = NBC_Sched_op(redbuf, redbuf, handle->tmpbuf, count, datatype, op, schedule);
          }
        }
        if (NBC_OK != res) { free(handle->tmpbuf); printf("Error in NBC_Sched_op() (%i)\n", res); return res; }
        /* this cannot be done until handle->tmpbuf is unused :-( */
        res = NBC_Sched_barrier(schedule);
        if (NBC_OK != res) { free(handle->tmpbuf); printf("Error in NBC_Sched_barrier() (%i)\n", res); return res; }
      }
    } else {
      /* we have to send this round */
      vpeer = vrank - (1<<(r-1));
      VRANK2RANK(peer, vpeer, root)
      if(firstred) {
        /* we did not reduce anything */
        res = NBC_Sched_send(sendbuf, count, datatype, peer, schedule);
      } else {
        /* we have to use the redbuf the root (which works in receivebuf) is never sending .. */
        res = NBC_Sched_send(redbuf, count, datatype, peer, schedule);
      }
      if (NBC_OK != res) { free(handle->tmpbuf); printf("Error in NBC_Sched_send() (%i)\n", res); return res; }
      /* leave the game */
      break;
    }
  }
  
  res = NBC_Sched_commit(schedule);
  if (NBC_OK != res) { free(handle->tmpbuf); printf("Error in NBC_Sched_commit() (%i)\n", res); return res; }
  
  res = NBC_Start(handle, comm, schedule);
  if (NBC_OK != res) { free(handle->tmpbuf); printf("Error in NBC_Start() (%i)\n", res); return res; }
  
  /* tmpbuf is freed with the handle */
  return NBC_OK;
}
