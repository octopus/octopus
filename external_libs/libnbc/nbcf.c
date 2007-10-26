#include "nbc.h"
#include "../../config.h"

void FC_FUNC_(nbcf_newhandle, NBCF_NEWHANDLE)
     (NBC_Handle **handle) {
  *handle = (NBC_Handle *)malloc(sizeof(NBC_Handle));
}

void FC_FUNC_(nbcf_freehandle, NBCF_FREEHANDLE)
  (NBC_Handle **handle) {
  free(*handle);
}

void FC_FUNC_(nbcf_ibarrier, NBCF_IBARRIER)
     (int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  ccomm = MPI_Comm_f2c(*comm);

  *ierr = NBC_Ibarrier(ccomm, *handle);
}

void FC_FUNC_(nbcf_ibcast, NBCF_IBCAST)
     (void *buffer, int *count, int *datatype, int *root,
      int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype cdatatype;

  ccomm = MPI_Comm_f2c(*comm);
  cdatatype = MPI_Type_f2c(*datatype);

  *ierr = NBC_Ibcast(buffer, *count, cdatatype, *root, ccomm, *handle);
}

void FC_FUNC_(nbcf_ialltoallv, NBCF_IALLTOALLV)
     (void *sendbuf, int *sendcounts, int *sdispls, int *sendtype,
      void *recvbuf, int *recvcounts, int *rdispls, int *recvtype,
      int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype csendtype, crecvtype;

  ccomm = MPI_Comm_f2c(*comm);
  csendtype = MPI_Type_f2c(*sendtype);
  crecvtype = MPI_Type_f2c(*recvtype);

  *ierr = NBC_Ialltoallv(sendbuf, sendcounts, sdispls, csendtype,
			 recvbuf, recvcounts, rdispls, crecvtype,
			 ccomm, *handle);
}

void FC_FUNC_(nbcf_igather, NBCF_IGATHER)
     (void *sendbuf, int *sendcount, int *sendtype,
      void *recvbuf, int *recvcount, int *recvtype,
      int *root, int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype csendtype, crecvtype;

  ccomm = MPI_Comm_f2c(*comm);
  csendtype = MPI_Type_f2c(*sendtype);
  crecvtype = MPI_Type_f2c(*recvtype);

  *ierr = NBC_Igather(sendbuf, *sendcount, csendtype,
		      recvbuf, *recvcount, crecvtype,
		      *root, ccomm, *handle);
}

void FC_FUNC_(nbcf_igatherv, NBCF_IGATHERV)
     (void *sendbuf, int *sendcount, int *sendtype,
      void *recvbuf, int *recvcounts, int *displs, int *recvtype,
      int *root, int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype csendtype, crecvtype;

  ccomm = MPI_Comm_f2c(*comm);
  csendtype = MPI_Type_f2c(*sendtype);
  crecvtype = MPI_Type_f2c(*recvtype);

  *ierr = NBC_Igatherv(sendbuf, *sendcount, csendtype,
		       recvbuf, recvcounts, displs, crecvtype,
		       *root, ccomm, *handle);
}

void FC_FUNC_(nbcf_iscatter, NBCF_ISCATTER)
     (void* sendbuf, int *sendcount, int *sendtype,
      void* recvbuf, int *recvcount, int *recvtype,
      int *root, int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype csendtype, crecvtype;

  ccomm = MPI_Comm_f2c(*comm);
  csendtype = MPI_Type_f2c(*sendtype);
  crecvtype = MPI_Type_f2c(*recvtype);

  *ierr = NBC_Iscatter(sendbuf, *sendcount, csendtype,
		       recvbuf, *recvcount, crecvtype,
		       *root, ccomm, *handle);
}

void FC_FUNC_(nbcf_iscatterv, NBCF_ISCATTERV)
     (void* sendbuf, int *sendcounts, int *displs, int *sendtype,
      void* recvbuf, int *recvcount, int *recvtype,
      int *root, int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype csendtype, crecvtype;

  ccomm = MPI_Comm_f2c(*comm);
  csendtype = MPI_Type_f2c(*sendtype);
  crecvtype = MPI_Type_f2c(*recvtype);

  *ierr = NBC_Iscatterv(sendbuf, sendcounts, displs, csendtype,
			recvbuf, *recvcount, crecvtype,
			*root, ccomm, *handle);
}

void FC_FUNC_(nbcf_iallgather, NBCF_IALLGATHER)
     (void *sendbuf, int *sendcount, int *sendtype,
      void *recvbuf, int *recvcount, int *recvtype,
      int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype csendtype, crecvtype;

  ccomm = MPI_Comm_f2c(*comm);
  csendtype = MPI_Type_f2c(*sendtype);
  crecvtype = MPI_Type_f2c(*recvtype);

  *ierr = NBC_Iallgather(sendbuf, *sendcount, csendtype,
			 recvbuf, *recvcount, crecvtype,
			 ccomm, *handle);
}

void FC_FUNC_(nbcf_iallgatherv, NBCF_IALLGATHERV)
     (void *sendbuf, int *sendcount, int *sendtype,
      void *recvbuf, int *recvcounts, int *displs, int *recvtype,
      int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype csendtype, crecvtype;

  ccomm = MPI_Comm_f2c(*comm);
  csendtype = MPI_Type_f2c(*sendtype);
  crecvtype = MPI_Type_f2c(*recvtype);

  *ierr = NBC_Iallgatherv(sendbuf, *sendcount, csendtype,
			  recvbuf, recvcounts, displs, crecvtype,
			  ccomm, *handle);
}

void FC_FUNC_(nbcf_ialltoall, NBCf_IALLTOALL)
     (void *sendbuf, int *sendcount, int *sendtype,
      void *recvbuf, int *recvcount, int *recvtype,
      int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Datatype csendtype, crecvtype;

  ccomm = MPI_Comm_f2c(*comm);
  csendtype = MPI_Type_f2c(*sendtype);
  crecvtype = MPI_Type_f2c(*recvtype);

  *ierr = NBC_Ialltoall(sendbuf, *sendcount, csendtype,
			recvbuf, *recvcount, crecvtype,
			ccomm, *handle);
}

void FC_FUNC_(nbcf_ireduce, NBCF_IREDUCE)
     (void *sendbuf, void *recvbuf, int *count, int *datatype,
      int *op, int *root, int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Op cop;
  MPI_Datatype cdatatype;

  ccomm = MPI_Comm_f2c(*comm);
  cop = MPI_Op_f2c(*op);
  cdatatype = MPI_Type_f2c(*datatype);

  *ierr = NBC_Ireduce(sendbuf, recvbuf, *count, cdatatype,
		      cop, *root, ccomm, *handle);
}

void FC_FUNC_(nbcf_iallreduce, NBCF_IALLREDUCE)
     (void *sendbuf, void *recvbuf, int *count, int *datatype,
      int *op, int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Op cop;
  MPI_Datatype cdatatype;

  ccomm = MPI_Comm_f2c(*comm);
  cop = MPI_Op_f2c(*op);
  cdatatype = MPI_Type_f2c(*datatype);

  *ierr = NBC_Iallreduce(sendbuf, recvbuf, *count, cdatatype,
			 cop, ccomm, *handle);
}

void FC_FUNC_(nbcf_iscan, NBCF_ISCAN)
     (void *sendbuf, void *recvbuf, int *count, int *datatype,
      int *op, int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Op cop;
  MPI_Datatype cdatatype;

  ccomm = MPI_Comm_f2c(*comm);
  cop = MPI_Op_f2c(*op);
  cdatatype = MPI_Type_f2c(*datatype);

  *ierr = NBC_Iscan(sendbuf, recvbuf, *count, cdatatype,
		    cop, ccomm, *handle);
}

void FC_FUNC_(nbcf_ireduce_scatter, NBCF_IREDUCE_SCATTER)
     (void *sendbuf, void *recvbuf, int *recvcounts, int *datatype,
      int *op, int *comm, NBC_Handle **handle, int *ierr) {
  MPI_Comm ccomm;
  MPI_Op cop;
  MPI_Datatype cdatatype;

  ccomm = MPI_Comm_f2c(*comm);
  cop = MPI_Op_f2c(*op);
  cdatatype = MPI_Type_f2c(*datatype);

  *ierr = NBC_Ireduce_scatter(sendbuf, recvbuf, recvcounts, cdatatype,
			      cop, ccomm, *handle);
}

void FC_FUNC_(nbcf_test, NBCF_TEST)
     (NBC_Handle **handle, int *ierr) {
  *ierr = NBC_Test(*handle);
}

void FC_FUNC_(nbcf_wait, NBCF_WAIT)
     (NBC_Handle **handle, int *ierr) {
  *ierr = NBC_Wait(*handle);
}

void FC_FUNC_(nbcf_wait_poll, NBCF_WAIT_POLL)
     (NBC_Handle **handle, int *ierr) {
  *ierr = NBC_Wait_poll(*handle);
}
