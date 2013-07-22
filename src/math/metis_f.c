/*
 Copyright (C) 2013 M. Oliveira

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

 $Id$
*/


#include <config.h>

#include <stdlib.h>

#if defined(HAVE_METIS)
#include <metis.h>
#endif
#if defined(HAVE_PARMETIS)
#include <parmetis.h>
#include <mpi.h>
#endif


#ifdef HAVE_METIS
 
void FC_FUNC_(oct_metis_setdefaultoptions, OCT_METIS_SETDEFAULTOPTIONS)
     (idx_t *options)
{
  METIS_SetDefaultOptions(options);
}


void FC_FUNC_(oct_metis_partgraphrecursive, OCT_METIS_PARTGRAPHRECURSIVE)
     (idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, idx_t *nparts, 
      real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, idx_t *part)
{

  METIS_PartGraphRecursive(nvtxs, ncon, xadj, adjncy, NULL, NULL, NULL, nparts, 
			   tpwgts, ubvec, options, objval, part);
}


void FC_FUNC_(oct_metis_partgraphkway, OCT_METIS_PARTGRAPHKWAY)
     (idx_t *nvtxs, idx_t *ncon, idx_t *xadj, idx_t *adjncy, idx_t *nparts, 
      real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *objval, idx_t *part)
{

  METIS_PartGraphKway(nvtxs, ncon, xadj, adjncy, NULL, NULL, NULL, nparts, 
		      tpwgts, ubvec, options, objval, part);
}

#endif



#ifdef HAVE_PARMETIS

void FC_FUNC_(oct_parmetis_v3_partkway, OCT_PARMETIS_PARTKWAY)
     (idx_t *vtxdist, idx_t *xadj, idx_t *adjncy, idx_t *ncon, 
      idx_t *nparts, real_t *tpwgts, real_t *ubvec, idx_t *options, 
      idx_t *edgecut, idx_t *part, MPI_Fint *fcomm)
{
  idx_t wgtflag = 0, numflag = 1;

  MPI_Comm comm;

#ifdef HAVE_MPI2
  comm = MPI_Comm_f2c(*fcomm);
#else
  comm = *fcomm;
#endif

  ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, NULL, NULL, &wgtflag, &numflag, 
		       ncon, nparts, tpwgts, ubvec, options, edgecut, part, &comm);
}

#endif
