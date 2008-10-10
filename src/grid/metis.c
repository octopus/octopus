/*
 Copyright (C) 2002-2008 M. Marques, A. Castro, A. Rubio, G. Bertsch

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
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.

 $Id$
*/

#include <config.h>

/* -------------------------- interface to METIS ----------------------------- */
#if defined(HAVE_METIS)
#include <metis.h>

void FC_FUNC_(oct_metis_part_mesh_nodal, OCT_METIS_PART_MESH_NODAL)
  (int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, 
   int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}

void FC_FUNC_(oct_metis_part_mesh_dual, OCT_METIS_PART_MESH_DUAL)
  (int *ne, int *nn, idxtype *elmnts, int *etype, int *numflag, int *nparts, 
   int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}

void FC_FUNC_(oct_metis_part_graph_recursive, OCT_METIS_PART_GRAPH_RECURSIVE)
  (int *n, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
   int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

void FC_FUNC_(oct_metis_part_graph_kway, OCT_METIS_PART_GRAPH_KWAY)
  (int *n, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
   int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}

void FC_FUNC_(oct_metis_part_graph_vkway, OCT_METIS_PART_GRAPH_VKWAY)
  (int *n, idxtype *xadj, idxtype *adjncy, idxtype *vwgt,
   idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts,
   int *options, int *edgecut, idxtype *part)
{
  METIS_PartGraphVKway(n, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
#endif
