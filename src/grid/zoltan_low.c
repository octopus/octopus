/*
  Copyright (C) 2008 X. Andrade

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

  $Id: recipes.c 2146 2006-05-23 17:36:00Z xavier $
*/
#include <assert.h>
#include <config.h>

#ifdef HAVE_MPI

#include <zoltan.h>

typedef struct {
  int np;
  int np_part_global;
  int dim;
  double * x_global;
  int * part;
  int ipart;
  int nedges;
  int estart;
  int * xedges;
  int * edges;
} mesh_t;

static mesh_t mesh;

/* callbacks */

int get_dimension(void * data, int *err)
{
  *err = 0;
  return mesh.dim;
}

int get_num_objects(void * data, int *err)
{
  int nn, i;
  *err = 0;

  nn = 0;

  for (i = 0; i < mesh.np; i++) {
    if(mesh.part[i] == mesh.ipart) nn++;
  }

  return nn;
}

void get_local_objects_list(void * data, int num_gids, int num_lids,
		     ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int wgt_dim, float *obj_wgts,
		     int *err)
{
/* ZOLTAN_OBJ_LIST_FN callback function.
** Returns list of objects owned by this processor.
*/
  int nn, i;
    
  *err = 0;
  nn = 0;

  for (i = 0; i < mesh.np; i++){
    if(mesh.part[i] == mesh.ipart){
      gids[nn] = i;
      nn++;
    }
  }
  
}

void get_objects_coords(void *data, int num_gids, int num_lids, int num_objs,
  ZOLTAN_ID_PTR gids, ZOLTAN_ID_PTR lids, int numDim, double *pts, int *err)
{ 
/* ZOLTAN_GEOM_MULTI_FN callback.
** Returns coordinates of objects listed in gids and lids.
*/
  int i, id;
  int idir;
  
  if (numDim != mesh.dim){
    *err = 1;
    return;
  }

  for (i=0; i < num_objs; i++){
    id = gids[i];
    
    if ((id < 0) || (id >= mesh.np))
      {
	*err = 1;
	return;
      }
    
    for(idir = 0; idir < mesh.dim; idir++){
      pts[i*mesh.dim + idir] = mesh.x_global[idir*mesh.np_part_global + id];
    }
  }

} 

/* graph functions */

void get_num_edges(void *data, int dim_gid, int dim_lid, int num_objs, 
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *num_edges, int *err){

  int iobj;
  int gid;

  for(iobj = 0; iobj < num_objs; iobj++) {
    gid = global_ids[iobj] - mesh.estart;
    assert(gid >= 0);
    num_edges[iobj] = mesh.xedges[gid + 1] - mesh.xedges[gid];
    assert(num_edges[iobj] > 0);
  }
  *err = 0;
}

void get_edges(void *data, int dim_gid, int dim_lid, ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, 
  ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs, int wgt_dim, float *ewgts, int *err){

  int iedge;
  int gid = global_id[0];
  int ii;

  assert(gid < mesh.np);
  
  gid -= mesh.estart;

  assert(gid >= 0);

  ii = 0;
  for(iedge = mesh.xedges[gid]; iedge < mesh.xedges[gid + 1]; iedge++){
    nbor_global_id[ii] = mesh.edges[iedge - 1] - 1;
    assert(nbor_global_id[ii] < mesh.np);
    nbor_procs[ii] = mesh.part[nbor_global_id[ii]] - 1; 
    ii++;
  }
  
  assert(ii == mesh.xedges[gid + 1] - mesh.xedges[gid]);
  
  *err = 0;

}

/* these values have to match with the ones defined in zoltan.F90 and
   mesh_init.F90*/

#define RCB         2
#define RIB         3
#define HSFC        4
#define REFTREE     5
#define GRAPH       6
#define HYPERGRAPH  7

void FC_FUNC_(zoltan_partition, ZOLTAN_PARTITION)(const int * method,
						  const int * sbdim, 
						  const int * np, 
						  const int * np_part_global, 
						  double * x_global,
						  int * estart,
						  int * xedges,
						  int * edges,
						  const int * ipart,
						  int * part,
						  int * fcomm
						  ){
  struct Zoltan_Struct *zz;
  int ii;
  char name[] = "octopus_mpi";
  char* argv = name; 
  float version;
  int rc;
  MPI_Comm comm;

  /* all these variables are for the partition function */
  int changes;
  int num_gid_entries;
  int num_lid_entries;
  int num_import;
  ZOLTAN_ID_PTR import_global_ids;
  ZOLTAN_ID_PTR import_local_ids;
  int *import_procs;
  int *import_to_part;
  int num_export;
  ZOLTAN_ID_PTR export_global_ids;
  ZOLTAN_ID_PTR export_local_ids;
  int *export_procs;
  int *export_to_part;

  mesh.dim = sbdim[0];
  mesh.np = np[0];
  mesh.np_part_global = np_part_global[0];
  mesh.x_global = x_global;
  mesh.xedges = xedges;
  mesh.estart = estart[0] - 1; /* convert to C index convention */
  mesh.edges = edges;
  mesh.ipart = ipart[0];
  mesh.part = part;

  Zoltan_Initialize(1, &argv, &version);

  comm = MPI_Comm_f2c(*fcomm);

  zz = Zoltan_Create(comm);

  rc = Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");

  /* the method to partition the grid */
  switch(*method){
  case(RCB):
    rc = Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
    break;
  case(RIB):
    rc = Zoltan_Set_Param(zz, "LB_METHOD", "RIB");
    break;
  case(HSFC):
    rc = Zoltan_Set_Param(zz, "LB_METHOD", "HSFC");
    break;
  case(REFTREE):
    rc = Zoltan_Set_Param(zz, "LB_METHOD", "REFTREE");
    break;
  case(GRAPH):
    rc = Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
    break;
  case(HYPERGRAPH):
    rc = Zoltan_Set_Param(zz, "LB_METHOD", "HYPERGRAPH");
    break;
  }

  assert(rc == ZOLTAN_OK);

  /* always required callbacks */
  rc = Zoltan_Set_Num_Obj_Fn(zz, get_num_objects, NULL);
  assert(rc == ZOLTAN_OK);

  rc = Zoltan_Set_Obj_List_Fn(zz, get_local_objects_list, NULL);
  assert(rc == ZOLTAN_OK);

  if(*method == RCB || *method == RIB || *method == HSFC || *method == REFTREE) {
    /* register geometry callbacks */
    rc = Zoltan_Set_Num_Geom_Fn(zz, get_dimension, NULL);
    assert(rc == ZOLTAN_OK);
    
    rc = Zoltan_Set_Geom_Multi_Fn(zz, get_objects_coords, NULL);
    assert(rc == ZOLTAN_OK);
  }

  if(*method == GRAPH || *method == HYPERGRAPH ) {
    /* since our original distribution is very bad, we tell zoltan to
       start the partition from scratch */
    rc = Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION");
    assert(rc == ZOLTAN_OK);

    /* register graph callbacks */
    rc = Zoltan_Set_Num_Edges_Multi_Fn(zz, get_num_edges, NULL);
    assert(rc == ZOLTAN_OK);
    
    rc = Zoltan_Set_Edge_List_Fn(zz, get_edges, NULL);
    assert(rc == ZOLTAN_OK);
  }

  /* tell zoltan that we want a list of the assignment of all local
     points as a result */
  rc = Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTITION ASSIGNMENTS"); 
  assert(rc == ZOLTAN_OK);

  Zoltan_LB_Partition (zz,
		       &changes, /* boolean indicating if the partition has changed */
		       &num_gid_entries, /* number of dimension of the global indices, 1 for us */
		       &num_lid_entries, /* number of dimension of the local indices, 1 for us */
		       /* since we use "RESULT_LISTS" = "PARTITION_ASSIGNMENTS" 
			  many arguments don't have values */
		       &num_import,        /* -1 */
		       &import_global_ids, /* NULL */
		       &import_local_ids,  /* NULL */
		       &import_procs,      /* NULL */
		       &import_to_part,    /* NULL */
		       &num_export,        /* number of points in this processor */
		       &export_global_ids, 
		       &export_local_ids,
		       &export_procs,      /* the partition to which each processor is assigned */
		       &export_to_part);

  /* safety measure */
  for(ii = 0; ii < mesh.np; ii++) part[ii] = -1; 

  /* convert to fortran index convention and place it in the
     correspondind region of the part array (the algathering will be
     performed by the caller) */
  for(ii = 0; ii < num_export; ii++) part[mesh.estart + ii] = export_procs[ii] + 1; 

  Zoltan_LB_Free_Part (&export_global_ids,
		       &export_local_ids,
		       &export_procs, 
		       &export_to_part);

  Zoltan_Destroy(&zz);

}

#endif
