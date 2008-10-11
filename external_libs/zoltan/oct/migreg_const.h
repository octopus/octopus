/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: migreg_const.h,v $
 *    $Author: kddevin $
 *    $Date: 2002/06/19 23:56:41 $
 *    Revision: 1.11 $
 ****************************************************************************/

#ifndef __MIGREG_CONST_H
#define __MIGREG_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


extern int Zoltan_Oct_migreg_migrate_orphans(ZZ *zz, pRegion RegionList, int nreg, 
				      int level, Map *array, int *c1, int *c2);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
