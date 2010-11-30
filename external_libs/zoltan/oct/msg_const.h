/*****************************************************************************
 * CVS File Information :
 *    $RCSfile: msg_const.h,v $
 *    $Author: kddevin $
 *    $Date: 2002/06/19 23:56:42 $
 *    Revision: 1.9 $
 ****************************************************************************/

#ifndef __MSG_CONST_H
#define __MSG_CONST_H

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


extern int    Zoltan_Oct_msg_int_scan(MPI_Comm, int, int value);
extern float  Zoltan_Oct_msg_float_scan(MPI_Comm, int, float value);

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif

#endif
