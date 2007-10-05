/*
 * Copyright (c) 2006 The Trustees of Indiana University and Indiana
 *                    University Research and Technology
 *                    Corporation.  All rights reserved.
 * Copyright (c) 2006 The Technical University of Chemnitz. All 
 *                    rights reserved.
 */
/* defines the Fortran 77 name mangling scheme */

#ifdef MANGLE_NBC_F77_0
#define NBC_F77_FUNC(name,NAME) name
#define NBC_F77_FOUND
#endif

#ifdef MANGLE_NBC_F77_1
#define NBC_F77_FUNC_(name,NAME) name ## _
#define NBC_F77_FOUND
#endif

#ifdef MANGLE_NBC_F77_2
#define NBC_F77_FUNC_(name,NAME) name ## __
#define NBC_F77_FOUND
#endif

#ifdef MANGLE_NBC_F77_3
#define NBC_F77_FUNC_(name,NAME) NAME
#define NBC_F77_FOUND
#endif

#ifdef MANGLE_NBC_F77_4
#define NBC_F77_FUNC_(name,NAME) NAME ## _
#define NBC_F77_FOUND
#endif

#ifdef MANGLE_NBC_F77_5
#define NBC_F77_FUNC_(name,NAME) NAME ## __
#define NBC_F77_FOUND
#endif

/* default, e.g., if no cmake could be used */
#ifndef NBC_F77_FOUND
#define NBC_F77_FUNC_(name,NAME) name ## _
#endif

