/*
 Copyright (C) 2016 X. Andrade

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

#ifdef HAVE_CUDA
#include <cuda.h>
#include <nvrtc.h>
#else
typedef int CUcontext;
typedef int CUdevice;
typedef int CUmodule;
typedef int CUfunction;
#endif

#include <stdlib.h> //we have to include this before cmath to workaround a bug in the PGI "compiler".
#include <cmath>

#include <iostream>

#include <fstream>

#include "string_f.h" /* fortran <-> c string compatibility issues */

#include <vector>
#include <sstream>
#include <iterator>

#include <fortran_types.h>

using namespace std;

extern "C" void FC_FUNC_(cuda_init, CUDA_INIT)(CUcontext ** context, CUdevice ** device){

#ifdef HAVE_CUDA
  CUresult err = cuInit(0);

  *context = new CUcontext;
  *device = new CUdevice;

  int ndevices;

  cuDeviceGetCount(&ndevices);
  
  cout << "Number of devices : " << ndevices << endl;

  if (ndevices == 0) {
    cerr << "Error: no CUDA devices" << std::endl;
    exit(1);
  }
  cuDeviceGet(*device, 0);

  char devicename[200];

  cuDeviceGetName(devicename, sizeof(devicename), **device);

  cout << "CUDA device name : " << devicename << endl;

  int major = 0, minor = 0;
  cuDeviceComputeCapability(&major, &minor, **device);
  cout << "CUDA compute capability : " <<  major << '.' <<  minor << endl;
  
  size_t mem;
  cuDeviceTotalMem(&mem, **device);
  cout << "CUDA device memory : " << mem/pow(1024.0, 3) << " GB" << endl;

  cuCtxCreate(*context, 0, **device);
#endif
}

extern "C" void FC_FUNC_(cuda_end, CUDA_END)(CUcontext ** context, CUdevice ** device){
#ifdef HAVE_CUDA

  cuCtxDestroy(**context);
  
  delete *context;
  delete *device;
#endif
}

extern "C" void FC_FUNC_(cuda_build_program, CUDA_BUILD_PROGRAM)(CUmodule ** module, STR_F_TYPE include_path, STR_F_TYPE const fname, STR_F_TYPE const flags STR_ARG3){
#ifdef HAVE_CUDA
  char *include_path_c;
  char *fname_c;
  char *flags_c;

  TO_C_STR1(include_path, include_path_c);
  TO_C_STR2(fname, fname_c);
  TO_C_STR3(flags, flags_c);
  
  // read the source
  
  string source;

  source = "#include \"" + string(fname_c) + "\"\n";

  // cout << source << "|" << endl;

  cout << fname_c << endl;

  nvrtcProgram prog;
  nvrtcCreateProgram(&prog, source.c_str(), "kernel_include.c", 0, NULL,NULL);
  free(fname_c);

  string all_flags = "-DCUDA -default-device " + string("-I") + include_path_c + string(" ") + string(flags_c);

  stringstream flags_stream(all_flags);
  istream_iterator<string> iter(flags_stream);
  istream_iterator<string> end;
  vector<string> tokens(iter, end);
  
  const char ** opts = new const char*[tokens.size()];
  for (unsigned ii = 0; ii < tokens.size(); ii++) opts[ii] = tokens[ii].c_str();

  cout << "OPTS " << all_flags << endl;

  nvrtcResult err = nvrtcCompileProgram(prog, tokens.size(), opts);

  free(flags_c);
  
  size_t logSize;
  nvrtcGetProgramLogSize(prog, &logSize);
  char *log = new char[logSize];
  nvrtcGetProgramLog(prog, log);  

  cout << "=====================================" << endl;

  cout << log << endl;

  if(CUDA_SUCCESS != err){
    cerr << "Error in compiling" << endl;
    exit(1);
  }
  
  delete [] log;

  // Obtain PTX from the program.
  size_t ptxSize;
  nvrtcGetPTXSize(prog, &ptxSize);
  char *ptx = new char[ptxSize];
  nvrtcGetPTX(prog, ptx);

  nvrtcDestroyProgram(&prog);

  *module = new CUmodule;

  cuModuleLoadDataEx(*module, ptx, 0, 0, 0);

  delete [] ptx;
#endif
}

extern "C" void FC_FUNC_(cuda_create_kernel, CUDA_CREATE_KERNEL)(CUfunction ** kernel, CUmodule ** module, STR_F_TYPE kernel_name STR_ARG1){
#ifdef HAVE_CUDA  
  char *kernel_name_c;

  TO_C_STR1(kernel_name, kernel_name_c);
  
  *kernel = new CUfunction;

  cuModuleGetFunction(*kernel, **module, kernel_name_c);

  free(kernel_name_c);
#endif
}

extern "C" void FC_FUNC_(cuda_release_module, CUDA_RELEASE_MODULE)(CUmodule ** module){
#ifdef HAVE_CUDA
  cuModuleUnload(**module); 
  delete *module;
#endif
}

extern "C" void FC_FUNC_(cuda_release_kernel, CUDA_RELEASE_KERNEL)(CUfunction ** kernel){
#ifdef HAVE_CUDA
  delete *kernel;
#endif
}

extern "C" void FC_FUNC_(cuda_device_max_threads_per_block, CUDA_DEVICE_MAX_THREADS_PER_BLOCK)(CUdevice ** device, fint * max_threads){
#ifdef HAVE_CUDA
  int value;
  cuDeviceGetAttribute (&value, CU_DEVICE_ATTRIBUTE_MAX_THREADS_PER_BLOCK, **device);
  *max_threads = value;
#endif
}

extern "C" void FC_FUNC_(cuda_mem_alloc, CUDA_MEM_ALLOC)(CUdeviceptr ** cuda_ptr, const fint8 * size){
#ifdef HAVE_CUDA
  *cuda_ptr = new CUdeviceptr;
  cuMemAlloc(*cuda_ptr, *size);
#endif
}

extern "C" void FC_FUNC_(cuda_mem_free, CUDA_MEM_FREE)(CUdeviceptr ** cuda_ptr){
#ifdef HAVE_CUDA
  cuMemFree(**cuda_ptr);
  delete *cuda_ptr;
#endif
}


