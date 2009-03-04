#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include "../src/spglib.h"

static PyObject * get_spacegroup(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {
  {"spacegroup", get_spacegroup, METH_VARARGS, "Space group finder"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initlibpyspg(void)
{
/*   PyObject* m = */
    Py_InitModule3("libpyspg", functions,
		   "C-extension for spglib\n\n...\n");
}

static PyObject * get_spacegroup(PyObject *self, PyObject *args)
{
  int i;
  int size;
  double symprec;
  char symbol[21];
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "iOOOd", &size, &lattice, &position, &atom_type, &symprec))
    return NULL;

  double* lat = (double*)lattice->data;
  double* pos = (double*)position->data;
  int* typat = (int*)atom_type->data;

  i = spg_get_international(symbol, lat, pos, typat, size, symprec);
  printf("%s (%d)\n", symbol, i);

  Py_RETURN_NONE;
}

