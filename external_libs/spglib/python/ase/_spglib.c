#include <Python.h>
#include <stdio.h>
#include <numpy/arrayobject.h>
#include <spglib.h>

static PyObject * get_spacegroup(PyObject *self, PyObject *args);
static PyObject * get_symmetry(PyObject *self, PyObject *args);
static PyObject * get_multiplicity(PyObject *self, PyObject *args);

static PyMethodDef functions[] = {
  {"spacegroup", get_spacegroup, METH_VARARGS, "International symbol"},
  {"symmetry", get_symmetry, METH_VARARGS, "Symmetry operations"},
  {"multiplicity", get_multiplicity, METH_VARARGS, "Number of symmetry operations"},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_spglib(void)
{
  Py_InitModule3("_spglib", functions, "C-extension for spglib\n\n...\n");
  return;
}

static PyObject * get_spacegroup(PyObject *self, PyObject *args)
{
  int i;
  double symprec;
  char symbol[26];
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOd", &lattice, &position, &atom_type, &symprec))
    return NULL;

  const double* lat = (double*)lattice->data;
  const double* pos = (double*)position->data;
  const int num_atom = position->dimensions[0];
  const long* typat_long = (long*)atom_type->data;

  int typat[num_atom];
  for (i = 0; i < num_atom; i++)
    typat[i] = (int)typat_long[i];

  const num_spg = spg_get_international(symbol, lat, pos, typat, num_atom, symprec);
  sprintf(symbol, "%s (%d)", symbol, num_spg);

  return PyString_FromString(symbol);
}

static PyObject * get_multiplicity(PyObject *self, PyObject *args)
{
  int i;
  double symprec;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOd", &lattice, &position, &atom_type, &symprec))
    return NULL;

  const double* lat = (double*)lattice->data;
  const double* pos = (double*)position->data;
  const int num_atom = position->dimensions[0];
  const long* types_long = (long*)atom_type->data;

  int types[num_atom];
  for (i = 0; i < num_atom; i++)
    types[i] = (int)types_long[i];
  
  const num_sym = spg_get_multiplicity(lat, pos, types, num_atom, symprec);

  return PyInt_FromLong((long) num_sym);
}

static PyObject * get_symmetry(PyObject *self, PyObject *args)
{
  int i, j, k;
  double symprec;
  PyArrayObject* lattice;
  PyArrayObject* position;
  PyArrayObject* rotation;
  PyArrayObject* translation;
  PyArrayObject* atom_type;
  if (!PyArg_ParseTuple(args, "OOOOOd", &rotation, &translation,
			&lattice, &position, &atom_type, &symprec))
    return NULL;

  const double* lat = (double*)lattice->data;
  const double* pos = (double*)position->data;
  const long* types_long = (long*)atom_type->data;
  const int num_atom = position->dimensions[0];
  long *rot_long = (long*)rotation->data;
  double* trans = (double*)translation->data;
  const int num_sym_from_array_size = rotation->dimensions[0];

  int rot[num_sym_from_array_size][3][3];
  
  int types[num_atom];
  for (i = 0; i < num_atom; i++)
    types[i] = (int)types_long[i];
  
  /* num_sym has to be larger than num_sym_from_array_size. */
  const int num_sym = spg_get_symmetry(rot, trans, num_sym_from_array_size,
				       lat, pos, types, num_atom, symprec);
  for (i = 0; i < num_sym; i++)
    for (j = 0; j < 3; j++)
      for (k = 0; k < 3; k++)
	rot_long[i*9+j*3+k] = (long)rot[i][j][k];

  return PyInt_FromLong((long) num_sym);
}


