.. spglib for ASE documentation master file, created by sphinx-quickstart on Mon Feb  2 14:53:46 2009.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

pyspglib for ASE
==========================================

This is document for Pyspglib for ASE (Atomic Simulation Environment). 
Pyspglib is the python module to use spglib library. 

How to build spglib python module
=================================
The C sources of spglib and interface for the python C/API are compiled. The development environment for python and gcc are required before starting to build.

1. Go to the :file:`python/ase` directory
2. Type the command::

    % python setup.py --home=<my-directory>

The :file:`{<my-directory>}` is possibly current directory, :file:`.`.

3. You will find :file:`_spglib.so` in the :file:`{my-directory}/lib/python`. Add the directory where you copy :file:`_spglib.so` and :file:`spglib.py` to :envvar:`$PYTHONPATH`.

Example
=============
To run this example, you need to copy :file:`_spglib.so` and :file:`spglib.py` to a current directory or the directory where the :envvar:`$PYTHONPATH` is set.

.. literalinclude:: example_ase.py

How to use it
=============
1. Import spglib::

    import spglib

2. Call the methods with ASE Atoms object.

Methods
=======

get_spacegroup
--------------
::

    get_spacegroup(atoms, precision=1e-5)

``atoms`` is the object of ASE Atoms class. ``precision`` is the float variable, which is used as tolerance in symmetry search. The unit is about fractional coordinates. 

International space group symbol and the number are obtained as a string.

get_symmetry
------------
::

    get_symmetry(Atoms, precition=1e-5)

``atoms`` is the object of ASE Atoms class. ``precision`` is the float variable, which is used as tolerance in symmetry search. The unit is about fractional coordinates. 

Symmetry operations are obtained as a dictionary. The key ``rotation`` contains a numpy array of integer, which is "number of symmetry operations" x "3x3 matrices". The key ``translation`` contains a numpy array of float, which is "number of symmetry operations" x "vectors". The orders of the rotation matrices and the translation vectors correspond with each other, e.g. , the second symmetry operation is organized by the second rotation matrix and second translation vector in the respective arrays. The operations are applied for the fractional coordinates (not for Cartesian coordinates).

The rotation matrix and translation vector are used as follows::

    new_vector[3x1] = rotation[3x3] * vector[3x1] + translation[3x1]

The three values in the vector are given for the a, b, and c axes, respectively.
