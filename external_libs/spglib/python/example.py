#!/usr/bin/evn python
from numpy import *
import libpyspg as spg

lattice = array([[4.,0.,0.], [0.,4.,0.], [0.,0.,3.]])
position = array([
	[0.0, 0.0, 0.0],
	[0.5, 0.5, 0.5],
	[0.3, 0.3, 0.0],
	[0.7, 0.7, 0.0],
	[0.2, 0.8, 0.5],
	[0.8, 0.2, 0.5]
	])
atom_type = array([1,1,2,2,2,2])

spg.spacegroup(atom_type.size, lattice, position, atom_type, 1e-5)
