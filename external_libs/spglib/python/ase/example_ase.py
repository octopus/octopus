#!/usr/bin/evn python

import sys
from ase import *
import spglib as spg
import numpy as npy

Si = Atoms(symbols='Si8',
             cell=[(4,0,0),
                   (0,4,0),
                   (0,0,4)],
             scaled_positions=[(0, 0, 0),
                        (0, 0.5, 0.5),
                        (0.5, 0, 0.5),
                        (0.5, 0.5, 0),
                        (0.25, 0.25, 0.25),
                        (0.25, 0.75, 0.75),
                        (0.75, 0.25, 0.75),
                        (0.75, 0.75, 0.25)],
             pbc=True)

SiO2 = Atoms(symbols='Si2O4',
             cell=[(4,0,0),
                   (0,4,0),
                   (0,0,3)],
             scaled_positions=[(0, 0, 0),
                               (0.5, 0.5, 0.5),
                               (0.3, 0.3, 0.0),
                               (0.7, 0.7, 0.0),
                               (0.2, 0.8, 0.5),
                               (0.8, 0.2, 0.5)],
             pbc=True)


# For VASP case
# import ase.io.vasp as vasp
# bulk = vasp.read_vasp(sys.argv[1])

print "Spacegroup of Si is ", spg.get_spacegroup(Si)
print ""
print "Spacegroup of SiO2 is ", spg.get_spacegroup(SiO2)
print ""
print "Symmetry operations of SiO2 unitcell is:"
print ""
symmetry = spg.get_symmetry(SiO2)
for i in range(symmetry['rotation'].shape[0]):
  print "--------------- %4d ---------------" % (i+1)
  rot = symmetry['rotation'][i]
  trans = symmetry['translation'][i]
  print "rotation:"
  for x in rot:
    print "   [%2d %2d %2d]" % (x[0], x[1], x[2])
  print "translation:"
  print "   (%8.5f %8.5f %8.5f)" % (trans[0], trans[1], trans[2])
