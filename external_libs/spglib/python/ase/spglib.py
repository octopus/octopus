"""
Spglib interface for ASE
"""

import _spglib as spg
import numpy as npy

def get_symmetry(bulk, symprec = 1e-5):
  """
  Return symmetry operations as hash.
  Hash key 'rotation' gives the numpy integer array
  of the rotation matrices for scaled positions
  Hash key 'translation' gives the numpy float64 array
  of the translation vectors in scaled positions
  """

  # Atomic positions have to be specified by scaled positions for spglib.
  positions = bulk.get_scaled_positions()
  cell = bulk.get_cell().transpose().copy()
  numbers = bulk.get_atomic_numbers()

  # Get number of symmetry operations and allocate symmetry operations
  multi = spg.multiplicity(cell, positions, numbers, symprec)
  rotation = npy.zeros((multi, 3, 3), dtype=int)
  translation = npy.zeros((multi, 3))

  # Get symmetry operations
  num_sym = spg.symmetry(rotation, translation, cell, positions, numbers, symprec)

  return {'rotation': rotation, 'translation': translation}

def get_spacegroup(bulk, symprec = 1e-5):
  """
  Return space group in international table symbol and number
  as a string.
  """
  # Atomic positions have to be specified by scaled positions for spglib.
  return spg.spacegroup(bulk.get_cell().transpose().copy(), bulk.get_scaled_positions(), bulk.get_atomic_numbers(), symprec)
