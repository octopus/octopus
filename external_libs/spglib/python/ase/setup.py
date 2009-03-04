from distutils.core import setup, Extension

extension = Extension('_spglib',
                      include_dirs = ['../../src'],
                      sources = ['_spglib.c',
                                 '../../src/bravais_art.c',
                                 '../../src/bravais.c',
                                 '../../src/cell.c',
                                 '../../src/debug.c',
                                 '../../src/mathfunc.c',
                                 '../../src/pointgroup.c',
                                 '../../src/primitive.c',
                                 '../../src/spacegroup.c',
                                 '../../src/spacegroup_database.c',
                                 '../../src/spacegroup_data.c',
                                 '../../src/spglib.c',
                                 '../../src/symmetry.c'])

setup (name = 'spglib',
       version = '0.5.0',
       description = 'This is the spglib module.',
       author = 'Atsushi Togo',
       author_email = 'atz.togo@gmail.com',
       url = 'http://spglib.sourceforge.net/',
       ext_modules = [extension])
