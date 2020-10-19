/*
Copyright (C) 2014 Adapted from CGAL example (Author: Pierre Alliez) by Vladimir Fuka
Copyright (C) 2020 Heiko Appel

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301, USA.
*/

#include <config.h>
// undef CC here because that is used in CGAL...
#undef CC

#ifdef HAVE_CGAL
#include <iostream>
#include <fstream>
#include <list>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/algorithm.h>
#include <CGAL/Side_of_triangle_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::Point_3 Point;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef CGAL::Side_of_triangle_mesh<Polyhedron, K> Point_inside;
typedef CGAL::Bbox_3 Bbox_3;

typedef struct {double x,y,z;} d3;

using std::cout;
using std::endl;

extern "C" int debuglevel;

extern "C" {

  void polyhedron_from_file (Polyhedron **poly, const char *fname, int verbose, int * const ierr){
    Polyhedron *polyhedron = new Polyhedron;

    std::ifstream in(fname);

    if (verbose) {cout << " Reading file " << fname << " " << endl;}

    try {
      in >> *polyhedron;
    }
    catch(...) {
      *ierr = 2;
      return;
    }

    if (verbose) {
      cout << " facets: " << polyhedron->size_of_facets() << endl;
      cout << " halfedges: " << polyhedron->size_of_halfedges() << endl;
      cout << " vertices: " << polyhedron->size_of_vertices() << endl;
    }

    if (polyhedron->size_of_facets()==0 ||
        polyhedron->size_of_halfedges()==0 ||
        polyhedron->size_of_vertices()==0){
          *ierr = 1;
          return;
        };

    *poly = polyhedron;
    *ierr = 0;
  }

  void polyhedron_build_AABB_tree(Tree **tree_out, Polyhedron **polyhedron){
    // Construct AABB tree with a KdTree
    Tree *tree = new Tree(faces(**polyhedron).first, faces(**polyhedron).second, **polyhedron);
    tree->accelerate_distance_queries();
    *tree_out = tree;
  }

  bool polyhedron_point_inside(Tree **tree, Point *query) {
    // Initialize the point-in-polyhedron tester
    Point_inside inside_tester(**tree);

    // Determine the side and return true if inside or on the boundary.
    return inside_tester(*query) == CGAL::ON_BOUNDED_SIDE || inside_tester(*query) == CGAL::ON_BOUNDARY;
  }

  void polyhedron_finalize_AABB_tree(Tree **tree){
    delete *tree; *tree = NULL;
  }

  void polyhedron_finalize_polyhedron(Polyhedron **polyhedron){
    delete *polyhedron; *polyhedron = NULL;
  }
}
#endif
