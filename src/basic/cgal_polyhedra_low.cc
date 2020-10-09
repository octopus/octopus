
// Adapted from CGAL example (Author: Pierre Alliez) by Vladimir Fuka.

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


typedef struct {Polyhedron *poly; Tree *tree;} Polytree;

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
  
  void polyhedron_closest (const Polytree *ptree, const d3 *query, d3 *near){
    Point query_point(query->x,query->y,query->z);
    
    Point closest = ptree->tree->closest_point(query_point);
    
    near->x = closest.x();
    near->y = closest.y();
    near->z = closest.z();
  }

  bool polyhedron_point_inside(Polyhedron *polyhedron, Point *query) {
    // Construct AABB tree with a KdTree
    Tree tree(faces(*polyhedron).first, faces(*polyhedron).second, *polyhedron);
    tree.accelerate_distance_queries();
    // Initialize the point-in-polyhedron tester
    Point_inside inside_tester(tree);

    // Determine the side and return true if outside!
    return inside_tester(*query) == CGAL::ON_BOUNDED_SIDE || inside_tester(*query) == CGAL::ON_BOUNDARY;
  }

  void polyhedron_bbox(const Polytree *ptree, d3 *const min, d3 *const max){
    Bbox_3 bbox = ptree->tree->bbox();
    *min = {bbox.xmin(), bbox.ymin(), bbox.zmin()};
    *max = {bbox.xmax(), bbox.ymax(), bbox.zmax()};
  }
  
  void polyhedron_finalize(Polytree **pptree){
    delete (*pptree)->tree; (*pptree)->tree = NULL;
    delete (*pptree)->poly; (*pptree)->poly = NULL;
    delete *pptree; *pptree = NULL;
  }
 
}
