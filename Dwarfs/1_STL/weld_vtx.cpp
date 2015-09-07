#include <iostream>
#include <vector>
#include <algorithm>
#include "timer.h"
#include "parse_arguments.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

using Kernel   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Delaunay = CGAL::Delaunay_triangulation_2<Kernel>;
using Point    = Kernel::Point_2;


using namespace std;

using Real = double;
using Vertex = tuple<Real,Real>;

template<class ForwardIterator,
         class InputIterator,
         class OutputIterator>
OutputIterator my_find_all_lower_bounds(
    ForwardIterator haystack_begin,
    ForwardIterator haystack_end,
    InputIterator   needles_begin,
    InputIterator   needles_end,
    OutputIterator  result)
{
  return transform(needles_begin, needles_end, result,
      [=](auto& needle)
  {
    auto iter = lower_bound(haystack_begin, haystack_end, needle);
    return distance(haystack_begin, iter);
  });
}
      

vector<size_t> weld_vertices(vector<Vertex> input)
{
  auto vertices = input;
  vector<size_t> indices(input.size());

  // sort vertices to bring duplicated together
  sort(vertices.begin(), vertices.end());

  // find unique vertcies and erase redundancies
  auto redundant_begin = unique(vertices.begin(), vertices.end());
  vertices.erase(redundant_begin, vertices.end());
  cout << "\tUnique vertices: " << vertices.size() << endl;

  // find index of each vertex in the list of unique vertices
  my_find_all_lower_bounds(vertices.begin(), vertices.end(),
      input.begin(), input.end(),
      indices.begin());
  return indices;
}

vector<Vertex> generate_triangulation(int n_vertices)
{
  vector<Point> points(n_vertices);
  for (int i = 0; i < n_vertices; ++i)
  {
    points[i] = Point(drand48(), drand48());
  }

  Delaunay triangulation;
  triangulation.insert(points.begin(), points.end());

  vector<Vertex> vertices;
  for (auto face = triangulation.finite_faces_begin(); face != triangulation.finite_faces_end(); ++face)
  {
    auto triangle = triangulation.triangle(face);
#if 0
    cout << "Triangle:\t" << triangle << endl;
    cout << "Vertex 0:\t" << triangle[0][1] << endl;
    cout << "Vertex 1:\t" << triangle[1] << endl;
    cout << "Vertex 2:\t" << triangle[2] << endl;
    cout << " ------\n";
#endif
    vertices.emplace_back(triangle[0][0],triangle[0][1]);
    vertices.emplace_back(triangle[1][0],triangle[1][1]);
    vertices.emplace_back(triangle[2][0],triangle[2][1]);
  }

  return vertices;
}

int main(int argc, char *argv[])
{
  srand48((size_t)main);
  using namespace parse_arguments;

  auto n_vertices = 1024;
  auto params = pack(argc, argv,
      param("number of vertices", n_vertices, "n", "")
      );

  cout << params.parse_all();

  cout << "Generating triangulation...\n";
  auto vertices = generate_triangulation(n_vertices);
 
  { 
    Timer _timer("weld_vertices",Timer::verbose_destructor{});
    cout << "Welding vertices... \n";
    cout << "\tNumber of vertices port triangulation: " << vertices.size() << endl;
    _timer.tbeg();
    auto indices = weld_vertices(vertices);
    _timer.tend();
    _timer.tbeg();
    indices = weld_vertices(vertices);
    _timer.tend();
    _timer.tbeg();
    indices = weld_vertices(vertices);
    _timer.tend();
    cout << "\tNumber of indices: " << indices.size() << endl;
  }

  return 0;
}
