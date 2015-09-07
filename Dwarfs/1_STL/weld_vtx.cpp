#include <iostream>
#include <vector>
#include <algorithm>
#include "parse_arguments.h"

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

  // find index of each vertex in the list of unique vertices
  auto result_end = my_find_all_lower_bounds(vertices.begin(), vertices.end(),
      input.begin(), input.end(),
      indices.begin());
  indices.erase(result_end, indices.end());
  return indices;
}

int main(int argc, char *argv[])
{
  using namespace parse_arguments;

  auto n_vertices = 1024;
  auto params = pack(argc, argv,
      param("number of vertices", n_vertices, "n", "")
      );

  cerr << params.parse_all();

  return 0;
}
