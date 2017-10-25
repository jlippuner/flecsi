#include <iostream>
#include "flecsi/topology/structured_mesh_topology.h"
#include "flecsi/data/data.h"

using namespace std;
using namespace flecsi;
using namespace topology;
using namespace flecsi::data;


class Vertex : public structured_mesh_entity_t<0, 1>{
public:
  Vertex(){}

  Vertex(structured_mesh_topology_base_t &){}
};

class Edge : public structured_mesh_entity_t<1, 1>{
public:
  Edge(){}

  Edge(structured_mesh_topology_base_t &){}
};

class Face : public structured_mesh_entity_t<2, 1>{
public:

  Face(){}

  Face(structured_mesh_topology_base_t &){}
};

class Cell : public structured_mesh_entity_t<3, 1>{
public:

  Cell(){}

  Cell(structured_mesh_topology_base_t &){}
};


class TestMesh2dType{
public:
  static constexpr size_t num_dimensions = 2;
  static constexpr size_t num_domains = 1;

  static constexpr std::array<size_t,num_dimensions> lower_bounds = {0,0};
  static constexpr std::array<size_t,num_dimensions> upper_bounds = {1,1};

  using entity_types = std::tuple<
  std::pair<domain_<0>, Vertex>,
  std::pair<domain_<0>, Edge>,
  std::pair<domain_<0>, Face>,
  std::pair<domain_<0>, Cell>>;

};

constexpr std::array<size_t,TestMesh2dType::num_dimensions> TestMesh2dType::lower_bounds;
constexpr std::array<size_t,TestMesh2dType::num_dimensions> TestMesh2dType::upper_bounds;

using id_vector_t = std::vector<size_t>;
using TestMesh = structured_mesh_topology_t<TestMesh2dType>;

int main(int, char**) {
  TestMesh mesh;

  printf("Num vertices: %lu\n", mesh.num_entities<0,0>());
  printf("Num edges: %lu\n", mesh.num_entities<1,0>());
  printf("Num faces: %lu\n", mesh.num_entities<2,0>());
  printf("Num cells: %lu\n", mesh.num_entities<3,0>());

  auto idxs = mesh.get_indices<1>(9);
  printf("Indices = [%lu, %lu]\n", idxs[0], idxs[1]);

  return EXIT_SUCCESS;
} // TEST
