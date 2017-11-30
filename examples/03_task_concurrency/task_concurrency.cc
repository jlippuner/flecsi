/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------------~~*/

#include "flecsi/execution/context.h"
#include "flecsi/execution/execution.h"
#include "flecsi/data/sparse_accessor.h"
#include "flecsi/supplemental/coloring/add_colorings.h"
#include "flecsi/supplemental/mesh/test_mesh_2d.h"

#include <thread>
#include <chrono>

namespace flecsi {
namespace execution {

//----------------------------------------------------------------------------//
// Type definitions
//----------------------------------------------------------------------------//

using point_t = flecsi::supplemental::point_t;
using vertex_t = flecsi::supplemental::vertex_t;
using cell_t = flecsi::supplemental::cell_t;
using mesh_t = flecsi::supplemental::test_mesh_2d_t;

using coloring_info_t = flecsi::coloring::coloring_info_t;
using adjacency_info_t = flecsi::coloring::adjacency_info_t;

template<
  size_t PS
>
using mesh = data_client_handle__<mesh_t, PS>;

template<
  size_t EP,
  size_t SP,
  size_t GP
>
using field = dense_accessor<double, EP, SP, GP>;

template<
  size_t EP,
  size_t SP,
  size_t GP
>
using sparse_field = sparse_accessor<double, EP, SP, GP>;

template<
  size_t PS
>
using scalar_int = global_accessor<int, PS>;


//----------------------------------------------------------------------------//
// Variable registration
//----------------------------------------------------------------------------//

flecsi_register_data_client(mesh_t, clients, m);

flecsi_register_field(mesh_t, data, A, double, dense, 1, index_spaces::cells);
flecsi_register_field(mesh_t, data, B, double, dense, 1, index_spaces::cells);
flecsi_register_field(mesh_t, data, C, double, dense, 1, index_spaces::cells);
flecsi_register_field(mesh_t, data, glob, int, global, 1);

//----------------------------------------------------------------------------//
// Initialization task
//----------------------------------------------------------------------------//

void initialize_mesh(mesh<wo> m) {

  auto & context { execution::context_t::instance() };

  auto & vertex_map { context.index_map(index_spaces::vertices) };
  auto & reverse_vertex_map
    { context.reverse_index_map(index_spaces::vertices) };
  auto & cell_map { context.index_map(index_spaces::cells) };

  std::vector<vertex_t *> vertices;

  const size_t width { 8 };
  const double dt { 1.0/width };

  for(auto & vm: vertex_map) {
    const size_t mid { vm.second };
    const size_t row { mid/(width+1) };
    const size_t column { mid%(width+1) };

    point_t p { column*dt, row*dt };
    vertices.push_back(m.make<vertex_t>(p));
  } // for

  size_t count{0};
  for(auto & cm: cell_map) {
    const size_t mid { cm.second };

    const size_t row { mid/width };
    const size_t column { mid%width };

    const size_t v0 { (column    ) + (row    ) * (width + 1) };
    const size_t v1 { (column + 1) + (row    ) * (width + 1) };
    const size_t v2 { (column + 1) + (row + 1) * (width + 1) };
    const size_t v3 { (column    ) + (row + 1) * (width + 1) };

    const size_t lv0 { reverse_vertex_map[v0] };
    const size_t lv1 { reverse_vertex_map[v1] };
    const size_t lv2 { reverse_vertex_map[v2] };
    const size_t lv3 { reverse_vertex_map[v3] };

    auto c { m.make<cell_t>() };
    m.init_cell<0>(c, { vertices[lv0], vertices[lv1],
      vertices[lv2], vertices[lv3] });
  } // for

  m.init<0>();

} // initialize_mesh

flecsi_register_task(initialize_mesh, loc, single);

//----------------------------------------------------------------------------//
// Initialize pressure
//----------------------------------------------------------------------------//

void compute_a(mesh<ro> m, field<rw, rw, ro> a) {
  int cnt = 0;
  for(auto c: m.cells(owned)) {
    std::this_thread::sleep_for(std::chrono::milliseconds(10));

    a(c) = (double)cnt;

    if (cnt % 50 == 0)
      printf("Computing A(%lu) = %f...\n", c, a(c));
    ++cnt;
  } // for

} // compute_a

flecsi_register_task(compute_a, loc, single);

void sparse(mesh<ro> m, sparse_field<rw, rw, ro> a) {

} // sparse

flecsi_register_task(sparse, loc, single);

void compute_b(mesh<ro> m, field<rw, rw, ro> b) {
  int cnt = 0;
  for(auto c: m.cells(owned)) {
    std::this_thread::sleep_for(std::chrono::milliseconds(5));

    b(c) = (double)(cnt * cnt);

    if (cnt % 50 == 0)
      printf("Computing B(%lu) = %f...\n", c, b(c));

    ++cnt;
  } // for

} // compute_b

flecsi_register_task(compute_b, loc, single);

void init_c(mesh<ro> m, field<rw, rw, ro> c) {
  int cnt = 0;
  for(auto cell: m.cells(owned)) {
    std::this_thread::sleep_for(std::chrono::milliseconds(1));

    if (cnt % 50 == 0)
      printf("Initializing C...\n");

    c(cell) = 0.0;
    ++cnt;
  } // for

} // init_c

flecsi_register_task(init_c, loc, single);

void compute_c(mesh<ro> m, field<ro, ro, ro> a, field<ro, ro, ro> b,
    field<rw, rw, ro> c) {
  int cnt = 0;
  for(auto cl: m.cells(owned)) {
    std::this_thread::sleep_for(std::chrono::milliseconds(2));

    c(cl) = a(cl) + b(cl);

    if (cnt % 50 == 0)
      printf("[%i] Computing C (%f = %f + %f)...\n", cnt, c(cl), a(cl), b(cl));
    ++cnt;
  } // for

} // compute_c

flecsi_register_task(compute_c, loc, single);

//----------------------------------------------------------------------------//
// Update pressure
//----------------------------------------------------------------------------//

void global_task(mesh<ro> m, scalar_int<ro> glob) {
//  printf("Executing global task (glob = %i)\n", glob);
  std::this_thread::sleep_for(std::chrono::milliseconds(50));
  //++glob.data();
//  printf("Global task done (glob = %i)\n", glob);
} // global_task

flecsi_register_task(global_task, loc, single);


//----------------------------------------------------------------------------//
// Top-Level Specialization Initialization
//----------------------------------------------------------------------------//

void specialization_tlt_init(int argc, char ** argv) {
  printf("start specialization_tlt_init\n");

  coloring_map_t map { index_spaces::vertices, index_spaces::cells };

  flecsi_execute_mpi_task(add_colorings, map);

  auto & context { execution::context_t::instance() };

  auto & vinfo { context.coloring_info(index_spaces::vertices) };
  auto & cinfo { context.coloring_info(index_spaces::cells) };

  adjacency_info_t ai;
  ai.index_space = index_spaces::cells_to_vertices;
  ai.from_index_space = index_spaces::cells;
  ai.to_index_space = index_spaces::vertices;
  ai.color_sizes.resize(cinfo.size());

  for(auto & itr : cinfo){
    size_t color{itr.first};
    const coloring_info_t & ci = itr.second;
    ai.color_sizes[color] = (ci.exclusive + ci.shared + ci.ghost) * 4;
  } // for

  context.add_adjacency(ai);

  printf("end specialization_tlt_init\n");
} // specialization_tlt_init

//----------------------------------------------------------------------------//
// SPMD Specialization Initialization
//----------------------------------------------------------------------------//

void specialization_spmd_init(int argc, char ** argv) {
  printf("start specialization_spmd_init\n");
  auto mh = flecsi_get_client_handle(mesh_t, clients, m);
  flecsi_execute_task(initialize_mesh, single, mh);
  printf("end specialization_spmd_init\n");
} // specialization_spmd_ini

//----------------------------------------------------------------------------//
// User driver.
//----------------------------------------------------------------------------//

void run_something() {
  auto mh2 = flecsi_get_client_handle(mesh_t, clients, m);
  auto ah2 = flecsi_get_handle(mh2, data, A, double, dense, 0);
  auto bh2 = flecsi_get_handle(mh2, data, B, double, dense, 0);
  auto ch2 = flecsi_get_handle(mh2, data, C, double, dense, 0);

  flecsi_execute_task(compute_a, single, mh2, ah2);
  flecsi_execute_task(compute_b, single, mh2, bh2);
  flecsi_execute_task(init_c, single, mh2, ch2);
  flecsi_execute_task(compute_c, single, mh2, ah2, bh2, ch2);
}

#define SIMPLE

#ifdef SIMPLE
void driver(int argc, char ** argv) {
  printf("start driver\n");

  auto mh = flecsi_get_client_handle(mesh_t, clients, m);
  auto ah = flecsi_get_handle(mh, data, A, double, dense, 0);
  auto bh = flecsi_get_handle(mh, data, B, double, dense, 0);
  auto ch = flecsi_get_handle(mh, data, C, double, dense, 0);
  auto globh = flecsi_get_handle(mh, data, glob, int, global, 0);

  flecsi_execute_task(compute_a, single, mh, ah);
  flecsi_execute_task(compute_b, single, mh, bh);

  flecsi_execute_task(init_c, single, mh, ch);
  flecsi_execute_task(compute_c, single, mh, ah, bh, ch);

  flecsi_execute_task(global_task, single, mh, globh);

  printf("end driver\n");
} // driver

#else

void driver(int argc, char ** argv) {
  printf("start driver\n");

  auto mh = flecsi_get_client_handle(mesh_t, clients, m);
  auto ah = flecsi_get_handle(mh, data, A, double, dense, 0);
  auto bh = flecsi_get_handle(mh, data, B, double, dense, 0);
  auto ch = flecsi_get_handle(mh, data, C, double, dense, 0);
  auto globh = flecsi_get_handle(mh, data, glob, int, global, 0);

  bool do_global = true;

  if (do_global)
    flecsi_execute_task(global_task, single, mh, globh);

  for (int i = 0; i < 3; ++i) {
    flecsi_execute_task(compute_a, single, mh, ah);
    flecsi_execute_task(compute_b, single, mh, bh);

    if (i == 2)
      flecsi_execute_task(init_c, single, mh, ch);
    else
      flecsi_execute_task(global_task, single, mh, globh);

    int j = 0;
    do {
      flecsi_execute_task(compute_c, single, mh, ah, bh, ch);
      ++j;
    } while (j < i);
  }

  run_something();

  int i = 0;
  while (i < 3) {
    flecsi_execute_task(init_c, single, mh, ch);
    run_something();
    ++i;
  }

  switch (i) {
  case 0:
    flecsi_execute_task(compute_a, single, mh, ah);
    break;
  case 1:
    flecsi_execute_task(compute_b, single, mh, bh);
//    break;
  case 2:
    if (do_global) {
      flecsi_execute_task(compute_c, single, mh, ah, bh, ch);
    } else if (i == 3) {
      flecsi_execute_task(init_c, single, mh, ch);
    } else {
      for (int k = 0; k < 3; ++k) {
        flecsi_execute_task(global_task, single, mh, globh);
      }
    }
    break;
  default:
    flecsi_execute_task(global_task, single, mh, globh);
  }

  if (do_global)
    int a = 0;
  else
    flecsi_execute_task(compute_c, single, mh, ah, bh, ch);

  if (do_global)
    run_something();

  flecsi_execute_task(compute_c, single, mh, ah, bh, ch);

  flecsi_execute_task(global_task, single, mh, globh);

  printf("end driver\n");
} // driver

#endif // SIMPLE


} // namespace execution
} // namespace flecsi

int main(int argc, char ** argv) {
  printf("start main\n");

#ifdef GASNET_CONDUIT_MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    // If you fail this assertion, then your version of MPI
    // does not support calls from multiple threads and you
    // cannot use the GASNet MPI conduit
    if (provided < MPI_THREAD_MULTIPLE)
      printf("ERROR: Your implementation of MPI does not support "
           "MPI_THREAD_MULTIPLE which is required for use of the "
           "GASNet MPI conduit with the Legion-MPI Interop!\n");
    assert(provided == MPI_THREAD_MULTIPLE);
#else
std::cout <<"MPI MPI  MPI"<< std::endl;
   MPI_Init(&argc, &argv);
#endif

  std::cout <<"task concurrency example"<<std::endl;

  // Call FleCSI runtime initialization
  auto retval = flecsi::execution::context_t::instance().initialize(argc, argv);

#ifndef GASNET_CONDUIT_MPI
  MPI_Finalize();
#endif

  printf("end main\n");
  return retval;
} // main

/*~------------------------------------------------------------------------~--*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
