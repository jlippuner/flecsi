/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------------~~*/

//----------------------------------------------------------------------------//
//! \file
//! \date Initial file creation: Jul 26, 2016
//----------------------------------------------------------------------------//

#include <cstddef>
#include <cstdint>

#include "flecsi/data/data.h"

clog_register_tag(runtime_driver);

namespace flecsi {
namespace execution {

//----------------------------------------------------------------------------//
// Implementation of FleCSI runtime driver task.
//----------------------------------------------------------------------------//

void
remap_shared_entities()
{
  // TODO: Is this superseded by index_map/reverse_index_map?
  auto& flecsi_context = context_t::instance();
  const int my_color = flecsi_context.color();

  for (auto coloring_info_pair : flecsi_context.coloring_info_map()) {
    auto index_space = coloring_info_pair.first;
    auto coloring_info = coloring_info_pair.second;

    auto &my_coloring_info = flecsi_context.coloring_info(index_space).at(
      my_color);
    auto index_coloring = flecsi_context.coloring(index_space);

    size_t index = 0;
    for (auto shared : index_coloring.shared) {
//      clog_rank(warn, 0) << "myrank: " << my_color
//                         << " shared id: " << shared.id
//                         << ", rank: " << shared.rank
//                         << ", offset: " << shared.offset
//                         << ", index: " << index << std::endl;
      for (auto peer : shared.shared) {
        MPI_Send(&index, 1, MPI_UNSIGNED_LONG_LONG, peer, 77, MPI_COMM_WORLD);
      }
      index++;
    }

    std::vector<size_t> ghost_index(index_coloring.ghost.size());

    MPI_Status status;
    std::set<flecsi::coloring::entity_info_t> new_ghost;

    for (auto ghost : index_coloring.ghost) {
      MPI_Recv(&index, 1, MPI_UNSIGNED_LONG_LONG,
               ghost.rank, 77, MPI_COMM_WORLD, &status);
      new_ghost.insert(
        flecsi::coloring::entity_info_t(ghost.id, ghost.rank, index, {}));
    }
//    for (auto ghost : index_coloring.ghost) {
//      clog_rank(warn, 1) << "myrank: " << my_color
//                         << " old ghost id: " << ghost.id
//                         << ", rank: " << ghost.rank
//                         << ", offset: " << ghost.offset
//                         << std::endl;
//    }
//    for (auto ghost : new_ghost) {
//      clog_rank(warn, 1) << "myrank: " << my_color
//                         << " new ghost id: " << ghost.id
//                         << ", rank: " << ghost.rank
//                         << ", offset: " << ghost.offset
//                         << std::endl;
//    }
    context_t::instance().coloring(index_space).ghost.swap(new_ghost);
  }
}

void
runtime_driver(
  int argc,
  char ** argv
)
{
  {
  clog_tag_guard(runtime_driver);
  clog(info) << "In MPI runtime driver" << std::endl;
  }

#if defined(FLECSI_ENABLE_SPECIALIZATION_TLT_INIT)
  {
  clog_tag_guard(runtime_driver);
  clog(info) << "Executing specialization tlt task" << std::endl;
  }

  // Execute the specialization driver.
  specialization_tlt_init(argc, argv);
#endif // FLECSI_ENABLE_SPECIALIZATION_TLT_INIT

  //--------------------------------------------------------------------------//
  // Invoke callbacks for entries in the client registry.
  //
  // NOTE: This needs to be called before the field registry below because
  //       The client callbacks register field callbacks with the field
  //       registry.
  //--------------------------------------------------------------------------//

  auto & client_registry =
    flecsi::data::storage_t::instance().client_registry();

  for(auto & c: client_registry) {
    for(auto & d: c.second) {
      d.second.second(d.second.first);
    } // for
  } // for

  //--------------------------------------------------------------------------//
  // Invoke callbacks for entries in the field registry.
  //--------------------------------------------------------------------------//

  auto & field_registry =
    flecsi::data::storage_t::instance().field_registry();

  for(auto & c: field_registry) {
    for(auto & f: c.second) {
      f.second.second(f.first, f.second.first);
    } // for
  } // for

  auto& flecsi_context = context_t::instance();
  for (auto fi : flecsi_context.registered_fields()) {
    flecsi_context.put_field_info(fi);
  }

  remap_shared_entities();

  // Setup maps from mesh to compacted (local) index space and vice versa
  //
  // This depends on the ordering of the BLIS data structure setup.
  // Currently, this is Exclusive - Shared - Ghost.

  for(auto is: flecsi_context.coloring_map()) {
    std::map<size_t, size_t> _map;
    size_t counter(0);

    for(auto index: is.second.exclusive) {
      _map[counter++] = index.id;
    } // for

    for(auto index: is.second.shared) {
      _map[counter++] = index.id;
    } // for

    for(auto index: is.second.ghost) {
      _map[counter++] = index.id;
    } // for

    flecsi_context.add_index_map(is.first, _map);
  } // for

  flecsi_context.advance_state();
  // Call the specialization color initialization function.
#if defined(FLECSI_ENABLE_SPECIALIZATION_SPMD_INIT)
  specialization_spmd_init(argc, argv);
#endif // FLECSI_ENABLE_SPECIALIZATION_SPMD_INIT
  
  flecsi_context.advance_state();

  // Execute the user driver.
  driver(argc, argv);

} // runtime_driver

} // namespace execution 
} // namespace flecsi

/*~------------------------------------------------------------------------~--*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
