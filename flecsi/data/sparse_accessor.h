/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_sparse_accessor_h
#define flecsi_sparse_accessor_h

#include "flecsi/data/sparse_data_handle.h"

//----------------------------------------------------------------------------//
//! @file
//! @date Initial file creation: Nov 14, 2017
//----------------------------------------------------------------------------//

namespace flecsi {

//----------------------------------------------------------------------------//
//! The sparse_accessor_base_t type provides an empty base type for 
//! compile-time
//! identification of data handle objects.
//!
//! @ingroup data
//----------------------------------------------------------------------------//

struct sparse_accessor_base_t {};

//----------------------------------------------------------------------------//
//! The sparse_accessor type captures information about permissions
//! and specifies a data policy.
//!
//! @tparam T                     The data type referenced by the handle.
//! @tparam EXCLUSIVE_PERMISSIONS The permissions required on the exclusive
//!                               indices of the index partition.
//! @tparam SHARED_PERMISSIONS    The permissions required on the shared
//!                               indices of the index partition.
//! @tparam GHOST_PERMISSIONS     The permissions required on the ghost
//!                               indices of the index partition.
//! @tparam DATA_POLICY           The data policy for this handle type.
//!
//! @ingroup data
//----------------------------------------------------------------------------//

template<
  typename T,
  size_t EXCLUSIVE_PERMISSIONS,
  size_t SHARED_PERMISSIONS,
  size_t GHOST_PERMISSIONS
>
struct sparse_accessor : public sparse_accessor_base_t {
  using handle_t = 
    sparse_data_handle__<
      T,
      EXCLUSIVE_PERMISSIONS,
      SHARED_PERMISSIONS,
      GHOST_PERMISSIONS
    >;

  using offset_t = typename handle_t::offset_t;
  using entry_value_t = typename handle_t::entry_value_t;

  //--------------------------------------------------------------------------//
  //! Copy constructor.
  //--------------------------------------------------------------------------//

  sparse_accessor(const sparse_data_handle__<T, 0, 0, 0>& h)
  : handle(reinterpret_cast<const handle_t&>(h)){

  }

  T &
  operator () (
    size_t index,
    size_t entry
  )
  {
    assert(index < handle.num_total_ && 
      "sparse accessor: index out of bounds");

    const offset_t& oi = handle.offsets[index];

    entry_value_t * start = handle.entries + oi.start();
    entry_value_t * end = start + oi.count();

    entry_value_t * itr =
      std::lower_bound(start, end, entry_value_t(entry),
        [](const entry_value_t & k1, const entry_value_t & k2) -> bool {
          return k1.entry < k2.entry;
        });

    assert(itr != end && 
      "sparse accessor: unmapped entry");

    return itr->value;
  } // operator ()

  void dump() const{
    for(size_t i = 0; i < handle.num_total_; ++i){
      const offset_t& offset = handle.offsets[i];
      std::cout << "index: " << i << std::endl;
      for(size_t j = 0; j < offset.count(); ++j){
        size_t k = offset.start() + j;
        std::cout << "  " << handle.entries[k].entry << 
          " = " << handle.entries[k].value << std::endl;
      }
    }
  }

  handle_t handle;  
};

} // namespace flecsi

#endif // flecsi_sparse_accessor_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
