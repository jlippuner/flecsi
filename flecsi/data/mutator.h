/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_mutator_h
#define flecsi_mutator_h

#include "flecsi/data/mutator_handle.h"

//----------------------------------------------------------------------------//
//! @file
//! @date Initial file creation: Nov 14, 2017
//----------------------------------------------------------------------------//

namespace flecsi {

//----------------------------------------------------------------------------//
//! The mutator_base_t type provides an empty base type for 
//! compile-time
//! identification of data handle objects.
//!
//! @ingroup data
//----------------------------------------------------------------------------//

struct mutator_base_t {};

//----------------------------------------------------------------------------//
//! The mutator type captures information about permissions
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
  typename T
>
struct mutator : public mutator_base_t {
  using handle_t = mutator_handle__<T>;
  using offset_t = typename handle_t::offset_t;
  using entry_value_t = typename handle_t::entry_value_t;
  using erase_set_t = typename handle_t::erase_set_t;

  //--------------------------------------------------------------------------//
  //! Copy constructor.
  //--------------------------------------------------------------------------//

  mutator(const mutator_handle__<T>& h)
  : h_(h){

  }

  T &
  operator () (
    size_t index,
    size_t entry
  )
  {
    assert(h_.offsets_ && "uninitialized mutator");
    assert(index < h_.num_entries_);

    offset_t& offset = h_.offsets_[index]; 

    size_t n = offset.count();

    if(n >= h_.num_slots_) {
      if(index < h_.num_exclusive_){
        (*h_.num_exclusive_insertions)++;
      }
      
      return h_.spare_map_->emplace(index,
        entry_value_t(entry))->second.value;
    } // if

    entry_value_t * start = h_.entries_ + index * h_.num_slots_;
    entry_value_t * end = start + n;

    entry_value_t * itr =
      std::lower_bound(start, end, entry_value_t(entry),
        [](const entry_value_t & e1, const entry_value_t & e2) -> bool{
          return e1.entry < e2.entry;
        });

    // if we are attempting to create an entry that already exists
    // just over-write the value and exit.  
    if ( itr != end && itr->entry == entry) {
      return itr->value;
    }

    while(end != itr) {
      *(end) = *(end - 1);
      --end;
    } // while

    itr->entry = entry;

    if(index < h_.num_exclusive_){
      (*h_.num_exclusive_insertions)++;
    }

    offset.set_count(n + 1);

    return itr->value;
  } // operator ()

  void dump(){
    for(size_t p = 0; p < 3; ++p){
      switch(p){
        case 0:
          std::cout << "exclusive: " << std::endl;
          break;
        case 1:
          std::cout << "shared: " << std::endl;
          break;
        case 2:
          std::cout << "ghost: " << std::endl;
          break;
        default:
          break;
      }

      size_t start = h_.pi_.start[p];
      size_t end = h_.pi_.end[p];

      for(size_t i = start; i < end; ++i){
        const offset_t& offset = h_.offsets_[i];
        std::cout << "  index: " << i << std::endl;
        for(size_t j = 0; j < offset.count(); ++j){
          std::cout << "    " << h_.entries_[i * h_.num_slots_ + j].entry << 
            " = " << h_.entries_[i * h_.num_slots_ + j].value << std::endl;
        }

        auto p = h_.spare_map_->equal_range(i);
        auto itr = p.first;
        while(itr != p.second){
          std::cout << "    +" << itr->second.entry << " = " << 
            itr->second.value << std::endl;
          ++itr;
        }
      }      
    }
  }

  void
  erase(
    size_t index,
    size_t entry
  )
  {
    if(!h_.erase_set_){
      h_.erase_set_ = new erase_set_t;
    }

    h_.erase_set_->emplace(std::make_pair(index, entry));
  }

  handle_t h_;  
};

} // namespace flecsi

#endif // flecsi_mutator_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
