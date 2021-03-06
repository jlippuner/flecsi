/*~--------------------------------------------------------------------------~*
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_topology_entity_storage_h
#define flecsi_topology_entity_storage_h

#include "flecsi/runtime/flecsi_runtime_entity_storage_policy.h"

//----------------------------------------------------------------------------//
//! @file
//! @date Initial file creation: Aug 01, 2016
//----------------------------------------------------------------------------//

namespace flecsi {
namespace topology {

template<typename T>
using entity_storage_t = FLECSI_RUNTIME_ENTITY_STORAGE_TYPE<T>;

using offset_storage_t = FLECSI_RUNTIME_OFFSET_STORAGE_TYPE;

} // namespace topology
} // namespace flecsi

#endif // flecsi_topology_entity_storage_h

/*~-------------------------------------------------------------------------~-*
*~-------------------------------------------------------------------------~-*/
