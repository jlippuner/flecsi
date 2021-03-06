/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------------~~*/

//----------------------------------------------------------------------------//
//! @file
//! @date Initial file creation: Jul 26, 2016
//----------------------------------------------------------------------------//

#include "flecsi/execution/mpi/context_policy.h"

namespace flecsi {
namespace execution {

//----------------------------------------------------------------------------//
// Implementation of mpi_context_policy_t::initialize.
//----------------------------------------------------------------------------//

int
mpi_context_policy_t::initialize(
  int argc,
  char ** argv
)
{
  MPI_Comm_rank(MPI_COMM_WORLD, &color_);

  runtime_driver(argc, argv);

  return 0;
} // mpi_context_policy_t::initialize

} // namespace execution 
} // namespace flecsi

/*~------------------------------------------------------------------------~--*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
