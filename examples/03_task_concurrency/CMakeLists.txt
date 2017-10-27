#~----------------------------------------------------------------------------~#
#  @@@@@@@@  @@           @@@@@@   @@@@@@@@ @@
# /@@/////  /@@          @@////@@ @@////// /@@
# /@@       /@@  @@@@@  @@    // /@@       /@@
# /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@
# /@@////   /@@/@@@@@@@/@@       ////////@@/@@
# /@@       /@@/@@//// //@@    @@       /@@/@@
# /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@
# //       ///  //////   //////  ////////  // 
# 
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

#------------------------------------------------------------------------------#
# Compile definition to include driver
#------------------------------------------------------------------------------#

if(FLECSI_RUNTIME_MODEL STREQUAL "legion")

  add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_TLT_INIT)
  add_definitions(-DFLECSI_ENABLE_SPECIALIZATION_SPMD_INIT)
  add_definitions(-DCINCH_OVERRIDE_DEFAULT_INITIALIZATION_DRIVER)

  set(exec_path ${PROJECT_SOURCE_DIR}/flecsi/execution)

  add_custom_command(
    OUTPUT task_concurrency_input
    COMMAND ${CMAKE_COMMAND} -E copy ${exec_path}/test/simple2d-16x16.msh .
  )

  #----------------------------------------------------------------------------#
  # Add source files for example application.
  #----------------------------------------------------------------------------#

  set (exmaple_src      task_concurrency.cc
                        ${exec_path}/../supplemental/coloring/add_colorings.cc
                        ${exec_path}/driver_initialization.cc
                        ${exec_path}/legion/runtime_driver.cc)

  #----------------------------------------------------------------------------#
  # Add a rule to build the executable
  #----------------------------------------------------------------------------#


  add_executable(task_concurrency ${exmaple_src} task_concurrency_input)

  #----------------------------------------------------------------------------#
  # Add link dependencies
  #----------------------------------------------------------------------------#

  target_link_libraries(task_concurrency
    flecsi
    ${FLECSI_RUNTIME_LIBRARIES}
  )

endif()
#~---------------------------------------------------------------------------~-#
# Formatting options
# vim: set tabstop=2 shiftwidth=2 expandtab :
#~---------------------------------------------------------------------------~-#
