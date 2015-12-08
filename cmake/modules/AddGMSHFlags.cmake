# module providing convenience methods for compiling
# binaries with GMSH support
#
# Provides the following functions:
#
# add_dune_gmsh_flags(target1 target2...)
#
# adds GMSH flags to the targets for compilation and linking
function(add_dune_gmsh_flags _targets)
  if(GMSH_FOUND)
    foreach(_target ${_targets})
      target_link_libraries(${_target} ${GMSH_DUNE_LIBRARIES})
      get_target_property(_props ${_target} COMPILE_FLAGS)
      string(REPLACE "_props-NOTFOUND" "" _props "${_props}")
      set_target_properties(${_target} PROPERTIES COMPILE_FLAGS
        "${_props} ${GMSH_DUNE_COMPILE_FLAGS}")
    endforeach()
  endif()
endfunction()
