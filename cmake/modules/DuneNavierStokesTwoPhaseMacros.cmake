find_package(Alberta REQUIRED)
include(AddAlbertaFlags)

find_package(GMSH REQUIRED)
include(AddGMSHFlags)

find_package(SuiteSparse OPTIONAL_COMPONENTS LDL SPQR UMFPACK REQUIRED)
include(AddSuiteSparseFlags)
