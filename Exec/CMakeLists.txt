set(erf_lib_name erf_srclib)
add_library(${erf_lib_name} OBJECT)
include(${CMAKE_SOURCE_DIR}/CMake/BuildERFExe.cmake)
build_erf_lib(${erf_lib_name})

#set(AMRW_SRC_DIR /projects/erf/bperry/amr-wind)
#add_subdirectory(${AMRW_SRC_DIR} amrwind)

if (ERF_ENABLE_MULTIBLOCK)
  add_subdirectory(MultiBlock)
  add_subdirectory(MultiBlockABL)
else ()
  add_subdirectory(ScalarAdvDiff)
endif()
