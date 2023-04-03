
namespace amrex {

const char* buildInfoGetBuildDate() {

  static const char BUILD_DATE[] = "2022-08-10 17:05:38.652554";
  return BUILD_DATE;
}

const char* buildInfoGetBuildDir() {

  static const char BUILD_DIR[] = "/lustre/eaglefs/projects/erf/bperry/ERFmb/Exec/MultiBlock";
  return BUILD_DIR;
}

const char* buildInfoGetBuildMachine() {

  static const char BUILD_MACHINE[] = "Linux r9i1n7 3.10.0-1062.9.1.el7.x86_64 #1 SMP Fri Dec 6 15:49:49 UTC 2019 x86_64 x86_64 x86_64 GNU/Linux";
  return BUILD_MACHINE;
}

const char* buildInfoGetAMReXDir() {

  static const char AMREX_DIR[] = "../../Submodules/AMReX";
  return AMREX_DIR;
}

const char* buildInfoGetComp() {

  static const char COMP[] = "gnu";
  return COMP;
}

const char* buildInfoGetCompVersion() {

  static const char COMP_VERSION[] = "8.4.0";
  return COMP_VERSION;
}

// deprecated
const char* buildInfoGetFcomp() {

  static const char FCOMP[] = "";
  return FCOMP;
}

// deprecated
const char* buildInfoGetFcompVersion() {

  static const char FCOMP_VERSION[] = "";
  return FCOMP_VERSION;
}

const char* buildInfoGetCXXName() {

  static const char CXX_comp_name[] = "mpicxx";
  return CXX_comp_name;
}

const char* buildInfoGetFName() {

  static const char F_comp_name[] = "";
  return F_comp_name;
}

const char* buildInfoGetCXXFlags() {

  static const char CXX_flags[] = " -Werror=return-type -g -O0 -ggdb -ftrapv -Wall -Wextra -Wpedantic -Wnull-dereference -Wfloat-conversion -Wshadow -Woverloaded-virtual  -pthread   -Werror=return-type -g -O0 -ggdb -ftrapv -Wall -Wextra -Wpedantic -Wnull-dereference -Wfloat-conversion -Wshadow -Woverloaded-virtual  -pthread   -DAMREX_DEBUG -DAMREX_TESTING -DBL_USE_MPI -DAMREX_USE_MPI -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -DBL_FORT_USE_UNDERSCORE -DAMREX_FORT_USE_UNDERSCORE -DBL_Linux -DAMREX_Linux -DBL_USE_ASSERTION -DAMREX_USE_ASSERTION -DOMPI_SKIP_MPICXX -DERF_USE_MULTIBLOCK -Itmp_build_dir/s/3d.gnu.DEBUG.MPI.EXE -I. -I../../Source -I../../Source/BoundaryConditions -I../../Source/SpatialStencils -I../../Source/TimeIntegration -I../../Source/IO -I../../Exec/MultiBlock -I. -I../../Submodules/amrex-tutorials/ExampleCodes/Amr/Advection_AmrCore//Source -I../../Submodules/amrex-tutorials/ExampleCodes/Amr/Advection_AmrCore//Source/Src_K -I../../Submodules/AMReX/Src/Base -I../../Submodules/AMReX/Src/Base/Parser -I../../Submodules/AMReX/Src/Base -I../../Submodules/AMReX/Src/Base/Parser -I../../Submodules/AMReX/Src/Boundary -I../../Submodules/AMReX/Src/AmrCore -I../../Submodules/AMReX/Tools/C_scripts ";
  return CXX_flags;
}

const char* buildInfoGetFFlags() {

  static const char F_flags[] = "";
  return F_flags;
}

const char* buildInfoGetLinkFlags() {

  static const char link_flags[] = " -Werror=return-type -g -O0 -ggdb -ftrapv -Wall -Wextra -Wpedantic -Wnull-dereference -Wfloat-conversion -Wshadow -Woverloaded-virtual  -pthread   -DAMREX_DEBUG -DAMREX_TESTING -DBL_USE_MPI -DAMREX_USE_MPI -DBL_SPACEDIM=3 -DAMREX_SPACEDIM=3 -DBL_FORT_USE_UNDERSCORE -DAMREX_FORT_USE_UNDERSCORE -DBL_Linux -DAMREX_Linux -DBL_USE_ASSERTION -DAMREX_USE_ASSERTION -DOMPI_SKIP_MPICXX -DERF_USE_MULTIBLOCK -Itmp_build_dir/s/3d.gnu.DEBUG.MPI.EXE -I. -I../../Source -I../../Source/BoundaryConditions -I../../Source/SpatialStencils -I../../Source/TimeIntegration -I../../Source/IO -I../../Exec/MultiBlock -I. -I../../Submodules/amrex-tutorials/ExampleCodes/Amr/Advection_AmrCore//Source -I../../Submodules/amrex-tutorials/ExampleCodes/Amr/Advection_AmrCore//Source/Src_K -I../../Submodules/AMReX/Src/Base -I../../Submodules/AMReX/Src/Base/Parser -I../../Submodules/AMReX/Src/Base -I../../Submodules/AMReX/Src/Base/Parser -I../../Submodules/AMReX/Src/Boundary -I../../Submodules/AMReX/Src/AmrCore -I../../Submodules/AMReX/Tools/C_scripts  -L. -L/nopt/nrel/apps/base/2020-05-12/spack/opt/spack/linux-centos7-x86_64/gcc-4.8.5/gcc-8.4.0-xrmnckz436rnkcwooooi2r7n35pojkov/lib64/../lib64/";
  return link_flags;
}

const char* buildInfoGetLibraries() {

  static const char libraries[] = " -fexceptions -pthread -I/nopt/nrel/apps/openmpi/4.0.4-centos77/lib -L/nopt/nrel/apps/base/2020-05-12/spack/opt/spack/linux-centos7-x86_64/gcc-8.4.0/hwloc-1.11.11-mb5lwdajmllvrdtwltwe3r732aca76ny/lib -L/nopt/slurm/current/lib -Wl,-rpath -Wl,/nopt/nrel/apps/base/2020-05-12/spack/opt/spack/linux-centos7-x86_64/gcc-8.4.0/hwloc-1.11.11-mb5lwdajmllvrdtwltwe3r732aca76ny/lib -Wl,-rpath -Wl,/nopt/slurm/current/lib -Wl,-rpath -Wl,/nopt/nrel/apps/openmpi/4.0.4-centos77/lib -Wl,--enable-new-dtags -L/nopt/nrel/apps/openmpi/4.0.4-centos77/lib -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lmpi_cxx -lmpi -lgfortran -lquadmath";
  return libraries;
}

const char* buildInfoGetAux(int i) {

  //static const char AUX1[] = "${AUX[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

int buildInfoGetNumModules() {
  // int num_modules = X;
  int num_modules = 0;

  return num_modules;
}

const char* buildInfoGetModuleName(int i) {

  //static const char MNAME1[] = "${MNAME[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetModuleVal(int i) {

  //static const char MVAL1[] = "${MVAL[1]}";

  static const char EMPT[] = "";

  switch(i)
  {

    default: return EMPT;
  }
}

const char* buildInfoGetGitHash(int i) {

  //static const char HASH1[] = "${GIT[1]}";
  static const char HASH1[] = "last-amrlevel-version-449-g7d786e56-dirty";
  static const char HASH2[] = "22.07-19-g48702b488";

  static const char EMPT[] = "";

  switch(i)
  {
    case 1: return HASH1;
    case 2: return HASH2;

    default: return EMPT;
  }
}

const char* buildInfoGetBuildGitHash() {

  //static const char HASH[] = "${GIT}";
  static const char HASH[] = "";


  return HASH;
}

const char* buildInfoGetBuildGitName() {

  //static const char NAME[] = "";
  static const char NAME[] = "";


  return NAME;
}

#ifdef AMREX_USE_CUDA
const char* buildInfoGetCUDAVersion() {

  static const char CUDA_VERSION[] = "";
  return CUDA_VERSION;
}
#endif

}
