function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

function(build_erf_lib erf_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
  set(BIN_DIR ${CMAKE_BINARY_DIR}/Source/${erf_lib_name})

  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_lib_name})

  set(ERF_EOS_DIR "${CMAKE_SOURCE_DIR}/Source")
  target_sources(${erf_lib_name} PRIVATE
                 ${ERF_EOS_DIR}/EOS.H)
  target_include_directories(${erf_lib_name} SYSTEM PUBLIC ${ERF_EOS_DIR})

  if(ERF_ENABLE_MOISTURE)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MOISTURE)
  endif()

  if(ERF_ENABLE_MULTIBLOCK)
    set(ADV_SRC_DIR ${CMAKE_SOURCE_DIR}/Submodules/amrex-tutorials/ExampleCodes/Amr/Advection_AmrCore/Source/)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/MultiBlockContainer.H
                   ${SRC_DIR}/MultiBlockContainer.cpp
                   ${ADV_SRC_DIR}/AdvancePhiAllLevels.cpp
                   ${ADV_SRC_DIR}/AdvancePhiAtLevel.cpp
                   ${ADV_SRC_DIR}/AmrCoreAdv.cpp
                   ${ADV_SRC_DIR}/AmrCoreAdv.H
                   ${ADV_SRC_DIR}/bc_fill.H
                   ${ADV_SRC_DIR}/DefineVelocity.cpp
                   ${ADV_SRC_DIR}/face_velocity.H
                   ${ADV_SRC_DIR}/Kernels.H
                   ${ADV_SRC_DIR}/Tagging.H
                   ${ADV_SRC_DIR}/Prob.H
                   ${ADV_SRC_DIR}/Src_K/Adv_K.H
                   ${ADV_SRC_DIR}/Src_K/compute_flux_2D_K.H
                   ${ADV_SRC_DIR}/Src_K/compute_flux_3D_K.H
                   ${ADV_SRC_DIR}/Src_K/slope_K.H)
    target_include_directories(${erf_lib_name} SYSTEM PUBLIC ${ADV_SRC_DIR} ${ADV_SRC_DIR}/Src_K )
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_MULTIBLOCK)
    endif()

    #set(AMRW_SRC_DIR /projects/erf/bperry/amr-wind/amr-wind)
    target_sources(${erf_lib_name} PRIVATE
${AMRW_SRC_DIR}/incflo_advance.cpp
${AMRW_SRC_DIR}/incflo.cpp
${AMRW_SRC_DIR}/incflo_compute_dt.cpp
${AMRW_SRC_DIR}/incflo_regrid.cpp
${AMRW_SRC_DIR}/CFDSim.cpp
${AMRW_SRC_DIR}/core/SimTime.cpp
${AMRW_SRC_DIR}/core/Field.cpp
${AMRW_SRC_DIR}/core/IntField.cpp
${AMRW_SRC_DIR}/core/FieldRepo.cpp
${AMRW_SRC_DIR}/core/ScratchField.cpp
${AMRW_SRC_DIR}/core/ViewField.cpp
${AMRW_SRC_DIR}/core/MLMGOptions.cpp
${AMRW_SRC_DIR}/core/MeshMap.cpp
${AMRW_SRC_DIR}/boundary_conditions/BCInterface.cpp
${AMRW_SRC_DIR}/boundary_conditions/FixedGradientBC.cpp
${AMRW_SRC_DIR}/convection/incflo_godunov_weno.cpp
${AMRW_SRC_DIR}/convection/incflo_godunov_ppm.cpp
${AMRW_SRC_DIR}/convection/incflo_godunov_plm.cpp
${AMRW_SRC_DIR}/convection/incflo_godunov_predict.cpp
${AMRW_SRC_DIR}/convection/incflo_godunov_advection.cpp
${AMRW_SRC_DIR}/convection/incflo_mol_fluxes.cpp
${AMRW_SRC_DIR}/derive/field_algebra.cpp
${AMRW_SRC_DIR}/diffusion/incflo_diffusion.cpp
${AMRW_SRC_DIR}/projection/incflo_apply_nodal_projection.cpp
${AMRW_SRC_DIR}/setup/init.cpp
${AMRW_SRC_DIR}/utilities/diagnostics.cpp
${AMRW_SRC_DIR}/utilities/io.cpp
${AMRW_SRC_DIR}/utilities/bc_ops.cpp
${AMRW_SRC_DIR}/utilities/console_io.cpp
${AMRW_SRC_DIR}/utilities/IOManager.cpp
${AMRW_SRC_DIR}/utilities/FieldPlaneAveraging.cpp
${AMRW_SRC_DIR}/utilities/SecondMomentAveraging.cpp
${AMRW_SRC_DIR}/utilities/ThirdMomentAveraging.cpp
${AMRW_SRC_DIR}/utilities/PostProcessing.cpp
${AMRW_SRC_DIR}/utilities/DerivedQuantity.cpp
${AMRW_SRC_DIR}/utilities/DerivedQtyDefs.cpp
${AMRW_SRC_DIR}/utilities/tagging/RefinementCriteria.cpp
${AMRW_SRC_DIR}/utilities/tagging/CartBoxRefinement.cpp
${AMRW_SRC_DIR}/utilities/tagging/FieldRefinement.cpp
${AMRW_SRC_DIR}/utilities/tagging/GradientMagRefinement.cpp
${AMRW_SRC_DIR}/utilities/tagging/CurvatureRefinement.cpp
${AMRW_SRC_DIR}/utilities/tagging/QCriterionRefinement.cpp
${AMRW_SRC_DIR}/utilities/tagging/VorticityMagRefinement.cpp
${AMRW_SRC_DIR}/utilities/tagging/OversetRefinement.cpp
${AMRW_SRC_DIR}/utilities/tagging/GeometryRefinement.cpp
${AMRW_SRC_DIR}/utilities/tagging/BoxRefiner.cpp
${AMRW_SRC_DIR}/utilities/tagging/CylinderRefiner.cpp
${AMRW_SRC_DIR}/utilities/sampling/Sampling.cpp
${AMRW_SRC_DIR}/utilities/sampling/SamplingContainer.cpp
${AMRW_SRC_DIR}/utilities/sampling/LineSampler.cpp
${AMRW_SRC_DIR}/utilities/sampling/LidarSampler.cpp
${AMRW_SRC_DIR}/utilities/sampling/DTUSpinnerSampler.cpp
${AMRW_SRC_DIR}/utilities/sampling/PlaneSampler.cpp
${AMRW_SRC_DIR}/utilities/sampling/ProbeSampler.cpp
${AMRW_SRC_DIR}/utilities/sampling/FieldNorms.cpp
${AMRW_SRC_DIR}/utilities/sampling/KineticEnergy.cpp
${AMRW_SRC_DIR}/utilities/sampling/Enstrophy.cpp
${AMRW_SRC_DIR}/utilities/sampling/FreeSurface.cpp
${AMRW_SRC_DIR}/utilities/sampling/WaveEnergy.cpp
${AMRW_SRC_DIR}/utilities/averaging/TimeAveraging.cpp
${AMRW_SRC_DIR}/utilities/averaging/ReAveraging.cpp
${AMRW_SRC_DIR}/utilities/averaging/ReynoldsStress.cpp
#${AMRW_SRC_DIR}/utilities/ncutils/nc_interface.cpp
#${AMRW_SRC_DIR}/utilities/ascent/ascent.H
#${AMRW_SRC_DIR}/utilities/ascent/ascent.cpp
${AMRW_SRC_DIR}/wind_energy/ABL.cpp
${AMRW_SRC_DIR}/wind_energy/ABLStats.cpp
${AMRW_SRC_DIR}/wind_energy/ABLFieldInit.cpp
${AMRW_SRC_DIR}/wind_energy/ABLWallFunction.cpp
${AMRW_SRC_DIR}/wind_energy/ABLFillInflow.cpp
${AMRW_SRC_DIR}/wind_energy/ABLBoundaryPlane.cpp
${AMRW_SRC_DIR}/wind_energy/MOData.cpp
${AMRW_SRC_DIR}/wind_energy/ABLFillMPL.cpp
${AMRW_SRC_DIR}/wind_energy/ABLModulatedPowerLaw.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/actuator_utils.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/Actuator.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/ActuatorContainer.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/aero/AirfoilTable.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/wing/wing_ops.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/wing/FlatPlate.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/wing/FixedWing.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/turbine/turbine_utils.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/turbine/fast/FastIface.cpp
#${AMRW_SRC_DIR}/wind_energy/actuator/turbine/fast/TurbineFast.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/disk/ActuatorDisk.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/disk/disk_ops.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/disk/Joukowsky_ops.cpp
${AMRW_SRC_DIR}/wind_energy/actuator/disk/uniform_ct_ops.cpp
${AMRW_SRC_DIR}/equation_systems/PDEBase.cpp
${AMRW_SRC_DIR}/equation_systems/DiffusionOps.cpp
${AMRW_SRC_DIR}/equation_systems/icns/icns_advection.cpp
${AMRW_SRC_DIR}/equation_systems/icns/icns.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/ABLForcing.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/GeostrophicForcing.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/BoussinesqBuoyancy.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/DensityBuoyancy.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/GravityForcing.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/CoriolisForcing.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/BodyForce.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/ABLMeanBoussinesq.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/ActuatorForcing.cpp
${AMRW_SRC_DIR}/equation_systems/icns/source_terms/SynthTurbForcing.cpp
${AMRW_SRC_DIR}/equation_systems/temperature/temperature.cpp
${AMRW_SRC_DIR}/equation_systems/density/density.cpp
${AMRW_SRC_DIR}/equation_systems/tke/TKE.cpp
${AMRW_SRC_DIR}/equation_systems/tke/source_terms/KsgsM84Src.cpp
${AMRW_SRC_DIR}/equation_systems/tke/source_terms/KwSSTSrc.cpp
${AMRW_SRC_DIR}/equation_systems/sdr/SDR.cpp
${AMRW_SRC_DIR}/equation_systems/sdr/source_terms/SDRSrc.cpp
${AMRW_SRC_DIR}/equation_systems/levelset/levelset.cpp
${AMRW_SRC_DIR}/equation_systems/vof/SplitAdvection.cpp
${AMRW_SRC_DIR}/equation_systems/vof/vof.cpp
${AMRW_SRC_DIR}/turbulence/LaminarModel.cpp
${AMRW_SRC_DIR}/turbulence/turb_utils.cpp
${AMRW_SRC_DIR}/turbulence/LES/Smagorinsky.cpp
${AMRW_SRC_DIR}/turbulence/LES/OneEqKsgs.cpp
${AMRW_SRC_DIR}/turbulence/RANS/KOmegaSST.cpp
${AMRW_SRC_DIR}/turbulence/RANS/KOmegaSSTIDDES.cpp
${AMRW_SRC_DIR}/physics/BoussinesqBubble.cpp
${AMRW_SRC_DIR}/physics/BoussinesqBubbleFieldInit.cpp
${AMRW_SRC_DIR}/physics/ChannelFlow.cpp
${AMRW_SRC_DIR}/physics/RayleighTaylor.cpp
${AMRW_SRC_DIR}/physics/RayleighTaylorFieldInit.cpp
${AMRW_SRC_DIR}/physics/TaylorGreenVortex.cpp
${AMRW_SRC_DIR}/physics/FreeStream.cpp
${AMRW_SRC_DIR}/physics/ConvectingTaylorVortex.cpp
${AMRW_SRC_DIR}/physics/EkmanSpiral.cpp
${AMRW_SRC_DIR}/physics/SyntheticTurbulence.cpp
${AMRW_SRC_DIR}/physics/HybridRANSLESABL.cpp
${AMRW_SRC_DIR}/physics/VortexRing.cpp
${AMRW_SRC_DIR}/physics/multiphase/MultiPhase.cpp
${AMRW_SRC_DIR}/physics/multiphase/VortexPatch.cpp
${AMRW_SRC_DIR}/physics/udfs/UDF.cpp
${AMRW_SRC_DIR}/physics/udfs/LinearProfile.cpp
${AMRW_SRC_DIR}/physics/udfs/PowerLawProfile.cpp
#${AMRW_SRC_DIR}/physics/mms/MMS.cpp
#${AMRW_SRC_DIR}/physics/mms/MMSForcing.cpp
${AMRW_SRC_DIR}/overset/TiogaInterface.cpp
${AMRW_SRC_DIR}/immersed_boundary/IB.cpp
${AMRW_SRC_DIR}/immersed_boundary/bluff_body/bluff_body_ops.cpp
${AMRW_SRC_DIR}/immersed_boundary/bluff_body/Box.cpp
${AMRW_SRC_DIR}/immersed_boundary/bluff_body/Cylinder.cpp
${AMRW_SRC_DIR}/immersed_boundary/bluff_body/Sphere.cpp
${AMRW_SRC_DIR}/mesh_mapping_models/ChannelFlowMap.cpp
${AMRW_SRC_DIR}/mesh_mapping_models/ConstantMap.cpp
    )

  target_include_directories(${erf_lib_name} SYSTEM PUBLIC
${AMRW_SRC_DIR}/../
  )

  if(ERF_ENABLE_NETCDF)
    target_sources(${erf_lib_name} PRIVATE
                   ${SRC_DIR}/IO/NCInterface.H
                   ${SRC_DIR}/IO/NCWpsFile.H
                   ${SRC_DIR}/IO/NCPlotFile.H
                   ${SRC_DIR}/IO/NCBuildFABs.cpp
                   ${SRC_DIR}/IO/NCInterface.cpp
                   ${SRC_DIR}/IO/NCPlotFile.cpp
                   ${SRC_DIR}/IO/NCCheckpoint.cpp
                   ${SRC_DIR}/IO/NCMultiFabFile.cpp
                   ${SRC_DIR}/IO/ReadFromWRFBdy.cpp
                   ${SRC_DIR}/IO/ReadFromWRFInput.cpp
                   ${SRC_DIR}/IO/NCColumnFile.cpp)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_NETCDF)
  endif()

  if(ERF_ENABLE_HDF5)
    target_compile_definitions(${erf_lib_name} PUBLIC ERF_USE_HDF5)
  endif()

  target_sources(${erf_lib_name}
     PRIVATE
       ${SRC_DIR}/DataStruct.H
       ${SRC_DIR}/ERF_Constants.H
       ${SRC_DIR}/Derive.H
       ${SRC_DIR}/Derive.cpp
       ${SRC_DIR}/IndexDefines.H
       ${SRC_DIR}/prob_common.H
       ${SRC_DIR}/ERF_Math.H
       ${SRC_DIR}/ERF.H
       ${SRC_DIR}/ERF.cpp
       ${SRC_DIR}/ERF_init.cpp
       ${SRC_DIR}/ERF_init1d.cpp
       ${SRC_DIR}/ERF_SumIQ.cpp
       ${SRC_DIR}/ERF_Tagging.cpp
       ${SRC_DIR}/BoundaryConditions/ABLMost.H
       ${SRC_DIR}/BoundaryConditions/ABLMost.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_cons.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_xvel.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_yvel.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_zvel.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_bndryreg.cpp
       ${SRC_DIR}/BoundaryConditions/BoundaryConditions_wrfbdy.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_FillPatch.cpp
       ${SRC_DIR}/BoundaryConditions/ERF_PhysBCFunct.cpp
       ${SRC_DIR}/BoundaryConditions/PlaneAverage.H
       ${SRC_DIR}/BoundaryConditions/VelPlaneAverage.H
       ${SRC_DIR}/BoundaryConditions/DirectionSelector.H
       ${SRC_DIR}/IO/Checkpoint.cpp
       ${SRC_DIR}/IO/ERF_ReadBndryPlanes.H
       ${SRC_DIR}/IO/ERF_ReadBndryPlanes.cpp
       ${SRC_DIR}/IO/ERF_WriteBndryPlanes.H
       ${SRC_DIR}/IO/ERF_WriteBndryPlanes.cpp
       ${SRC_DIR}/IO/Plotfile.cpp
       ${SRC_DIR}/IO/writeJobInfo.cpp
       ${SRC_DIR}/SpatialStencils/Advection.cpp
       ${SRC_DIR}/SpatialStencils/Diffusion.cpp
       ${SRC_DIR}/SpatialStencils/EddyViscosity.H
       ${SRC_DIR}/SpatialStencils/ExpansionRate.H
       ${SRC_DIR}/SpatialStencils/StrainRate.H
       ${SRC_DIR}/SpatialStencils/StressTerm.H
       ${SRC_DIR}/SpatialStencils/Interpolation.H
       ${SRC_DIR}/SpatialStencils/ComputeTurbulentViscosity.cpp
       ${SRC_DIR}/SpatialStencils/MomentumToVelocity.cpp
       ${SRC_DIR}/SpatialStencils/VelocityToMomentum.cpp
       ${SRC_DIR}/SpatialStencils/TerrainMetrics.H
       ${SRC_DIR}/SpatialStencils/TerrainMetrics.cpp
       ${SRC_DIR}/TimeIntegration/ERF_ComputeTimestep.cpp
       ${SRC_DIR}/TimeIntegration/ERF_TimeStepping.cpp
       ${SRC_DIR}/TimeIntegration/ERF_MRI.H
       ${SRC_DIR}/TimeIntegration/ERF_SRI.H
       ${SRC_DIR}/TimeIntegration/TimeIntegration.H
       ${SRC_DIR}/TimeIntegration/TimeIntegration_driver.cpp
       ${SRC_DIR}/TimeIntegration/ERF_slow_rhs.cpp
       ${SRC_DIR}/TimeIntegration/ERF_fast_rhs.cpp
  )

  if(NOT "${erf_exe_name}" STREQUAL "erf_unit_tests")
    target_sources(${erf_lib_name}
       PRIVATE
         ${SRC_DIR}/main.cpp
    )
  endif()

  include(AMReXBuildInfo)
  generate_buildinfo(${erf_lib_name} ${CMAKE_SOURCE_DIR})
  target_include_directories(${erf_lib_name} PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

  if(ERF_ENABLE_NETCDF)
    if(NETCDF_FOUND)
      #Link our executable to the NETCDF libraries, etc
      target_link_libraries(${erf_lib_name} PUBLIC ${NETCDF_LINK_LIBRARIES})
      target_include_directories(${erf_lib_name} PUBLIC ${NETCDF_INCLUDE_DIRS})
    endif()
  endif()

  if(ERF_ENABLE_MPI)
    target_link_libraries(${erf_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
  endif()

  #ERF include directories
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR})
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR}/BoundaryConditions)
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR}/SpatialStencils)
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR}/TimeIntegration)
  target_include_directories(${erf_lib_name} PUBLIC ${SRC_DIR}/IO)
  target_include_directories(${erf_lib_name} PUBLIC ${CMAKE_BINARY_DIR})

  if (AMReX_SUNDIALS)
    target_link_libraries(${erf_lib_name} PUBLIC SUNDIALS::cvode)
    target_link_libraries(${erf_lib_name} PUBLIC SUNDIALS::arkode)
    target_link_libraries(${erf_lib_name} PUBLIC SUNDIALS::nvecmanyvector)
  endif ()

  #Link to amrex library
  target_link_libraries_system(${erf_lib_name} PUBLIC amrex)
  target_link_libraries_system(${erf_lib_name} PUBLIC AMReX-Hydro::amrex_hydro_api)
  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ERF_SOURCES}")
    endforeach()
    set_target_properties(
    ${erf_lib_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  #Define what we want to be installed during a make install
  install(TARGETS ${erf_lib_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction(build_erf_lib)

function(build_erf_exe erf_exe_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  target_link_libraries(${erf_exe_name} PRIVATE ${erf_lib_name} AMReX-Hydro::amrex_hydro_api)
  target_link_libraries(${erf_exe_name}  PUBLIC ${erf_lib_name})
  include(${CMAKE_SOURCE_DIR}/CMake/SetERFCompileFlags.cmake)
  set_erf_compile_flags(${erf_exe_name})

  target_sources(${erf_exe_name}
     PRIVATE
       ${SRC_DIR}/ERF_init_bcs.cpp

  )
  if(ERF_ENABLE_CUDA)
    set(pctargets "${erf_exe_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(ERF_SOURCES ${tgt} SOURCES)
      list(FILTER ERF_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${ERF_SOURCES} PROPERTIES LANGUAGE CUDA)
      message(STATUS "setting cuda for ${ERF_SOURCES}")
    endforeach()
    set_target_properties(
    ${erf_exe_name} PROPERTIES
    LANGUAGE CUDA
    CUDA_SEPARABLE_COMPILATION ON
    CUDA_RESOLVE_DEVICE_SYMBOLS ON)
  endif()

  install(TARGETS ${erf_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
