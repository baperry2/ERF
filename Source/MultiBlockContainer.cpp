#include <MultiBlockContainer.H>
#include <AMReX_NonLocalBC.H>
#include <ERF.H>

// Vector input constructor
MultiBlockContainer::MultiBlockContainer(const std::vector<amrex::RealBox>& rb_v,
                                         std::vector<int> max_level_in_v,
                                         const std::vector<amrex::Vector<int>>& n_cell_in_v,
                                         std::vector<int> coord_v,
                                         const std::vector<amrex::Vector<amrex::IntVect>>& ref_ratios_v,
                                         const std::vector<amrex::Array<int,AMREX_SPACEDIM>>& is_per_v,
                                         std::vector<std::string> prefix_v,
                                         int max_step)
: m_max_step(max_step),
  erf1(rb_v[0],max_level_in_v[0],n_cell_in_v[0],coord_v[0],ref_ratios_v[0],is_per_v[0],prefix_v[0])
{
    // Store ptr to container to call member functions
    amrwind.SetMultiBlockPointer(this);
    erf1.SetMultiBlockPointer(this);

    // Set the permutation/sign of dtos
    dtos_efroma.permutation = amrex::IntVect{AMREX_D_DECL(   0,   1,   2)};
    dtos_efroma.sign        = amrex::IntVect{AMREX_D_DECL(   1,   1,   1)};
    dtos_afrome.permutation = amrex::IntVect{AMREX_D_DECL(   0,   1,   2)};
    dtos_afrome.sign        = amrex::IntVect{AMREX_D_DECL(   1,   1,   1)};

    // Set offset of dtos (NOTE: i_dst =  i_src - i_off -> [0] - [1])
    {
      const amrex::Geometry geom = amrwind.repo().mesh().Geom(0);
      amrex::Real dx = geom.CellSize(0);
      amrex::Real dy = geom.CellSize(1);
      amrex::Real dz = geom.CellSize(2);
      int offx = amrex::Math::floor((geom.ProbLo(0) - rb_v[0].lo(0)) / dx);
      int offy = amrex::Math::floor((geom.ProbLo(1) - rb_v[0].lo(1)) / dy);
      int offz = amrex::Math::floor((geom.ProbLo(2) - rb_v[0].lo(2)) / dz);
      dtos_efroma.offset = amrex::IntVect{AMREX_D_DECL(offx, offy, offz)};
      dtos_afrome.offset = amrex::IntVect{AMREX_D_DECL(-offx, -offy, -offz)};
    }
}

// Destructor
MultiBlockContainer::~MultiBlockContainer()
{
}

// Initialize block data
void
MultiBlockContainer::InitializeBlocks()
{
  amrex::Print() << "    STARTING INITIALIZATION : \n";
  amrex::Print() << "===================================="  << "\n";
  amrex::Print() << "         ERF1 INITIALIZATION        "  << "\n";
  amrex::Print() << "------------------------------------"  << "\n";
  erf1.InitData();
  amrex::Print() << '\n';
  amrex::Print() << "------------------------------------"  << "\n";
  amrex::Print() << "       AMRWIND INITIALIZATION       "  << "\n";
  amrex::Print() << "------------------------------------"  << "\n";
  amrwind.InitData();
  amrex::Print() << '\n';
  amrex::Print() << "------------------------------------"  << "\n";
  amrex::Print() << "     MultiBlock Intitialization     "  << "\n";
  amrex::Print() << "------------------------------------"  << "\n";
  SetBoxLists();
  SetBlockCommMetaData();
  amrex::ParmParse pp("mbc");
  do_two_way_coupling = false;
  pp.query("do_two_way_coupling", do_two_way_coupling);
  FillPatchBlocksAE();
  amrex::Print() << "------------------------------------"  << "\n";
  amrex::Print() << '\n';

}

// Set up BoxList vector for use with Communication Meta Data
void
MultiBlockContainer::SetBoxLists()
{
  // when copying from amr-wind to erf, data from the entire amr-wind domain is used
    amrex::Box awbox = amrwind.repo().mesh().Geom(0).Domain();
    abox_efroma = amrex::Box(awbox.smallEnd(), awbox.bigEnd());
    amrex::Dim3 se = dtos_efroma.Inverse(amrex::lbound(abox_efroma));
    amrex::Dim3 be = dtos_efroma.Inverse(amrex::ubound(abox_efroma));
    ebox_efroma = amrex::Box({se.x, se.y, se.z}, {be.x, be.y, be.z});

    // amrex::Print() << "A-W BOX in A-W Coords: " << abox_efroma << std::endl;
    // amrex::Print() << "A-W BOX in ERF Coords: " << ebox_efroma << std::endl;

    // when copying from erf to amr-wind, we only copy data at the boundaries of the amr-wind domain
    bool ok_to_continue = true;
    for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
      auto ori = oit();
      amrex::IntVect sev = awbox.smallEnd();
      amrex::IntVect bev = awbox.bigEnd();
      if (ori.faceDir() == 1) {
        sev[ori.coordDir()] = bev[ori.coordDir()];
        bev[ori.coordDir()] += 1;
      }
      else {
        bev[ori.coordDir()] = sev[ori.coordDir()];
        sev[ori.coordDir()] -= 1;
      }
      amrex::Box abx(sev,bev);
      aboxvec_afrome.push_back(abx);
      amrex::Dim3 ese = dtos_afrome(amrex::lbound(abx));
      amrex::Dim3 ebe = dtos_afrome(amrex::ubound(abx));
      amrex::Box ebx({ese.x, ese.y, ese.z}, {ebe.x, ebe.y, ebe.z});
      eboxvec_afrome.push_back(ebx);
      
      // amrex::Print() << "Ori: " << ori << " A-W BNDRY in A-W Coords: " << abx << std::endl;
      // amrex::Print() << "Ori: " << ori << " A-W BNDRY in ERF Coords: " << ebx << std::endl;

      // Do some error checking to ensure that the AMR-Wind is fully within the ERF domain (not touching any boundaries)
      // Note: in theory this error check could trip for a non-bopundary plane mass_inflow in AMR-Wind
      bool need_bndry = amrwind.repo().get_field("velocity").bc_type()[ori] == BC::mass_inflow;
      if ( !(erf1.domain_p[0].contains(ebx)) and need_bndry ) {
        amrex::Print() << "ERF domain must fully contain the AMR-Wind boundary planes, does not in direction "
                       << ori.coordDir() << " on face " << ori.faceDir() << std::endl;
        ok_to_continue = false;
      }
    }
    if (!ok_to_continue) amrex::Abort("Invalid domains for ERF/AMR-Wind coupling");
    
}

// Set up MB Communication Meta Data
void
MultiBlockContainer::SetBlockCommMetaData()
{
    // Hard-coded bounds for now
    int nvars  = erf1.vars_new[0].size(); // Destination MF
    int ndirs  = AMREX_SPACEDIM;

    amrex::IntVect nghost(0);
    amrex::NonLocalBC::MultiBlockCommMetaData *cmd_efroma_tmp =
      new amrex::NonLocalBC::MultiBlockCommMetaData(erf1.vars_new[0][Vars::cons], ebox_efroma,
                                                    amrwind.repo().get_field("temperature")(0), nghost, dtos_efroma);
    cmd_efroma.push_back(cmd_efroma_tmp);
}

// Advance blocks
void
MultiBlockContainer::AdvanceBlocks()
{
    amrex::Print() << "STARTING MAIN DRIVER FOR: " << m_max_step << " STEPS" << "\n";
    amrex::Print() << "\n";

    for (int step(1); step <= m_max_step; ++step) {
        amrex::Print() << "    STARTING ADVANCE DRIVER: " << step << "\n";
        amrex::Print() << "===================================="  << "\n";
        amrex::Print() << "           ERF BLOCK STARTS         "  << "\n";
        amrex::Print() << "------------------------------------"  << "\n";
        erf1.Evolve_MB(step,1);
        amrex::Print() << '\n';
        amrex::Print() << "------------------------------------"  << "\n";
        amrex::Print() << "        AMR-WIND BLOCK STARTS       "  << "\n";
        amrex::Print() << "------------------------------------"  << "\n";
        amrwind.Evolve_MB(step,1);
        if (do_two_way_coupling) {
          amrex::Print() << '\n';
          amrex::Print() << "------------------------------------"  << "\n";
          amrex::Print() << "           FILLPATCH A->E           "  << "\n";
          amrex::Print() << "------------------------------------"  << "\n";
          FillPatchBlocksAE ();
        }
        amrex::Print() << '\n';
        amrex::Print() << "------------------------------------"  << "\n";
        amrex::Print() << "          COMPLETE                  "  << "\n";
        amrex::Print() << "------------------------------------"  << "\n";
        amrex::Print() << "\n";
    }
}


// Fill AMR-Wind Boundary Regsiter from ERF1
void MultiBlockContainer::CopyERFtoAMRWindBoundaryReg (amrex::BndryRegister& receive_br,
                                                       amrex::Orientation ori,
                                                       amrex::Real time,
                                                       const std::string &field) {

  // Need a ghost cell in case AMR-Wind boundary to be filled coincides with ERF boundary
  amrex::IntVect nghost(1);

  // ERF level where we are getting the data
  const int erf_source_level = 0;

  /*
  // WARNING: for this to work properly we need to make sure the new state data is FillPatched
  //          old data is FillPatched at beginning of timestep and should be good
  amrex::Vector<amrex::MultiFab>* erf_data;
  if (time == erf1.get_t_old() ) {
    erf_data = &erf1.vars_old[erf_source_level];
  } else if (time == erf1.get_t_new() ) {
    erf_data = &erf1.vars_new[erf_source_level];
  } else {
    amrex::Abort("MultiBlockContainer::CopyERFtoAMRWindBoundaryReg invalid time requested");
  }
  
  if (field == "temperature") {
    amrex::NonLocalBC::MultiBlockCommMetaData cmd =
      amrex::NonLocalBC::MultiBlockCommMetaData(receive_br[ori].multiFab(), bv_afrome[ori],
                                                (*erf_data)[Vars::cons], nghost, dtos_afrome);

    // RHOTheta -> Theta (FIXME bad, should only do for relevant cells, etc)
    amrex::MultiFab::Divide((*erf_data)[Vars::cons], (*erf_data)[Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);

    
    amrex::Abort("I AM JUST STOPPING HERE FOR FUN");
    
    // Copy data
    amrex::NonLocalBC::ParallelCopy(receive_br[ori].multiFab(), (*erf_data)[Vars::cons],
                                    cmd, Cons::RhoScalar, 0, 1, dtos_afrome);

    // Theta -> RhoTheta (FIXME bad, should only do for relevant cells, etc)
    amrex::MultiFab::Multiply((*erf_data)[Vars::cons], (*erf_data)[Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);
  }
  */
  
  // WARNING: for this to work properly we need to make sure the new state data is FillPatched
  //          old data is FillPatched at beginning of timestep and should be good
  //  amrex::Vector<amrex::MultiFab>& erf_data;
  bool on_old_time{time == erf1.get_t_old()};
  bool on_new_time{time == erf1.get_t_new()};
  AMREX_ALWAYS_ASSERT(on_new_time || on_old_time);
  amrex::Vector<amrex::MultiFab>& erf_data = on_old_time ? erf1.vars_old[erf_source_level] : erf1.vars_new[erf_source_level];
  if (field == "temperature") {
    if (on_old_time) {
      amrex::Print() << std::endl << "FILLPATCHING FROM ERF TO AMR WIND ON OLD TIME, ORIENTATION " << ori << std::endl << std::endl;
    }
    else {
      amrex::Print() << std::endl << "FILLPATCHING FROM ERF TO AMR WIND ON NEW TIME, ORIENTATION " << ori << std::endl << std::endl;
    }

    // JUNK for selecting only subdomain

    std::map<int,int> mfmap;
    amrex::Vector<int> new_dm;
    amrex::BoxList new_bl;
    const amrex::BoxArray& old_ba = erf_data[Vars::cons].boxArray();
    const amrex::DistributionMapping& old_dm = erf_data[Vars::cons].DistributionMap();
    for (int i = 0; i < old_ba.size(); i++) {
      amrex::Box isect = old_ba[i] & eboxvec_afrome[ori];
      // ISSUES ISSUES TISSUES
      // when to multiply/divide by density
      
      amrex::Print() << "SETUP " << amrex::ParallelDescriptor::MyProc() << " " << isect  << isect.ok() << std::endl;
      if (isect.ok()) {
        new_bl.push_back(isect);
        new_dm.push_back(old_dm[i]);
        mfmap.insert({new_dm.size()-1, i});
      } 
    }
    amrex::BoxArray new_ba(new_bl);
    amrex::DistributionMapping new_new_dm(new_dm);
    amrex::MultiFab newmf(new_ba, new_new_dm, 1, 0);

    amrex::Print() << "OLD DM" << old_dm << std::endl;
    amrex::Print() << "NEW DM" << new_new_dm << std::endl;
    
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
      for ( amrex::MFIter mfi(newmf,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Array4<const amrex::Real> old_data = erf_data[Vars::cons].const_array(mfmap[mfi.index()]);
        std::cout << "MFITER " << " " << amrex::ParallelDescriptor::MyProc() << " " << mfi.tilebox() << std::endl;
      }
    }
    
    // END JUNK
    
    amrex::NonLocalBC::MultiBlockCommMetaData cmd =
      amrex::NonLocalBC::MultiBlockCommMetaData(receive_br[ori].multiFab(), aboxvec_afrome[ori],
                                                erf_data[Vars::cons], nghost, dtos_afrome);

    // RHOTheta -> Theta (FIXME bad, should only do for relevant cells, etc)
    amrex::MultiFab::Divide(erf_data[Vars::cons], erf_data[Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);

    // amrex::Abort("I AM JUST STOPPING HERE FOR FUN");
    
    // Copy data
    amrex::NonLocalBC::ParallelCopy(receive_br[ori].multiFab(), erf_data[Vars::cons],
                                    cmd, Cons::RhoScalar, 0, 1, dtos_afrome );

    // Theta -> RhoTheta (FIXME bad, should only do for relevant cells, etc)
    amrex::MultiFab::Multiply(erf_data[Vars::cons], erf_data[Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);
  }
  
}

void
MultiBlockContainer::PopulateErfTimesteps (amrex::Real* tsteps) {
  tsteps[0] = erf1.get_t_old();
  tsteps[1] = erf1.get_t_new();
}

// Wrapper for ParallelCopy between ERF and AMRWIND
void
MultiBlockContainer::FillPatchBlocksAE()
{
  // Create temporary multifabs to store data from AMR-Wind 
  // Note: AMR-Wind data is cell centered for all variables, and we only care about the valid data
  //       so no ghost cells for Temp or Dens. But we will be interpolating velocity to faces so we
  //       need one ghost cell for velocity
  // TODO: make these only correspond to the needed region (box_efroma) rather than the full ERF domain
  const amrex::BoxArray& ba            = erf1.vars_new[0][Vars::cons].boxArray();
  const amrex::DistributionMapping& dm = erf1.vars_new[0][Vars::cons].DistributionMap();
  amrex::MultiFab Temp_AW{ba, dm, 1, 0};
  amrex::MultiFab Dens_AW{ba, dm, 1, 0};  
  amrex::MultiFab Vel_AW{ba, dm, 3, 1};

  // Bring AMR-Wind data to temporary multifabs
  amrex::NonLocalBC::ParallelCopy(Temp_AW, amrwind.repo().get_field("temperature")(0),
                                  *(cmd_efroma[0]), 0, 0, 1, dtos_efroma);
  amrex::NonLocalBC::ParallelCopy(Dens_AW, amrwind.repo().get_field("density")(0),
                                  *(cmd_efroma[0]), 0, 0, 1, dtos_efroma);
  amrex::NonLocalBC::ParallelCopy(Vel_AW, amrwind.repo().get_field("velocity")(0),
                                  *(cmd_efroma[0]), 0, 0, 3, dtos_efroma);

  // Compute ERF variables from AMR-Wind variables and store in ERF data structures

  // Cell centered quantities
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(erf1.vars_new[0][Vars::cons]); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
    // compute cell-centered scalar and velocities
    auto cons_arr = erf1.vars_new[0][Vars::cons][mfi].array();
    auto vel_aw_arr = Vel_AW[mfi].array();
    auto dens_aw_arr = Dens_AW[mfi].array();
    auto temp_aw_arr = Temp_AW[mfi].array();
    amrex::Box ibox = box & ebox_efroma; // intersection of boxes
    amrex::ParallelFor(ibox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                // Save AMR-Wind Scalar into ERF data
                                cons_arr(i,j,k,Cons::RhoScalar) = cons_arr(i,j,k,Cons::Rho)*temp_aw_arr(i,j,k);
                                // Cell-centered velocity using energy preserving correction from Sprague & Satkauskas 2015
                                amrex::Real dens_correction = std::sqrt(dens_aw_arr(i,j,k)/cons_arr(i,j,k,Cons::Rho));
                                vel_aw_arr(i,j,k,0) *= dens_correction;
                                vel_aw_arr(i,j,k,1) *= dens_correction;
                                vel_aw_arr(i,j,k,2) *= dens_correction;
                            });
  }

  // We need to fill the interior boundary cells of velocity so we can interpolate to faces
  // Notes: Because the AMR-Wind domain is assumed to not exceed the extent of the ERF domain
  //         and we will only be filling velocity at face centers within (rather than on the
  //         boundary of) the AMR-Wind domain, we can get away with just doing a FillBoundary
  //        This fills all boundary cells but in theory we only need to fill boundary cells
  //         surrounding ebox_efroma
  Vel_AW.FillBoundary(erf1.Geom(0).periodicity());

  // Move cell centered velocity to faces and store for ERF
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(erf1.vars_new[0][Vars::cons]); mfi.isValid(); ++mfi) {
    auto vel_aw_arr = Vel_AW[mfi].array();
    // interpolate corrected velocities to face centers (only on interior)
    auto velx_arr = erf1.vars_new[0][Vars::xvel][mfi].array();
    amrex::Box xbox = mfi.nodaltilebox(0) & surroundingNodes(ebox_efroma,0).grow(0,-1);
    amrex::ParallelFor(xbox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                               velx_arr(i,j,k) = 0.5*(vel_aw_arr(i,j,k,0) + vel_aw_arr(i-1,j,k,0));
                             });
    
    auto vely_arr = erf1.vars_new[0][Vars::yvel][mfi].array();
    amrex::Box ybox =  mfi.nodaltilebox(1) & surroundingNodes(ebox_efroma,1).grow(1,-1);
    amrex::ParallelFor(ybox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                               vely_arr(i,j,k) = 0.5*(vel_aw_arr(i,j,k,1) + vel_aw_arr(i,j-1,k,1));
                             });
    
    auto velz_arr = erf1.vars_new[0][Vars::zvel][mfi].array();
    amrex::Box zbox =  mfi.nodaltilebox(2) & surroundingNodes(ebox_efroma,2).grow(2,-1);
    amrex::ParallelFor(zbox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                               velz_arr(i,j,k) = 0.5*(vel_aw_arr(i,j,k,2) + vel_aw_arr(i,j,k-1,2));
                             });
  }
}
