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
    dtos.permutation = amrex::IntVect{AMREX_D_DECL(   0,   1,   2)};
    dtos.sign        = amrex::IntVect{AMREX_D_DECL(   1,   1,   1)};

    // Set offset of dtos (NOTE: i_dst =  i_src - i_off -> [0] - [1])
    {
      amrex::Real dx = ( rb_v[0].hi(0) - rb_v[0].lo(0) ) / n_cell_in_v[0][0];
      amrex::Real dy = ( rb_v[0].hi(1) - rb_v[0].lo(1) ) / n_cell_in_v[0][1];
      amrex::Real dz = ( rb_v[0].hi(2) - rb_v[0].lo(2) ) / n_cell_in_v[0][2];
      int offx = amrex::Math::floor(( rb_v[0].lo(0) - rb_v[1].lo(0) ) / dx);
      int offy = amrex::Math::floor(( rb_v[0].lo(1) - rb_v[1].lo(1) ) / dy);
      int offz = amrex::Math::floor(( rb_v[0].lo(2) - rb_v[1].lo(2) ) / dz);
      // DEBUG
      //offx=0; offy=0; offz=0;
      dtos.offset = amrex::IntVect{AMREX_D_DECL(offx, offy, offz)};
    }

    // Set the permutation/sign of dtos
    dtos_etoa.permutation = amrex::IntVect{AMREX_D_DECL(   0,   1,   2)};
    dtos_etoa.sign        = amrex::IntVect{AMREX_D_DECL(   1,   1,   1)};
    dtos_atoe.permutation = amrex::IntVect{AMREX_D_DECL(   0,   1,   2)};
    dtos_atoe.sign        = amrex::IntVect{AMREX_D_DECL(   1,   1,   1)};

    // Set offset of dtos (NOTE: i_dst =  i_src - i_off -> [0] - [1])
    {
      const amrex::Geometry geom = amrwind.repo().mesh().Geom(0);
      amrex::Real dx = geom.CellSize(0);
      amrex::Real dy = geom.CellSize(1);
      amrex::Real dz = geom.CellSize(2);
      int offx = amrex::Math::floor((geom.ProbLo(0) - rb_v[0].lo(0)) / dx);
      int offy = amrex::Math::floor((geom.ProbLo(1) - rb_v[0].lo(1)) / dy);
      int offz = amrex::Math::floor((geom.ProbLo(2) - rb_v[0].lo(2)) / dz);
      dtos_etoa.offset = amrex::IntVect{AMREX_D_DECL(offx, offy, offz)};
      dtos_atoe.offset = amrex::IntVect{AMREX_D_DECL(-offx, -offy, -offz)};
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
    // Hard-coded bounds for now
    int nvars  = erf1.vars_new[0].size();
    int ndirs  = AMREX_SPACEDIM;

    amrex::Box awbox = amrwind.repo().mesh().Geom(0).Domain();
    amrex::Dim3 send = dtos_etoa.Inverse(amrex::lbound(awbox));
    amrex::Dim3 bend = dtos_etoa.Inverse(amrex::ubound(awbox));
    box_etoa  = amrex::Box({send.x, send.y, send.z}, {bend.x, bend.y, bend.z});

    for (int i(0); i<nvars; ++i) {
      for (amrex::OrientationIter oit; oit != nullptr; ++oit) {
        auto ori = oit();
        amrex::IntVect se = awbox.smallEnd();
        amrex::IntVect be = awbox.bigEnd();
        if (ori.faceDir() == 1) {
          se[ori.coordDir()] = be[ori.coordDir()];
          be[ori.coordDir()] += 1;
        }
        else {
          be[ori.coordDir()] = se[ori.coordDir()];
          se[ori.coordDir()] -= 1;
        }
        amrex::Box bx(se,be);
        blv_atoe.push_back(bx);
      }
    }
    
    /*
    // DEBUG BOX LIST
    for (int i(0); i<nvars; ++i) {
        amrex::Print() << "DOM: " << erf2.domain_p[i] << "\n";
        amrex::Print() << "BA: " << erf2.vars_new[0][i].boxArray() << "\n";
        for (int j(0); j<6; ++j)
            amrex::Print() << (blv[i].data())[j] << "\n";

        amrex::Print() << "\n";
    }
    exit(0);
    */

}

// Set up MB Communication Meta Data
void
MultiBlockContainer::SetBlockCommMetaData()
{
    // Hard-coded bounds for now
    int nvars  = erf1.vars_new[0].size(); // Destination MF
    int ndirs  = AMREX_SPACEDIM;

    amrex::IntVect nghost(0);
    amrex::NonLocalBC::MultiBlockCommMetaData *cmd_etoa_tmp =
      new amrex::NonLocalBC::MultiBlockCommMetaData(erf1.vars_new[0][Vars::cons], box_etoa,
                                                    amrwind.repo().get_field("temperature")(0), nghost, dtos_etoa);
    cmd_etoa.push_back(cmd_etoa_tmp);
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
void MultiBlockContainer::CopyToBoundaryRegister (amrex::BndryRegister& receive_br_old, amrex::BndryRegister& receive_br_new, amrex::Orientation ori) {

  // Need a ghost cell in case AMR-Wind boundary to be filled coincides with ERF boundary
  amrex::IntVect nghost(1);

  // WARNING: for this to work properly we need to make sure the new state data is FillPatched
  //          old data is FillPatched at beginning of timestep and should be good
  amrex::Vector<amrex::MultiFab>* erf_data;
  erf_data = &erf1.vars_old[0];
  std::cout << (*erf_data)[Vars::cons].min(Cons::Rho);
  
  amrex::NonLocalBC::MultiBlockCommMetaData *cmd_old =
    new amrex::NonLocalBC::MultiBlockCommMetaData(receive_br_old[ori].multiFab(), blv_atoe[ori],
                                                  erf1.vars_old[0][Vars::cons], nghost, dtos_atoe);

  // RHOTheta -> Theta (FIXME bad, should only do for relevant cells, etc)
  amrex::MultiFab::Divide(erf1.vars_old[0][Vars::cons], erf1.vars_old[0][Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);

  // Copy data
  amrex::NonLocalBC::ParallelCopy(receive_br_old[ori].multiFab(), erf1.vars_old[0][Vars::cons],
                                  *cmd_old, Cons::RhoScalar, 0, 1, dtos_atoe);

  // Theta -> RhoTheta (FIXME bad, should only do for relevant cells, etc)
  amrex::MultiFab::Multiply(erf1.vars_old[0][Vars::cons], erf1.vars_old[0][Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);
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
  //       so no ghost cells for Temp of Dens. But we will be interpolating velocity to faces so we
  //       need one ghost cell for velocity
  // TODO: make these only correspond to the needed region (box_etoa) rather than the full ERF domain
  const amrex::BoxArray& ba            = erf1.vars_new[0][Vars::cons].boxArray();
  const amrex::DistributionMapping& dm = erf1.vars_new[0][Vars::cons].DistributionMap();
  amrex::MultiFab Temp_AW{ba, dm, 1, 0};
  amrex::MultiFab Dens_AW{ba, dm, 1, 0};  
  amrex::MultiFab Vel_AW{ba, dm, 3, 1};

  // Bring AMR-Wind data to temporary multifabs
  amrex::NonLocalBC::ParallelCopy(Temp_AW, amrwind.repo().get_field("temperature")(0),
                                  *(cmd_etoa[0]), 0, 0, 1, dtos_etoa);
  amrex::NonLocalBC::ParallelCopy(Dens_AW, amrwind.repo().get_field("density")(0),
                                  *(cmd_etoa[0]), 0, 0, 1, dtos_etoa);
  amrex::NonLocalBC::ParallelCopy(Vel_AW, amrwind.repo().get_field("velocity")(0),
                                  *(cmd_etoa[0]), 0, 0, 3, dtos_etoa);

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
    amrex::Box ibox = box & box_etoa; // intersection of boxes
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
  //         surrounding box_etoa
  Vel_AW.FillBoundary(erf1.Geom(0).periodicity());

  // Move cell centered velocity to faces and store for ERF
#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(erf1.vars_new[0][Vars::cons]); mfi.isValid(); ++mfi) {
    auto vel_aw_arr = Vel_AW[mfi].array();
    // interpolate corrected velocities to face centers (only on interior)
    auto velx_arr = erf1.vars_new[0][Vars::xvel][mfi].array();
    amrex::Box xbox = mfi.nodaltilebox(0) & surroundingNodes(box_etoa,0).grow(0,-1);
    amrex::ParallelFor(xbox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                               velx_arr(i,j,k) = 0.5*(vel_aw_arr(i,j,k,0) + vel_aw_arr(i-1,j,k,0));
                             });
    
    auto vely_arr = erf1.vars_new[0][Vars::yvel][mfi].array();
    amrex::Box ybox =  mfi.nodaltilebox(1) & surroundingNodes(box_etoa,1).grow(1,-1);
    amrex::ParallelFor(ybox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                               vely_arr(i,j,k) = 0.5*(vel_aw_arr(i,j,k,1) + vel_aw_arr(i,j-1,k,1));
                             });
    
    auto velz_arr = erf1.vars_new[0][Vars::zvel][mfi].array();
    amrex::Box zbox =  mfi.nodaltilebox(2) & surroundingNodes(box_etoa,2).grow(2,-1);
    amrex::ParallelFor(zbox, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                               velz_arr(i,j,k) = 0.5*(vel_aw_arr(i,j,k,2) + vel_aw_arr(i,j,k-1,2));
                             });
  }
}
