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
  erf1(rb_v[0],max_level_in_v[0],n_cell_in_v[0],coord_v[0],ref_ratios_v[0],is_per_v[0],prefix_v[0]),
  erf2(rb_v[1],max_level_in_v[1],n_cell_in_v[1],coord_v[1],ref_ratios_v[1],is_per_v[1],prefix_v[1])
{
    // Store ptr to container to call member functions
    amrwind.SetMultiBlockPointer(this);
    erf1.SetMultiBlockPointer(this);
    erf2.SetMultiBlockPointer(this);

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
  erf1.InitData();
  amrex::Print() << "       AMRWIND INITIALIZATION       "  << "\n";
  amrex::Print() << "------------------------------------"  << "\n";
  amrwind.InitData();
  amrex::Print() << '\n';
  amrex::Print() << "         ERF2 INITIALIZATION        "  << "\n";
  amrex::Print() << "------------------------------------"  << "\n";
  erf2.InitData();
  amrex::Print() << '\n';
  amrex::Print() << "COMPLETE" << "\n";
  amrex::Print() << "\n";
  amrex::Print() << "------------------------------------"  << "\n";
}

void
MultiBlockContainer::SetBoxListsAtoE()
{
  const int nvars = 4;
  const int ndirs  = AMREX_SPACEDIM;
}


// Set up BoxList vector for use with Communication Meta Data
void
MultiBlockContainer::SetBoxLists()
{
    // Hard-coded bounds for now
    int nvars  = erf2.vars_new[0].size();
    int ndirs  = AMREX_SPACEDIM;

    for (int i(0); i<nvars; ++i) {
        // Get ghost cells, domain & grown box
        amrex::IntVect nghost = erf2.vars_new[0][i].nGrowVect();
        amrex::Box dom = erf2.domain_p[i];
        amrex::Box gbx = grow(dom,nghost);
        // Tmp BoxList
        amrex::BoxList bl;
        bl.clear();
        bl.set(gbx.ixType());
        for (int j(0); j<ndirs; ++j) {
            // Local box copies
            amrex::Box lgbx(gbx);
            amrex::Box ugbx(gbx);
            // Get lower & upper bound
            int se = dom.smallEnd(j) - 1;
            int be = dom.bigEnd(j)   + 1;
            // Modify bounds for nodal vars
            if (gbx.ixType().nodeCentered(j)) {
                se += 1;
                be -= 1;
            }
            // Populate BoxList
            bl.push_back( lgbx.setBig(j,se) );
            bl.push_back( ugbx.setSmall(j,be) );
        }
        blv.push_back(bl);
        blv_full.push_back(gbx);

        amrex::Box awbox = amrwind.repo().mesh().Geom(0).Domain();
        amrex::Dim3 send = dtos_etoa.Inverse(amrex::lbound(awbox));
        amrex::Dim3 bend = dtos_etoa.Inverse(amrex::ubound(awbox));
        amrex::Box awbox_erf({send.x, send.y, send.z}, {bend.x, bend.y, bend.z});
        blv_etoa.push_back(awbox_erf);

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
    int nvars  = erf2.vars_new[0].size(); // Destination MF
    int ndirs  = AMREX_SPACEDIM;

    // Loop over num_vars to set communicator
    for (int i(0); i<nvars; ++i) {
        // Make space
        cmd.push_back(std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*>());
        // Get ghost cell vector for multifab growth
        amrex::IntVect nghost = erf2.vars_new[0][i].nGrowVect();
        for (int j(0); j<2*ndirs; ++j) {

            // Store temp ptr to communicator for i^th variable
            amrex::NonLocalBC::MultiBlockCommMetaData *cmd_tmp =
                new amrex::NonLocalBC::MultiBlockCommMetaData(erf2.vars_new[0][i], (blv[i].data())[j],
                                                              erf1.vars_new[0][i], nghost, dtos);
            // Populate cmd vector
            cmd[i].push_back(cmd_tmp);
        }
        amrex::NonLocalBC::MultiBlockCommMetaData *cmd_full_tmp =
                new amrex::NonLocalBC::MultiBlockCommMetaData(erf2.vars_new[0][i], blv_full[i],
                                                              erf1.vars_old[0][i], nghost, dtos);
        cmd_full.push_back(cmd_full_tmp);
    }

    amrex::IntVect nghost(0);
    amrex::NonLocalBC::MultiBlockCommMetaData *cmd_etoa_tmp =
      new amrex::NonLocalBC::MultiBlockCommMetaData(erf1.vars_new[0][Vars::cons], blv_etoa[0],
                                                    amrwind.repo().get_field("temperature")(0), nghost, dtos_etoa);
    cmd_etoa.push_back(cmd_etoa_tmp);
}

// Advance blocks
void
MultiBlockContainer::AdvanceBlocks()
{
    amrex::Print() << "STARTING MAIN DRIVER FOR: " << m_max_step << " STEPS" << "\n";
    amrex::Print() << "\n";

    FillPatchBlocksAE();
    FillPatchBlocksFull();

    for (int step(1); step <= m_max_step; ++step) {
        amrex::Print() << "    STARTING ADVANCE DRIVER: " << step << "\n";
        amrex::Print() << "===================================="  << "\n";
        erf1.Evolve_MB(step,1);
        amrex::Print() << '\n';
        amrex::Print() << "        SECOND BLOCK STARTS         "  << "\n";
        amrex::Print() << "------------------------------------"  << "\n";
        erf2.Evolve_MB(step,1);
        amrex::Print() << '\n';
        amrex::Print() << "        AMRWIND BLOCK STARTS        "  << "\n";
        amrex::Print() << "------------------------------------"  << "\n";
        amrwind.Evolve_MB(step,1);
        amrex::Print() << '\n';
        amrex::Print() << "           FILLPATCH A->E           "  << "\n";
        amrex::Print() << "------------------------------------"  << "\n";
        FillPatchBlocksAE ();
        amrex::Print() << "COMPLETE" << "\n";
        amrex::Print() << "\n";

    }
}


// Fill AMR-Wind Boundary Regsiter from ERF1
void MultiBlockContainer::CopyToBoundaryRegister (amrex::BndryRegister& receive_br_old, amrex::BndryRegister& receive_br_new, amrex::Orientation ori) {

  // Need a ghost cell in case AMR-Wind boundary to be filled coincides with ERF boundary
  amrex::IntVect nghost(1);

  std::cout << "TEST TEST TEST" << std::endl
            << " xvel " <<  erf1.vars_old[0][Vars::xvel].min(0) << " " << erf1.vars_old[0][Vars::xvel].max(0)
            << " yvel " <<  erf1.vars_old[0][Vars::yvel].min(0) << " " << erf1.vars_old[0][Vars::yvel].max(0)
            << " zvel " <<  erf1.vars_old[0][Vars::zvel].min(0) << " " << erf1.vars_old[0][Vars::zvel].max(0)
            << std::endl;

  amrex::NonLocalBC::MultiBlockCommMetaData *cmd_old =
    new amrex::NonLocalBC::MultiBlockCommMetaData(receive_br_old[ori].multiFab(), blv_atoe[ori],
                                                  erf1.vars_old[0][Vars::cons], nghost, dtos_atoe);
  amrex::NonLocalBC::MultiBlockCommMetaData *cmd_new =
    new amrex::NonLocalBC::MultiBlockCommMetaData(receive_br_new[ori].multiFab(), blv_atoe[ori],
                                                  erf1.vars_new[0][Vars::cons], nghost, dtos_atoe);

  // RHOTheta -> Theta (FIXME bad, should only do for relevant cells, etc)
  amrex::MultiFab::Divide(erf1.vars_old[0][Vars::cons], erf1.vars_old[0][Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);
  amrex::MultiFab::Divide(erf1.vars_new[0][Vars::cons], erf1.vars_new[0][Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);

  // Copy data
  amrex::NonLocalBC::ParallelCopy(receive_br_old[ori].multiFab(), erf1.vars_old[0][Vars::cons],
                                  *cmd_old, Cons::RhoScalar, 0, 1, dtos_atoe);
  amrex::NonLocalBC::ParallelCopy(receive_br_new[ori].multiFab(), erf1.vars_new[0][Vars::cons],
  *cmd_new, Cons::RhoScalar, 0, 1, dtos_atoe);

  // RHOTheta -> Theta (FIXME bad, should only do for relevant cells, etc)
  amrex::MultiFab::Multiply(erf1.vars_old[0][Vars::cons], erf1.vars_old[0][Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);
  amrex::MultiFab::Multiply(erf1.vars_new[0][Vars::cons], erf1.vars_new[0][Vars::cons], Cons::Rho, Cons::RhoScalar, 1, nghost);
}

void
MultiBlockContainer::PopulateErfTimesteps (amrex::Real* tsteps) {
  tsteps[0] = erf1.get_t_old();
  tsteps[1] = erf1.get_t_new();
}


// Wrapper for ParallelCopy between classes
void
MultiBlockContainer::FillPatchBlocks(int src_ind, int dst_ind)
{
    // Hard-coded bounds for now
    int ndirs  = AMREX_SPACEDIM;

    // Loop faces of box to perform ParallelCopy
    // NOTE - cmd built with ERF2 so uses dst_ind
    for (int j(0); j<2*ndirs; ++j)
        amrex::NonLocalBC::ParallelCopy(erf2.vars_new[0][dst_ind], erf1.vars_new[0][src_ind],
                                        *(cmd[dst_ind][j]), 0, 0, 1, dtos);
}

// Wrapper for ParallelCopy between classes
void
MultiBlockContainer::FillPatchBlocksFull()
{
  amrex::NonLocalBC::ParallelCopy(erf2.vars_new[0][Vars::cons], erf1.vars_old[0][Vars::cons],
                                *(cmd_full[Vars::cons]), 0, 0, Cons::NumVars, dtos);

  amrex::NonLocalBC::ParallelCopy(erf2.vars_new[0][Vars::xvel], erf1.vars_old[0][Vars::xvel],
                                    *(cmd_full[Vars::xvel]), 0, 0, 1, dtos);
  amrex::NonLocalBC::ParallelCopy(erf2.vars_new[0][Vars::yvel], erf1.vars_old[0][Vars::yvel],
                                    *(cmd_full[Vars::yvel]), 0, 0, 1, dtos);
  amrex::NonLocalBC::ParallelCopy(erf2.vars_new[0][Vars::zvel], erf1.vars_old[0][Vars::zvel],
                                    *(cmd_full[Vars::zvel]), 0, 0, 1, dtos);
}

// Wrapper for ParallelCopy between ERF and AMRWIND
void
MultiBlockContainer::FillPatchBlocksAE()
{
  std::cout << "STARTING AE FILLPATCH" << std::endl;
  std::cout << "VALUES "
    //<< erf1.vars_new[0][Vars::cons].min(Cons::RhoTheta) << " "
    //      << erf1.vars_new[0][Vars::cons].min(Cons::RhoTheta) << " "
            << erf1.vars_new[0][Vars::cons].max(Cons::RhoScalar) << " "
            << erf1.vars_new[0][Vars::cons].max(Cons::RhoScalar) << " "
            << amrwind.repo().get_field("temperature")(0).min(0) << " "
            << amrwind.repo().get_field("temperature")(0).max(0) << " "
            << std::endl;
  amrex::NonLocalBC::ParallelCopy(erf1.vars_new[0][Vars::cons], amrwind.repo().get_field("temperature")(0),
                                  *(cmd_etoa[Vars::cons]), 0, Cons::RhoScalar, 1, dtos_etoa);
  //  *(cmd_etoa[Vars::cons]), 0, Cons::RhoTheta, 1, dtos_etoa);

#ifdef _OPENMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(erf1.vars_new[0][Vars::cons]); mfi.isValid(); ++mfi) {
    const amrex::Box& box = mfi.validbox();
        auto cons_arr = erf1.vars_new[0][Vars::cons][mfi].array();

        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                                  if (blv_etoa[0].contains(i,j,k)) {
                                    // TODO: FIXME: HACK
                                    /*
                                    amrex::Real rhotheta0= 300.0*1.161440171;
                                    cons_arr(i,j,k,Cons::Rho) = rhotheta0/cons_arr(i,j,k,Cons::RhoTheta);
                                    cons_arr(i,j,k,Cons::RhoTheta) = rhotheta0; */
                                    // cons_arr(i,j,k,Cons::RhoTheta) *= cons_arr(i,j,k,Cons::Rho);
                                    cons_arr(i,j,k,Cons::RhoScalar) *= cons_arr(i,j,k,Cons::Rho);
                                    }
                                });
    }

  std::cout << "FINISHED AE FILLPATCH" << std::endl;
  std::cout << "VALUES "
    //<< erf1.vars_new[0][Vars::cons].min(Cons::RhoTheta) << " "
    //        << erf1.vars_new[0][Vars::cons].max(Cons::RhoTheta) << " "
    << erf1.vars_new[0][Vars::cons].min(Cons::RhoScalar) << " "
     << erf1.vars_new[0][Vars::cons].max(Cons::RhoScalar) << " "
            << amrwind.repo().get_field("temperature")(0).min(0) << " "
            << amrwind.repo().get_field("temperature")(0).max(0) << " "
            << std::endl;

}
/*
void
MultiBlockContainer::PopulateErfTimesteps (amrex::Vector<amrex::Real> tsteps) {
  erf1.CoarseTimeState(tsteps);
}
*/
