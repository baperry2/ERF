/**
 * \file ERF_init1d.cpp
 */
#include <ERF_EOS.H>
#include <ERF.H>
#include <ERF_TileNoZ.H>
#include <ERF_prob_common.H>
#include <ERF_ParFunctions.H>
#include <ERF_Utils.H>

#include <ERF_Interpolation_1D.H>

using namespace amrex;

/**
 * Initialization function for host and device vectors
 * used to store averaged quantities when calculating
 * the effects of Rayleigh Damping.
 */
void
ERF::initRayleigh ()
{
    const int khi = geom[0].Domain().bigEnd(2);
    solverChoice.rayleigh_ztop = (solverChoice.use_terrain) ? zlevels_stag[0][khi+1] : geom[0].ProbHi(2);

    h_rayleigh_ptrs.resize(max_level+1);
    d_rayleigh_ptrs.resize(max_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        // These have 4 components: ubar, vbar, wbar, thetabar
        h_rayleigh_ptrs[lev].resize(Rayleigh::nvars);
        d_rayleigh_ptrs[lev].resize(Rayleigh::nvars);

        const int zlen_rayleigh = geom[lev].Domain().length(2);

        // Allocate space for these 1D vectors
        for (int n = 0; n < Rayleigh::nvars; n++) {
            h_rayleigh_ptrs[lev][n].resize(zlen_rayleigh, 0.0_rt);
            d_rayleigh_ptrs[lev][n].resize(zlen_rayleigh, 0.0_rt);
        }

        // Init the host vectors
        prob->erf_init_rayleigh(h_rayleigh_ptrs[lev], geom[lev], z_phys_nd[lev], solverChoice.rayleigh_zdamp);

        // Copy from host vectors to device vectors
        for (int n = 0; n < Rayleigh::nvars; n++) {
            Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][n].begin(), h_rayleigh_ptrs[lev][n].end(),
                      d_rayleigh_ptrs[lev][n].begin());
        }
    }
}

/**
 * Sets the Rayleigh Damping averaged quantities from an
 * externally supplied input sounding data file.
 *
 * @param[in] restarting Boolean parameter that indicates whether
                         we are currently restarting from a checkpoint file.
 */
void
ERF::setRayleighRefFromSounding (bool restarting)
{
    // If we are restarting then we haven't read the input_sounding file yet
    //    so we need to read it here
    // TODO: should we store this information in the checkpoint file instead?
    if (restarting) {
        input_sounding_data.resize_arrays();
        for (int n = 0; n < input_sounding_data.n_sounding_files; n++) {
            input_sounding_data.read_from_file(geom[0], zlevels_stag[0], n);
        }
    }

    const Real* z_inp_sound     = input_sounding_data.z_inp_sound[0].dataPtr();
    const Real* U_inp_sound     = input_sounding_data.U_inp_sound[0].dataPtr();
    const Real* V_inp_sound     = input_sounding_data.V_inp_sound[0].dataPtr();
    const Real* theta_inp_sound = input_sounding_data.theta_inp_sound[0].dataPtr();
    const int   inp_sound_size  = input_sounding_data.size(0);

    int refine_fac{1};
    for (int lev = 0; lev <= finest_level; lev++)
    {
        const int klo = geom[lev].Domain().smallEnd(2);
        const int khi = geom[lev].Domain().bigEnd(2);
        const int Nz = khi - klo + 1;

        Vector<Real> zcc(Nz);
        Vector<Real> zlevels_sub(zlevels_stag[0].begin()+klo/refine_fac,
                                 zlevels_stag[0].begin()+khi/refine_fac+2);
        expand_and_interpolate_1d(zcc, zlevels_sub, refine_fac, true);
#if 0
        amrex::AllPrint() << "lev="<<lev<<" : (refine_fac="<<refine_fac<<",klo="<<klo<<",khi="<<khi<<") ";
        for (int k = 0; k < zlevels_sub.size(); k++) { amrex::AllPrint() << zlevels_sub[k] << " "; }
        amrex::AllPrint() << " --> ";
        for (int k = 0; k < Nz; k++) { amrex::AllPrint() << zcc[k] << " "; }
        amrex::AllPrint() << std::endl;
#endif

        for (int k = 0; k < Nz; k++)
        {
            h_rayleigh_ptrs[lev][Rayleigh::ubar][k]     = interpolate_1d(z_inp_sound, U_inp_sound, zcc[k], inp_sound_size);
            h_rayleigh_ptrs[lev][Rayleigh::vbar][k]     = interpolate_1d(z_inp_sound, V_inp_sound, zcc[k], inp_sound_size);
            h_rayleigh_ptrs[lev][Rayleigh::wbar][k]     = Real(0.0);
            h_rayleigh_ptrs[lev][Rayleigh::thetabar][k] = interpolate_1d(z_inp_sound, theta_inp_sound, zcc[k], inp_sound_size);
        }

        // Copy from host version to device version
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][Rayleigh::ubar].begin(), h_rayleigh_ptrs[lev][Rayleigh::ubar].end(),
                         d_rayleigh_ptrs[lev][Rayleigh::ubar].begin());
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][Rayleigh::vbar].begin(), h_rayleigh_ptrs[lev][Rayleigh::vbar].end(),
                         d_rayleigh_ptrs[lev][Rayleigh::vbar].begin());
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][Rayleigh::wbar].begin(), h_rayleigh_ptrs[lev][Rayleigh::wbar].end(),
                         d_rayleigh_ptrs[lev][Rayleigh::wbar].begin());
        Gpu::copy(Gpu::hostToDevice, h_rayleigh_ptrs[lev][Rayleigh::thetabar].begin(), h_rayleigh_ptrs[lev][Rayleigh::thetabar].end(),
                         d_rayleigh_ptrs[lev][Rayleigh::thetabar].begin());

        refine_fac *= ref_ratio[lev][2];
    }
}

/**
 * Initialization function for host and device vectors
 * used to store the effects of sponge Damping.
 */
void
ERF::initSponge ()
{
    h_sponge_ptrs.resize(max_level+1);
    d_sponge_ptrs.resize(max_level+1);

    for (int lev = 0; lev <= finest_level; lev++)
    {
        // These have 2 components: ubar, vbar
        h_sponge_ptrs[lev].resize(Sponge::nvars_sponge);
        d_sponge_ptrs[lev].resize(Sponge::nvars_sponge);

        const int zlen_sponge = geom[lev].Domain().length(2);

        // Allocate space for these 1D vectors
        for (int n = 0; n < Sponge::nvars_sponge; n++) {
            h_sponge_ptrs[lev][n].resize(zlen_sponge, 0.0_rt);
            d_sponge_ptrs[lev][n].resize(zlen_sponge, 0.0_rt);
        }

    }
}

/**
 * Sets the sponge damping averaged quantities from an
 * externally supplied input sponge data file.
 *
 * @param[in] restarting Boolean parameter that indicates whether
                         we are currently restarting from a checkpoint file.
 */
void
ERF::setSpongeRefFromSounding (bool restarting)
{
    // If we are restarting then we haven't read the input_sponge file yet
    //    so we need to read it here
    // TODO: should we store this information in the checkpoint file instead?
    if (restarting) {
        input_sponge_data.read_from_file(geom[0], zlevels_stag[0]);
    }

    const Real* z_inp_sponge     = input_sponge_data.z_inp_sponge.dataPtr();
    const Real* U_inp_sponge     = input_sponge_data.U_inp_sponge.dataPtr();
    const Real* V_inp_sponge     = input_sponge_data.V_inp_sponge.dataPtr();
    const int   inp_sponge_size  = input_sponge_data.size();

    for (int lev = 0; lev <= finest_level; lev++)
    {
        const int khi = geom[lev].Domain().bigEnd()[2];
        Vector<Real> zcc(khi+1);

        if (z_phys_cc[lev]) {
            // use_terrain=1
            // calculate the damping strength based on the max height at each k
            reduce_to_max_per_height(zcc, z_phys_cc[lev]);
        } else {
            const auto *const prob_lo = geom[lev].ProbLo();
            const auto *const dx = geom[lev].CellSize();
            for (int k = 0; k <= khi; k++)
            {
                zcc[k] = prob_lo[2] + (k+0.5) * dx[2];
            }
        }

        for (int k = 0; k <= khi; k++)
        {
            h_sponge_ptrs[lev][Sponge::ubar_sponge][k] = interpolate_1d(z_inp_sponge, U_inp_sponge, zcc[k], inp_sponge_size);
            h_sponge_ptrs[lev][Sponge::vbar_sponge][k] = interpolate_1d(z_inp_sponge, V_inp_sponge, zcc[k], inp_sponge_size);
        }

        // Copy from host version to device version
        Gpu::copy(Gpu::hostToDevice, h_sponge_ptrs[lev][Sponge::ubar_sponge].begin(), h_sponge_ptrs[lev][Sponge::ubar_sponge].end(),
                         d_sponge_ptrs[lev][Sponge::ubar_sponge].begin());
        Gpu::copy(Gpu::hostToDevice, h_sponge_ptrs[lev][Sponge::vbar_sponge].begin(), h_sponge_ptrs[lev][Sponge::vbar_sponge].end(),
                         d_sponge_ptrs[lev][Sponge::vbar_sponge].begin());
    }
}

/**
 * Initialize density and pressure base state in
 * hydrostatic equilibrium.
 */
void
ERF::initHSE (int lev)
{
    // This integrates up through column to update p_hse, pi_hse;
    // r_hse is not const b/c FillBoundary is called at the end for r_hse and p_hse

    MultiFab r_hse (base_state[lev], make_alias, 0, 1); // r_0  is first  component
    MultiFab p_hse (base_state[lev], make_alias, 1, 1); // p_0  is second component
    MultiFab pi_hse(base_state[lev], make_alias, 2, 1); // pi_0 is third  component

    bool all_boxes_touch_bottom = true;
    Box domain(geom[lev].Domain());

    if (lev == 0) {
        BoxArray ba(base_state[lev].boxArray());
        for (int i = 0; i < ba.size(); i++) {
            if (ba[i].smallEnd(2) != domain.smallEnd(2)) {
                all_boxes_touch_bottom = false;
            }
        }
    }

    if (all_boxes_touch_bottom || lev > 0) {

        // Initial r_hse may or may not be in HSE -- defined in ERF_prob.cpp
        if (solverChoice.use_moist_background){
            prob->erf_init_dens_hse_moist(r_hse, z_phys_nd[lev], geom[lev]);
        } else {
            prob->erf_init_dens_hse(r_hse, z_phys_nd[lev], z_phys_cc[lev], geom[lev]);
        }

        erf_enforce_hse(lev, r_hse, p_hse, pi_hse, z_phys_cc[lev]);

    } else {

        BoxArray ba_new(domain);

        ChopGrids2D(ba_new, domain, ParallelDescriptor::NProcs());

        DistributionMapping dm_new(ba_new);

        MultiFab new_base_state(ba_new, dm_new, 3, 1);
        new_base_state.ParallelCopy(base_state[lev],0,0,3,1,1);

        MultiFab new_r_hse (new_base_state, make_alias, 0, 1); // r_0  is first  component
        MultiFab new_p_hse (new_base_state, make_alias, 1, 1); // p_0  is second component
        MultiFab new_pi_hse(new_base_state, make_alias, 2, 1); // pi_0 is third  component

        std::unique_ptr<MultiFab> new_z_phys_cc;
        std::unique_ptr<MultiFab> new_z_phys_nd;
        if (solverChoice.use_terrain) {
            new_z_phys_cc = std::make_unique<MultiFab>(ba_new,dm_new,1,1);
            new_z_phys_cc->ParallelCopy(*z_phys_cc[lev],0,0,1,1,1);

            BoxArray ba_new_nd(ba_new);
            ba_new_nd.surroundingNodes();
            new_z_phys_nd = std::make_unique<MultiFab>(ba_new_nd,dm_new,1,1);
            new_z_phys_nd->ParallelCopy(*z_phys_nd[lev],0,0,1,1,1);
        }

        // Initial r_hse may or may not be in HSE -- defined in ERF_prob.cpp
        if (solverChoice.use_moist_background){
            prob->erf_init_dens_hse_moist(new_r_hse, new_z_phys_nd, geom[lev]);
        } else {
            prob->erf_init_dens_hse(new_r_hse, new_z_phys_nd, new_z_phys_cc, geom[lev]);
        }

        erf_enforce_hse(lev, new_r_hse, new_p_hse, new_pi_hse, new_z_phys_cc);

        // Now copy back into the original arrays
        base_state[lev].ParallelCopy(new_base_state,0,0,3,1,1);
    }

    // Make sure to fill the ghost cells of the base state
    if (lev > 0) {
        // Interp all three components: rho, p, pi
        int  icomp = 0; int bccomp = 0; int ncomp = 3;

        PhysBCFunctNoOp null_bc;
        Interpolater* mapper = &cell_cons_interp;

        Real time = 0.;

        Vector<MultiFab*> fmf = {&base_state[lev  ], &base_state[lev  ]};
        Vector<MultiFab*> cmf = {&base_state[lev-1], &base_state[lev-1]};
        Vector<Real> ftime    = {time, time};
        Vector<Real> ctime    = {time, time};
        FillPatchTwoLevels(base_state[lev], time,
                           cmf, ctime, fmf, ftime,
                           icomp, icomp, ncomp, geom[lev-1], geom[lev],
                           null_bc, 0, null_bc, 0, refRatio(lev-1),
                           mapper, domain_bcs_type, bccomp);
    }
}

void
ERF::initHSE ()
{
    AMREX_ALWAYS_ASSERT(!init_sounding_ideal);
    for (int lev = 0; lev <= finest_level; lev++)
    {
        initHSE(lev);
    }
}

/**
 * Enforces hydrostatic equilibrium when using terrain.
 *
 * @param[in]  lev  Integer specifying the current level
 * @param[out] dens MultiFab storing base state density
 * @param[out] pres MultiFab storing base state pressure
 * @param[out] pi   MultiFab storing base state Exner function
 * @param[in]  z_cc Pointer to MultiFab storing cell centered z-coordinates
 */
void
ERF::erf_enforce_hse (int lev,
                      MultiFab& dens, MultiFab& pres, MultiFab& pi,
                      std::unique_ptr<MultiFab>& z_cc)
{
    Real l_gravity = solverChoice.gravity;
    bool l_use_terrain = solverChoice.use_terrain;

    const auto geomdata = geom[lev].data();
    const Real dz = geomdata.CellSize(2);

    const Box& domain = geom[lev].Domain();

    for ( MFIter mfi(dens, TileNoZ()); mfi.isValid(); ++mfi )
    {
        // Create a flat box with same horizontal extent but only one cell in vertical
        const Box& tbz = mfi.nodaltilebox(2);
        int klo = tbz.smallEnd(2);
        int khi = tbz.bigEnd(2);

        Box b2d = tbz; // Copy constructor
        b2d.grow(0,1); b2d.grow(1,1); // Grow by one in the lateral directions
        b2d.setRange(2,0);

        // We integrate to the first cell (and below) by using rho in this cell
        // If gravity == 0 this is constant pressure
        // If gravity != 0, hence this is a wall, this gives gp0 = dens[0] * gravity
        // (dens_hse*gravity would also be dens[0]*gravity because we use foextrap for rho at k = -1)
        // Note ng_pres_hse = 1

       // We start by assuming pressure on the ground is p_0 (in ERF_Constants.H)
       // Note that gravity is positive

        Array4<Real>  rho_arr = dens.array(mfi);
        Array4<Real> pres_arr = pres.array(mfi);
        Array4<Real>   pi_arr =   pi.array(mfi);
        Array4<Real> zcc_arr;
        if (l_use_terrain) {
           zcc_arr = z_cc->array(mfi);
        }

        const Real rdOcp = solverChoice.rdOcp;

        ParallelFor(b2d, [=] AMREX_GPU_DEVICE (int i, int j, int)
        {
            // Set value at surface from Newton iteration for rho
            if (klo == 0) {
                // Physical height of the terrain at cell center
                Real hz;
                if (l_use_terrain) {
                    hz = zcc_arr(i,j,klo);
                } else {
                    hz = 0.5*dz;
                }

                pres_arr(i,j,klo  ) = p_0 - hz * rho_arr(i,j,klo) * l_gravity;
                pi_arr(i,j,klo  ) = getExnergivenP(pres_arr(i,j,klo  ), rdOcp);

                // Set ghost cell with dz and rho at boundary
                pres_arr(i,j,klo-1) = p_0 + hz * rho_arr(i,j,klo) * l_gravity;
                pi_arr(i,j,klo-1) = getExnergivenP(pres_arr(i,j,klo-1), rdOcp);

            } else {
                // If level > 0 and klo > 0, we need to use the value of pres_arr(i,j,klo-1) which was
                //    filled from FillPatch-ing it.
                Real dz_loc;
                if (l_use_terrain) {
                    dz_loc = (zcc_arr(i,j,klo) - zcc_arr(i,j,klo-1));
                } else {
                    dz_loc = dz;
                }

                Real dens_interp = 0.5*(rho_arr(i,j,klo) + rho_arr(i,j,klo-1));
                pres_arr(i,j,klo) = pres_arr(i,j,klo-1) - dz_loc * dens_interp * l_gravity;

                pi_arr(i,j,klo  ) = getExnergivenP(pres_arr(i,j,klo  ), rdOcp);
                pi_arr(i,j,klo-1) = getExnergivenP(pres_arr(i,j,klo-1), rdOcp);
            }

            Real dens_interp;
            if (l_use_terrain) {
                for (int k = klo+1; k <= khi; k++) {
                    Real dz_loc = (zcc_arr(i,j,k) - zcc_arr(i,j,k-1));
                    dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                    pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz_loc * dens_interp * l_gravity;
                    pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
                }
            } else {
                for (int k = klo+1; k <= khi; k++) {
                    dens_interp = 0.5*(rho_arr(i,j,k) + rho_arr(i,j,k-1));
                    pres_arr(i,j,k) = pres_arr(i,j,k-1) - dz * dens_interp * l_gravity;
                    pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
                }
            }
        });

        int domlo_x = domain.smallEnd(0); int domhi_x = domain.bigEnd(0);
        int domlo_y = domain.smallEnd(1); int domhi_y = domain.bigEnd(1);

        if (pres[mfi].box().smallEnd(0) < domlo_x)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(0,domlo_x-1);
            bx.setBig(0,domlo_x-1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(domlo_x,j,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().bigEnd(0) > domhi_x)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(0,domhi_x+1);
            bx.setBig(0,domhi_x+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(domhi_x,j,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().smallEnd(1) < domlo_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domlo_y-1);
            bx.setBig(1,domlo_y-1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domlo_y,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }

        if (pres[mfi].box().bigEnd(1) > domhi_y)
        {
            Box bx = mfi.nodaltilebox(2);
            bx.setSmall(1,domhi_y+1);
            bx.setBig(1,domhi_y+1);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                pres_arr(i,j,k) = pres_arr(i,domhi_y,k);
                pi_arr(i,j,k) = getExnergivenP(pres_arr(i,j,k), rdOcp);
            });
        }
    }

    dens.FillBoundary(geom[lev].periodicity());
    pres.FillBoundary(geom[lev].periodicity());
}

void ERF::init_geo_wind_profile(const std::string input_file,
                                Vector<Real>& u_geos,
                                Gpu::DeviceVector<Real>& u_geos_d,
                                Vector<Real>& v_geos,
                                Gpu::DeviceVector<Real>& v_geos_d,
                                const Geometry& lgeom,
                                const Vector<Real>& zlev_stag)
{
    const int klo = 0;
    const int khi = lgeom.Domain().bigEnd()[AMREX_SPACEDIM-1];
    const amrex::Real dz = lgeom.CellSize()[AMREX_SPACEDIM-1];

    const bool grid_stretch = (zlev_stag.size() > 0);
    const Real zbot = (grid_stretch) ? zlev_stag[klo]   : lgeom.ProbLo(AMREX_SPACEDIM-1);
    const Real ztop = (grid_stretch) ? zlev_stag[khi+1] : lgeom.ProbHi(AMREX_SPACEDIM-1);

    amrex::Print() << "Reading geostrophic wind profile from " << input_file << std::endl;
    std::ifstream profile_reader(input_file);
    if(!profile_reader.is_open()) {
        amrex::Error("Error opening the abl_geo_wind_table\n");
    }

    // First, read the input data into temp vectors
    std::string line;
    Vector<Real> z_inp, Ug_inp, Vg_inp;
    Real z, Ug, Vg;
    amrex::Print() << "z  Ug  Vg" << std::endl;
    while(std::getline(profile_reader, line)) {
        std::istringstream iss(line);
        iss >> z >> Ug >> Vg;
        amrex::Print() << z << " " << Ug << " " << Vg << std::endl;
        z_inp.push_back(z);
        Ug_inp.push_back(Ug);
        Vg_inp.push_back(Vg);
        if (z >= ztop) break;
    }

    const int Ninp = z_inp.size();
    AMREX_ALWAYS_ASSERT(z_inp[0] <= zbot);
    AMREX_ALWAYS_ASSERT(z_inp[Ninp-1] >= ztop);

    // Now, interpolate vectors to the cell centers
    for (int k = 0; k <= khi; k++) {
        z = (grid_stretch) ? 0.5 * (zlev_stag[k] + zlev_stag[k+1])
                           : zbot + (k + 0.5) * dz;
        u_geos[k] = interpolate_1d(z_inp.dataPtr(), Ug_inp.dataPtr(), z, Ninp);
        v_geos[k] = interpolate_1d(z_inp.dataPtr(), Vg_inp.dataPtr(), z, Ninp);
    }

    // Copy from host version to device version
    Gpu::copy(Gpu::hostToDevice, u_geos.begin(), u_geos.end(), u_geos_d.begin());
    Gpu::copy(Gpu::hostToDevice, v_geos.begin(), v_geos.end(), v_geos_d.begin());

    profile_reader.close();
}
