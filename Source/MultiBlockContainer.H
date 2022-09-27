#ifndef MultiBlockContainer_H_
#define MultiBlockContainer_H_

#include <AMReX_NonLocalBC.H>
#include <ERF.H>
#include <amr-wind/incflo.H>

class MultiBlockContainer
{
public:

    // Constructor with vector of parameters
    MultiBlockContainer (const std::vector<amrex::RealBox>& rb_v,
                         std::vector<int> max_level_in_v,
                         const std::vector<amrex::Vector<int>>& n_cell_in_v,
                         std::vector<int> coord_v,
                         const std::vector<amrex::Vector<amrex::IntVect>>& ref_ratios_v,
                         const std::vector<amrex::Array<int,AMREX_SPACEDIM>>& is_per_v,
                         std::vector<std::string> prefix_v,
                         int max_step);

    // Destructor
    ~MultiBlockContainer ();

    // Initialize block data
    void InitializeBlocks ();

    // Set up BoxList vector for use with Communication Meta Data
    void SetBoxLists ();

    void SetBoxListsAtoE();

    // Set up MB Communication Meta Data
    void SetBlockCommMetaData ();

    // Advance blocks
    void AdvanceBlocks ();

    // Fill MFs from other block(s)
    void FillPatchBlocks (int src_ind, int dst_ind);

    // Fill MFs from other block(s)
    void FillPatchBlocksFull ();

    // Fill MFs from other block(s)
    void FillPatchBlocksAE ();

    // Fill AMR-Wind Boundary Regsiter from ERF1
  void CopyToBoundaryRegister (amrex::BndryRegister& receive_br_old, amrex::BndryRegister& receive_br_new, amrex::Orientation ori);

  //void PopulateErfTimesteps (amrex::Vector<amrex::Real> tsteps) const {erf1.CoarseTimeState(/*tsteps*/);}
  void PopulateErfTimesteps(amrex::Real* tsteps);

  std::vector<amrex::Box> blv_atoe;
private:
    int m_max_step;
  // TODO: why does amrwind have to come before erf?
    incflo amrwind;
    ERF erf1;
    ERF erf2;

    amrex::NonLocalBC::MultiBlockIndexMapping dtos;
    amrex::NonLocalBC::MultiBlockIndexMapping dtos_etoa;
    amrex::NonLocalBC::MultiBlockIndexMapping dtos_atoe;
    std::vector<std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*>> cmd;
    std::vector<amrex::BoxList> blv;
    std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*> cmd_full;
    std::vector<amrex::Box> blv_full;
    std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*> cmd_etoa;
  //std::vector<std::vector<amrex::NonLocalBC::MultiBlockCommMetaData*>> cmd_atoe;
    std::vector<amrex::Box> blv_etoa;
};

#endif
