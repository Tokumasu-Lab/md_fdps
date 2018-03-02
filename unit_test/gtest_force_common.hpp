//=======================================================================================
//  This is common routine for Force test.
//=======================================================================================

#include <cassert>

#include "gtest_common.hpp"


void test_init(const PS::S64 n_step){
    if(PS::Comm::getRank() != 0) return;

    //--- system parameters
    System::profile.coef_ema      = 0.3;
    System::profile.theta         = 0.5;
    System::profile.n_leaf_limit  = 8;
    System::profile.n_group_limit = 64;

    System::profile.cut_off_LJ    = 12.0;
    System::profile.cut_off_intra =  9.0;

    //--- set loop condition
    System::profile.dt       = 1.0;
    System::profile.istep    = 0;
    System::profile.nstep_st = 0;
    System::profile.nstep_ed = n_step;

    //--- set domain size
    Normalize::setBoxSize( PS::F32vec{ 40.0,
                                       40.0,
                                       40.0 } );
}

template <class Tptcl, class Tdinfo, class Tforce,
          class Tdata>
void execute_force_calc(Tptcl              &atom,
                        Tdinfo             &dinfo,
                        Tforce             &force,
                        std::vector<Tdata> &force_log,
                        std::vector<Tdata> &force_ref ){

    //--- sync settings
    System::broadcast_profile(0);
    MODEL::coef_table.broadcast(0);
    System::InitDinfo(dinfo);

    //--- split domain & particle
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    //--- initialize force calculator
    const PS::S64 n_total = atom.getNumberOfParticleGlobal();
    force.init(n_total);

    //--- calculate force
    PS::S32 record_count = 0;
    while( System::isLoopContinue() ){
        //--- exchange particle
        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        //--- calculate intermolecular force in FDPS
        force.update_intra_pair_list(atom, dinfo, MODEL::coef_table.mask_scaling);
        force.update_force(atom, dinfo);

        //--- recording
        test_record(atom, record_count, force_log);
        ++record_count;

        //--- move
        test_move(atom);
        atom.adjustPositionIntoRootDomain(dinfo);

        //--- nest step
        System::StepNext();
    }
}
