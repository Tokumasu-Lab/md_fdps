//***************************************************************************************
//  This is unit test of loading system conditions.
//***************************************************************************************

#include <iostream>

#include <particle_simulator.hpp>

#include "md_defs.hpp"
#include "ext_sys_control.hpp"
#include "md_loading_condition.hpp"

int main(int argc, char *argv[]) {
    PS::Initialize(argc, argv);

    PS::DomainInfo dinfo;

    //--- make ext_sys controller object
    EXT_SYS::Sequence   ext_sys_sequence;
    EXT_SYS::Controller ext_sys_controller;

    if(PS::Comm::getRank() == 0) {

        //--- display total threads for FDPS
        std::ostringstream oss;
        oss << "Number of processes          : " << PS::Comm::getNumberOfProc()   << "\n"
            << "Number of threads per process: " << PS::Comm::getNumberOfThread() << "\n";
        std::cout << oss.str() << std::flush;

        //--- initialize
        oss.str("");
        oss << "\n"
            << "--- settings loaded in rank: " << PS::Comm::getRank() << "\n"
            << "\n";
        std::cout << oss.str() << std::flush;

        System::loading_sequence_condition(MD_DEFS::condition_sequence_file,
                                           ext_sys_sequence,
                                           ext_sys_controller );
        System::loading_molecular_condition(MD_DEFS::condition_molecule_file);
    }


    System::broadcast_profile();
    ext_sys_sequence.broadcast(0);
    ext_sys_controller.broadcast(0);

    {
        std::ostringstream oss;
        oss << "--- sync in MPI broadcast ---" << "\n";
        std::cout << oss.str() << std::flush;
    }

    System::InitDinfo(dinfo);

    //--- display settings
    if(PS::Comm::getRank() == PS::Comm::getNumberOfProc()-1) {
        std::cout << "\n"
                  << "--- settings displayed in rank: " << PS::Comm::getRank() << "\n"
                  << std::endl;

        System::print_profile();
        System::print_initializer_setting();

        ext_sys_controller.print();
        ext_sys_sequence.print();
    }


    PS::Finalize();
    return 0;
}
