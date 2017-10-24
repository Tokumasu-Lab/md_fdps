//***************************************************************************************
//  This is unit test of loading model parameter function.
//***************************************************************************************

#include <fstream>
#include <iostream>
#include <vector>
#include <random>

#include <particle_simulator.hpp>

#include "md_intra_pair_table.hpp"


int main(int argc, char *argv[]) {
    PS::Initialize(argc, argv);

    //--- display total threads for FDPS
    if(PS::Comm::getRank() == 0){
        fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }


    std::cout << "\n"
              << " --- typical structure test ---" << std::endl;
    //--- define test connection structure
    //
    //    6
    //    |
    //  0-1-2-3-4-5
    //      |
    //      7-8
    //

    PS::S64 n_atom = 9;
    MODEL::intra_pair_manager.setAtomNumber(n_atom);

    MODEL::intra_pair_manager.addBond(0, 1);
    MODEL::intra_pair_manager.addBond(6, 1);
    MODEL::intra_pair_manager.addBond(2, 1);

    MODEL::intra_pair_manager.addBond(2, 7);
    MODEL::intra_pair_manager.addBond(2, 3);

    MODEL::intra_pair_manager.addBond(3, 4);
    MODEL::intra_pair_manager.addBond(4, 5);

    MODEL::intra_pair_manager.addBond(7, 8);


    //--- construct pair list
    MODEL::intra_pair_manager.makeIntraList();


    //--- show result
    for(PS::S64 i=0; i<n_atom; ++i){
        MODEL::print_connection(i);
    }

    std::cout << "\n\n"
              << " --- circulation structure test ---" << std::endl;
    //--- circulation structure test
    //
    //  0-3
    //  | |
    //  1-2-4-7
    //      | |
    //      5-6
    //

    n_atom = 8;
    MODEL::intra_pair_manager.clear();
    MODEL::intra_pair_manager.setAtomNumber(n_atom);

    MODEL::intra_pair_manager.addBond(0, 1);
    MODEL::intra_pair_manager.addBond(1, 2);
    MODEL::intra_pair_manager.addBond(2, 3);
    MODEL::intra_pair_manager.addBond(3, 0);

    MODEL::intra_pair_manager.addBond(2, 4);

    MODEL::intra_pair_manager.addBond(4, 7);
    MODEL::intra_pair_manager.addBond(7, 6);
    MODEL::intra_pair_manager.addBond(6, 5);
    MODEL::intra_pair_manager.addBond(5, 4);

    //--- construct pair list
    MODEL::intra_pair_manager.makeIntraList();


    //--- show result
    for(PS::S64 i=0; i<n_atom; ++i){
        MODEL::print_connection(i);
    }


    PS::Finalize();
    return 0;
}
