//***************************************************************************************
//  This is unit test of loading model parameter function.
//***************************************************************************************

#include <cmath>
#include <iostream>
#include <vector>

#include <particle_simulator.hpp>

#include "md_fdps_atom_class.hpp"
#include "md_fdps_coef_table.hpp"
#include "md_fdps_ext_sys_control.hpp"
#include "md_fdps_fileIO.hpp"
#include "md_fdps_initialize.hpp"


int main(int argc, char *argv[]) {
    PS::Initialize(argc, argv);

    //--- display total threads for FDPS
    if(PS::Comm::getRank() == 0){
        fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }

    //--- make particle system object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;

    atom.initialize();
    atom.setNumberOfParticleLocal(0);

    //--- make ext_sys controller object
    EXT_SYS::Sequence   ext_sys_sequence;
    EXT_SYS::Controller ext_sys_controller;

    //--- load settings
    if(PS::Comm::getRank() == 0){
        //--- display normalizing units
        Unit::print_unit();

        System::loading_sequence_condition("condition_sequence.imp",
                                           ext_sys_sequence,
                                           ext_sys_controller );
        System::loading_molecular_condition("condition_molecule.imp");

        for(size_t i=0; i<System::model_list.size(); ++i){
            MODEL::loading_model_parameter(ENUM::whatis(System::model_list.at(i).first),
                                           System::model_template.at(i),
                                           System::bond_template.at(i),
                                           MODEL::coefTable_elem,
                                           MODEL::coefTable_bond,
                                           MODEL::coefTable_angle,
                                           MODEL::coefTable_torsion);
        }
        System::print_setting();
        ext_sys_controller.print();
        ext_sys_sequence.print();
        PS::F64 init_temperature = ext_sys_sequence.getSetting(0).temperature;
        Initialize::InitParticle(atom, init_temperature);
    }

    //--- send settings to all MPI processes
    System::broadcast_setting(0);
    MODEL::broadcast_coefTable(0);
    MODEL::intra_pair_manager.broadcast(0);
    ext_sys_sequence.broadcast(0);
    ext_sys_controller.broadcast(0);

    //--- devide atom particle in MPI processes
    System::InitDinfo(dinfo);
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    std::cout << "proc = " << PS::Comm::getRank() << " / atoms = " << atom.getNumberOfParticleLocal() << std::endl;

    //--- output result
    FILE_IO::Init( System::setting );
    FILE_IO::pdb_start    = 0;
    FILE_IO::pdb_interval = 1;
    FILE_IO::recordVMD(atom, 0);
    FILE_IO::recordPos(atom, 0);
    FILE_IO::recordResume(atom, 0);

    PS::S64 n_atom_sample = atom.getNumberOfParticleLocal();
    n_atom_sample = std::min( n_atom_sample, PS::S64(10) );
    if(PS::Comm::getRank() == 1){
        std::cout << "\nnumber of sampling particle = " << n_atom_sample << std::endl;
        for(PS::S64 i=0; i<n_atom_sample; ++i){
            std::cout << "\n"
                      << " ID = " << atom[i].getAtomID() << "\n"
                      << MODEL::print_atom(atom[i]) << std::endl;
            MODEL::print_connection( atom[i].getAtomID() );
        }
    }

    PS::Finalize();
    return 0;
}
