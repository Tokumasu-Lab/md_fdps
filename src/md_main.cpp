/**************************************************************************************************/
/**
* @file  md_main.cpp
* @brief main routine of "md_fdps" code.
*          This code is using the Framework for Developing Particle Simulator (FDPS).
*          https://github.com/FDPS
*/
/**************************************************************************************************/

//--- standard library header
#include <chrono>

//--- FDPS header
#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

//--- external library for MD
#include <molecular_dynamics_ext.hpp>

//--- MD function headers
//------ definition of data set
#include "unit.hpp"
#include "md_enum.hpp"
#include "atom_class.hpp"
#include "md_coef_table.hpp"
#include "md_intra_pair_table.hpp"
//------ calculate interaction
#include "ff_intra_force.hpp"
#include "ff_inter_force.hpp"
//------ system observer
#include "md_observer.hpp"
//------ kick & drift
//#include "md_atom_move.hpp"
//------ external system control
#include "md_ext_sys_control.hpp"
//------ file I/O
#include "md_fileIO.hpp"
//------ initialize
#include "md_initialize.hpp"


/**
* @brief force calculater interface class.
*/
class CalcForce {
private:
    PS::TreeForForceShort<Force_FP, Atom_EP, Atom_EP>::Scatter tree_atom;
    PS::PM::ParticleMesh pm;

public:
    void init(const PS::S64 &n_total){
        tree_atom.initialize(n_total,
                             System::setting.theta,
                             System::setting.n_leaf_limit,
                             System::setting.n_group_limit);
    }

    /**
    * @brief update cutoff length in normalized space.
    */
    void setRcut(){
        Atom_EP::setRcut_LJ(      Normalize::normCutOff( System::get_cutoff_LJ() ) );
        Atom_EP::setRcut_coulomb( Normalize::normCutOff_PM() );

        //--- check cut off length
        assert(Atom_EP::getRcut_LJ() < 0.5);
    }

    /**
    * @brief update force on atom.
    */
    template <class Tpsys, class Tdinfo>
    void update(Tpsys  &atom,
                Tdinfo &dinfo){

        //--- clear force
        PS::S64 n_local = atom.getNumberOfParticleLocal();
        for(PS::S64 i=0; i<n_local; ++i){
            atom[i].clear();
        }

        this->setRcut();

        //=================
        //* PM part
        //=================
        pm.setDomainInfoParticleMesh(dinfo);
        pm.setParticleParticleMesh(atom, true);   // clear previous charge information
        pm.calcMeshForceOnly();

        //--- get potential and field
        for(PS::S64 i=0; i<n_local; ++i){
            PS::F64vec pos = atom[i].getPos();

            atom[i].addFieldCoulomb( Normalize::realPMForce(     -pm.getForce(pos)     ) );
            atom[i].addPotCoulomb(   Normalize::realPMPotential( -pm.getPotential(pos) ) );
        }

        //=================
        //* PP part
        //=================
        this->tree_atom.calcForceAll(calcForceShort<Force_FP, Atom_EP, Atom_EP>,
                                     atom, dinfo);
        for(PS::S64 i=0; i<n_local; ++i){
            Force_FP result = tree_atom.getForce(i);
            atom[i].addFieldCoulomb( result.getFieldCoulomb() );
            atom[i].addPotCoulomb(   result.getPotCoulomb() );
            atom[i].addForceLJ(  result.getForceLJ()  );
            atom[i].addPotLJ(    result.getPotLJ()    );
            atom[i].addVirialLJ( result.getVirialLJ() );
        }

        //=================
        //* Intra force part
        //=================
        calcForceIntra<Atom_EP, ForceIntra>(this->tree_atom, atom);
    }
};


/**
* @brief main routine
*/
int main(int argc, char* argv[]){
    PS::Initialize(argc, argv);

    //--- record real start time
    auto real_start_time = std::chrono::system_clock::now();

    //--- display total threads for FDPS
    if(PS::Comm::getRank() == 0){
        std::cerr << "Number of processes          : " << PS::Comm::getNumberOfProc()   << "\n"
                  << "Number of threads per process: " << PS::Comm::getNumberOfThread() << std::endl;
    }

    //--- make particle system object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- initialize
    atom.initialize();
    atom.setNumberOfParticleLocal(0);

    //--- make ext_sys controller object
    EXT_SYS::Sequence   ext_sys_sequence;
    EXT_SYS::Controller ext_sys_controller;

    //------ load settings
    if(PS::Comm::getRank() == 0){
        //--- load setting files
        System::loading_sequence_condition("condition_sequence.imp",
                                           ext_sys_sequence,
                                           ext_sys_controller );
        System::loading_molecular_condition("condition_molecule.imp");

        //--- load model parameters
        for(size_t i=0; i<System::model_list.size(); ++i){
            MODEL::loading_model_parameter(ENUM::what(System::model_list.at(i).first),
                                           System::model_template.at(i),
                                           System::bond_template.at(i),
                                           MODEL::coefTable_atom,
                                           MODEL::coefTable_residue,
                                           MODEL::coefTable_bond,
                                           MODEL::coefTable_angle,
                                           MODEL::coefTable_torsion);
        }
        //--- display settings
        Unit::print_unit();
        System::print_setting();
        ext_sys_controller.print();
        ext_sys_sequence.print();

        //--- insert molecule
        PS::F64 init_temperature = ext_sys_sequence.getSetting(0).temperature;
        Initialize::InitParticle(atom, init_temperature);
    }

    //--- send settings to all MPI processes
    System::broadcast_setting(0);
    MODEL::broadcast_coefTable(0);
    MODEL::intra_pair_manager.broadcast(0);
    ext_sys_sequence.broadcast(0);
    ext_sys_controller.broadcast(0);

    System::InitDinfo(dinfo);
    FILE_IO::Init( System::setting );

    //--- split domain & particle
    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    //--- initialize force culculator
    PS::S64 n_local = atom.getNumberOfParticleLocal();
    PS::S64 n_total = atom.getNumberOfParticleGlobal();
    force.init(n_total);

    //--- initialize observer
    Observer::Energy   eng;
    Observer::Property prop;

    Observer::MovingAve<Observer::Energy>   eng_ave;
    Observer::MovingAve<Observer::Property> prop_ave;

    eng.file_init("energy_raw.dat");
    prop.file_init("property_raw.dat");
    eng_ave.file_init( "eng_ave.dat",  System::get_eng_start(),  System::get_eng_interval() );
    prop_ave.file_init("prop_ave.dat", System::get_prop_start(), System::get_prop_interval() );


    force.update(atom, dinfo);


    //--- main loop
    if(PS::Comm::getRank() == 0) std::cout << "\n --- main loop start! ---\n" << std::endl;
    while( System::isLoopContinue() ){

        //--- get system property
        eng.getEnergy(atom);
        prop.getProperty(eng);

        //--- affect external system controller
        PS::S64 n_rigid = 0;
        ext_sys_controller.apply(n_rigid,
                                 ext_sys_sequence.getSetting( System::get_istep() ),
                                 System::get_dt(),
                                 atom,
                                 eng);

        //--- output record
        FILE_IO::recordPos(   atom, System::get_istep() );
        FILE_IO::recordVMD(   atom, System::get_istep() );
        FILE_IO::recordResume(atom, System::get_istep() );

        //--- energy log
        //------ raw data
        eng.record(  System::get_dt(), System::get_istep() );
        prop.record( System::get_dt(), System::get_istep() );
        //------ moving average
        eng_ave.record(  System::get_dt(), System::get_istep(), eng );
        prop_ave.record( System::get_dt(), System::get_istep(), prop );

        //--- kick
        //kick(0.5*System::get_dt(), atom);
        ext_sys_controller.kick(0.5*System::get_dt(), atom);

        //--- drift
        //drift(System::get_dt(), atom);
        ext_sys_controller.drift(System::get_dt(), atom);

        //--- exchange particle
        if( System::isDinfoUpdate() ){ dinfo.decomposeDomainAll(atom); }
        atom.exchangeParticle(dinfo);
        n_local = atom.getNumberOfParticleLocal();

        //--- calculate intermolecular force in FDPS
        force.update(atom, dinfo);

        //--- kick
        //kick(0.5*System::get_dt(), atom);
        ext_sys_controller.kick(0.5*System::get_dt(), atom);

        //--- nest step
        System::StepNext();
    }
    if(PS::Comm::getRank() == 0) std::cout << "\n --- main loop ends! ---\n" << std::endl;

    //--- show elapsed time
    auto real_end_time = std::chrono::system_clock::now();
    if(PS::Comm::getRank() == 0){
        std::cout << "  elapsed time = "
                  << chrono_str::to_str_h_m_s_ms(real_end_time - real_start_time)
                  << "\n" << std::endl;
    }

    //--- finalize FDPS
    PS::Finalize();
    return 0;
}
