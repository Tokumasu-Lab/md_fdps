//***************************************************************************************
//  This program is unit test of dihedral torsion interaction.
//***************************************************************************************

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

//--- external library for MD
#include <molecular_dynamics_ext.hpp>

//--- user defined headers
//------ definition of data set
#include "unit.hpp"
#include "md_enum.hpp"
#include "atom_class.hpp"
#include "md_coef_table.hpp"
//------ calculate interaction
#include "ff_intra_force.hpp"
#include "ff_inter_force.hpp"
//------ kick & drift
#include "md_atom_move.hpp"
//------ system observer
#include "md_observer.hpp"
//------ external system control
//#include "md_ext_sys_control.hpp"
//------ file I/O
#include "md_fileIO.hpp"
//------ initialize
#include "md_initialize.hpp"


const PS::S64 n_atom = 4;
const PS::S64 n_loop = 1000;

template <class Tpsys>
void test_init(Tpsys &psys){
    if(PS::Comm::getRank() != 0) return;

    //--- allocate arrays
    psys.setNumberOfParticleLocal(n_atom);
    MODEL::intra_pair_manager.setAtomNumber(n_atom);

    //--- set domain size
    Normalize::setBoxSize( PS::F32vec{ 40.0,
                                       40.0,
                                       40.0 } );

    //--- set initial position & parameters
    //------ for angle test
    for(PS::S32 i=0; i<n_atom; ++i){
        psys[i].setAtomID(i);
        psys[i].setMolID(i);
        psys[i].setCharge( 0.0 );
        psys[i].setVDW_R( 3.0 );
        psys[i].setVDW_D( 0.0 );
        psys[i].clear();
    }
    psys[0].setAtomType( AtomName::CH3 );
    psys[0].setMolType(  MolName::UA_propan_1_ol );
    psys[1].setAtomType( AtomName::CH2 );
    psys[1].setMolType(  MolName::UA_propan_1_ol );
    psys[2].setAtomType( AtomName::CH2 );
    psys[2].setMolType(  MolName::UA_propan_1_ol );
    psys[3].setAtomType( AtomName::O );
    psys[3].setMolType(  MolName::UA_propan_1_ol );

    MODEL::intra_pair_manager.addBond(0, 1);
    MODEL::intra_pair_manager.addBond(1, 2);
    MODEL::intra_pair_manager.addBond(2, 3);

    psys[0].setPos( PS::F64vec(0.5) + Normalize::normPos(PS::F64vec( 0.5, 1.0, 0.0)) );
    psys[1].setPos( PS::F64vec(0.5) + Normalize::normPos(PS::F64vec( 0.5, 0.0, 0.0)) );
    psys[2].setPos( PS::F64vec(0.5) + Normalize::normPos(PS::F64vec(-0.5, 0.0, 0.0)) );
    psys[3].setPos( PS::F64vec(0.5) + Normalize::normPos(PS::F64vec(-0.5, 1.0, 0.0)) );
    //psys[3].setPos( PS::F64vec(0.5) + VEC_EXT::rot_x( Normalize::normPos(PS::F64vec(-0.5, 1.0, 0.0)), -Unit::pi*(65.1/180.0) ) );

    //--- preview start point
    std::cout << "\n" << "initial position:" << "\n"
              << std::setw(8) << psys[0].getAtomType() << "  " << Normalize::realPos( psys[0].getPos() ) << "\n"
              << std::setw(8) << psys[1].getAtomType() << "  " << Normalize::realPos( psys[1].getPos() ) << "\n"
              << std::setw(8) << psys[2].getAtomType() << "  " << Normalize::realPos( psys[2].getPos() ) << "\n"
              << std::setw(8) << psys[3].getAtomType() << "  " << Normalize::realPos( psys[3].getPos() ) << std::endl;

    std::cout << "\n" << "torsion parameter:" << std::endl;
    MODEL::print_coef_table(MODEL::coefTable_torsion);

    //--- disable all bond and angle potential
    for(auto &c_bond : MODEL::coefTable_bond){
        c_bond.second.form = IntraFuncForm::none;
    }
    for(auto &c_angle : MODEL::coefTable_angle){
        c_angle.second.form = IntraFuncForm::none;
    }

    //--- set loop condition
    System::setting.nstep_st = 0;
    System::setting.nstep_ed = n_loop;
}

template <class Tpsys>
void test_move(Tpsys &psys){

    for(PS::S64 i=0; i<psys.getNumberOfParticleLocal(); ++i){
        PS::F64vec pos_tmp = Normalize::realPos( psys[i].getPos() );

        if(psys[i].getAtomID() == 3){
            //--- move for angle test
            PS::F64vec pos_local = pos_tmp - Normalize::realPos( PS::F64vec(0.5) );
                       pos_local = VEC_EXT::rot_x(pos_local, 2.0*Unit::pi/PS::F64(n_loop));
            pos_tmp = Normalize::realPos( PS::F64vec(0.5) ) + pos_local;
        }

        pos_tmp = Normalize::normPos(pos_tmp);
        pos_tmp = Normalize::periodicAdjustNorm(pos_tmp);
        psys[i].setPos(pos_tmp);
    }
}

static PS::S32 record_count = 0;

template <class Tpsys, class TPrinter>
void test_record(const Tpsys    &psys,
                       TPrinter &printer){

    //--- get target atom property
    Atom_FP buf;
    PS::S32 data_proc = -1;
    PS::S32 tgt_id    = 3;

    for(PS::S32 i=0; i<psys.getNumberOfParticleLocal(); ++i){
        if(psys[i].getAtomID() == tgt_id){
            buf = psys[i];
            data_proc = PS::Comm::getRank();
        }
    }
    data_proc = PS::Comm::getMaxValue(data_proc);
//    std::cout << "  data_proc = " << data_proc << std::endl;
    COMM_TOOL::broadcast(buf, data_proc);

    if(PS::Comm::getRank() != 0) return;

    //--- write file
    std::ostringstream oss;
    PS::F64vec pos_tmp = buf.getPos();
               pos_tmp = Normalize::realPos( pos_tmp - PS::F64vec(0.5) );

    oss << std::setw(10) << record_count;
    //--- for Intra force
    oss << std::setw(15) << std::setprecision(7) << PS::F64(record_count)*360.0/PS::F64(n_loop);  // 1 rotate
    oss << std::setw(15) << std::setprecision(7) << buf.getPotTorsion();
    oss << std::setw(15) << std::setprecision(7) << buf.getForceIntra().x;
    oss << std::setw(15) << std::setprecision(7) << buf.getForceIntra().y;
    oss << std::setw(15) << std::setprecision(7) << buf.getForceIntra().z;

    oss << std::setw(15) << std::setprecision(7) << buf.getPos().x;
    oss << std::setw(15) << std::setprecision(7) << buf.getPos().y;
    oss << std::setw(15) << std::setprecision(7) << buf.getPos().z;
    oss << "\n";

    printer.print(oss.str());
    ++record_count;
}

template <class TPrinter>
void test_header(TPrinter &printer){
    std::ostringstream oss;

    printer.print(oss.str());
}


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

    void setRcut(){
        Atom_EP::setRcut_LJ(      Normalize::normCutOff( System::get_cutoff_LJ() ) );
        Atom_EP::setRcut_coulomb( Normalize::normCutOff_PM() );
    }

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


int main(int argc, char* argv[]) {
        PS::Initialize(argc, argv);

        //--- display total threads for FDPS
        if(PS::Comm::getRank() == 0){
            fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
            fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
        }

        //--- make particle system object
        PS::DomainInfo              dinfo;
        PS::ParticleSystem<Atom_FP> atom;
        CalcForce                   force;

        //--- make ext_sys controller object
        EXT_SYS::Sequence   ext_sys_sequence;
        EXT_SYS::Controller ext_sys_controller;

        //--- data logger
        Observer::Energy      eng;
                              eng.file_init("test_force_dihedral_eng.dat");
        Observer::FilePrinter printer;
                              printer.file_init("test_force_dihedral_log.dat", 0);

        //--- initialize
        atom.initialize();
        atom.setNumberOfParticleLocal(0);

        //------ load settings
        if(PS::Comm::getRank() == 0){
            System::loading_sequence_condition("condition_sequence.imp",
                                               ext_sys_sequence,
                                               ext_sys_controller );
            System::loading_molecular_condition("condition_molecule.imp");

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
        //    System::print_setting();
        //    Initialize::InitParticle(atom);
            test_init(atom);
        }
        //------ sync settings
        System::broadcast_setting(0);
        MODEL::broadcast_coefTable(0);
        MODEL::intra_pair_manager.broadcast(0);

        System::InitDinfo(dinfo);
        FILE_IO::Init( System::setting );

        //--- split domain & particle
        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        //--- initialize force culculator
        PS::S64 n_local = atom.getNumberOfParticleLocal();
        PS::S64 n_total = atom.getNumberOfParticleGlobal();
        force.init(n_total);

        //--- main loop
        if(PS::Comm::getRank() == 0) std::cout << "test loop start!" << std::endl;
        while( System::isLoopContinue() ){

            //--- exchange particle
            dinfo.decomposeDomainAll(atom);
            atom.exchangeParticle(dinfo);

            //--- calculate intermolecular force in FDPS
            force.update(atom, dinfo);

            //--- recording
            eng.getEnergy(atom);
            eng.record( System::get_dt(), System::get_istep() );
            test_record(atom, printer);

            //--- move
            test_move(atom);

            //--- nest step
            System::StepNext();
        }
        if(PS::Comm::getRank() == 0) std::cout << "test loop ends!" << std::endl;

        //--- finalize FDPS
        PS::Finalize();
        return 0;
}
