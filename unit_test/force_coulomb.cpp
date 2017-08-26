//***************************************************************************************
//  This program is unit test of coulomb interaction.
//***************************************************************************************

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

//--- external library for MD
#include <molecular_dynamics_ext.hpp>

//--- user defined headers
//------ definition of data set
#include "md_fdps_unit.hpp"
#include "md_fdps_enum_model.hpp"
#include "md_fdps_atom_class.hpp"
#include "md_fdps_coef_table.hpp"
//------ calculate interaction
#include "ff_intra_force.hpp"
#include "ff_inter_force.hpp"
//------ kick & drift
#include "md_fdps_atom_move.hpp"
//------ system observer
#include "md_fdps_observer.hpp"
//------ external system control
//#include "md_fdps_ext_sys_control.hpp"
//------ file I/O
#include "md_fdps_fileIO.hpp"
//------ initialize
#include "md_fdps_initialize.hpp"


const PS::S64 n_atom = 2;
const PS::S64 n_loop = 1000;

template <class Tpsys>
void test_init(Tpsys &psys){
    if(PS::Comm::getRank() != 0) return;

    //--- allocate arrays
    psys.setNumberOfParticleLocal(n_atom);
    MODEL::intra_pair_manager.setAtomNumber(n_atom);

    //--- set initial position & parameters
    //------ for bond test
    for(PS::S32 i=0; i<n_atom; ++i){
        psys[i].setAtomID(i);
        psys[i].setAtomType( AtomName::Ar );
        psys[i].setMolType(  MolName::AA_Ar );
        psys[i].setMolID(i);
        psys[i].setCharge( 1.0*Unit::coef_coulomb );
        psys[i].setVDW_R( 3.0 );
        psys[i].setVDW_D( 0.0 );
        psys[i].clear();
    }

    psys[0].setPos( PS::F64vec(0.5) );
    psys[1].setPos( PS::F64vec(0.5) + PS::F64vec(0.001, 0.0, 0.0) );

    //--- set loop
    System::setting.nstep_st = 0;
    System::setting.nstep_ed = n_loop;
}

template <class Tpsys>
void test_move(Tpsys &psys){

    for(PS::S64 i=0; i<psys.getNumberOfParticleLocal(); ++i){
        PS::F64vec pos_tmp = psys[i].getPos();

        if(psys[i].getAtomID() == 1){
            //--- move for 2-body potential test (LJ, coulomb, bond)
            pos_tmp.x += 0.998/PS::F64(n_loop);         //  setting test distance range (normalized)
        }

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
    PS::S32 tgt_id    = 1;

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
    oss << std::setw(15) << std::setprecision(7) << pos_tmp.x;
    oss << std::setw(15) << std::setprecision(7) << pos_tmp.y;
    oss << std::setw(15) << std::setprecision(7) << pos_tmp.z;
    oss << std::setw(15) << std::setprecision(7) << buf.getPotCoulomb()*buf.getCharge();
    oss << std::setw(15) << std::setprecision(7) << ( buf.getFieldCoulomb()*buf.getCharge() ).x;
    oss << std::setw(15) << std::setprecision(7) << ( buf.getFieldCoulomb()*buf.getCharge() ).y;
    oss << std::setw(15) << std::setprecision(7) << ( buf.getFieldCoulomb()*buf.getCharge() ).z;
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
                              eng.file_init("test_force_coulomb_eng.dat");
        Observer::FilePrinter printer;
                              printer.file_init("test_force_coulomb_log.dat", 0);

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
                MODEL::loading_model_parameter(ENUM::whatis(System::model_list.at(i).first),
                                               System::model_template.at(i),
                                               System::bond_template.at(i),
                                               MODEL::coefTable_elem,
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
