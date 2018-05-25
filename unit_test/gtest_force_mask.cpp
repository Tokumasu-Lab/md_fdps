//***************************************************************************************
//  This program is unit test of intramolecular mask for coulomb & LJ interaction.
//
//    notice: sometimes deadlock occurs at 4 or more MPI process.
//            the cause is taht number of particle is < that of MPI proc? (did not be confirmed)
//***************************************************************************************

#include <gtest/gtest.h>

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

//--- external library for MD
#include <molecular_dynamics_ext.hpp>

//--- user defined headers
//------ definition of data set
#include "unit.hpp"
#include "md_enum.hpp"
#include "md_defs.hpp"
#include "atom_class.hpp"
#include "md_coef_table.hpp"
//------ calculate interaction
#include "md_force.hpp"
//------ loading model parameter
#include "md_loading_model.hpp"

//--- common tool of unit test for force
#include "gtest_force_common.hpp"


namespace TEST_DEFS {
    const PS::S64 n_loop = 200;
    const PS::F64 range  = 10.0;

    const MD_DEFS::ID_type id_tgt = 1;

    const PS::F32 eps_abs = 1.e-4;
    const PS::F32 eps_rel = 1.e-4;

    const std::string mask_100_naive_log_file{"test_bin/force_mask_100_naive_log.dat"};
    const std::string mask_100_opt_log_file{  "test_bin/force_mask_100_opt_log.dat"};
    const std::string mask_100_ref_file{      "unit_test/ref/force_mask_100_ref.dat"};

    const std::string mask_075_naive_log_file{"test_bin/force_mask_075_naive_log.dat"};
    const std::string mask_075_opt_log_file{  "test_bin/force_mask_075_opt_log.dat"};
    const std::string mask_075_ref_file{      "unit_test/ref/force_mask_075_ref.dat"};

    const std::string mask_050_naive_log_file{"test_bin/force_mask_050_naive_log.dat"};
    const std::string mask_050_opt_log_file{  "test_bin/force_mask_050_opt_log.dat"};
    const std::string mask_050_ref_file{      "unit_test/ref/force_mask_050_ref.dat"};

    const std::string mask_025_naive_log_file{"test_bin/force_mask_025_naive_log.dat"};
    const std::string mask_025_opt_log_file{  "test_bin/force_mask_025_opt_log.dat"};
    const std::string mask_025_ref_file{      "unit_test/ref/force_mask_025_ref.dat"};

    const std::string mask_000_naive_log_file{"test_bin/force_mask_000_naive_log.dat"};
    const std::string mask_000_opt_log_file{  "test_bin/force_mask_000_opt_log.dat"};
    const std::string mask_000_ref_file{      "unit_test/ref/force_mask_000_ref.dat"};
};

template <class Tpsys>
void test_move(Tpsys &atom){

    for(PS::S64 i=0; i<atom.getNumberOfParticleLocal(); ++i){
        PS::F64vec pos_tmp = Normalize::realPos( atom[i].getPos() );

        if(atom[i].getAtomID() == 1){
            //--- move for 2-body potential test (LJ, coulomb, bond)
            pos_tmp.x += TEST_DEFS::range/PS::F64(TEST_DEFS::n_loop);
        }

        pos_tmp = Normalize::normPos(pos_tmp);
        atom[i].setPos(pos_tmp);
    }
}

class ForceData {
public:
    PS::S32    count;
    PS::F32vec pos;
    PS::F32    pot_LJ;
    PS::F32vec force_LJ;
    PS::F32    pot_coulomb;
    PS::F32vec field_coulomb;

    bool read_line(std::string line){
        STR_TOOL::removeCR(line);
        const auto str_list = STR_TOOL::split(line, " ");

        const size_t data_field_len = 12;
        if(str_list.size() < data_field_len)     return false;
        if( !STR_TOOL::isInteger(str_list[0]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[1]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[2]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[3]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[4]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[5]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[6]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[7]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[8]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[9]) )  return false;
        if( !STR_TOOL::isNumeric(str_list[10]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[11]) ) return false;

        this->count = std::stoi(str_list[0]);
        this->pos   = PS::F32vec{ std::stof(str_list[1]),
                                  std::stof(str_list[2]),
                                  std::stof(str_list[3]) };
        this->pot_LJ   = std::stof(str_list[4]);
        this->force_LJ = PS::F32vec{ std::stof(str_list[5]),
                                     std::stof(str_list[6]),
                                     std::stof(str_list[7]) };
        this->pot_coulomb   = std::stof(str_list[8]);
        this->field_coulomb = PS::F32vec{ std::stof(str_list[9]),
                                          std::stof(str_list[10]),
                                          std::stof(str_list[11]) };

        return true;
    }
};
std::ostream& operator << (std::ostream &s, const ForceData &d){
    s << std::setw(10) << d.count << " "
      << std::scientific
      << std::setw(15) << std::setprecision(8) << d.pos.x << " "
      << std::setw(15) << std::setprecision(8) << d.pos.y << " "
      << std::setw(15) << std::setprecision(8) << d.pos.z << " "
      << std::setw(15) << std::setprecision(8) << d.pot_LJ << " "
      << std::setw(15) << std::setprecision(8) << d.force_LJ.x << " "
      << std::setw(15) << std::setprecision(8) << d.force_LJ.y << " "
      << std::setw(15) << std::setprecision(8) << d.force_LJ.z << " "
      << std::setw(15) << std::setprecision(8) << d.pot_coulomb << " "
      << std::setw(15) << std::setprecision(8) << d.field_coulomb.x << " "
      << std::setw(15) << std::setprecision(8) << d.field_coulomb.y << " "
      << std::setw(15) << std::setprecision(8) << d.field_coulomb.z << "\n";
    return s;
}

void check_result(const std::vector<ForceData> &result,
                  const std::vector<ForceData> &ref    ){

    EXPECT_EQ(result.size(), ref.size());
    const PS::S32 n = std::min(result.size(), ref.size());

    if(PS::Comm::getRank() == 0){
        for(PS::S32 i=0; i<n; ++i){
            EXPECT_EQ(result[i].count, ref[i].count) << " i= " << i;
            EXPECT_FLOAT_EQ(result[i].pos.x    , ref[i].pos.x    ) << " i= " << i;
            EXPECT_FLOAT_EQ(result[i].pos.y    , ref[i].pos.y    ) << " i= " << i;
            EXPECT_FLOAT_EQ(result[i].pos.z    , ref[i].pos.z    ) << " i= " << i;
            EXPECT_TRUE( float_relative_eq(result[i].pot_LJ    , ref[i].pot_LJ    , TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
            EXPECT_TRUE( float_relative_eq(result[i].force_LJ.x, ref[i].force_LJ.x, TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
            EXPECT_TRUE( float_relative_eq(result[i].force_LJ.y, ref[i].force_LJ.y, TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
            EXPECT_TRUE( float_relative_eq(result[i].force_LJ.z, ref[i].force_LJ.z, TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
            EXPECT_TRUE( float_relative_eq(result[i].pot_coulomb    , ref[i].pot_coulomb    , TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
            EXPECT_TRUE( float_relative_eq(result[i].field_coulomb.x, ref[i].field_coulomb.x, TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
            EXPECT_TRUE( float_relative_eq(result[i].field_coulomb.y, ref[i].field_coulomb.y, TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
            EXPECT_TRUE( float_relative_eq(result[i].field_coulomb.z, ref[i].field_coulomb.z, TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i= " << i;
        }
    }
    COMM_TOOL::barrier();
}

template <class Tptcl, class Tdata>
void test_record(const Tptcl              &atom,
                 const PS::S32             count,
                       std::vector<Tdata> &logger){

    Atom_FP buf;
    PS::S32 data_proc = -1;

    for(PS::S32 i=0; i<atom.getNumberOfParticleLocal(); ++i){
        if(atom[i].getAtomID() == TEST_DEFS::id_tgt){
            buf       = atom[i];
            data_proc = PS::Comm::getRank();
        }
    }
    data_proc = PS::Comm::getMaxValue(data_proc);
    COMM_TOOL::broadcast(buf, data_proc);

    logger.push_back( ForceData{count,
                                Normalize::realPos( buf.getPos() - PS::F32vec{0.5, 0.5, 0.5}),
                                buf.getPotLJ(),
                                buf.getForceLJ(),
                                buf.getCharge()*buf.getPotCoulomb(),
                                buf.getCharge()*buf.getFieldCoulomb()} );
}

template <class Tptcl, class Tdinfo, class Tforce,
          class Tdata>
void execute_force_calc_naive(Tptcl              &atom,
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
        force.update_force_naive(atom, dinfo);

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

template <class Tdinfo, class Tpsys>
void setup_test_data(Tdinfo &dinfo,
                     Tpsys  &atom){

    test_init(TEST_DEFS::n_loop);

    atom.initialize();
    atom.setNumberOfParticleLocal(0);

    //--- atom settings
    if(PS::Comm::getRank() == 0){
        const PS::S64 n_atom = 2;
        atom.setNumberOfParticleLocal(n_atom);
        for(PS::S32 i=0; i<n_atom; ++i){
            atom[i].setAtomID(i);
            atom[i].setMolID(i);
            atom[i].setAtomType( AtomName::Ow );
            atom[i].setMolType( MolName::AA_wat_SPC_Fw );
            atom[i].setCharge( -1.0*Unit::coef_coulomb );
            atom[i].setVDW_R( 0.5*3.165492 );
            atom[i].setVDW_D( std::sqrt(0.1554253) );
            atom[i].clear();
        }

        atom[0].bond.add(1);
        atom[1].bond.add(0);

        atom[0].setPos( PS::F64vec(0.5) );
        atom[1].setPos( PS::F64vec(0.5) + Normalize::normPos( PS::F64vec(0.5, 0.0, 0.0) ) );
    }
}


TEST(TestForceMask, Naive100){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "1.0"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "1.0"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc_naive(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_100_naive_log_file, force_log);
    load_log_file( TEST_DEFS::mask_100_ref_file      , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Naive075){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "0.75"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "0.75"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc_naive(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_075_naive_log_file, force_log);
    load_log_file( TEST_DEFS::mask_075_ref_file      , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Naive050){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "0.5"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "0.5"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc_naive(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_050_naive_log_file, force_log);
    load_log_file( TEST_DEFS::mask_050_ref_file      , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Naive025){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "0.25"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "0.25"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc_naive(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_025_naive_log_file, force_log);
    load_log_file( TEST_DEFS::mask_025_ref_file      , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Naive000){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "0.0"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "0.0"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc_naive(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_000_naive_log_file, force_log);
    load_log_file( TEST_DEFS::mask_000_ref_file      , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Opt100){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "1.0"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "1.0"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_100_opt_log_file, force_log);
    load_log_file( TEST_DEFS::mask_100_ref_file    , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Opt075){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "0.75"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "0.75"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_075_opt_log_file, force_log);
    load_log_file( TEST_DEFS::mask_075_ref_file    , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Opt050){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "0.5"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "0.5"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_050_opt_log_file, force_log);
    load_log_file( TEST_DEFS::mask_050_ref_file    , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Opt025){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "0.25"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "0.25"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_025_opt_log_file, force_log);
    load_log_file( TEST_DEFS::mask_025_ref_file    , force_ref);

    check_result(force_log, force_ref);
}

TEST(TestForceMask, Opt000){
    //--- ParticleSystem object
    PS::DomainInfo              dinfo;
    PS::ParticleSystem<Atom_FP> atom;
    CalcForce                   force;

    //--- data logger
    std::vector<ForceData> force_log;
    std::vector<ForceData> force_ref;

    //--- initialize
    setup_test_data(dinfo, atom);
    force_log.clear();
    force_ref.clear();

    if(PS::Comm::getRank() == 0){
        //--- model parameters
        MODEL::coef_table.clear();
        std::vector<std::string> param_line_bond = {"Ow", "Ow", "none"};
        MODEL::loading_param_bond(ENUM::what(MolName::AA_wat_SPC_Fw),
                                  param_line_bond,
                                  MODEL::coef_table.bond);
        std::vector<std::string> mask_line_LJ      = {"scaling_LJ"     , "0.0"};
        std::vector<std::string> mask_line_coulomb = {"scaling_coulomb", "0.0"};
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_LJ,
                                     MODEL::coef_table.mask_scaling);
        MODEL::loading_param_scaling(ENUM::what(MolName::AA_wat_SPC_Fw),
                                     mask_line_coulomb,
                                     MODEL::coef_table.mask_scaling);
    }

    execute_force_calc(atom, dinfo, force, force_log, force_ref);

    write_log_file(TEST_DEFS::mask_000_opt_log_file, force_log);
    load_log_file( TEST_DEFS::mask_000_ref_file    , force_ref);

    check_result(force_log, force_ref);
}


#include "gtest_main_mpi.hpp"
