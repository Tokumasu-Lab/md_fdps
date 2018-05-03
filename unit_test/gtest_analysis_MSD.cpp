//***************************************************************************************
//  This program is unit test of MSD analyzer.
//***************************************************************************************

#include <gtest/gtest.h>
#include <random>

#include <particle_simulator.hpp>

//--- external library for MD
#include <molecular_dynamics_ext.hpp>

//--- user defined headers
#include "unit.hpp"
#include "md_defs.hpp"

#include "enum_tgt_state.hpp"

#include "analysis_MSD.hpp"
#include "analysis_target.hpp"

#include "gtest_common.hpp"


namespace TEST_DEFS {
    const PS::S64 mt_seed  = 1234567;

    const PS::S64 n_A_atom = 80;
    const PS::S64 n_B_atom = 40;
    const PS::S64 n_atom   = n_A_atom + n_B_atom;

    const PS::S64 n_sample = 1000;

    const PS::F32vec vec_A{1.0, 1.0, 1.0};
    const PS::F32vec vec_B{2.0, 1.0, 2.0};

    const PS::F32 eps_abs = 1.e-4;
    const PS::F32 eps_rel = 1.e-3;

    const PS::S32 sample_interval = 100;
    const PS::F64 dt              = Unit::to_real_time(1.0*Unit::femto_second/Unit::norm_time);
    const PS::F64 time_interval   = static_cast<PS::F64>(sample_interval)*dt;

    const std::string msd_A_log_file{"test_bin/msd_A_log.tsv"};
    const std::string msd_A_raw_file{"test_bin/msd_A_raw.tsv"};
    const std::string msd_A_log_ref{ "unit_test/ref/msd_A_log.tsv"};
    const std::string msd_A_raw_ref{ "unit_test/ref/msd_A_raw.tsv"};

    const std::string msd_B_log_file{"test_bin/msd_B_log.tsv"};
    const std::string msd_B_raw_file{"test_bin/msd_B_raw.tsv"};
    const std::string msd_B_log_ref{ "unit_test/ref/msd_B_log.tsv"};
    const std::string msd_B_raw_ref{ "unit_test/ref/msd_B_raw.tsv"};

    const std::string msd_AB_log_file{"test_bin/msd_AB_log.tsv"};
    const std::string msd_AB_raw_file{"test_bin/msd_AB_raw.tsv"};
    const std::string msd_AB_log_ref{ "unit_test/ref/msd_AB_log.tsv"};
    const std::string msd_AB_raw_ref{ "unit_test/ref/msd_AB_raw.tsv"};

    const PS::F32vec box_size{1.0, 1.0, 1.0};  // no meaning for MSD analysis. used in initialize only.
};


class MatchTarget {
public:
    Analysis::TargetValue<TGT_STATE> state;

    MatchTarget() = default;
    MatchTarget(const TGT_STATE s){
        this->state = s;
    }
    ~MatchTarget() = default;

    template <class Tptcl>
    bool isMatch(const Tptcl &ptcl) const {
        return this->state.isMatch(ptcl.state);
    }
};

class AtomMock {
public:
    MD_DEFS::ID_type id{-1};
    PS::F32vec       pos{0.0, 0.0, 0.0};
    PS::F32vec       trj{0.0, 0.0, 0.0};
    TGT_STATE        state;

    MD_DEFS::ID_type getId()  const { return this->id;  }
    PS::F32vec       getPos() const { return this->pos; }
    PS::F32vec       getTrj() const { return this->trj; }

    void setId(    const MD_DEFS::ID_type  id     ) { this->id    = id;      }
    void setPos(   const PS::F32vec       &pos_new) { this->pos   = pos_new; }
    void set_state(const TGT_STATE         state  ) { this->state = state;   }

    void addTrj(const PS::F32vec &move   ) { this->trj += move; }
    void clearTrj()                        { this->trj  = 0.0;  }

    template <class Tptcl>
    void copyFromFP(const Tptcl &ptcl){
        this->id    = ptcl.getId();
        this->pos   = ptcl.getPos();
        this->trj   = ptcl.getTrj();
        this->state = ptcl.state;
    }
    template <class Tforce>
    void copyFromForce(const Tforce &force){}
    void clear(){}
};

class ProfileMock {
public:
    PS::F64 trj_time = 0.0;

    void    set_trj_time(const PS::F64 t)       { this->trj_time = t;    }
    PS::F64 get_trj_time()                const { return this->trj_time; }
};

template <class Tpsys, class Tstate>
void make_test_particle(      Tpsys           &psys,
                              std::mt19937_64 &rand,
                        const PS::S64         &n_ptcl,
                        const PS::F64vec      &vec_ptcl,
                        const Tstate          &state    ){

    assert(n_ptcl >= 0);

    const PS::S64 n_start = psys.getNumberOfParticleLocal();
    const PS::S64 n_end   = n_start + n_ptcl;
    psys.setNumberOfParticleLocal(n_end);

    std::uniform_real_distribution<PS::F64> dist_f64(0.0, 1.0);
    std::uniform_int_distribution<PS::S32>  dist_int(-1, 1);

    for(PS::S64 i=n_start; i<n_end; ++i){
        const PS::S64 id_ptcl = i;
        psys[id_ptcl].setId(id_ptcl);
        psys[id_ptcl].setPos( PS::F64vec{dist_f64(rand),
                                         dist_f64(rand),
                                         dist_f64(rand) } );
        psys[id_ptcl].clearTrj();
        psys[id_ptcl].addTrj( PS::F64vec{dist_int(rand)*vec_ptcl.x,
                                         dist_int(rand)*vec_ptcl.y,
                                         dist_int(rand)*vec_ptcl.z } );
        psys[id_ptcl].set_state(state);
    }
}

template <class Tdinfo, class Tpsys>
void init_system(      Tdinfo  &dinfo,
                       Tpsys   &psys,
                 const PS::S64  n_total){

    dinfo.initialize(0.3);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain( PS::F32vec{0.0, 0.0, 0.0},
                            PS::F32vec{1.0, 1.0, 1.0} );

    psys.initialize();
    psys.setNumberOfParticleLocal(0);

    Normalize::setBoxSize(TEST_DEFS::box_size);
}

struct MSD_Data {
    PS::F32 t;
    PS::F32 r2_3d;
    PS::F32 r2_x;
    PS::F32 r2_y;
    PS::F32 r2_z;

    bool read_line(std::string line){
        STR_TOOL::removeCR(line);
        const auto str_list = STR_TOOL::split(line, "\t");

        if(str_list.size() < 5)                 return false;
        if( !STR_TOOL::isNumeric(str_list[0]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[1]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[2]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[3]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[4]) ) return false;

        this->t     = std::stof(str_list[0]);
        this->r2_3d = std::stof(str_list[1]);
        this->r2_x  = std::stof(str_list[2]);
        this->r2_y  = std::stof(str_list[3]);
        this->r2_z  = std::stof(str_list[4]);

        return true;
    }
};

void msd_result_compare(const std::string &file_log,
                        const std::string &file_ref){

    const PS::S32 rank = 0;
    if(PS::Comm::getRank() != rank) return;

    std::vector<MSD_Data> log_data;
    std::vector<MSD_Data> ref_data;

    FS_TOOL::file_load(file_log, log_data, rank);
    FS_TOOL::file_load(file_ref, ref_data, rank);

    ASSERT_EQ(log_data.size(), ref_data.size());
    for(size_t i=0; i<log_data.size(); ++i){
        EXPECT_TRUE( float_relative_eq(log_data[i].t,
                                       ref_data[i].t,
                                       TEST_DEFS::eps_abs,
                                       TEST_DEFS::eps_rel) ) << " i = " << i;
        EXPECT_TRUE( float_relative_eq(log_data[i].r2_3d,
                                       ref_data[i].r2_3d,
                                       TEST_DEFS::eps_abs,
                                       TEST_DEFS::eps_rel) ) << " i = " << i;
        EXPECT_TRUE( float_relative_eq(log_data[i].r2_x,
                                       ref_data[i].r2_x,
                                       TEST_DEFS::eps_abs,
                                       TEST_DEFS::eps_rel) ) << " i = " << i;
        EXPECT_TRUE( float_relative_eq(log_data[i].r2_y,
                                       ref_data[i].r2_y,
                                       TEST_DEFS::eps_abs,
                                       TEST_DEFS::eps_rel) ) << " i = " << i;
        EXPECT_TRUE( float_relative_eq(log_data[i].r2_z,
                                       ref_data[i].r2_z,
                                       TEST_DEFS::eps_abs,
                                       TEST_DEFS::eps_rel) ) << " i = " << i;
    }
}

TEST(MSD, TracerFuncWhichProc){
    Analysis::MSD_Tracer<MD_DEFS::ID_type> msd_tracer;
    msd_tracer.init(TEST_DEFS::n_atom);

    for(PS::S64 i=0; i<TEST_DEFS::n_atom; ++i){
        EXPECT_EQ(msd_tracer.which_proc(i), msd_tracer.which_proc_naive(i)) << " i = " << i;
    }

    EXPECT_THROW(msd_tracer.which_proc_naive(-1)               , std::out_of_range);
    EXPECT_THROW(msd_tracer.which_proc_naive(TEST_DEFS::n_atom), std::out_of_range);

    EXPECT_THROW(msd_tracer.which_proc(-1)               , std::out_of_range);
    EXPECT_THROW(msd_tracer.which_proc(TEST_DEFS::n_atom), std::out_of_range);
}

TEST(MSD, TracerFuncAddTrace){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, TEST_DEFS::n_atom);

    Analysis::MSD_Tracer<MD_DEFS::ID_type, PS::F32vec> msd_tracer;
    msd_tracer.init(TEST_DEFS::n_atom);

    for(PS::S32 i=0; i<10; ++i){
        atom.setNumberOfParticleLocal(0);
        if(PS::Comm::getRank() == 0){
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_A_atom,
                               TEST_DEFS::vec_A*TEST_DEFS::sample_interval,
                               TGT_STATE::A        );
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_B_atom,
                               TEST_DEFS::vec_B*TEST_DEFS::sample_interval,
                               TGT_STATE::B        );
            atom.adjustPositionIntoRootDomain(dinfo);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        ProfileMock time_prof;
        time_prof.set_trj_time(TEST_DEFS::time_interval);  // regular interval

        msd_tracer.add_trace(time_prof, atom);
    }

    atom.setNumberOfParticleLocal(0);
    if(PS::Comm::getRank() == 0){
        make_test_particle(atom,
                           mt,
                           TEST_DEFS::n_A_atom,
                           TEST_DEFS::vec_A*TEST_DEFS::sample_interval,
                           TGT_STATE::A        );
        make_test_particle(atom,
                           mt,
                           TEST_DEFS::n_B_atom,
                           TEST_DEFS::vec_B*TEST_DEFS::sample_interval,
                           TGT_STATE::B        );
        atom.adjustPositionIntoRootDomain(dinfo);

        //--- missing 1 particle at last
        const PS::U64 n_total_false = atom.getNumberOfParticleLocal()-1;
        atom.setNumberOfParticleLocal(n_total_false);
    }


    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    ProfileMock time_prof;
    time_prof.set_trj_time(TEST_DEFS::time_interval);

    if(PS::Comm::getRank() == PS::Comm::getNumberOfProc()-1){
        EXPECT_THROW(msd_tracer.add_trace(time_prof, atom),
                     std::invalid_argument                 );
    } else {
        EXPECT_NO_THROW(msd_tracer.add_trace(time_prof, atom));
    }
}

TEST(MSD, TracerFuncIntervalChecker){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, TEST_DEFS::n_atom);

    Analysis::MSD_Tracer<MD_DEFS::ID_type, PS::F32vec> msd_tracer;
    msd_tracer.init(TEST_DEFS::n_atom);

    for(PS::S32 i=0; i<10; ++i){
        atom.setNumberOfParticleLocal(0);
        if(PS::Comm::getRank() == 0){
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_A_atom,
                               TEST_DEFS::vec_A*TEST_DEFS::sample_interval,
                               TGT_STATE::A                                );
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_B_atom,
                               TEST_DEFS::vec_B*TEST_DEFS::sample_interval,
                               TGT_STATE::B                                );
            atom.adjustPositionIntoRootDomain(dinfo);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        ProfileMock time_prof;
        time_prof.set_trj_time(TEST_DEFS::time_interval);  // regular interval

        msd_tracer.add_trace(time_prof, atom);
    }

    EXPECT_TRUE(msd_tracer.isRegularInterval());

    atom.setNumberOfParticleLocal(0);
    if(PS::Comm::getRank() == 0){
        make_test_particle(atom,
                           mt,
                           TEST_DEFS::n_A_atom,
                           TEST_DEFS::vec_A*TEST_DEFS::sample_interval,
                           TGT_STATE::A                                );
        make_test_particle(atom,
                           mt,
                           TEST_DEFS::n_B_atom,
                           TEST_DEFS::vec_B*TEST_DEFS::sample_interval,
                           TGT_STATE::B                                );
        atom.adjustPositionIntoRootDomain(dinfo);
    }

    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    ProfileMock time_prof;
    time_prof.set_trj_time(TEST_DEFS::time_interval*1.01);  // irregular interval

    msd_tracer.add_trace(time_prof, atom);

    EXPECT_FALSE(msd_tracer.isRegularInterval());
}

TEST(MSD, ptclA){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, TEST_DEFS::n_atom);

    Analysis::MSD_Tracer<MD_DEFS::ID_type, PS::F32vec> msd_tracer;
    msd_tracer.init(TEST_DEFS::n_atom);

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "  generate sample..." << "\n";
        std::cout << oss.str() << std::flush;
    }

    for(PS::S32 i_sample=0; i_sample<TEST_DEFS::n_sample; ++i_sample){
        atom.setNumberOfParticleLocal(0);
        if(PS::Comm::getRank() == 0){
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_A_atom,
                               TEST_DEFS::vec_A*TEST_DEFS::sample_interval,
                               TGT_STATE::A        );
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_B_atom,
                               TEST_DEFS::vec_B*TEST_DEFS::sample_interval,
                               TGT_STATE::B        );
            atom.adjustPositionIntoRootDomain(dinfo);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        ProfileMock time_prof;
        time_prof.set_trj_time(TEST_DEFS::time_interval);

        msd_tracer.add_trace(time_prof, atom);
    }

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "     ...done." << "\n";
        std::cout << oss.str() << std::flush;
    }

    //--- set up sampler
    std::vector<MatchTarget> tgt_list = { MatchTarget{TGT_STATE::A} };

    Analysis::MSD_Sampler<MatchTarget, MD_DEFS::ID_type> msd_sampler;
    msd_sampler.input_setting(TEST_DEFS::msd_A_log_file,
                              tgt_list,
                              Analysis::MSD_SAMPLING_MODE::resampling);
    msd_sampler.broadcast_setting();
    msd_sampler.make_id_table(atom, msd_tracer);
    msd_sampler.output(msd_tracer);

    msd_result_compare(TEST_DEFS::msd_A_log_file,
                       TEST_DEFS::msd_A_log_ref  );

    Analysis::MSD_Sampler<MatchTarget, MD_DEFS::ID_type> msd_sampler_raw;
    msd_sampler_raw.input_setting(TEST_DEFS::msd_A_raw_file,
                                  tgt_list,
                                  Analysis::MSD_SAMPLING_MODE::raw);
    msd_sampler_raw.broadcast_setting();
    msd_sampler_raw.make_id_table(atom, msd_tracer);
    msd_sampler_raw.output(msd_tracer);

    msd_result_compare(TEST_DEFS::msd_A_raw_file,
                       TEST_DEFS::msd_A_raw_ref  );
}

TEST(MSD, ptclB){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, TEST_DEFS::n_atom);

    Analysis::MSD_Tracer<MD_DEFS::ID_type, PS::F32vec> msd_tracer;
    msd_tracer.init(TEST_DEFS::n_atom);

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "  generate sample..." << "\n";
        std::cout << oss.str() << std::flush;
    }

    for(PS::S32 i_sample=0; i_sample<TEST_DEFS::n_sample; ++i_sample){
        atom.setNumberOfParticleLocal(0);
        if(PS::Comm::getRank() == 0){
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_A_atom,
                               TEST_DEFS::vec_A*TEST_DEFS::sample_interval,
                               TGT_STATE::A                                );
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_B_atom,
                               TEST_DEFS::vec_B*TEST_DEFS::sample_interval,
                               TGT_STATE::B                                );
            atom.adjustPositionIntoRootDomain(dinfo);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        ProfileMock time_prof;
        time_prof.set_trj_time(TEST_DEFS::time_interval);

        msd_tracer.add_trace(time_prof, atom);
    }

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "     ...done." << "\n";
        std::cout << oss.str() << std::flush;
    }

    //--- set up sampler
    std::vector<MatchTarget> tgt_list = { MatchTarget{TGT_STATE::B} };

    Analysis::MSD_Sampler<MatchTarget, MD_DEFS::ID_type> msd_sampler;
    msd_sampler.input_setting(TEST_DEFS::msd_B_log_file,
                              tgt_list,
                              Analysis::MSD_SAMPLING_MODE::resampling);
    msd_sampler.broadcast_setting();
    msd_sampler.make_id_table(atom, msd_tracer);
    msd_sampler.output(msd_tracer);

    msd_result_compare(TEST_DEFS::msd_B_log_file,
                       TEST_DEFS::msd_B_log_ref  );

    Analysis::MSD_Sampler<MatchTarget, MD_DEFS::ID_type> msd_sampler_raw;
    msd_sampler_raw.input_setting(TEST_DEFS::msd_B_raw_file,
                                  tgt_list,
                                  Analysis::MSD_SAMPLING_MODE::raw);
    msd_sampler_raw.broadcast_setting();
    msd_sampler_raw.make_id_table(atom, msd_tracer);
    msd_sampler_raw.output(msd_tracer);

    msd_result_compare(TEST_DEFS::msd_B_raw_file,
                       TEST_DEFS::msd_B_raw_ref  );
}

TEST(MSD, ptclAB){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, TEST_DEFS::n_atom);

    Analysis::MSD_Tracer<MD_DEFS::ID_type, PS::F32vec> msd_tracer;
    msd_tracer.init(TEST_DEFS::n_atom);

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "  generate sample..." << "\n";
        std::cout << oss.str() << std::flush;
    }

    for(PS::S32 i_sample=0; i_sample<TEST_DEFS::n_sample; ++i_sample){
        atom.setNumberOfParticleLocal(0);
        if(PS::Comm::getRank() == 0){
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_A_atom,
                               TEST_DEFS::vec_A*TEST_DEFS::sample_interval,
                               TGT_STATE::A                                );
            make_test_particle(atom,
                               mt,
                               TEST_DEFS::n_B_atom,
                               TEST_DEFS::vec_B*TEST_DEFS::sample_interval,
                               TGT_STATE::B                                );
            atom.adjustPositionIntoRootDomain(dinfo);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        ProfileMock time_prof;
        time_prof.set_trj_time(TEST_DEFS::time_interval);

        msd_tracer.add_trace(time_prof, atom);
    }

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "     ...done." << "\n";
        std::cout << oss.str() << std::flush;
    }

    //--- set up sampler
    std::vector<MatchTarget> tgt_list = { MatchTarget{TGT_STATE::A},
                                          MatchTarget{TGT_STATE::B} };

    Analysis::MSD_Sampler<MatchTarget, MD_DEFS::ID_type> msd_sampler;
    msd_sampler.input_setting(TEST_DEFS::msd_AB_log_file,
                              tgt_list,
                              Analysis::MSD_SAMPLING_MODE::resampling);
    msd_sampler.broadcast_setting();
    msd_sampler.make_id_table(atom, msd_tracer);
    msd_sampler.output(msd_tracer);

    msd_result_compare(TEST_DEFS::msd_AB_log_file,
                       TEST_DEFS::msd_AB_log_ref  );

    Analysis::MSD_Sampler<MatchTarget, MD_DEFS::ID_type> msd_sampler_raw;
    msd_sampler_raw.input_setting(TEST_DEFS::msd_AB_raw_file,
                                  tgt_list,
                                  Analysis::MSD_SAMPLING_MODE::raw);
    msd_sampler_raw.broadcast_setting();
    msd_sampler_raw.make_id_table(atom, msd_tracer);
    msd_sampler_raw.output(msd_tracer);

    msd_result_compare(TEST_DEFS::msd_AB_raw_file,
                       TEST_DEFS::msd_AB_raw_ref  );
}

TEST(MSD, exception){
    std::vector<MatchTarget>                    tgt_list;
    Analysis::MSD_Sampler<MatchTarget, PS::S64> msd_sampler;

    EXPECT_THROW(msd_sampler.input_setting("",
                                           tgt_list,
                                           Analysis::MSD_SAMPLING_MODE::raw,
                                           -1),
                 std::invalid_argument                                      ) << " all" ;

    EXPECT_THROW(msd_sampler.input_setting("a",
                                           tgt_list,
                                           Analysis::MSD_SAMPLING_MODE::raw,
                                           -1),
                 std::invalid_argument                                      ) << " tgt, rank" ;

    tgt_list = { MatchTarget{TGT_STATE::A} };
    EXPECT_THROW(msd_sampler.input_setting("a",
                                           tgt_list,
                                           Analysis::MSD_SAMPLING_MODE::raw,
                                           -1),
                 std::out_of_range                                          ) << " rank" ;

    EXPECT_NO_THROW(msd_sampler.input_setting("a",
                                              tgt_list,
                                              Analysis::MSD_SAMPLING_MODE::raw,
                                              0)                               );

    EXPECT_THROW(msd_sampler.input_setting("a",
                                           tgt_list,
                                           Analysis::MSD_SAMPLING_MODE::raw,
                                           0),
                 std::logic_error                                           ) << " already set";
}


#include "gtest_main_mpi.hpp"
