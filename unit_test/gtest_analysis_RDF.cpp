//***************************************************************************************
//  This program is unit test of RDF analyzer.
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

#include "analysis_RDF.hpp"
#include "analysis_target.hpp"

#include "gtest_common.hpp"


namespace TEST_DEFS {
    const PS::S64 mt_seed  = 19937;

    const PS::S64 n_atom   = 3000;
    const PS::S64 n_sample = 100;

    const PS::F32    ex_r     = 1.0;
    const PS::F32vec box_size = {30.0, 30.0, 30.0};
    const PS::S32    n_stripe = 30;

    const PS::F32 range_min  = 0.0;
    const PS::F32 range_max  = 6.0;
    const PS::S32 resolution = 100;

    const PS::F32 eps_abs = 1.e-4;
    const PS::F32 eps_rel = 1.e-3;

    const std::string rdf_AtoA_log_file{"test_bin/rdf_AtoA_log.tsv"};
    const std::string rdf_AtoA_ref_file{"unit_test/ref/rdf_AtoA_ref.tsv"};

    const std::string rdf_AtoB_log_file{"test_bin/rdf_AtoB_log.tsv"};
    const std::string rdf_AtoB_ref_file{"unit_test/ref/rdf_AtoB_ref.tsv"};

    const std::string rdf_ABtoAB_log_file{"test_bin/rdf_ABtoAB_log.tsv"};
    const std::string rdf_ABtoAB_ref_file{"unit_test/ref/rdf_ABtoAB_ref.tsv"};

    const PS::F32 range_shift = 3.0;

    const std::string rdf_ABtoABshift_log_file{"test_bin/rdf_ABtoABshift_log.tsv"};
    const std::string rdf_ABtoABshift_ref_file{"unit_test/ref/rdf_ABtoABshift_ref.tsv"};
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
    TGT_STATE        state;

    MD_DEFS::ID_type getId()  const { return this->id;  }
    PS::F32vec       getPos() const { return this->pos; }

    void setId( const MD_DEFS::ID_type  id     ) { this->id  = id;      }
    void setPos(const PS::F32vec       &pos_new) { this->pos = pos_new; }

    template <class Tptcl>
    void copyFromFP(const Tptcl &ptcl){
        this->id    = ptcl.getId();
        this->pos   = ptcl.getPos();
        this->state = ptcl.state;
    }
    template <class Tforce>
    void copyFromForce(const Tforce &force){}
    void clear(){}

    static PS::F32 r_cut;
    static PS::F32 getRSearch(){ return AtomMock::r_cut; }
};
PS::F32 AtomMock::r_cut = 0.0;


template <class Tpsys, class Tindex>
typename std::vector<Tindex>::const_iterator check_excluded_vol(const PS::F64vec           pos_new,
                                                                const Tpsys               &psys,
                                                                const PS::F64              ex_r_norm,
                                                                const std::vector<Tindex> &index_list){

    PS::F32 r2_ex = ex_r_norm*ex_r_norm;
    for(auto itr = index_list.begin(); itr != index_list.end(); ++itr){
        PS::F32vec r_vec = psys[*itr].getPos() - pos_new;
                   r_vec = Normalize::relativePosAdjustNorm(r_vec);
        PS::F32    r2    = r_vec*r_vec;

        if(r2 <= r2_ex){
            return itr;
        }
    }

    return index_list.end();
}

template <class Tpsys>
void make_test_particle(      Tpsys           &psys,
                              std::mt19937_64 &rand,
                        const PS::S64          n,
                        const PS::F32          ex_r_norm,
                        const PS::S32          try_limit){

    assert(0.0 < ex_r_norm && ex_r_norm < 1.0);
    assert(n         > 0);
    assert(try_limit > 0);

    psys.setNumberOfParticleLocal(n);

    std::uniform_real_distribution<PS::F32> dist_32(0.0, 1.0);
    MD_EXT::CellIndex<PS::S64>              cell_index;

    cell_index.init( PS::F32vec{0.0, 0.0, 0.0},
                     PS::F32vec{1.0, 1.0, 1.0},
                     ex_r_norm                 );

    for(PS::S64 i=0; i<n; ++i){
        PS::S32 try_count = 0;
        while(true){
            ++try_count;
            if(try_count >= try_limit){
                std::ostringstream oss;
                oss << "n or ex_r_norm are too large." << "\n"
                    << "    i    = " << i << " / " << n << "\n"
                    << "    ex_r = " << ex_r_norm << "\n";
                throw std::invalid_argument(oss.str());
            }

            PS::F32vec pos_new = PS::F32vec{dist_32(rand), dist_32(rand), dist_32(rand)};

            const auto index_list = cell_index.get_data_list(pos_new, ex_r_norm);
            auto       itr        = check_excluded_vol(pos_new, psys, ex_r_norm, index_list);

            if(itr == index_list.end()){
                psys[i].setId(i);
                psys[i].setPos(pos_new);
                cell_index.add(pos_new, i);
                break;
            }
        }
    }
}

template <class Tpsys>
void make_stripe(      Tpsys   &psys,
                 const PS::S32  n_stripe){

    const PS::F32 stripe_size_inv = n_stripe/1.0; // normalized box size: PS::F32vec{1.0, 1.0, 1.0}

    const PS::S64 n_atom = psys.getNumberOfParticleLocal();
    for(PS::S64 i=0; i<n_atom; ++i){
        if( (static_cast<int>(psys[i].getPos().x*stripe_size_inv) % 2) == 0 ){
            psys[i].state = TGT_STATE::A;
        } else {
            psys[i].state = TGT_STATE::B;
        }
    }
}

template <class Tdinfo, class Tpsys, class Ttree>
void init_system(      Tdinfo  &dinfo,
                       Tpsys   &psys,
                       Ttree   &tree,
                 const PS::S64  n_total){

    dinfo.initialize(0.3);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain( PS::F32vec{0.0, 0.0, 0.0},
                            PS::F32vec{1.0, 1.0, 1.0} );

    psys.initialize();
    psys.setNumberOfParticleLocal(0);

    tree.initialize(n_total, 0.5, 8, 64);

    Normalize::setBoxSize(TEST_DEFS::box_size);
}

struct RDF_Data {
    PS::F32 r;
    PS::F32 rdf;
    PS::F32 coordinate;

    bool read_line(std::string line){
        STR_TOOL::removeCR(line);
        const auto str_list = STR_TOOL::split(line, "\t");

        if(str_list.size() < 3)                 return false;
        if( !STR_TOOL::isNumeric(str_list[0]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[1]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[2]) ) return false;

        this->r          = std::stof(str_list[0]);
        this->rdf        = std::stof(str_list[1]);
        this->coordinate = std::stof(str_list[2]);

        return true;
    }
};

void rdf_result_compare(const std::string &file_log,
                        const std::string &file_ref){

    const PS::S32 rank = 0;
    if(PS::Comm::getRank() != rank) return;

    std::vector<RDF_Data> log_data;
    std::vector<RDF_Data> ref_data;

    FS_TOOL::file_load(file_log, log_data, rank);
    FS_TOOL::file_load(file_ref, ref_data, rank);

    ASSERT_EQ(log_data.size(), ref_data.size());
    for(size_t i=0; i<log_data.size(); ++i){
        EXPECT_TRUE( float_relative_eq(log_data[i].r,
                                       ref_data[i].r,
                                       TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i = " << i;
        EXPECT_TRUE( float_relative_eq(log_data[i].rdf,
                                       ref_data[i].rdf,
                                       TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i = " << i;
        EXPECT_TRUE( float_relative_eq(log_data[i].coordinate,
                                       ref_data[i].coordinate,
                                       TEST_DEFS::eps_abs, TEST_DEFS::eps_rel) ) << " i = " << i;
    }
}


TEST(RDF, AtoA){
    PS::DomainInfo                           dinfo;
    PS::ParticleSystem<AtomMock>             atom;
    PS::TreeForForceShort<AtomMock,
                          AtomMock,
                          AtomMock>::Scatter tree;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, tree, TEST_DEFS::n_atom);

    std::vector<MatchTarget> sight_list;
    std::vector<MatchTarget> tgt_list;

    sight_list = { MatchTarget{TGT_STATE::A} };
    tgt_list   = { MatchTarget{TGT_STATE::A} };

    Analysis::RDF_Sampler<MatchTarget> rdf_sampler;
    rdf_sampler.input_setting(TEST_DEFS::rdf_AtoA_log_file,
                              sight_list,
                              tgt_list,
                              TEST_DEFS::range_min,
                              TEST_DEFS::range_max,
                              TEST_DEFS::resolution       );

    rdf_sampler.broadcast_setting(0);

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
                               TEST_DEFS::n_atom,
                               Normalize::normCutOff(TEST_DEFS::ex_r),
                               5000);
            atom.adjustPositionIntoRootDomain(dinfo);
            make_stripe(atom, TEST_DEFS::n_stripe);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        atom[0].r_cut = Normalize::normCutOff( TEST_DEFS::range_max );
        rdf_sampler.sampling(atom, tree, dinfo);
        rdf_sampler.convert_count_to_raw( Normalize::getVol() );
    }

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "     ...done." << "\n";
        std::cout << oss.str() << std::flush;
    }

    rdf_sampler.output();
    rdf_result_compare(TEST_DEFS::rdf_AtoA_log_file,
                       TEST_DEFS::rdf_AtoA_ref_file);
}


TEST(RDF, AtoB){
    PS::DomainInfo                           dinfo;
    PS::ParticleSystem<AtomMock>             atom;
    PS::TreeForForceShort<AtomMock,
                          AtomMock,
                          AtomMock>::Scatter tree;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, tree, TEST_DEFS::n_atom);

    std::vector<MatchTarget> sight_list;
    std::vector<MatchTarget> tgt_list;

    sight_list = { MatchTarget{TGT_STATE::A} };
    tgt_list   = { MatchTarget{TGT_STATE::B} };

    Analysis::RDF_Sampler<MatchTarget> rdf_sampler;
    rdf_sampler.input_setting(TEST_DEFS::rdf_AtoB_log_file,
                              sight_list,
                              tgt_list,
                              TEST_DEFS::range_min,
                              TEST_DEFS::range_max,
                              TEST_DEFS::resolution       );

    rdf_sampler.broadcast_setting(0);

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
                               TEST_DEFS::n_atom,
                               Normalize::normCutOff(TEST_DEFS::ex_r),
                               5000);
            atom.adjustPositionIntoRootDomain(dinfo);
            make_stripe(atom, TEST_DEFS::n_stripe);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        atom[0].r_cut = Normalize::normCutOff( TEST_DEFS::range_max );
        rdf_sampler.sampling(atom, tree, dinfo);
        rdf_sampler.convert_count_to_raw( Normalize::getVol() );
    }

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "     ...done." << "\n";
        std::cout << oss.str() << std::flush;
    }

    rdf_sampler.output();
    rdf_result_compare(TEST_DEFS::rdf_AtoB_log_file,
                       TEST_DEFS::rdf_AtoB_ref_file);
}

TEST(RDF, ABtoAB){
    PS::DomainInfo                           dinfo;
    PS::ParticleSystem<AtomMock>             atom;
    PS::TreeForForceShort<AtomMock,
                          AtomMock,
                          AtomMock>::Scatter tree;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, tree, TEST_DEFS::n_atom);

    std::vector<MatchTarget> sight_list;
    std::vector<MatchTarget> tgt_list;

    sight_list = { MatchTarget{TGT_STATE::A}, MatchTarget{TGT_STATE::B} };
    tgt_list   = { MatchTarget{TGT_STATE::A}, MatchTarget{TGT_STATE::B} };

    Analysis::RDF_Sampler<MatchTarget> rdf_sampler;
    rdf_sampler.input_setting(TEST_DEFS::rdf_ABtoAB_log_file,
                              sight_list,
                              tgt_list,
                              TEST_DEFS::range_min,
                              TEST_DEFS::range_max,
                              TEST_DEFS::resolution       );

    rdf_sampler.broadcast_setting(0);

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
                               TEST_DEFS::n_atom,
                               Normalize::normCutOff(TEST_DEFS::ex_r),
                               5000);
            atom.adjustPositionIntoRootDomain(dinfo);
            make_stripe(atom, TEST_DEFS::n_stripe);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        atom[0].r_cut = Normalize::normCutOff( TEST_DEFS::range_max );
        rdf_sampler.sampling(atom, tree, dinfo);
        rdf_sampler.convert_count_to_raw( Normalize::getVol() );
    }

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "     ...done." << "\n";
        std::cout << oss.str() << std::flush;
    }

    rdf_sampler.output();
    rdf_result_compare(TEST_DEFS::rdf_ABtoAB_log_file,
                       TEST_DEFS::rdf_ABtoAB_ref_file);
}

TEST(RDF, ABtoABshift){
    PS::DomainInfo                           dinfo;
    PS::ParticleSystem<AtomMock>             atom;
    PS::TreeForForceShort<AtomMock,
                          AtomMock,
                          AtomMock>::Scatter tree;

    std::mt19937_64 mt;
    mt.seed( TEST_DEFS::mt_seed );

    init_system(dinfo, atom, tree, TEST_DEFS::n_atom);

    std::vector<MatchTarget> sight_list;
    std::vector<MatchTarget> tgt_list;

    sight_list = { MatchTarget{TGT_STATE::A}, MatchTarget{TGT_STATE::B} };
    tgt_list   = { MatchTarget{TGT_STATE::A}, MatchTarget{TGT_STATE::B} };

    Analysis::RDF_Sampler<MatchTarget> rdf_sampler;
    rdf_sampler.input_setting(TEST_DEFS::rdf_ABtoABshift_log_file,
                              sight_list,
                              tgt_list,
                              TEST_DEFS::range_min + TEST_DEFS::range_shift,
                              TEST_DEFS::range_max + TEST_DEFS::range_shift,
                              TEST_DEFS::resolution       );

    rdf_sampler.broadcast_setting(0);

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
                               TEST_DEFS::n_atom,
                               Normalize::normCutOff(TEST_DEFS::ex_r),
                               5000);
            atom.adjustPositionIntoRootDomain(dinfo);
            make_stripe(atom, TEST_DEFS::n_stripe);
        }

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        atom[0].r_cut = Normalize::normCutOff( TEST_DEFS::range_max + TEST_DEFS::range_shift );
        rdf_sampler.sampling(atom, tree, dinfo);
        rdf_sampler.convert_count_to_raw( Normalize::getVol() );
    }

    if(PS::Comm::getRank() == 0){
        std::ostringstream oss;
        oss << "     ...done." << "\n";
        std::cout << oss.str() << std::flush;
    }

    rdf_sampler.output();
    rdf_result_compare(TEST_DEFS::rdf_ABtoABshift_log_file,
                       TEST_DEFS::rdf_ABtoABshift_ref_file);
}


TEST(RDF, exception){
    std::vector<MatchTarget>           sight_list;
    std::vector<MatchTarget>           tgt_list;
    Analysis::RDF_Sampler<MatchTarget> rdf_sampler;

    EXPECT_THROW(rdf_sampler.input_setting("",
                                           sight_list,
                                           tgt_list,
                                           -1.0,
                                           -2.0,
                                           0,
                                           -1        ),
                 std::invalid_argument                 ) << " all";

    EXPECT_THROW(rdf_sampler.input_setting("a",
                                           sight_list,
                                           tgt_list,
                                           -1.0,
                                           -2.0,
                                           0,
                                           -1        ),
                 std::invalid_argument                 ) << " sight, tgt, range, resolution, rank";

    sight_list = { MatchTarget{TGT_STATE::A} };
    EXPECT_THROW(rdf_sampler.input_setting("a",
                                           sight_list,
                                           tgt_list,
                                           -1.0,
                                           -2.0,
                                           0,
                                           -1        ),
                 std::invalid_argument                 ) << " tgt, range, resolution, rank";

    tgt_list = { MatchTarget{TGT_STATE::A} };
    EXPECT_THROW(rdf_sampler.input_setting("a",
                                           sight_list,
                                           tgt_list,
                                           -1.0,
                                           -2.0,
                                           0,
                                           -1        ),
                 std::invalid_argument                 ) << " range, resolution, rank";

    EXPECT_THROW(rdf_sampler.input_setting("a",
                                           sight_list,
                                           tgt_list,
                                           0.0,
                                           -2.0,
                                           0,
                                           -1        ),
                 std::invalid_argument                 ) << " range_max(1), resolution, rank";

    EXPECT_THROW(rdf_sampler.input_setting("a",
                                           sight_list,
                                           tgt_list,
                                           0.0,
                                           0.0,
                                           0,
                                           -1        ),
                 std::invalid_argument                 ) << " range_max(2), resolution, rank";

    EXPECT_THROW(rdf_sampler.input_setting("a",
                                           sight_list,
                                           tgt_list,
                                           0.0,
                                           1.0,
                                           0,
                                           -1        ),
                 std::invalid_argument                 ) << " resolution, rank";

    EXPECT_THROW(rdf_sampler.input_setting("a",
                                           sight_list,
                                           tgt_list,
                                           0.0,
                                           1.0,
                                           2,
                                           -1        ),
                 std::out_of_range                     ) << " rank";

    EXPECT_NO_THROW(rdf_sampler.input_setting("a",
                                              sight_list,
                                              tgt_list,
                                              0.0,
                                              1.0,
                                              2,
                                              0           ) );
}

#include "gtest_main_mpi.hpp"
