//***************************************************************************************
//  This program is unit test of Clustering analyzer.
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

#include "analysis_clustering.hpp"
#include "analysis_target.hpp"

#include "gtest_common.hpp"


namespace TEST_DEFS {
    //--- setting for grid structure test
    const std::vector<std::pair<PS::F32, TGT_STATE>> edge = { {0.5, TGT_STATE::A},
                                                              {1.0, TGT_STATE::A},
                                                              {1.5, TGT_STATE::A},
                                                              {1.5, TGT_STATE::A},
                                                              {1.5, TGT_STATE::A},
                                                              {2.0, TGT_STATE::B},
                                                              {1.5, TGT_STATE::A},
                                                              {0.5, TGT_STATE::A},
                                                              {3.0, TGT_STATE::B},
                                                              {0.5, TGT_STATE::C} };  // terminator

    const PS::S64 n_grid = 4;

    const PS::F64 range       = 2.1;
    const size_t  n_pin_10x   = 6;
    const PS::S64 cycle_limit = 1000;

    const PS::F32 eps_abs = 1.e-4;
    const PS::F32 eps_rel = 1.e-3;

    const PS::S32 n_sample = 8;
    const PS::S32 n_blank  = 2;

    const std::vector<std::pair<PS::F32, TGT_STATE>> blank_edge = { {1.0, TGT_STATE::C},
                                                                    {1.0, TGT_STATE::C},
                                                                    {1.0, TGT_STATE::C},
                                                                    {1.0, TGT_STATE::C},
                                                                    {1.0, TGT_STATE::C},
                                                                    {1.0, TGT_STATE::C}, };  // state C is not target in this test.

    //--- for grid test data
    const std::string cls_naive_A_statistics_file{      "test_bin/cls_naive_A_statistics.tsv"};
    const std::string cls_naive_A_histogram_file{       "test_bin/cls_naive_A_histogram.tsv"};
    const std::string cls_naive_A_particle_file{        "test_bin/cls_naive_A_particle.pdb"};
    const std::string cls_naive_A_particle_extract_file{"test_bin/cls_naive_A_particle_extract.pdb"};

    const std::string cls_A_statistics_file{      "test_bin/cls_A_statistics.tsv"};
    const std::string cls_A_histogram_file{       "test_bin/cls_A_histogram.tsv"};
    const std::string cls_A_particle_file{        "test_bin/cls_A_particle.pdb"};
    const std::string cls_A_particle_extract_file{"test_bin/cls_A_particle_extract.pdb"};

    const std::string cls_A_wb_statistics_file{"test_bin/cls_A_wb_statistics.tsv"};
    const std::string cls_A_wb_histogram_file{ "test_bin/cls_A_wb_histogram.tsv" };

    const std::string cls_naive_AB_statistics_file{      "test_bin/cls_naive_AB_statistics.tsv"};
    const std::string cls_naive_AB_histogram_file{       "test_bin/cls_naive_AB_histogram.tsv"};
    const std::string cls_naive_AB_particle_file{        "test_bin/cls_naive_AB_particle.pdb"};
    const std::string cls_naive_AB_particle_extract_file{"test_bin/cls_naive_AB_particle_extract.pdb"};

    const std::string cls_AB_statistics_file{      "test_bin/cls_AB_statistics.tsv"};
    const std::string cls_AB_histogram_file{       "test_bin/cls_AB_histogram.tsv"};
    const std::string cls_AB_particle_file{        "test_bin/cls_AB_particle.pdb"};
    const std::string cls_AB_particle_extract_file{"test_bin/cls_AB_particle_extract.pdb"};

    const std::string cls_A_statistics_ref{      "unit_test/ref/cls_A_statistics.tsv"};
    const std::string cls_A_histogram_ref{       "unit_test/ref/cls_A_histogram.tsv"};
    const std::string cls_A_particle_ref{        "unit_test/ref/cls_A_particle.pdb"};
    const std::string cls_A_particle_extract_ref{"unit_test/ref/cls_A_particle_extract.pdb"};

    const std::string cls_A_wb_statistics_ref{"unit_test/ref/cls_A_wb_statistics.tsv"};
    const std::string cls_A_wb_histogram_ref{ "unit_test/ref/cls_A_wb_histogram.tsv" };

    const std::string cls_AB_statistics_ref{      "unit_test/ref/cls_AB_statistics.tsv"};
    const std::string cls_AB_histogram_ref{       "unit_test/ref/cls_AB_histogram.tsv"};
    const std::string cls_AB_particle_ref{        "unit_test/ref/cls_AB_particle.pdb"};
    const std::string cls_AB_particle_extract_ref{"unit_test/ref/cls_AB_particle_extract.pdb"};
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
        return this->state.isMatch(ptcl.getAtomType());
    }
};

class AtomMock {
public:
    MD_DEFS::ID_type id{-1};
    PS::F32vec       pos{0.0, 0.0, 0.0};
    TGT_STATE        type;

    MD_DEFS::ID_type getId()       const { return this->id;   }
    MD_DEFS::ID_type getAtomID()   const { return this->id;   }
    TGT_STATE        getAtomType() const { return this->type; }
    PS::F32vec       getPos()      const { return this->pos;  }

    void setId(      const MD_DEFS::ID_type  id     ) { this->id   = id;      }
    void setPos(     const PS::F32vec       &pos_new) { this->pos  = pos_new; }
    void setAtomType(const TGT_STATE         type   ) { this->type = type;    }

    template <class Tptcl>
    void copyFromFP(const Tptcl &ptcl){
        this->id   = ptcl.getId();
        this->pos  = ptcl.getPos();
        this->type = ptcl.getAtomType();
    }
    template <class Tforce>
    void copyFromForce(const Tforce &force){}
    void clear(){}
};

class ProfileMock {
public:
    PS::F64 time = 0.0;

    PS::F64 theta         = 0.5;
    PS::S32 n_leaf_limit  = 8;
    PS::S32 n_group_limit = 64;

    void    set_time(const PS::F64 t)       { this->time = t;    }
    PS::F64 get_time()                const { return this->time; }
};

template <class Tpsys, class Tedge>
void make_test_particle_grid(      Tpsys   &psys,
                             const Tedge   &edge_disp,
                             const PS::S32  n         ){

    assert(n                > 0);
    assert(edge_disp.size() > 0);

    //--- make edge ( displacement -> position (at real space) )
    auto edge = edge_disp;
    for(size_t i=1; i<edge.size(); ++i){
        edge.at(i).first += edge.at(i-1).first;
    }

    const PS::F32 len_edge = edge.back().first;

    //--- make cell template
    PS::S64 id = 0;
    std::vector<AtomMock> cell_template;

    //------ center of cell
    cell_template.resize(1);
    cell_template[0].setId(id);
    cell_template[0].setPos( PS::F32vec{ static_cast<PS::F32>(0.5*len_edge),
                                         static_cast<PS::F32>(0.5*len_edge),
                                         static_cast<PS::F32>(0.5*len_edge) } );
    cell_template[0].setAtomType( TGT_STATE::A );

    //------ edge for x,y,z
    for(const auto& p : edge){
        for(PS::S32 dir=0; dir<3; ++dir){
            ++id;

            PS::F32vec pos_tmp{0.0, 0.0, 0.0};
            pos_tmp[dir] = p.first;

            AtomMock atom_tmp;
            atom_tmp.setId(id);
            atom_tmp.setPos(pos_tmp);
            atom_tmp.setAtomType(p.second);

            cell_template.push_back(atom_tmp);
        }
    }
    //--- check sum
    assert(static_cast<PS::S64>(cell_template.size())-1 == id);

/*
    std::ostringstream oss;
    oss << "cell_template:" << "\n";
    for(const auto& p : cell_template){
        oss << "  id  = " << p.getId() << "\n"
            << "    pos  = " << p.getPos() << "\n"
            << "    type = " << p.getAtomType() << "\n";
    }
    std::cout << oss.str() << std::flush;
*/

    //--- make grid
    const PS::S64 n_total = cell_template.size()*n*n*n;
    psys.setNumberOfParticleLocal(n_total);

    Normalize::setBoxSize( PS::F64vec{len_edge*n, len_edge*n, len_edge*n} );

    PS::S64    id_cell = 0;
    PS::F32vec pos_cell{0.0, 0.0, 0.0};
    for(PS::S32 iz=0; iz<n; ++iz){
        for(PS::S32 iy=0; iy<n; ++iy){
            for(PS::S32 ix=0; ix<n; ++ix){
                pos_cell = { len_edge*ix, len_edge*iy, len_edge*iz };

                for(PS::S64 i=0; i<static_cast<PS::S64>(cell_template.size()); ++i){
                    const PS::S64    index   = id_cell + i;
                    const PS::F32vec pos_tmp = Normalize::normPos( pos_cell + cell_template.at(i).getPos() );
                    assert(index < n_total);

                    psys[index].setId( id_cell + cell_template.at(i).getId() );
                    psys[index].setPos( pos_tmp );
                    psys[index].setAtomType( cell_template.at(i).getAtomType() );
                }

                id_cell += cell_template.size();
            }
        }
    }
    //--- check sum
    assert(n_total == id_cell);
/*
    oss.str("");
    oss << "psys:" << "\n";
    for(PS::S64 i=0; i<n_total; ++i){
        oss << "  id  = "    << psys[i].getId() << "\n"
            << "    pos  = " << psys[i].getPos() << "\n"
            << "    type = " << psys[i].getAtomType() << "\n";
    }
    std::cout << oss.str() << std::flush;
*/
}

template <class Tdinfo, class Tpsys>
void init_system(      Tdinfo  &dinfo,
                       Tpsys   &psys  ){

    dinfo.initialize(0.3);
    dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
    dinfo.setPosRootDomain( PS::F32vec{0.0, 0.0, 0.0},
                            PS::F32vec{1.0, 1.0, 1.0} );

    psys.initialize();
    psys.setNumberOfParticleLocal(0);
}

struct StatisticsData{
    PS::F64 time;
    PS::S64 n_cls;
    PS::S64 max_size;
    PS::F64 ave_size;
    PS::F64 ave_size_woMax;

    bool read_line(const std::string &line){
        auto line_tmp = line;
        STR_TOOL::removeCR(line_tmp);
        const auto str_list = STR_TOOL::split(line_tmp, "\t");

        if(str_list.size() < 5)                 return false;
        if( !STR_TOOL::isNumeric(str_list[0]) ) return false;
        if( !STR_TOOL::isInteger(str_list[1]) ) return false;
        if( !STR_TOOL::isInteger(str_list[2]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[3]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[4]) ) return false;

        this->time           = std::stod(str_list[0]);
        this->n_cls          = std::stoi(str_list[1]);
        this->max_size       = std::stoi(str_list[2]);
        this->ave_size       = std::stod(str_list[3]);
        this->ave_size_woMax = std::stod(str_list[4]);

        return true;
    }
};
void clustering_statistics_compare(const std::string &file_log,
                                   const std::string &file_ref){

    const PS::S32 rank = 0;
    if(PS::Comm::getRank() != rank) return;

    std::vector<StatisticsData> log_data;
    std::vector<StatisticsData> ref_data;

    FS_TOOL::file_load(file_log, log_data, rank);
    FS_TOOL::file_load(file_ref, ref_data, rank);

    ASSERT_EQ(log_data.size(), ref_data.size());
    for(size_t i=0; i<log_data.size(); ++i){
        EXPECT_FLOAT_EQ(log_data[i].time          , ref_data[i].time          ) << " i = " << i;
        EXPECT_EQ(      log_data[i].n_cls         , ref_data[i].n_cls         ) << " i = " << i;
        EXPECT_EQ(      log_data[i].max_size      , ref_data[i].max_size      ) << " i = " << i;
        EXPECT_FLOAT_EQ(log_data[i].ave_size      , ref_data[i].ave_size      ) << " i = " << i;
        EXPECT_FLOAT_EQ(log_data[i].ave_size_woMax, ref_data[i].ave_size_woMax) << " i = " << i;
    }
}

struct HistogramData{
    PS::F64 cls_size;
    PS::F64 probability;

    bool read_line(const std::string &line){
        auto line_tmp = line;
        STR_TOOL::removeCR(line_tmp);
        const auto str_list = STR_TOOL::split(line_tmp, "\t");

        if(str_list.size() < 2)                 return false;
        if( !STR_TOOL::isNumeric(str_list[0]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[1]) ) return false;

        this->cls_size    = std::stod(str_list[0]);
        this->probability = std::stod(str_list[1]);

        return true;
    }
};
void clustering_histogram_compare(const std::string &file_log,
                                  const std::string &file_ref){

    const PS::S32 rank = 0;
    if(PS::Comm::getRank() != rank) return;

    std::vector<HistogramData> log_data;
    std::vector<HistogramData> ref_data;

    FS_TOOL::file_load(file_log, log_data, rank);
    FS_TOOL::file_load(file_ref, ref_data, rank);

    ASSERT_EQ(log_data.size(), ref_data.size());
    for(size_t i=0; i<log_data.size(); ++i){
        EXPECT_FLOAT_EQ(log_data[i].cls_size   , ref_data[i].cls_size   ) << " i = " << i;
        EXPECT_FLOAT_EQ(log_data[i].probability, ref_data[i].probability) << " i = " << i;
    }
}

struct PtclData {
    PS::S64     id;
    TGT_STATE   state;
    std::string res;
    PS::S64     cls_id;
    PS::F32vec  pos;

    bool read_line(const std::string &line){
        auto line_tmp = line;
        STR_TOOL::removeCR(line_tmp);
        const auto str_list = STR_TOOL::split(line_tmp, "\t");

        if(str_list.size() < 7)                 return false;
        if( !STR_TOOL::isInteger(str_list[1]) ) return false;
        if( !STR_TOOL::isInteger(str_list[4]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[5]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[6]) ) return false;
        if( !STR_TOOL::isNumeric(str_list[7]) ) return false;

        this->id     = std::stoi(str_list[1]);
        this->state  = ENUM::which_TGT_STATE(str_list[2]);
        this->res    = str_list[3];
        this->cls_id = std::stoi(str_list[4]);
        this->pos    = PS::F32vec{ std::stof(str_list[5]),
                                   std::stof(str_list[6]),
                                   std::stof(str_list[7]) };

        return true;
    }
};
//--- for std::sort()
bool operator < (const PtclData &lhs, const PtclData &rhs){
    return lhs.id < rhs.id;
}

void clustering_ptcl_compare(const std::string &file_log,
                             const std::string &file_ref){

    const PS::S32 rank = 0;
    if(PS::Comm::getRank() != rank) return;

    std::vector<PtclData> log_data;
    std::vector<PtclData> ref_data;

    FS_TOOL::file_load(file_log, log_data, rank);
    FS_TOOL::file_load(file_ref, ref_data, rank);

    std::sort(log_data.begin(), log_data.end());
    std::sort(ref_data.begin(), ref_data.end());

    ASSERT_EQ(log_data.size(), ref_data.size());
    for(size_t i=0; i<log_data.size(); ++i){
        EXPECT_EQ(log_data[i].id    , ref_data[i].id    ) << " i = " << i;
        EXPECT_EQ(log_data[i].state , ref_data[i].state ) << " i = " << i;
        EXPECT_EQ(log_data[i].res   , ref_data[i].res   ) << " i = " << i;
        EXPECT_EQ(log_data[i].cls_id, ref_data[i].cls_id) << " i = " << i;
        EXPECT_FLOAT_EQ(log_data[i].pos.x, ref_data[i].pos.x) << " i = " << i;
        EXPECT_FLOAT_EQ(log_data[i].pos.y, ref_data[i].pos.y) << " i = " << i;
        EXPECT_FLOAT_EQ(log_data[i].pos.z, ref_data[i].pos.z) << " i = " << i;
    }
}


TEST(Clustering, NaivePtclA){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    init_system(dinfo, atom);

    std::vector<MatchTarget> tgt_list = { MatchTarget{TGT_STATE::A} };

    Analysis::Clustering<MatchTarget, AtomMock> clustering_analyzer;
    clustering_analyzer.input_setting(TEST_DEFS::cls_naive_A_statistics_file,
                                      TEST_DEFS::cls_naive_A_histogram_file,
                                      tgt_list,
                                      TEST_DEFS::range,
                                      TEST_DEFS::n_pin_10x,
                                      TEST_DEFS::cycle_limit);
    clustering_analyzer.broadcast_setting();

    //--- make test data
    atom.setNumberOfParticleLocal(0);
    if(PS::Comm::getRank() == 0){
        make_test_particle_grid(atom, TEST_DEFS::edge, TEST_DEFS::n_grid);
    }
    Normalize::broadcast_boxSize(0);
    atom.adjustPositionIntoRootDomain(dinfo);

    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    ProfileMock time_prof;
    time_prof.set_time(0.0);

    //--- sampling
    clustering_analyzer.add_sample_naive(time_prof, atom, dinfo);

    //--- output
    clustering_analyzer.output_histogram();
    clustering_analyzer.output_particle(TEST_DEFS::cls_naive_A_particle_file);
    clustering_analyzer.output_particle(TEST_DEFS::cls_naive_A_particle_extract_file, true);

    //--- check result
    clustering_statistics_compare(TEST_DEFS::cls_naive_A_statistics_file,
                                  TEST_DEFS::cls_A_statistics_ref        );
    clustering_histogram_compare(TEST_DEFS::cls_naive_A_histogram_file,
                                 TEST_DEFS::cls_A_histogram_ref        );
    clustering_ptcl_compare(TEST_DEFS::cls_naive_A_particle_file,
                            TEST_DEFS::cls_A_particle_ref        );
    clustering_ptcl_compare(TEST_DEFS::cls_naive_A_particle_extract_file,
                            TEST_DEFS::cls_A_particle_extract_ref        );
}
TEST(Clustering, PtclA){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    init_system(dinfo, atom);

    std::vector<MatchTarget> tgt_list = { MatchTarget{TGT_STATE::A} };

    Analysis::Clustering<MatchTarget, AtomMock> clustering_analyzer;
    clustering_analyzer.input_setting(TEST_DEFS::cls_A_statistics_file,
                                      TEST_DEFS::cls_A_histogram_file,
                                      tgt_list,
                                      TEST_DEFS::range);
    clustering_analyzer.broadcast_setting();

    //--- make test data
    atom.setNumberOfParticleLocal(0);
    if(PS::Comm::getRank() == 0){
        make_test_particle_grid(atom, TEST_DEFS::edge, TEST_DEFS::n_grid);
    }
    Normalize::broadcast_boxSize(0);
    atom.adjustPositionIntoRootDomain(dinfo);

    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    ProfileMock time_prof;
    time_prof.set_time(0.0);

    //--- sampling
    clustering_analyzer.add_sample(time_prof, atom, dinfo);

    //--- output
    clustering_analyzer.output_histogram();
    clustering_analyzer.output_particle(TEST_DEFS::cls_A_particle_file);
    clustering_analyzer.output_particle(TEST_DEFS::cls_A_particle_extract_file, true);

    //--- check result
    clustering_statistics_compare(TEST_DEFS::cls_A_statistics_file,
                                  TEST_DEFS::cls_A_statistics_ref  );
    clustering_histogram_compare(TEST_DEFS::cls_A_histogram_file,
                                 TEST_DEFS::cls_A_histogram_ref  );
    clustering_ptcl_compare(TEST_DEFS::cls_A_particle_file,
                            TEST_DEFS::cls_A_particle_ref  );
    clustering_ptcl_compare(TEST_DEFS::cls_A_particle_extract_file,
                            TEST_DEFS::cls_A_particle_extract_ref  );
}
TEST(Clustering, PtclAwithBlank){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    init_system(dinfo, atom);

    std::vector<MatchTarget> tgt_list = { MatchTarget{TGT_STATE::A} };

    Analysis::Clustering<MatchTarget, AtomMock> clustering_analyzer;
    clustering_analyzer.input_setting(TEST_DEFS::cls_A_wb_statistics_file,
                                      TEST_DEFS::cls_A_wb_histogram_file,
                                      tgt_list,
                                      TEST_DEFS::range);
    clustering_analyzer.broadcast_setting();

    //--- make test data
    for(PS::S32 i_sample=0; i_sample<TEST_DEFS::n_sample; ++i_sample){
        atom.setNumberOfParticleLocal(0);
        if(PS::Comm::getRank() == 0){
            make_test_particle_grid(atom, TEST_DEFS::edge, TEST_DEFS::n_grid);
        }
        Normalize::broadcast_boxSize(0);
        atom.adjustPositionIntoRootDomain(dinfo);

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        ProfileMock time_prof;
        time_prof.set_time(i_sample);

        //--- sampling
        clustering_analyzer.add_sample_naive(time_prof, atom, dinfo);
    }
    //--- adding blank sample
    for(PS::S32 i_sample=0; i_sample<TEST_DEFS::n_blank; ++i_sample){
        atom.setNumberOfParticleLocal(0);
        if(PS::Comm::getRank() == 0){
            std::cout << "-- input blank data --" << std::endl;
            make_test_particle_grid(atom, TEST_DEFS::blank_edge, TEST_DEFS::n_grid);
        }
        Normalize::broadcast_boxSize(0);
        atom.adjustPositionIntoRootDomain(dinfo);

        dinfo.decomposeDomainAll(atom);
        atom.exchangeParticle(dinfo);

        ProfileMock time_prof;
        time_prof.set_time(i_sample + TEST_DEFS::n_sample);

        //--- sampling
        clustering_analyzer.add_sample_naive(time_prof, atom, dinfo);
    }

    //--- output
    clustering_analyzer.output_histogram();

    //--- check result
    clustering_statistics_compare(TEST_DEFS::cls_A_wb_statistics_file,
                                  TEST_DEFS::cls_A_wb_statistics_ref  );
    clustering_histogram_compare(TEST_DEFS::cls_A_wb_histogram_file,
                                 TEST_DEFS::cls_A_wb_histogram_ref  );
}

TEST(Clustering, NaivePtclAB){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    init_system(dinfo, atom);

    std::vector<MatchTarget> tgt_list = { MatchTarget{TGT_STATE::A},
                                          MatchTarget{TGT_STATE::B} };

    Analysis::Clustering<MatchTarget, AtomMock> clustering_analyzer;
    clustering_analyzer.input_setting(TEST_DEFS::cls_naive_AB_statistics_file,
                                      TEST_DEFS::cls_naive_AB_histogram_file,
                                      tgt_list,
                                      TEST_DEFS::range);
    clustering_analyzer.broadcast_setting();

    //--- make test data
    atom.setNumberOfParticleLocal(0);
    if(PS::Comm::getRank() == 0){
        make_test_particle_grid(atom, TEST_DEFS::edge, TEST_DEFS::n_grid);
    }
    Normalize::broadcast_boxSize(0);
    atom.adjustPositionIntoRootDomain(dinfo);

    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    ProfileMock time_prof;
    time_prof.set_time(0.0);

    //--- sampling
    clustering_analyzer.add_sample_naive(time_prof, atom, dinfo);

    //--- output
    clustering_analyzer.output_histogram();
    clustering_analyzer.output_particle(TEST_DEFS::cls_naive_AB_particle_file);
    clustering_analyzer.output_particle(TEST_DEFS::cls_naive_AB_particle_extract_file, true);

    //--- check result
    clustering_statistics_compare(TEST_DEFS::cls_naive_AB_statistics_file,
                                  TEST_DEFS::cls_AB_statistics_ref        );
    clustering_histogram_compare(TEST_DEFS::cls_naive_AB_histogram_file,
                                 TEST_DEFS::cls_AB_histogram_ref        );
    clustering_ptcl_compare(TEST_DEFS::cls_naive_AB_particle_file,
                            TEST_DEFS::cls_AB_particle_ref        );
    clustering_ptcl_compare(TEST_DEFS::cls_naive_AB_particle_extract_file,
                            TEST_DEFS::cls_AB_particle_extract_ref        );
}
TEST(Clustering, PtclAB){
    PS::DomainInfo               dinfo;
    PS::ParticleSystem<AtomMock> atom;

    init_system(dinfo, atom);

    std::vector<MatchTarget> tgt_list = { MatchTarget{TGT_STATE::A},
                                          MatchTarget{TGT_STATE::B} };

    Analysis::Clustering<MatchTarget, AtomMock> clustering_analyzer;
    clustering_analyzer.input_setting(TEST_DEFS::cls_AB_statistics_file,
                                      TEST_DEFS::cls_AB_histogram_file,
                                      tgt_list,
                                      TEST_DEFS::range);
    clustering_analyzer.broadcast_setting();

    //--- make test data
    atom.setNumberOfParticleLocal(0);
    if(PS::Comm::getRank() == 0){
        make_test_particle_grid(atom, TEST_DEFS::edge, TEST_DEFS::n_grid);
    }
    Normalize::broadcast_boxSize(0);
    atom.adjustPositionIntoRootDomain(dinfo);

    dinfo.decomposeDomainAll(atom);
    atom.exchangeParticle(dinfo);

    ProfileMock time_prof;
    time_prof.set_time(0.0);

    //--- sampling
    clustering_analyzer.add_sample(time_prof, atom, dinfo);

    //--- output
    clustering_analyzer.output_histogram();
    clustering_analyzer.output_particle(TEST_DEFS::cls_AB_particle_file);
    clustering_analyzer.output_particle(TEST_DEFS::cls_AB_particle_extract_file, true);

    //--- check result
    clustering_statistics_compare(TEST_DEFS::cls_AB_statistics_file,
                                  TEST_DEFS::cls_AB_statistics_ref  );
    clustering_histogram_compare(TEST_DEFS::cls_AB_histogram_file,
                                 TEST_DEFS::cls_AB_histogram_ref  );
    clustering_ptcl_compare(TEST_DEFS::cls_AB_particle_file,
                            TEST_DEFS::cls_AB_particle_ref  );
    clustering_ptcl_compare(TEST_DEFS::cls_AB_particle_extract_file,
                            TEST_DEFS::cls_AB_particle_extract_ref  );
}


#include "gtest_main_mpi.hpp"
