//=======================================================================================
//  This is unit test of additional wrapper of PS::Comm::.
//     provides the interface for some STL container.
//     module location: ./generic_ext/comm_tool_allToAll.hpp
//=======================================================================================

#include <gtest/gtest.h>
#include <particle_simulator.hpp>
#include "comm_tool_allToAll.hpp"

#include <random>


namespace TEST_DEFS {
    const PS::S64 mt_seed = 7654321;
    const PS::S64 n_data  = 10000;
}

//==========================================
// MPI allToAll
//==========================================
struct DataBasic {
    std::vector<std::vector<int>>                   vec_vi;
    std::vector<std::vector<std::vector<int>>>      vec_vvi;
    std::vector<std::string>                        vec_s;
    std::vector<std::vector<std::string>>           vec_vs;
    std::vector<std::vector<std::pair<int, float>>> vec_vp_if;

    DataBasic() = default;
    ~DataBasic() = default;

    template <class Trand>
    void generate(const size_t  n_proc,
                  const size_t  n_data,
                        Trand  &mt     ){
        std::uniform_int_distribution<>  dist_int(0,10);
        std::uniform_real_distribution<> dist_real(-99.9,99.9);

        this->vec_vi.resize(n_proc);
        this->vec_vvi.resize(n_proc);
        this->vec_s.resize(n_proc);
        this->vec_vs.resize(n_proc);
        this->vec_vp_if.resize(n_proc);

        std::vector<int> buff_vi;
        std::string      buff_str;

        for(size_t i_proc=0; i_proc<n_proc; ++i_proc){
            auto& vec_int     = this->vec_vi[i_proc];
            auto& vec_vec_int = this->vec_vvi[i_proc];
            auto& str         = this->vec_s[i_proc];
            auto& vec_str     = this->vec_vs[i_proc];
            auto& vec_pair_if = this->vec_vp_if[i_proc];

            vec_int.clear();
            vec_vec_int.clear();
            str.clear();
            vec_str.clear();
            vec_pair_if.clear();

            buff_vi.clear();
            buff_str.clear();

            //--- for basic pattern
            for(size_t i=0; i<n_data; ++i){
                int   data_int  = dist_int(mt);
                float data_real = dist_real(mt);

                vec_int.push_back(data_int);

                if(data_int == 10){
                    vec_vec_int.push_back(buff_vi);
                    vec_str.push_back(buff_str);
                    buff_vi.clear();
                    buff_str = "";
                } else {
                    buff_vi.push_back(data_int);
                    buff_str += std::to_string(data_int);

                    str += std::to_string(data_int);
                }

                vec_pair_if.push_back( std::make_pair( data_int, data_real ) );
            }
        }

    }
};

class AllToAllBasic :
    public ::testing::Test{
    protected:
        //--- for basic data pattern
        std::vector<DataBasic> data;

        virtual void SetUp(){
            const int n_proc  = PS::Comm::getNumberOfProc();

            std::mt19937 mt;
            mt.seed(TEST_DEFS::mt_seed);

            this->data.resize(n_proc);
            for(int i=0; i<n_proc; ++i){
                this->data[i].generate(n_proc, TEST_DEFS::n_data, mt);
            }
        }
};

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(AllToAllBasic, VecInt){
    const int n_proc = PS::Comm::getNumberOfProc();
    const int rank   = PS::Comm::getRank();
    auto recv_vec_vi = COMM_TOOL::allToAll(data[rank].vec_vi);

    static_assert(std::is_same<decltype(recv_vec_vi), std::vector<std::vector<int>> >::value == true, "auto type check");

    for(int i=0; i<n_proc; ++i){
        EXPECT_EQ(recv_vec_vi.at(i), data.at(i).vec_vi.at(rank) ) << "from = " << i << ", to = " << rank;
    }
}

TEST_F(AllToAllBasic, VecVecInt){
    const int n_proc  = PS::Comm::getNumberOfProc();
    const int rank    = PS::Comm::getRank();
    auto recv_vec_vvi = COMM_TOOL::allToAll(data[rank].vec_vvi);

    static_assert(std::is_same<decltype(recv_vec_vvi),
                               std::vector<std::vector<std::vector<int>>> >::value == true, "auto type check");

    for(int i=0; i<n_proc; ++i){
        EXPECT_EQ(recv_vec_vvi.at(i), data.at(i).vec_vvi.at(rank) ) << "from = " << i << ", to = " << rank;
    }
}

TEST_F(AllToAllBasic, Str){
    const int n_proc = PS::Comm::getNumberOfProc();
    const int rank   = PS::Comm::getRank();

    auto recv_vs = COMM_TOOL::allToAll(data[rank].vec_s);

    static_assert(std::is_same<decltype(recv_vs),
                               std::vector<std::string> >::value == true, "auto type check");

    for(int i=0; i<n_proc; ++i){
        EXPECT_EQ(recv_vs.at(i), data.at(i).vec_s.at(rank) ) << "from = " << i << ", to = " << rank;
    }
}

TEST_F(AllToAllBasic, VecStr){
    const int n_proc = PS::Comm::getNumberOfProc();
    const int rank   = PS::Comm::getRank();
    auto recv_vec_vs = COMM_TOOL::allToAll(data[rank].vec_vs);

    static_assert(std::is_same<decltype(recv_vec_vs),
                               std::vector<std::vector<std::string>> >::value == true, "auto type check");

    for(int i=0; i<n_proc; ++i){
        EXPECT_EQ(recv_vec_vs.at(i), data.at(i).vec_vs.at(rank) ) << "from = " << i << ", to = " << rank;
    }
}

TEST_F(AllToAllBasic, VecPairIntFloat){
    const int n_proc    = PS::Comm::getNumberOfProc();
    const int rank      = PS::Comm::getRank();
    auto recv_vec_vp_if = COMM_TOOL::allToAll(data[rank].vec_vp_if);

    static_assert(std::is_same<decltype(recv_vec_vp_if),
                               std::vector<std::vector<std::pair<int,
                                                                 float>>> >::value == true, "auto type check");

    for(int i=0; i<n_proc; ++i){
        EXPECT_EQ(recv_vec_vp_if.at(i), data.at(i).vec_vp_if.at(rank) ) << "from = " << i << ", to = " << rank;
    }
}


struct DataRecursive {
    std::vector<std::vector<std::vector<std::vector<int>>>>            vec_vvvi;
    std::vector<std::vector<std::vector<std::string>>>                 vec_vvs;
    std::vector<std::vector<std::pair<std::string, std::vector<int>>>> vec_vp_s_vi;

    DataRecursive() = default;
    ~DataRecursive() = default;

    template<class Trand>
    void generate(const size_t  n_proc,
                  const size_t  n_data,
                        Trand  &mt     ){
        std::uniform_int_distribution<>  dist_int(0,11);

        this->vec_vvvi.resize(n_proc);
        this->vec_vvs.resize(n_proc);
        this->vec_vp_s_vi.resize(n_proc);

        std::vector<std::vector<int>> buff_v_vi;
        std::vector<std::string>      buff_v_str;
        std::vector<int> buff_vi;
        std::string      buff_str;

        for(size_t i_proc=0; i_proc<n_proc; ++i_proc){
            auto& vec_vec_vec_int = this->vec_vvvi[i_proc];
            auto& vec_vec_str     = this->vec_vvs[i_proc];
            auto& vec_pair_s_vi   = this->vec_vp_s_vi[i_proc];

            vec_vec_vec_int.clear();
            vec_vec_str.clear();
            vec_pair_s_vi.clear();

            buff_v_vi.clear();
            buff_v_str.clear();
            buff_vi.clear();
            buff_str.clear();

            for(size_t i=0; i<n_data; ++i){
                const int data = dist_int(mt);

                if(       data == 11){
                    vec_vec_vec_int.push_back(buff_v_vi);
                    vec_vec_str.push_back(buff_v_str);

                    buff_v_vi.clear();
                    buff_v_str.clear();

                } else if(data == 10){
                    buff_v_vi.push_back(buff_vi);
                    buff_v_str.push_back(buff_str);

                    vec_pair_s_vi.push_back( std::make_pair( buff_str,
                                                             buff_vi ) );

                    buff_vi.clear();
                    buff_str = "";

                } else {
                    buff_vi.push_back(data);
                    buff_str += std::to_string(data);
                }
            }
        }
    }
};

class AllToAllRecursive :
    public ::testing::Test{
    protected:
        //--- for basic data pattern
        std::vector<DataRecursive> data;

        virtual void SetUp(){
            const int n_proc = PS::Comm::getNumberOfProc();

            std::mt19937 mt;
            mt.seed(TEST_DEFS::mt_seed);

            this->data.resize(n_proc);
            for(int i=0; i<n_proc; ++i){;
                this->data[i].generate(n_proc, TEST_DEFS::n_data, mt);
            }
        }
};

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(AllToAllRecursive, VecVecVecInt){
    const int n_proc   = PS::Comm::getNumberOfProc();
    const int rank     = PS::Comm::getRank();
    auto recv_vec_vvvi = COMM_TOOL::allToAll(data[rank].vec_vvvi);

    static_assert(std::is_same<decltype(recv_vec_vvvi),
                               std::vector<std::vector<std::vector<std::vector<int>>>> >::value == true, "auto type check");

    for(int i=0; i<n_proc; ++i){
        EXPECT_EQ(recv_vec_vvvi.at(i), data.at(i).vec_vvvi.at(rank)) <<  "from = " << i << ", to = " << rank;
    }
}

TEST_F(AllToAllRecursive, VecVecStr){
    const int n_proc  = PS::Comm::getNumberOfProc();
    const int rank    = PS::Comm::getRank();
    auto recv_vec_vvs = COMM_TOOL::allToAll(data[rank].vec_vvs);

    static_assert(std::is_same<decltype(recv_vec_vvs),
                               std::vector<std::vector<std::vector<std::string>>> >::value == true, "auto type check");

    for(int i=0; i<n_proc; ++i){
        EXPECT_EQ(recv_vec_vvs.at(i), data.at(i).vec_vvs.at(rank)) <<  "from = " << i << ", to = " << rank;
    }
}

TEST_F(AllToAllRecursive, VecPairStrVecInt){
    const int n_proc      = PS::Comm::getNumberOfProc();
    const int rank        = PS::Comm::getRank();
    auto recv_vec_vp_s_vi = COMM_TOOL::allToAll(data[rank].vec_vp_s_vi);

    static_assert(std::is_same<decltype(recv_vec_vp_s_vi),
                               std::vector<std::vector<std::pair<std::string,
                                                                 std::vector<int>>>> >::value == true, "auto type check");

    for(int i=0; i<n_proc; ++i){
        EXPECT_EQ(recv_vec_vp_s_vi.at(i), data.at(i).vec_vp_s_vi.at(rank)) <<  "from = " << i << ", to = " << rank;
    }
}


#include "gtest_main_mpi.hpp"
