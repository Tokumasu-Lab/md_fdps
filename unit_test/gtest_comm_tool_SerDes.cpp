//=======================================================================================
//  This is unit test of internal tools for COMM_TOOL::.
//     provides serialize/deserialize interface for some STL container.
//     module location: ./generic_ext/comm_tool_SerDes.hpp
//=======================================================================================

#include <gtest/gtest.h>
#include "comm_tool_SerDes.hpp"

#include <random>


namespace TEST_DEFS {
    const PS::S64 mt_seed = 7654321;
    const PS::S64 n_data  = 100000;
}

//==========================================
// for Ser/Des std::vector<std::vector<T>>
//==========================================
class SerDes :
    public ::testing::Test{
    protected:
        std::vector<std::vector<int>>       vec_vi;
        std::vector<std::pair<int, double>> vec_p_i_f;
        std::vector<std::string>            vec_str;

        virtual void SetUp(){
            std::mt19937 mt;
            std::uniform_int_distribution<>  dist_int(0,10);
            std::uniform_real_distribution<> dist_real(-99.0, 99.0);
            mt.seed(TEST_DEFS::mt_seed);

            this->vec_vi.clear();
            this->vec_p_i_f.clear();
            this->vec_str.clear();

            std::vector<int> buff_vi;
            std::string      buff_str;
            buff_vi.clear();
            buff_str.clear();
            for(int64_t i=0; i<TEST_DEFS::n_data; ++i){
                int data = dist_int(mt);

                if(data == 10){
                    this->vec_vi.push_back(buff_vi);
                    this->vec_str.push_back(buff_str);
                    buff_vi.clear();
                    buff_str.clear();
                } else {
                    buff_vi.push_back(data);
                    buff_str += std::to_string(data);
                }

                this->vec_p_i_f.push_back( std::make_pair( dist_int(mt),
                                                           dist_real(mt) ) );
            }
        }
};

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST_F(SerDes, VecVec){
    std::vector<int>    serial_data;
    std::vector<size_t> serial_index;

    std::vector<std::vector<int>> result;

    COMM_TOOL::serialize_vector_vector(vec_vi,
                                       serial_data,
                                       serial_index);

    COMM_TOOL::deserialize_vector_vector(serial_data,
                                         serial_index,
                                         result       );

    EXPECT_EQ(result, vec_vi);
}

TEST_F(SerDes, VecPair){
    std::vector<int>    serial_Ta;
    std::vector<double> serial_Tb;

    std::vector<std::pair<int, double>> result;

    COMM_TOOL::split_vector_pair(vec_p_i_f,
                                 serial_Ta,
                                 serial_Tb);

    COMM_TOOL::combine_vector_pair(serial_Ta,
                                   serial_Tb,
                                   result    );

    EXPECT_EQ(result, vec_p_i_f);
}

TEST_F(SerDes, VecStr){
    std::vector<char>        serial_data;
    std::vector<std::string> result;

    COMM_TOOL::serialize_vector_string(vec_str,
                                       serial_data);

    COMM_TOOL::deserialize_vector_string(serial_data,
                                         result      );

    EXPECT_EQ(result, vec_str);
}

TEST_F(SerDes, Str){
    std::vector<char> serial_data;
    std::string       result;

    for(const auto& str : vec_str){
        COMM_TOOL::serialize_string(str,
                                    serial_data);

        COMM_TOOL::deserialize_string(serial_data,
                                      result      );

        EXPECT_EQ(result, str);
    }
}


#include "gtest_main.hpp"
