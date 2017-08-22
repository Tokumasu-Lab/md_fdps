//***************************************************************************************
//  This is unit test of "tple_array.hpp".
//***************************************************************************************

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>

#include <particle_simulator.hpp>

#include "comm_tool.hpp"


int main(int argc, char *argv[]) {
    PS::Initialize(argc, argv);

    //--- test data
    //------ for basic pattern
    std::vector<int>              vec_int;
    std::vector<std::vector<int>> vec_vec_int;
    std::vector<std::string>      vec_str;
    std::vector<std::pair<int, float>>  vec_pair_i_f;
    std::unordered_map<int, float>      map_i_f;
    std::unordered_multimap<int, float> m_map_i_f;

    //------ for recursive pattern sample
    std::vector<std::vector<std::vector<int>>> vec_vec_vec_int;
    std::vector<std::vector<std::string>>      vec_vec_str;
    std::vector<std::pair<std::string, std::vector<int>>>  vec_pair_s_vi;
    std::unordered_map<std::string, std::vector<int>>      map_s_vi;
    std::unordered_multimap<std::string, std::vector<int>> m_map_s_vi;

    if(PS::Comm::getRank() == 0) {
        std::cout << "TEST: COMM_TOOL wrapper for PS::Comm::broadcast()." << std::endl;
        fprintf(stderr, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stderr, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());

        //--- initialize
        std::cout << std::endl;
        std::cout << "--- settings test value in the process " << PS::Comm::getRank() << std::endl;
        std::cout << std::endl;

        std::cout << "  @ defined template functions:" << std::endl;
        std::cout << "   --- std::vector<T>"           << std::endl;
        for(int i=1; i<=10; ++i) { vec_int.emplace_back(i*2); }

        std::cout << "   --- std::vector<std::vector<T>>" << std::endl;
        for(int i=1; i<=3; ++i) {
            std::vector<int> v_tmp;
            v_tmp.clear();
            for(int j=1; j<=5; ++j){
                v_tmp.push_back(i*100 + j);
            }
            vec_vec_int.push_back(v_tmp);
        }

        std::cout << "   --- std::vector<std::string>"    << std::endl;
        vec_str.push_back("this");
        vec_str.push_back("is");
        vec_str.push_back("vector");
        vec_str.push_back("of");
        vec_str.push_back("string.");

        std::cout << "   --- std::vector<std::pair<Ta, Tb>>" << std::endl;
        for(int i=0; i<=4; ++i) { vec_pair_i_f.push_back(std::make_pair(i*10, float(i)*0.1)); }

        std::cout << "   --- std::unordered_map<Ta, Tb>"      << std::endl;
        for(int i=0; i<=4; ++i) { map_i_f[i*10] = -float(i)*0.1; }

        std::cout << "   --- std::unordered_multimap<Ta, Tb>" << std::endl;
        for(int i=0; i<=4; ++i) { m_map_i_f.insert( std::make_pair(i*10, -float(i)*0.1) ); }
        for(int i=2; i<=6; ++i) { m_map_i_f.insert( std::make_pair(i*10,  float(i)*0.1) ); }


        std::cout << "  @ recursive pattern sample:" << std::endl;
        std::cout << "   --- std::vector<std::vector<std::vector<T>>>" << std::endl;
        for(int i=1; i<=4; ++i){
            std::vector<std::vector<int>> tmp_vvi;
            for(int j=1; j<=i; ++j){
                std::vector<int> tmp_vi;
                for(int k=1; k<=j; ++k){
                    tmp_vi.push_back(k);
                }
                tmp_vvi.push_back(tmp_vi);
            }
            vec_vec_vec_int.push_back(tmp_vvi);
        }

        std::cout << "   --- std::vector<std::vector<std::string>>" << std::endl;
        std::vector<std::string> tmp_v_str;
        tmp_v_str.push_back("this");
        tmp_v_str.push_back("is");
        tmp_v_str.push_back("vector");
        tmp_v_str.push_back("of");
        tmp_v_str.push_back("vector");
        tmp_v_str.push_back("of");
        tmp_v_str.push_back("string.");

        for(size_t i=1; i<=tmp_v_str.size(); ++i){
            std::vector<std::string> tmp_v;
            for(size_t j=0; j<i; ++j){
                tmp_v.push_back(tmp_v_str.at(j));
            }
            vec_vec_str.push_back(tmp_v);
        }

        std::cout << "   --- std::vector<std::pair<Ta, Tb>>"        << std::endl;
        std::cout << "                Ta= std::string,  Tb= std::vector<T>" << std::endl;
        tmp_v_str.clear();
        tmp_v_str.push_back("this");
        tmp_v_str.push_back("is");
        tmp_v_str.push_back("vector");
        tmp_v_str.push_back("of");
        tmp_v_str.push_back("pair");
        tmp_v_str.push_back("of");
        tmp_v_str.push_back("string & vector<int>.");
        for(auto str : tmp_v_str){
            std::vector<int> v_tmp;
            for(size_t i=0; i<str.size(); ++i){
                v_tmp.push_back(static_cast<int>(str[i]));
            }
            vec_pair_s_vi.push_back(std::make_pair(str, v_tmp));
        }

        std::cout << "   --- std::unordered_map<Ta, Tb>"      << std::endl;
        std::cout << "                Ta= std::string,  Tb= std::vector<T>" << std::endl;
        for(auto str : tmp_v_str){
            std::vector<int> v_tmp;
            for(size_t i=0; i<str.size(); ++i){
                v_tmp.push_back(static_cast<int>(str[i]));
            }
            map_s_vi[str] = v_tmp;
        }

        std::cout << "   --- std::unordered_multimap<Ta, Tb>" << std::endl;
        std::cout << "                Ta= std::string,  Tb= std::vector<T>" << std::endl;
        for(auto str : tmp_v_str){
            std::vector<int> v_tmp;
            for(size_t i=0; i<str.size(); ++i){
                v_tmp.push_back(static_cast<int>(str[i]));
            }
            m_map_s_vi.insert(std::make_pair(str, v_tmp));
        }
    }

    COMM_TOOL::broadcast(vec_int, 0);
    COMM_TOOL::broadcast(vec_vec_int, 0);
    COMM_TOOL::broadcast(vec_str, 0);
    COMM_TOOL::broadcast(vec_pair_i_f, 0);
    COMM_TOOL::broadcast(map_i_f, 0);
    COMM_TOOL::broadcast(m_map_i_f, 0);

    COMM_TOOL::broadcast(vec_vec_vec_int, 0);
    COMM_TOOL::broadcast(vec_vec_str, 0);
    COMM_TOOL::broadcast(vec_pair_s_vi, 0);
    COMM_TOOL::broadcast(map_s_vi, 0);
    COMM_TOOL::broadcast(m_map_s_vi, 0);
    std::cout << "--- sync in MPI broadcast ---" << std::endl;

    if(PS::Comm::getRank() == 1) {

        //--- result
        std::cout << std::endl;
        std::cout << "--- display test value in the process " << PS::Comm::getRank() << std::endl;
        std::cout << std::endl;

        std::cout << "  @ defined template functions:" << std::endl;
        std::cout << "   --- std::vector<T>"           << std::endl;
        std::cout << "       [";
        for(auto e : vec_int) { std::cout << " " << e; }
        std::cout << " ]" << std::endl;
        std::cout << std::endl;

        std::cout << "   --- std::vector<std::vector<T>>" << std::endl;
        size_t count = 0;
        std::cout << "      v1 len = " << vec_vec_int.size() << std::endl;
        for(auto v : vec_vec_int){
            std::cout << "       v1_index = " << count << " [";
            for(auto e : v) { std::cout << " " << e; }
            std::cout << " ]" << std::endl;
            ++count;
        }
        std::cout << std::endl;

        std::cout << "   --- std::vector<std::string>"    << std::endl;
        std::cout << "       [";
        for(auto e : vec_str) { std::cout << " " << e; }
        std::cout << " ]" << std::endl;
        std::cout << std::endl;

        std::cout << "   --- std::vector<std::pair<Ta, Tb>>" << std::endl;
        std::cout << "       [";
        for(auto e : vec_pair_i_f) { std::cout << " (" << e.first << ", " << e.second << ")"; }
        std::cout << " ]" << std::endl;
        std::cout << std::endl;

        std::cout << "   --- std::unordered_map<Ta, Tb>"      << std::endl;
        std::cout << "       [";
        for(auto e : map_i_f) { std::cout << " (" << e.first << ", " << e.second << ")"; }
        std::cout << " ]" << std::endl;
        std::cout << std::endl;

        std::cout << "   --- std::unordered_multimap<Ta, Tb>" << std::endl;
        std::cout << "       [";
        for(auto e : m_map_i_f) { std::cout << " (" << e.first << ", " << e.second << ")"; }
        std::cout << " ]" << std::endl;
        std::cout << std::endl;



        std::cout << "  @ recursive pattern sample:" << std::endl;
        std::cout << "   --- std::vector<std::vector<std::vector<T>>>" << std::endl;
        size_t c_v1 = 0;
        for(auto v1 : vec_vec_vec_int){
            std::cout << "      v1_index = " << c_v1 << " / " << vec_vec_vec_int.size() << std::endl;
            size_t c_v2 = 0;
            for(auto v2 : v1) {
                std::cout << "          v2_index = " << c_v2 << " [";
                for(auto e : v2) { std::cout << " " << e; }
                ++c_v2;
                std::cout << " ]" << std::endl;
            }
            ++c_v1;
        }
        std::cout << std::endl;

        std::cout << "   --- std::vector<std::vector<std::string>>" << std::endl;
        count = 0;
        std::cout << "      v1 len = " << vec_vec_str.size() << std::endl;
        for(auto v : vec_vec_str){
            std::cout << "       v1_index = " << count << " [";
            for(auto e : v) { std::cout << " " << e; }
            std::cout << " ]" << std::endl;
            ++count;
        }
        std::cout << std::endl;

        std::cout << "   --- std::vector<std::pair<Ta, Tb>>"        << std::endl;
        std::cout << "                Ta= std::string,  Tb= std::vector<T>" << std::endl;
        for(auto ve : vec_pair_s_vi) {
            std::cout << "     (" << ve.first << ": ";
            for(auto e : ve.second) { std::cout << " " << e; }
            std::cout << ")" << std::endl;
        }
        std::cout << std::endl;

        std::cout << "   --- std::unordered_map<Ta, Tb>"      << std::endl;
        std::cout << "                Ta= std::string,  Tb= std::vector<T>" << std::endl;
        for(auto ve : map_s_vi) {
            std::cout << "     (" << ve.first << ": ";
            for(auto e : ve.second) { std::cout << " " << e; }
            std::cout << ")" << std::endl;
        }
        std::cout << std::endl;

        std::cout << "   --- std::unordered_multimap<Ta, Tb>" << std::endl;
        std::cout << "                Ta= std::string,  Tb= std::vector<T>" << std::endl;
        for(auto ve : m_map_s_vi) {
            std::cout << "     (" << ve.first << ": ";
            for(auto e : ve.second) { std::cout << " " << e; }
            std::cout << ")" << std::endl;
        }
        std::cout << std::endl;


        std::cout << "    the test succeeded." << std::endl;
    }

    PS::Finalize();
    return 0;
}
