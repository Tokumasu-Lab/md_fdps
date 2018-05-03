//=======================================================================================
//  This is unit test of MD_EXT::logspace_array<>.
//     module location: ./generic_ext/logspace_array.hpp
//=======================================================================================

#undef NDEBUG

#include <cmath>
#include <random>

#include <gtest/gtest.h>

#include <particle_simulator.hpp>

#include "logspace_array.hpp"

#include "str_tool.hpp"


namespace TEST_DEFS {
    const std::string logspace_lin_file = "./test_bin/logspace_count_lin.tsv";
    const std::string logspace_lin_ref  = "./unit_test/ref/logspace_count_lin.tsv";

    const std::string logspace_exp_file = "./test_bin/logspace_count_exp.tsv";
    const std::string logspace_exp_ref  = "./unit_test/ref/logspace_count_exp.tsv";
}

//--- unit test definition, CANNOT use "_" in test/test_case name.
TEST(LogspaceArray, init){
    MD_EXT::logspace_array<int> arr;

    //--- not initialized
    EXPECT_THROW(arr.fill(1)         , std::logic_error);
    EXPECT_THROW(arr.resize(100.0)   , std::logic_error);
    EXPECT_THROW(arr.resize_array(10), std::logic_error);
    EXPECT_THROW(arr.at(0)           , std::logic_error);

    EXPECT_TRUE(arr.empty());

    //--- invalid initialize parameter
    EXPECT_THROW(arr.init(-1.0, 1.0, 6), std::invalid_argument) << "negative value";
    EXPECT_THROW(arr.init( 1.0, 1.0, 6), std::invalid_argument) << "no range";
    EXPECT_THROW(arr.init( 1.0, 2.0, 0), std::invalid_argument) << "not enough array size";

    arr.init(1.0, 10.0, 6);

    EXPECT_EQ(arr.size(), 6);
    auto range = arr.range();
    EXPECT_FLOAT_EQ(range.first , 1.0);
    EXPECT_FLOAT_EQ(range.second, 10.0);
    EXPECT_FLOAT_EQ(arr.ratio() , std::pow(10.0/1.0, 1.0/6.0));

    //--- assign
    MD_EXT::logspace_array<int> arr_2;
    arr_2 = arr;

    EXPECT_EQ(arr_2.size(), 6);
    EXPECT_FLOAT_EQ(arr_2.range().first , 1.0);
    EXPECT_FLOAT_EQ(arr_2.range().second, 10.0);
    EXPECT_FLOAT_EQ(arr.ratio() , std::pow(10.0/1.0, 1.0/6.0));
}

TEST(LogspaceArray, access){
    MD_EXT::logspace_array<int> arr{1.0, 10.0, 6};

    //--- iterator access
    arr.fill(11);
    for(const auto& elem : arr){
        EXPECT_EQ(elem.second, 11);
    }
    arr.fill(1);
    for(const auto& elem : arr){
        EXPECT_EQ(elem.second, 1);
    }

    //--- logspace distance
    const double log_begin  = std::log(1.0);
    const double log_end    = std::log(10.0);
    const double log_ratio  = (log_end - log_begin)/6.0;
    const double real_ratio = std::exp(log_ratio);

    EXPECT_FLOAT_EQ(arr.ratio(), real_ratio);

    size_t count = 0;
    for(const auto& elem : arr){
        EXPECT_FLOAT_EQ(elem.first, std::exp(log_begin + double(count)*log_ratio));
        ++count;
    }

    //--- logspace range
    for(double f=1.0; f<10.0; f += 0.99){
        const auto range = arr.range(f);
        EXPECT_FLOAT_EQ(range.second, range.first*real_ratio);
    }

    EXPECT_THROW(arr.at(0.9) , std::out_of_range);
    EXPECT_THROW(arr.at(11.0), std::out_of_range);

    EXPECT_NO_THROW(arr.at(1.0));
    EXPECT_NO_THROW(arr.at(9.9));
}

TEST(LogspaceArray, resize){
    MD_EXT::logspace_array<int> arr{1.0, 10.0, 6};

    const double ratio    = arr.ratio();
    const size_t n_backet = arr.size();

    //--- resize by array size
    arr.resize_array(n_backet*2);
    EXPECT_FLOAT_EQ(arr.ratio(), ratio);
    EXPECT_FLOAT_EQ(arr.range().first , 1.0);
    EXPECT_FLOAT_EQ(arr.range().second, arr.range().first*std::pow(ratio, n_backet*2));
    EXPECT_FLOAT_EQ(arr.range().second, 100.0);

    for(auto itr = arr.begin(); itr != arr.end()-1; ++itr){
        const auto i_ratio = ((itr+1)->first)/((itr)->first);
        EXPECT_FLOAT_EQ(i_ratio, ratio);
    }

    //--- resize by max value
    const double eps = 1.e-6;
    arr.resize(1000.0 - eps);
    EXPECT_FLOAT_EQ(arr.ratio(), ratio);
    EXPECT_FLOAT_EQ(arr.range().first , 1.0);
    EXPECT_FLOAT_EQ(arr.range().second, 1000.0);
    EXPECT_EQ(arr.size(), 18);  //  max_range = 10.0^3, size = n_backet*3.

    for(auto itr = arr.begin(); itr != arr.end()-1; ++itr){
        const auto i_ratio = ((itr+1)->first)/((itr)->first);
        EXPECT_FLOAT_EQ(i_ratio, ratio);
    }

    //--- init as different array
    arr.init(10.0, 100.0, 10);
    EXPECT_FLOAT_EQ(arr.range().first , 10.0);
    EXPECT_FLOAT_EQ(arr.range().second, 100.0);
    EXPECT_EQ(arr.size(), 10);

    //--- throw
    EXPECT_THROW(arr.resize(0.5), std::length_error);
}

class TestData {
public:
    double x;
    int    y;

    bool read_line(const std::string &line_in){
        std::string line = line_in;
        STR_TOOL::removeCR(line);
        const auto str_list = STR_TOOL::split(line, "\t");

        if(str_list.size() < 2) return false;

        this->x = std::stod(str_list[0]);
        this->y = std::stoi(str_list[1]);

        return true;
    }
};
std::ostream& operator << (std::ostream& s, const TestData &d){
    s << std::hexfloat;
    s << d.x << "\t"
      << d.y << "\n";
    s << std::defaultfloat;

    return s;
}

template <class Tdata>
void write_log(const std::string        &file_name,
               const std::vector<Tdata> &data_list){

    std::ofstream file{file_name};
    if(!file.is_open()){
        throw std::ios_base::failure("failed to open the file: " + file_name);
    }

    for(const auto& data : data_list){
        file << data;
    }
    file.close();
}
template <class Tdata>
void read_log(const std::string        &file_name,
                    std::vector<Tdata> &data_list){

     data_list.clear();

     std::ifstream ifs{file_name};
     if(ifs.fail()){
         throw std::ios_base::failure("failed to open the file: " + file_name);
     }

     std::string line;
     while( getline(ifs, line) ){
         Tdata tmp;
         bool flag = tmp.read_line(line);
         if(flag){
             data_list.push_back(tmp);
         }
     }
}

TEST(LogspaceArray, countingLin){
    const double exp_range = 20.0;
    const double range     = std::exp(exp_range);

    const size_t n_sample = 100000;

    std::mt19937_64 mt;

    std::uniform_real_distribution<double> dist_lin(1.0, range);

    MD_EXT::logspace_array<int> arr_lin;

    arr_lin.init(1.0, range, 40);

    arr_lin.fill(0);
    mt.seed(987654321);
    for(size_t i=0; i<n_sample; ++i){
        const auto f = dist_lin(mt);
        ++(arr_lin.at(f));
    }

    std::vector<TestData> result;
    for(const auto& elem : arr_lin){
        result.push_back( TestData{elem.first, elem.second} );
    }

    write_log(TEST_DEFS::logspace_lin_file, result);

    std::vector<TestData> ref;
    read_log(TEST_DEFS::logspace_lin_ref, ref);

    EXPECT_EQ(result.size(), ref.size());
    for(size_t i=0; i<result.size(); ++i){
        EXPECT_FLOAT_EQ(result.at(i).x, ref.at(i).x);
        EXPECT_EQ(      result.at(i).y, ref.at(i).y);
    }
}

TEST(LogspaceArray, countingExp){
    const double exp_range = 20.0;
    const double range     = std::exp(exp_range);

    const size_t n_sample = 100000;

    std::mt19937_64 mt;

    std::uniform_real_distribution<double> dist_exp(0.0, exp_range);

    MD_EXT::logspace_array<int> arr_exp;

    arr_exp.init(1.0, range, 40);

    arr_exp.fill(0);
    mt.seed(987654321);
    for(size_t i=0; i<n_sample; ++i){
        const auto f = std::exp(dist_exp(mt));
        ++(arr_exp.at(f));
    }

    std::vector<TestData> result;
    for(const auto& elem : arr_exp){
        result.push_back( TestData{elem.first, elem.second} );
    }

    write_log(TEST_DEFS::logspace_exp_file, result);

    std::vector<TestData> ref;
    read_log(TEST_DEFS::logspace_exp_ref, ref);

    EXPECT_EQ(result.size(), ref.size());
    for(size_t i=0; i<result.size(); ++i){
        EXPECT_FLOAT_EQ(result.at(i).x, ref.at(i).x);
        EXPECT_EQ(      result.at(i).y, ref.at(i).y);
    }
}


#include "gtest_main.hpp"
