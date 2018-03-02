//=======================================================================================
//  This is common routine for unit test based on the googte test.
//=======================================================================================

#include <cassert>

::testing::AssertionResult float_relative_eq(const double &lhs,
                                             const double &rhs,
                                             const double &abs_err,
                                             const double &relative_err){

    assert(abs_err      >= 0.0);
    assert(relative_err >= 0.0);

    std::ostringstream oss;
    oss << "\n"
        << "  lhs = " << lhs << "\n"
        << "  rhs = " << rhs << "\n";

    const double diff = lhs - rhs;
    if(std::abs(diff) > abs_err){
        oss << "    diff = " << diff << " > " << "absolute error = " << abs_err << "\n";
    } else {
        return ::testing::AssertionSuccess();
    }

    double r_diff = 0.0;
    if(lhs == 0.0){
        r_diff = diff/rhs;
    } else {
        r_diff = diff/lhs;
    }

    if(std::abs(r_diff) > relative_err ){
        oss << "    relative diff =" << r_diff << " > " << "relative error = " << relative_err << "\n";
    } else {
        return ::testing::AssertionSuccess();
    }

    return ::testing::AssertionFailure() << oss.str();
}

::testing::AssertionResult float_relative_eq(const double &lhs,
                                             const double &rhs,
                                             const double &relative_err){
    return float_relative_eq(lhs, rhs, relative_err, relative_err);
}

template <class Tlog>
void write_log_file(const std::string &file_name,
                    const Tlog        &data_list){
    if(PS::Comm::getRank() != 0) return;

    std::ofstream file{file_name};
    for(const auto& data : data_list){
        file << data;
    }
    file.close();
}

template <class Tlog>
void load_log_file(const std::string &file_name,
                         Tlog        &data_list){
    data_list.clear();

    FS_TOOL::file_load(file_name, data_list, 0);
    COMM_TOOL::broadcast(data_list, 0);
}
