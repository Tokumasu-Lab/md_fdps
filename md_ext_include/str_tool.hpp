//***************************************************************************************
//  This program is loading model parameter function.
//***************************************************************************************
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>

namespace STR_TOOL {

    inline void removeCR(std::string &str){
        if(!str.empty() && str.back() == static_cast<char>(13)){
            str.pop_back();
        }
    }

    //--- split function for std::string
    inline std::vector<std::string> split(const std::string &str,
                                          const std::string &delim){

        assert(delim.size() == 1);  // use const std::string instead of const char
        std::istringstream ss{str};

        std::string item;
        std::vector<std::string> result;
        while (std::getline(ss, item, delim.at(0))) {
            if(!item.empty()){
                result.push_back(item);
            }
        }
        return result;
    }

    inline bool isInteger(const std::string &str){
        if(str.find_first_not_of("-0123456789 \t") != std::string::npos) {
            return false;
        } else {
            return true;
        }
    }

    inline bool isNumeric(const std::string &str){
        if(str.find_first_not_of("-0123456789.Ee \t") != std::string::npos) {
            return false;
        } else {
            return true;
        }
    }
}
