/**************************************************************************************************/
/**
* @file  str_tool.hpp
* @brief tools for reading parameter files.
*/
/**************************************************************************************************/
#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <stdexcept>


/**
* @brief tools for reading parameter files.
*/
namespace STR_TOOL {

    /**
    * @brief remove "CR" code from back of std::string.
    * @param[in, out] str std::string value after devided by "LF" code.
    * @details countermeasure for the difference of "\n" in Linux (LF) and Windows (CR+LF).
    */
    inline void removeCR(std::string &str){
        if(!str.empty() && str.back() == static_cast<char>(13)){
            str.pop_back();
        }
    }

    /**
    * @brief split function for std::string.
    * @param[in] str std::string value.
    * @param[in] delim delimiter.
    * @return std::vector<std::string> result.
    */
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
    inline std::vector<std::string> split(const std::string &str,
                                          const char        *delim){
        return split(str, std::string{delim});
    }
    inline std::vector<std::string> split(const std::string &str,
                                          const char         delim){
        return split(str, std::string{delim});
    }

    /**
    * @brief checking str is integer of not.
    * @param[in] str std::string value.
    * @return bool result. True: str is integer. False: str is not integer.
    */
    inline bool isInteger(const std::string &str){
        if(str.find_first_not_of("-+0123456789 \t") != std::string::npos) {
            return false;
        } else {
            return true;
        }
    }

    /**
    * @brief checking str is float of not.
    * @param[in] str std::string value.
    * @return bool result. True: str is float. False: str is not float.
    */
    inline bool isNumeric(const std::string &str){
        if(str.find_first_not_of("-+0123456789.Ee \t") != std::string::npos) {
            return false;
        } else {
            return true;
        }
    }

    /**
    * @brief justify plan indicator for "std::string shape_str_vec_2D()".
    */
    enum class STR_POS {
        left,
        right,
        internal,
        center,
    };

    namespace _Impl {
        void str_center(      std::ostream &s,
                        const std::string  &str,
                        const size_t        w   ){

            const size_t space   = w - str.size();
            const size_t s_back  = space/2;
            const size_t s_front = space - s_back;

            s << std::left;

            if(s_front > 0){
                s << std::setw(s_front) << " ";
            }
            s << str;
            if(s_back > 0){
                s << std::setw(s_back) << " ";
            }
        }
    }

    using StringList   =             std::vector<std::pair<std::string, STR_POS>>;
    using StringList2D = std::vector<std::vector<std::pair<std::string, STR_POS>>>;

    /**
    * @brief justifing string as 2D grid
    * @detail [line [word (string, pos)]]
    */
    std::string shape_str_vec_2d(const StringList2D &str_list_2D){

        //--- get word_length
        std::vector<size_t> w_len_list;
        for(const auto& line : str_list_2D){
            //--- allocate list
            if(line.size() > w_len_list.size()){
                size_t diff = line.size()-w_len_list.size();
                for(size_t i=0; i<diff; ++i){
                    w_len_list.push_back(0);
                }
            }

            //--- get length
            for(size_t i=0; i<line.size(); ++i){
                w_len_list[i] = std::max(w_len_list[i], line[i].first.size());
            }
        }

        //--- make output
        std::ostringstream oss;
        for(const auto& line : str_list_2D){
            for(size_t i=0; i<line.size(); ++i){
                oss << std::setw(w_len_list[i]);

                switch (line[i].second){
                    case STR_POS::left:
                        oss << std::setw(w_len_list[i]);
                        oss << std::left << line[i].first;
                    break;

                    case STR_POS::right:
                        oss << std::setw(w_len_list[i]);
                        oss << std::right << line[i].first;
                    break;

                    case STR_POS::internal:
                        oss << std::setw(w_len_list[i]);
                        oss << std::internal << line[i].first;
                    break;

                    case STR_POS::center:
                        _Impl::str_center(oss, line[i].first, w_len_list[i]);
                    break;

                    default:
                        throw std::invalid_argument("format error.");
                }
            }
            oss << "\n";
        }
        return oss.str();
    }
}
