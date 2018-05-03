//***************************************************************************************
/**
* @file  fstream_tool.hpp
* @brief simple file I/O tools.
*/
//***************************************************************************************
#pragma once

#include <string>
#include <fstream>
#include <iomanip>

#include <particle_simulator.hpp>
#include "comm_tool.hpp"

#include <sys/stat.h>

/*
*  @brief simple file I/O tools.
*/
namespace FS_TOOL {

    /*
    *  @brief directory maker.
    *  @param[in] dir_name target file.
    *  @param[in] rank     process rank. default value is 0.
    */
    void make_directory(const std::string &dir_name,
                        const PS::S32      rank = 0 ){

        COMM_TOOL::check_proc_rank(rank);
        if(PS::Comm::getRank() != rank) return;

        struct stat st;
        if(stat(dir_name.c_str(), &st) != 0) {
            int err = mkdir(dir_name.c_str(), 0777);

            if(err == -1){
                throw std::ios_base::failure("make the Directory: " + dir_name + " was failed.");
            } else {
                std::ostringstream oss;
                oss << "the Directory " << dir_name << " is successfully made." << "\n";
                std::cout << oss.str() << std::flush;
            }
        }

    }

    /*
    *  @brief simple file printer.
    */
    class FilePrinter {
    private:
        PS::S32       rank = 0;
        std::string   name;
        std::ofstream ofs;

    public:
        /*
         *  @param[in] name target file.
         *  @param[in] rank process rank. default value is 0.
        */
        void file_init(const std::string &file_name,
                       const PS::S32      rank = 0  ){

           COMM_TOOL::check_proc_rank(rank);

            if( this->ofs.is_open() ){ this->ofs.close(); }
            this->rank = rank;
            this->name = file_name;

            if(PS::Comm::getRank() != this->rank) return;

            this->ofs.open(file_name, std::ios::trunc);

            if(!ofs){
                throw std::ios_base::failure("failed to open the file: " + file_name);
            }
        }

        FilePrinter() = default;
        /*
         *  @param[in] name target file.
         *  @param[in] rank process rank. default value is 0.
        */
        FilePrinter(const std::string &name,
                    const PS::S32      rank = 0){
            this->file_init(name, rank);
        }
        ~FilePrinter() = default;

        void print(const std::string &str){
            if(PS::Comm::getRank() != this->rank) return;
            if( !this->ofs.is_open() ){
                throw std::ios_base::failure("the FilePrinter is not initialized.");
            }

            this->ofs << str;
            this->ofs.flush();   // flush write buffer at every call
        }

        template <class T>
        void operator << (const T &out){
            if(PS::Comm::getRank() != this->rank) return;
            if( !this->ofs.is_open() ){
                throw std::ios_base::failure("the FilePrinter is not initialized.");
            }

            this->ofs << out;
        }

        bool        is_open()   const { return this->ofs.is_open(); }
        std::string file_name() const { return this->name; }
        PS::S32     getRank()   const { return this->rank; }
    };


    /*
    *  @brief simple file loader for std::string.
    *  @param[in]  file_name target file.
    *  @param[out] line_list std::vector<> of line.
    *  @param[in]  rank      process rank. default value is 0.
    */
    void file_load(const std::string              &file_name,
                         std::vector<std::string> &line_list,
                   const PS::S32                   rank = 0  ){

        COMM_TOOL::check_proc_rank(rank);
        line_list.clear();

        if(PS::Comm::getRank() == rank){
            std::ifstream ifs{file_name};
            if(ifs.fail()){
                throw std::ios_base::failure("failed to open the file: " + file_name);
            }

            std::string line;
            while( getline(ifs, line) ){
                line_list.push_back(line);
            }
        }
    }

    /*
    *  @brief simple file loader for generic data type.
    *  @param[in]  file_name target file.
    *  @param[out] data_list std::vector<> of the data read from each lines.
    *  @param[in]  rank      process rank. default value is 0.
    *  @details the class of Tdata MUST have bool read_line(std::string &str) function
    *           to convert std::string of line into data.
    */
    template <class Tdata>
    void file_load(const std::string        &file_name,
                         std::vector<Tdata> &data_list,
                   const PS::S32             rank = 0  ){

       COMM_TOOL::check_proc_rank(rank);

        data_list.clear();

        if(PS::Comm::getRank() == rank){
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
    }

    /*
    *  @brief returning interface for simple file loader.
    *  @param[in]  file_name target file.
    *  @param[in]  rank      process rank. default value is 0.
    *  @return               std::vector<> of the data read from each lines.
    */
    template <class T = std::string>
    std::vector<T> file_load(const std::string &file_name,
                             const PS::S32      rank = 0  ){

        std::vector<T> buff;
        file_load(file_name, buff, rank);
        return buff;
    }

}
