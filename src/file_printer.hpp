//***************************************************************************************
//  This is observer functions.
//***************************************************************************************
#pragma once

#include <string>
#include <fstream>
#include <iomanip>

#include <particle_simulator.hpp>


namespace Observer {

    class FilePrinter {
    private:
        PS::S32       rank = 0;
        std::string   name;
        std::ofstream ofs;

    public:
        void file_init(const std::string &name, const PS::S32 &rank){
            if( this->ofs.is_open() ){ this->ofs.close(); }
            this->rank = rank;
            this->name = name;

            if(PS::Comm::getRank() != this->rank) return;
            this->ofs.open(name, std::ios::trunc);

            if(!ofs){
                std::ostringstream oss;
                oss << " failed to open the file: " << this->name << " \n";
                throw std::ios_base::failure(oss.str());
            }
        }

        void print(const std::string &str){
            if(PS::Comm::getRank() != this->rank) return;
            this->ofs << str;
            this->ofs.flush();   // flush write buffer at every call
        }

        std::string file_name() const { return this->name; }
        PS::S32     getRank()   const { return this->rank; }
    };

}
