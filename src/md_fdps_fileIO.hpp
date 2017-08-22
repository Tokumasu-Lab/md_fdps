//***************************************************************************************
//  This is record file I/O routine.
//    This code is used by "md_fdps_main.cpp"
//***************************************************************************************
#pragma once

#include <climits>
#include <cstdio>
#include <string>
#include <fstream>
#include <sys/stat.h>

#include <particle_simulator.hpp>

#include "md_fdps_atom_class.hpp"


//--- file I/O mode
enum class IO_MODE : int {
    pos,
    resume,
    VMD,
};

//--- std::string converter for enum
namespace ENUM {

    std::string whatis(const IO_MODE &e){
        switch (e) {
            case IO_MODE::pos:
                return "pos";
            break;

            case IO_MODE::resume:
                return "resume";
            break;

            case IO_MODE::VMD:
                return "VMD";
            break;

            default:
                throw std::out_of_range("undefined enum value.");
        }
    }

    IO_MODE which_IO_MODE(const std::string &str){
        if(        str == "pos" ){
            return IO_MODE::pos;
        } else if( str == "resume" ){
            return IO_MODE::resume;
        } else if( str == "VMD" ){
            return IO_MODE::VMD;
        } else {
            std::cerr << "  IO_MODE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in IO_MODE.");
        }
    }
}

inline std::ostream& operator << (std::ostream& s, const IO_MODE &e){
    s << ENUM::whatis(e);
    return s;
}

namespace FILE_IO {

    //--- prototype
    void makeOutputDirectory(char * dir_name);

    template <class Tset>
    void Init(const Tset &setting);

    template<class Tpsys>
    void recordPos(Tpsys & system,
                   const PS::S64 i_step);

    template<class Tpsys>
    void recordResume(Tpsys & system,
                      const PS::S64 i_step);

    template<class Tpsys>
    void recordVMD(Tpsys & system,
                   const PS::S64 i_step);




    //--- I/O mode
    static IO_MODE io_mode = IO_MODE::pos;

    //--- output cycle setting
    static PS::S32 pos_interval    = std::numeric_limits<PS::S32>::max();
    static PS::S32 pos_start       = std::numeric_limits<PS::S32>::max();
    static PS::S32 resume_interval = std::numeric_limits<PS::S32>::max();
    static PS::S32 resume_start    = std::numeric_limits<PS::S32>::max();
    static PS::S32 pdb_interval    = std::numeric_limits<PS::S32>::max();
    static PS::S32 pdb_start       = std::numeric_limits<PS::S32>::max();

    //--- output directry setting
    char dir_name_pos[1024];
    char dir_name_resume[1024];
    char dir_name_pdb[1024];

    //--- header class for file I/O
    class FileHeader {
      public:
        PS::S64 n_atom;
        PS::S64 i_step;
        char    buf[512];

        PS::S64 readAscii(FILE *fp){
            std::fscanf(fp, "%lld\n", &i_step);
            std::fscanf(fp, "%lld\n", &n_atom);
            std::fgets(buf, sizeof(buf), fp);
            return n_atom;
        }
        void writeAscii(FILE *fp) const {
            std::fprintf(fp, "%lld\n", i_step);
            std::fprintf(fp, "%lld\n", n_atom);
            std::fprintf(fp, buf);
        }
    };

    //--- directory manager
    void makeOutputDirectory(char *dir_name){
        struct stat st;
        if(stat(dir_name, &st) != 0) {
            PS::S32 ret_loc = 0;
            PS::S32 ret     = 0;
            if(PS::Comm::getRank() == 0){
                ret_loc = mkdir(dir_name, 0777);
            }
            PS::Comm::broadcast(&ret_loc, ret);
            if(ret == 0) {
                if(PS::Comm::getRank() == 0) std::fprintf(stderr, "Directory \"%s\" is successfully made.\n", dir_name);
            } else {
                std::fprintf(stderr, "Directory %s fails to be made.\n", dir_name);
                PS::Abort();
            }
        }
    }

    //--- initialize
    template <class Tsetting>
    void Init(const Tsetting &setting){
        std::sprintf(dir_name_pos,    "./posdata");
        std::sprintf(dir_name_resume, "./resume");
        std::sprintf(dir_name_pdb,    "./pdb");
        makeOutputDirectory(dir_name_pos);
        makeOutputDirectory(dir_name_resume);
        makeOutputDirectory(dir_name_pdb);

        pos_interval    = setting.pos_interval;
        pos_start       = setting.pos_start;
        resume_interval = setting.resume_interval;
        resume_start    = setting.resume_start;
        pdb_interval    = setting.pdb_interval;
        pdb_start       = setting.pdb_start;
    }

    //--- pos data output
    template<class Tpsys>
    void recordPos(Tpsys & system,
                   const PS::S64 i_step){

        //--- output cycle
        if(  i_step < pos_start ||
            (i_step % pos_interval) != 0 ) return;

        if(PS::Comm::getRank() == 0) std::cerr << "  output pos " << i_step << std::endl;

        char filename[256];
        std::sprintf(filename, "%s/pos%010lld.dat", dir_name_pos, i_step);
        FileHeader header;
        header.i_step = i_step;
        header.n_atom = system.getNumberOfParticleGlobal();
        std::sprintf(header.buf, "AtomId\tMolID\tAtomType\tMolType\tpos_x\tpos_y\tpos_z\n");

        //--- write data
        io_mode = IO_MODE::pos;
        system.writeParticleAscii(filename, header);
    }

    //--- resume data output
    template<class Tpsys>
    void recordResume(Tpsys & system,
                      const PS::S64 i_step){

        //--- output cycle
        if(  i_step < resume_start ||
            (i_step % resume_interval) != 0 ) return;

        if(PS::Comm::getRank() == 0) std::cerr << "  output resume " << i_step << std::endl;

        //--- common file name
        char filename[256];
      //  std::sprintf(filename, "%s/resume%010d", dir_name_resume, i_step);
        std::sprintf(filename, "%s/resume%010lld.dat", dir_name_resume, i_step);
        FileHeader header;
        header.i_step = i_step;
        header.n_atom = system.getNumberOfParticleLocal();
        std::sprintf(header.buf, "under developping now");

      //  //--- distributed file format
      //  char format[256];
      //  std::sprintf(format, "$s_%03d_%03d.dat");  // filename, total MPI Rank, process MPI Rank.

        //--- write data
        io_mode = IO_MODE::resume;
      //  system.writeParticleAscii(filename, format);
        system.writeParticleAscii(filename);
    }

    //--- VMD data output
    template<class Tpsys>
    void recordVMD(Tpsys & system,
                   const PS::S64 i_step){

        //--- output cycle
        if(  i_step < pdb_start ||
            (i_step % pdb_interval) != 0 ) return;

        if(PS::Comm::getRank() == 0) std::cerr << "  output pdb " << i_step << std::endl;

        char filename[256];
        std::sprintf(filename, "%s/vmd%010lld.pdb", dir_name_pdb, i_step);

        //--- write data
        io_mode = IO_MODE::VMD;
        system.writeParticleAscii(filename);
    }
}




//--- file I/O interface in FP class through FDPS library
void Atom_FP::writeAscii(FILE* fp) const {
    //--- value buffer
    PS::F64vec pos_real = Normalize::realPos( this->getPos() );

    std::string residue = ENUM::whatis(this->getMolType());
                residue = residue.substr(3);
    if(residue.size() < 3) { for(size_t i=0; i<3-residue.size(); ++i) residue += " "; }  // add space

    switch (FILE_IO::io_mode){
        case IO_MODE::pos:
            //--- position & ID data output
            std::fprintf(fp, "%8lld\t%8lld\t%8s\t%8s\t%8.3e\t%8.3e\t%8.3e\n",
                         this->getAtomID(),
                         this->getMolID(),
                         ENUM::whatis(this->getAtomType()).c_str(),
                         ENUM::whatis(this->getMolType()).c_str(),
                         pos_real.x,
                         pos_real.y,
                         pos_real.z);
        break;

        case IO_MODE::resume:
            //--- resume data output
        break;

        case IO_MODE::VMD:
            //--- VMD output
            char buff[128];
            char add_buff[32];

            std::sprintf(buff, "ATOM  ");                                               //  1~6  "ATOM  "
            std::sprintf(add_buff, "%5lld", this->getAtomID() );                        //  7~11 atom id
            std::strcat(buff, add_buff );
            std::strcat(buff, " ");                                                     // 12    space
            std::sprintf(add_buff, "%4s", ENUM::whatis(this->getAtomType()).c_str() );  // 13~16 atom name
            std::strcat(buff, add_buff);
            std::strcat(buff, " ");                                                     // 17    alternate location indicator
            std::sprintf(add_buff, "%.3s", residue.c_str() );                           // 18~20 residue name
            std::strcat(buff, add_buff);
            std::strcat(buff, " ");                                                     // 21    chain identifier
            std::sprintf(add_buff, "%4lld", this->getMolID() );                         // 22~25 residue sequence number
            std::strcat(buff, add_buff );
            std::strcat(buff, " ");                                                     // 26    code for insertion of residues
            std::strcat(buff, "   ");                                                   // 27~29 space

            std::sprintf(add_buff, " %7.3f", pos_real.x );                              // 30~37 X position [angstrom]
            std::strcat(buff, add_buff );
            std::sprintf(add_buff, " %7.3f", pos_real.y );                              // 38~45 Y position [angstrom]
            std::strcat(buff, add_buff );
            std::sprintf(add_buff, " %7.3f", pos_real.z );                              // 46~53 Z position [angstrom]
            std::strcat(buff, add_buff );

            std::strcat(buff, "\n");    // end of line

            std::fprintf(fp, buff);
        break;

        default:
            throw std::invalid_argument("undefined I/O mode.");
    }
}

void Atom_FP::readAscii(FILE* fp){

    //--- pos & ID data input  (for analysis program)
    switch (FILE_IO::io_mode){
        case IO_MODE::pos:
            //--- pos data input
        break;

        case IO_MODE::resume:
            //--- resume data input
        break;

        default:
            throw std::invalid_argument("undefined I/O mode.");
    }
}
