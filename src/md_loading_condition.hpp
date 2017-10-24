//***************************************************************************************
//  This program is loading model parameter function.
//***************************************************************************************
#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <tuple>
#include <cassert>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "atom_class.hpp"
#include "md_coef_table.hpp"
#include "md_ext_sys_control.hpp"


//--- file loading mode
enum class CONDITION_LOAD_MODE : int {
    time_step,
    tree,
    cut_off,
    ext_sys,
    ext_sys_sequence,
    record,

    molecule,
    box,
    ex_radius,
};

//--- std::string converter for enum
namespace ENUM {

    static const std::map<std::string, CONDITION_LOAD_MODE> table_str_CONDITION_LOAD_MODE{
        {"time_step"       , CONDITION_LOAD_MODE::time_step        },
        {"tree"            , CONDITION_LOAD_MODE::tree             },
        {"cut_off"         , CONDITION_LOAD_MODE::cut_off          },
        {"ext_sys"         , CONDITION_LOAD_MODE::ext_sys          },
        {"ext_sys_sequence", CONDITION_LOAD_MODE::ext_sys_sequence },
        {"record"          , CONDITION_LOAD_MODE::record           },

        {"molecule"        , CONDITION_LOAD_MODE::molecule         },
        {"box"             , CONDITION_LOAD_MODE::box              },
        {"ex_radius"       , CONDITION_LOAD_MODE::ex_radius        },
    };

    static const std::map<CONDITION_LOAD_MODE, std::string> table_CONDITION_LOAD_MODE_str{
        {CONDITION_LOAD_MODE::time_step       , "time_step"        },
        {CONDITION_LOAD_MODE::tree            , "tree"             },
        {CONDITION_LOAD_MODE::cut_off         , "cut_off"          },
        {CONDITION_LOAD_MODE::ext_sys         , "ext_sys"          },
        {CONDITION_LOAD_MODE::ext_sys_sequence, "ext_sys_sequence" },
        {CONDITION_LOAD_MODE::record          , "record"           },

        {CONDITION_LOAD_MODE::molecule        , "molecule"         },
        {CONDITION_LOAD_MODE::box             , "box"              },
        {CONDITION_LOAD_MODE::ex_radius       , "ex_radius"        },
    };

    CONDITION_LOAD_MODE which_CONDITION_LOAD_MODE(const std::string &str){
        if(table_str_CONDITION_LOAD_MODE.find(str) != table_str_CONDITION_LOAD_MODE.end()){
            return table_str_CONDITION_LOAD_MODE.at(str);
        } else {
            std::cerr << "  CONDITION_LOAD_MODE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in CONDITION_LOAD_MODE.");
        }
    }

    std::string what(const CONDITION_LOAD_MODE &e){
        if(table_CONDITION_LOAD_MODE_str.find(e) != table_CONDITION_LOAD_MODE_str.end()){
            return table_CONDITION_LOAD_MODE_str.at(e);
        } else {
            using type_base = typename std::underlying_type<CONDITION_LOAD_MODE>::type;
            std::cerr << "  CONDITION_LOAD_MODE: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in CONDITION_LOAD_MODE.");
        }
    }

}

inline std::ostream& operator << (std::ostream& s, const CONDITION_LOAD_MODE &e){
    s << ENUM::what(e);
    return s;
}


namespace System {

    //--- setting data class: DO NOT contain pointer or container.
    class setting_parameter{
    public:
        //--- for time step
        PS::S64 istep    = -1;
        PS::S64 nstep_st = -1;
        PS::S64 nstep_ed = -1;
        PS::F64 dt       = 0.0;

        //--- for Tree
        PS::S32 n_leaf_limit  = -1;
        PS::F64 coef_ema      = -1.0;
        PS::F64 theta         = -1.0;
        PS::S32 n_group_limit = -1;
        PS::S32 cycle_dinfo   = -1;

        //--- for cut_off radius
        PS::F64 cut_off_LJ    = -1.0;

        //--- for installing molecule at initialize
        PS::F64 ex_radius = -1.0;
        PS::S32 try_limit = -1;

        //--- for recording data
        PS::S64 pos_interval    = -1;
        PS::S64 pos_start       = -1;
        PS::S64 resume_interval = -1;
        PS::S64 resume_start    = -1;
        PS::S64 pdb_interval    = -1;
        PS::S64 pdb_start       = -1;

        PS::S64 eng_interval  = -1;
        PS::S64 eng_start     = -1;
        PS::S64 prop_interval = -1;
        PS::S64 prop_start    = -1;
    };
    setting_parameter setting;

    //--- for initialize particle
    std::vector<std::pair<MolName, PS::S64>>        model_list;
    std::vector<std::vector<Atom_FP>>               model_template;
    std::vector<std::vector<std::vector<PS::S64>>>  bond_template;

    //--- sysc settings in MPI processes
    void broadcast_setting(const PS::S32 root = 0){
        COMM_TOOL::broadcast(setting, root);

        Normalize::broadcast_boxSize(root);

        COMM_TOOL::broadcast(model_list,     root);
        COMM_TOOL::broadcast(model_template, root);
        COMM_TOOL::broadcast(bond_template,  root);
    }

    //--- load settings from file
    void loading_sequence_condition(const std::string         &file_name,
                                          EXT_SYS::Sequence   &sequence,
                                          EXT_SYS::Controller &controller){
        std::ifstream file_sys{file_name};
        std::string line;

        if(file_sys.fail()) throw std::ios_base::failure("file: " + file_name + " was not found.");

        //--- temporary value
        PS::S64 n_chain, n_rep, n_nys;
        PS::F64 NVT_freq, NPT_freq;

        CONDITION_LOAD_MODE mode = CONDITION_LOAD_MODE::time_step;
        while ( getline(file_sys, line) ) {
            STR_TOOL::removeCR(line);
            std::vector<std::string> str_list = STR_TOOL::split(line, " ");

            //--- skip empty or comment line(start as "!" or "//")
            if(str_list.empty()) continue;
            if(str_list[0].empty()) continue;
            if(str_list[0].substr(0,1) == "!" ||
               str_list[0].substr(0,2) == "//") continue;

            //--- header information
            if(       str_list[0] == "@<CONDITION>TIMESTEP"){
                mode = CONDITION_LOAD_MODE::time_step;
            } else if(str_list[0] == "@<CONDITION>TREE"){
                mode = CONDITION_LOAD_MODE::tree;
            } else if(str_list[0] == "@<CONDITION>CUT_OFF"){
                mode = CONDITION_LOAD_MODE::cut_off;
            } else if(str_list[0] == "@<CONDITION>EXT_SYS"){
                mode = CONDITION_LOAD_MODE::ext_sys;
            } else if(str_list[0] == "@<CONDITION>EXT_SYS_SEQUENCE"){
                mode = CONDITION_LOAD_MODE::ext_sys_sequence;
            } else if(str_list[0] == "@<CONDITION>RECORD"){
                mode = CONDITION_LOAD_MODE::record;
            }

            //--- loading data
            switch (mode) {
                case CONDITION_LOAD_MODE::time_step:
                    if( str_list.size() < 2) continue;

                    if( str_list[0] == "t_start"){
                        System::setting.nstep_st = std::stoi(str_list[1]);
                        System::setting.istep    = System::setting.nstep_st;
                    }
                    if( str_list[0] == "t_end")   System::setting.nstep_ed = std::stoi(str_list[1]);
                    if( str_list[0] == "dt")      System::setting.dt       = std::stof(str_list[1])
                                                                            *Unit::femto_second
                                                                            /Unit::norm_time;
                break;

                case CONDITION_LOAD_MODE::tree:
                    if( str_list.size() < 2) continue;

                    if( str_list[0] == "n_leaf_limit")  System::setting.n_leaf_limit  = std::stoi(str_list[1]);
                    if( str_list[0] == "coef_ema")      System::setting.coef_ema      = std::stof(str_list[1]);
                    if( str_list[0] == "theta")         System::setting.theta         = std::stof(str_list[1]);
                    if( str_list[0] == "n_group_limit") System::setting.n_group_limit = std::stoi(str_list[1]);
                    if( str_list[0] == "cycle_dinfo")   System::setting.cycle_dinfo   = std::stoi(str_list[1]);
                break;

                case CONDITION_LOAD_MODE::cut_off:
                    if( str_list.size() < 2) continue;

                    if( str_list[0] == "LJ")    System::setting.cut_off_LJ    = std::stof(str_list[1]);
                break;

                case CONDITION_LOAD_MODE::ext_sys:
                    if( str_list.size() < 2) continue;

                    //--- load controller settings
                    if( str_list[0] == "n_chain"  ) n_chain  = std::stoi(str_list[1]);
                    if( str_list[0] == "n_rep"    ) n_rep    = std::stoi(str_list[1]);
                    if( str_list[0] == "n_nys"    ) n_nys    = std::stoi(str_list[1]);
                    if( str_list[0] == "NVT_freq" ) NVT_freq = std::stof(str_list[1]);
                    if( str_list[0] == "NPT_freq" ) NPT_freq = std::stof(str_list[1]);

                    //--- load default control setting
                    if( str_list[0] == "default" ){
                        sequence.setDefault( EXT_SYS::Setting( ENUM::which_EXT_SYS_MODE( str_list[1] ),
                                                               0,   // period. no meaning for default setting.
                                                               std::stof(str_list[2]),
                                                               std::stof(str_list[3]) ) );
                    }
                break;

                case CONDITION_LOAD_MODE::ext_sys_sequence:
                    if( str_list.size() < 4) continue;

                    //--- load control sequence
                    sequence.addSetting( EXT_SYS::Setting( ENUM::which_EXT_SYS_MODE( str_list[0] ),
                                                           std::stoi(str_list[1]),
                                                           std::stof(str_list[2]),
                                                           std::stof(str_list[3]) ) );
                break;

                case CONDITION_LOAD_MODE::record:
                    if( str_list.size() < 2) continue;

                    if( str_list[0] == "pos_interval")    System::setting.pos_interval    = stoi(str_list[1]);
                    if( str_list[0] == "pos_start")       System::setting.pos_start       = stoi(str_list[1]);
                    if( str_list[0] == "resume_interval") System::setting.resume_interval = stoi(str_list[1]);
                    if( str_list[0] == "resume_start")    System::setting.resume_start    = stoi(str_list[1]);
                    if( str_list[0] == "pdb_interval")    System::setting.pdb_interval    = stoi(str_list[1]);
                    if( str_list[0] == "pdb_start")       System::setting.pdb_start       = stoi(str_list[1]);

                    if( str_list[0] == "eng_interval")  System::setting.eng_interval  = stoi(str_list[1]);
                    if( str_list[0] == "eng_start")     System::setting.eng_start     = stoi(str_list[1]);
                    if( str_list[0] == "prop_interval") System::setting.prop_interval = stoi(str_list[1]);
                    if( str_list[0] == "prop_start")    System::setting.prop_start    = stoi(str_list[1]);
                break;

                default:
                    std::cerr << "  file: " << file_name << std::endl;
                    throw std::invalid_argument("undefined loading mode.");
            }
        }

        //--- initialize ext_sys controller
        controller.init(n_chain,
                        n_rep,
                        n_nys,
                        NVT_freq,
                        NPT_freq);
    }

    void loading_molecular_condition(const std::string &file_name){
        std::ifstream file_sys{file_name};
        std::string line;

        if(file_sys.fail()) throw std::ios_base::failure("file: " + file_name + " was not found.");

        CONDITION_LOAD_MODE mode = CONDITION_LOAD_MODE::molecule;
        while ( getline(file_sys, line) ) {
            STR_TOOL::removeCR(line);
            std::vector<std::string> str_list = STR_TOOL::split(line, " ");

            //--- skip empty or comment line(start as "!" or "//")
            if(str_list.empty()) continue;
            if(str_list[0].empty()) continue;
            if(str_list[0].substr(0,1) == "!" ||
               str_list[0].substr(0,2) == "//") continue;

            //--- header information
            if(       str_list[0] == "@<CONDITION>MOLECULE"){
                mode = CONDITION_LOAD_MODE::molecule;
            } else if(str_list[0] == "@<CONDITION>BOX"){
                mode = CONDITION_LOAD_MODE::box;
            } else if(str_list[0] == "@<CONDITION>EX_RADIUS"){
                mode = CONDITION_LOAD_MODE::ex_radius;
            }

            //--- loading data
            switch (mode) {
                case CONDITION_LOAD_MODE::molecule:
                    if( str_list.size() < 2) continue;

                    System::model_list.push_back( std::make_pair( ENUM::which_MolName(str_list[0]),
                                                                  std::stoi(str_list[1])           ) );
                break;

                case CONDITION_LOAD_MODE::box:
                    if( str_list.size() < 3) continue;

                    Normalize::setBoxSize( PS::F32vec{ std::stof(str_list[0]),
                                                       std::stof(str_list[1]),
                                                       std::stof(str_list[2]) } );
                break;

                case CONDITION_LOAD_MODE::ex_radius:
                    if( str_list.size() < 2) continue;

                    System::setting.ex_radius = std::stof(str_list[0]);
                    System::setting.try_limit = std::stoi(str_list[1]);
                break;

                default:
                    std::cerr << "  file: " << file_name << std::endl;
                    throw std::invalid_argument("undefined loading mode.");
            }
        }
        //--- reach EOF

        //--- allocate model template vector
        model_template.resize(model_list.size());
        bond_template.resize(model_list.size());
    }

    //--- display system settings
    void print_setting(){
        std::ostringstream oss;

        oss << "system conditions:\n";
        oss << "  Time steps:\n";
        oss << "    istep    = " << setting.istep    << "\n";
        oss << "    nstep_st = " << setting.nstep_st << "\n";
        oss << "    nstep_ed = " << setting.nstep_ed << "\n";
        oss << "    dt       = " << setting.dt       << "\n";
        oss << "             = " << setting.dt/Unit::femto_second*Unit::norm_time << " [fs]\n";
        oss << "\n";

        oss << "  Tree settings:\n";
        oss << "    n_leaf_limit  = " << setting.n_leaf_limit  << "\n";
        oss << "    coef_ema      = " << setting.coef_ema      << "\n";
        oss << "    theta         = " << setting.theta         << "\n";
        oss << "    n_group_limit = " << setting.n_group_limit << "\n";
        oss << "    cycle_dinfo   = " << setting.cycle_dinfo   << "\n";
        oss << "\n";

        oss << "  Cut_off settings:\n";
        oss << "    cut_off_LJ      = " << std::setw(9) << std::setprecision(7) << setting.cut_off_LJ    << " [angstrom]\n";
        oss << "    cut_off_coulomb :  fixed as " << Normalize::normCutOff_PM() << " at normalized space.";
        oss << "\n";

        oss << "  Init molecule settings:\n";
        oss << "    box       = (" << Normalize::getBoxSize().x
            << ", "                << Normalize::getBoxSize().y
            << ", "                << Normalize::getBoxSize().z << ") [angstrom]\n";
        oss << "    ex_radius = "  << std::setw(9) << std::setprecision(7) << setting.ex_radius << " [angstrom]\n";
        oss << "    try_limit = "  << std::setw(9) <<                         setting.try_limit << " [times]\n";
        oss << "\n";

        oss << "model conditions:\n";
        if(System::model_list.size() != 0){
            for(auto m : System::model_list){
                oss << "    model_name = " << m.first  << "\n";
                oss << "    n_molecule = " << m.second << "\n";
                oss << "\n";
            }
        } else {
            oss << "    free.\n";
            oss << "\n";
        }

        oss << "recording conditions:\n";
        oss << "              " << std::setw(15) << "start    " << " | "
                                << std::setw(15) << "interval  \n";
        oss << "    pos     : " << std::setw(15) << setting.pos_start       << " | "
                                << std::setw(15) << setting.pos_interval    << "\n";
        oss << "    resume  : " << std::setw(15) << setting.resume_start    << " | "
                                << std::setw(15) << setting.resume_interval << "\n";
        oss << "    pdb     : " << std::setw(15) << setting.pdb_start       << " | "
                                << std::setw(15) << setting.pdb_interval    << "\n";
        oss << "\n";
        oss << "    energy  : " << std::setw(15) << setting.eng_start     << " | "
                                << std::setw(15) << setting.eng_interval  << "\n";
        oss << "    property: " << std::setw(15) << setting.prop_start    << " | "
                                << std::setw(15) << setting.prop_interval << "\n";
        oss << "\n";

        //--- output
        std::cout << oss.str() << std::flush;
    }



    //--- input setting to dinfo
    template<class Tdinfo>
    void InitDinfo(Tdinfo &dinfo){
        dinfo.initialize(setting.coef_ema);
        dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
        dinfo.setPosRootDomain(PS::F32vec{0.0, 0.0, 0.0}, PS::F32vec{1.0, 1.0, 1.0});  // fixed size for using PS::ParticleMesh
    }

    //--- input setting to tree
    template<class Ttree>
    void InitTree(const size_t &n_ptcl,
                        Ttree   &tree){
        tree.initialize(n_ptcl,
                        setting.theta,
                        setting.n_leaf_limit,
                        setting.n_group_limit);
    }

    //--- time step interface
    void StepNext(){ setting.istep++; }
    PS::S64 getIstep(){ return setting.istep; }

    //--- dinfo timing
    bool isDinfoUpdate(){
        assert(setting.cycle_dinfo > 0);
        return ( ((setting.istep - setting.nstep_st) % setting.cycle_dinfo) == 0 );
    }

    //--- judge continue lopp or not
    bool isLoopContinue(){
        assert(setting.istep    >= 0);
        assert(setting.nstep_ed >= 0);
        return ( setting.istep <= setting.nstep_ed );
    }

    //--- getter
    PS::S64 get_istep(){
        assert(setting.istep >= 0);
        return setting.istep;
    }
    PS::F64 get_dt() {
        assert(setting.dt > 0.0);
        return setting.dt;
    }
    PS::F64 get_cutoff_LJ()    { return setting.cut_off_LJ;    }

    PS::S64 get_eng_start()    { return setting.eng_start;     }
    PS::S64 get_prop_start()   { return setting.prop_start;    }
    PS::S64 get_eng_interval() { return setting.eng_interval;  }
    PS::S64 get_prop_interval(){ return setting.prop_interval; }
}
