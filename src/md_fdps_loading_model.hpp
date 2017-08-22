//***************************************************************************************
//  This program is loading model parameter function.
//***************************************************************************************
#pragma once

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include "str_tool.hpp"

#include "md_fdps_coef_table.hpp"
#include "md_fdps_atom_class.hpp"


//--- file loading mode
enum class MOL2_LOAD_MODE : int {
    header,
    atom,
    bond,
    angle,
    torsion,
};

//--- std::string converter for enum
namespace ENUM {

    std::string whatis(const MOL2_LOAD_MODE &e){
        switch (e) {
            case MOL2_LOAD_MODE::header:
                return "header";
            break;

            case MOL2_LOAD_MODE::atom:
                return "atom";
            break;

            case MOL2_LOAD_MODE::bond:
                return "bond";
            break;

            case MOL2_LOAD_MODE::angle:
                return "angle";
            break;

            case MOL2_LOAD_MODE::torsion:
                return "torsion";
            break;

            default:
                throw std::out_of_range("undefined enum value.");
        }
    }
}

inline std::ostream& operator << (std::ostream& s, const MOL2_LOAD_MODE &e){
    s << ENUM::whatis(e);
    return s;
}


namespace MODEL {

    const std::string model_dir{"./model/"};

    template<class Tptcl, class Tid>
    void loading_mol2_file(const std::string                   &model_name,
                                 std::vector<Tptcl>            &atom_list,
                                 std::vector<std::vector<Tid>> &bond_list){

        atom_list.clear();
        bond_list.clear();

        std::string file_name;
        if(model_name.find_first_of("/") != std::string::npos){
            file_name = model_name + ".mol2";
        } else {
            file_name = model_dir + model_name + ".mol2";
        }
        std::ifstream file_mol2{file_name};
        std::string line;

        if(file_mol2.fail()) throw std::ios_base::failure("file: " + model_name + " was not found.");

        //--- loading mol2 file
        PS::S32 line_count = 0;
        PS::S32 atom_count = 0;
        PS::S32 bond_count = 0;

        PS::S32 line_index = -1;  // illigal value
        PS::S32 n_atom     = -1;
        PS::S32 n_bond     = -1;
        MOL2_LOAD_MODE mode = MOL2_LOAD_MODE::header;
        while ( getline(file_mol2, line) ){
            STR_TOOL::removeCR(line);
        //    std::cout << "line = " << line << std::endl;
        //    for(int i=0; i<line.size(); i++){
        //        std::cout << static_cast<int>(line[i]) << std::endl;
        //    }

            line_count++;
            std::vector<std::string> str_list = STR_TOOL::split(line, " ");

            //--- skip empty or comment line(start as "!" or "//")
            if(str_list.empty()) continue;
            if(str_list[0].empty()) continue;
            if(str_list[0].substr(0,1) == "!" ||
               str_list[0].substr(0,2) == "//") continue;

            //--- header information
            if(       str_list[0] == "@<TRIPOS>MOLECULE"){
                line_index = line_count + 2;
                if(mode != MOL2_LOAD_MODE::header){
                    std::cerr << "  file: " << file_name << std::endl;
                    throw std::length_error("*.mol2 file must contain single molecule only.");
                }
            } else if(str_list[0] == "@<TRIPOS>ATOM"){
                    mode = MOL2_LOAD_MODE::atom;
            } else if(str_list[0] == "@<TRIPOS>BOND"){
                    mode = MOL2_LOAD_MODE::bond;

                    //--- initialize bond list
                    bond_list.resize(atom_count);
            }

            //--- loading data
            Atom_FP atom_tmp;
            PS::F64 id_i, id_j;
            switch (mode) {
                case MOL2_LOAD_MODE::header:
                    if(line_count == line_index){
                        n_atom = std::stoi(str_list[0]);
                        n_bond = std::stoi(str_list[1]);
                    }
                break;

                case MOL2_LOAD_MODE::atom:
                    //--- check format: 9 parameters, str_list[0] must be digit.
                    if( str_list.size() < 9) continue;
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;

                    atom_count++;
                    if( atom_count != std::stoi(str_list[0]) ){
                        std::cerr << "  file: " << file_name << std::endl;
                        std::cerr << "  atom_count = " << atom_count << std::endl;
                        std::cerr << "  id in file = " << str_list[0] << std::endl;
                        throw std::invalid_argument("atom id in *.mol2 file must be consective noumbers to start with 1.");
                    }

                    atom_tmp.setAtomID(std::stoi(str_list[0]));
                    atom_tmp.setMolID(-1);   // this is parameter table. not real atom.
                    atom_tmp.setAtomType(str_list[1]);
                    atom_tmp.setMolType(model_name);
                    atom_tmp.setPos( PS::F32vec{std::stof(str_list[2]),
                                                std::stof(str_list[3]),
                                                std::stof(str_list[4]) } );
                    atom_tmp.setCharge(std::stof(str_list[8]));
                    atom_list.push_back(atom_tmp);
                break;

                case MOL2_LOAD_MODE::bond:
                    //--- check format: 3 parameters, str_list[0], [1], and [2] must be integer.
                    if( str_list.size() < 3) continue;
                    if( !(STR_TOOL::isInteger(str_list[0]) &&
                          STR_TOOL::isInteger(str_list[1]) &&
                          STR_TOOL::isInteger(str_list[2]) ) ) continue;

                    bond_count++;
                    if( bond_count != std::stoi(str_list[0]) ){
                        std::cerr << "  file: " << file_name << std::endl;
                        std::cerr << "  bond_count = " << bond_count << std::endl;
                        std::cerr << "  id in file = " << str_list[0] << std::endl;
                        throw std::invalid_argument("bond id in *.mol2 file must be consective noumbers to start with 1.");
                    }

                    id_i = std::stoi(str_list[1]) - 1;
                    id_j = std::stoi(str_list[2]) - 1;

                    bond_list.at(id_i).push_back(id_j);
                    bond_list.at(id_j).push_back(id_i);
                break;

                default:
                    std::cerr << "  file: " << file_name << std::endl;
                    throw std::invalid_argument("undefined loading mode.");
            }
        }

        //--- check file consistency
        if(n_atom != atom_count ||
           n_bond != bond_count ){
            std::cerr << "  file: " << file_name << std::endl;
            std::cerr << "  n_atom    = " << n_atom    << " / atom_count    = " << atom_count    << std::endl;
            std::cerr << "  n_bond    = " << n_bond    << " / bond_count    = " << bond_count    << std::endl;
            throw std::invalid_argument("number of defined parameters are not match with header.");
        }
    }


    template<class Ttable_inter,class Ttable_bond, class Ttable_angle, class Ttable_torsion>
    void loading_param_file(const std::string    &model_name,
                                  Ttable_inter   &inter_table,
                                  Ttable_bond    &bond_table,
                                  Ttable_angle   &angle_table,
                                  Ttable_torsion &torsion_table){

        std::string file_name;
        if(model_name.find_first_of("/") != std::string::npos){
            file_name = model_name + ".param";
        } else {
            file_name = model_dir + model_name + ".param";
        }
        std::ifstream file_para{file_name};
        std::string line;

        if(file_para.fail()) throw std::ios_base::failure("file: " + model_name + ".para was not found.");

        //--- loading para file
        PS::S32 elem_count    = 0;
        PS::S32 bond_count    = 0;
        PS::S32 angle_count   = 0;
        PS::S32 torsion_count = 0;

        PS::S32 n_elem    = -1;
        PS::S32 n_bond    = -1;
        PS::S32 n_angle   = -1;
        PS::S32 n_torsion = -1;
        MOL2_LOAD_MODE mode = MOL2_LOAD_MODE::atom;
        while (getline(file_para, line)){
            STR_TOOL::removeCR(line);
            std::vector<std::string> str_list = STR_TOOL::split(line, " ");

            //--- skip empty or comment line(start as "!" or "//")
            if(str_list.empty()) continue;
            if(str_list[0].empty()) continue;
            if(str_list[0].substr(0,1) == "!" ||
               str_list[0].substr(0,2) == "//") continue;

            //--- header information
            if(       str_list[0] == "@<PARA>ELEMENT"){
                mode = MOL2_LOAD_MODE::atom;
                n_elem = std::stoi(str_list.at(1));
            } else if(str_list[0] == "@<PARA>BOND"){
                mode = MOL2_LOAD_MODE::bond;
                n_bond = std::stoi(str_list.at(1));
            } else if(str_list[0] == "@<PARA>ANGLE"){
                mode = MOL2_LOAD_MODE::angle;
                n_angle = std::stoi(str_list.at(1));
            } else if(str_list[0] == "@<PARA>TORSION"){
                mode = MOL2_LOAD_MODE::torsion;
                n_torsion = std::stoi(str_list.at(1));
            }

            //--- loading data
            CoefElement coef_elem;
            CoefBond    coef_bond;
            CoefAngle   coef_angle;
            CoefTorsion coef_torsion;
            MODEL::KeyElem    key_elem;
            MODEL::KeyBond    key_bond;
            MODEL::KeyAngle   key_angle;
            MODEL::KeyTorsion key_torsion;
            switch (mode) {
                case MOL2_LOAD_MODE::atom:
                    //--- check format: 5 parameters, str_list[0] must be integer.
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;
                    if( str_list.size() < 5) continue;

                    coef_elem.mass   = 1.e-3*std::stof(str_list[2]);  // convert [g/mol] -> [kg/mol]
                    coef_elem.charge = 0.0;                           // load from *.mol2 file.
                    coef_elem.vdw_d  = std::stof(str_list[3]);
                    coef_elem.vdw_r  = std::stof(str_list[4]);

                    key_elem = std::make_tuple( ENUM::which_MolName(model_name),
                                                ENUM::which_AtomName(str_list[1]) );
                    //--- check duplication
                    if( inter_table.find(key_elem) != inter_table.cend() ){
                        std::cerr << "WARNIG: 'coef_elem' of " + ENUM::whatis(key_elem) + " is overloaded." << std::endl;
                    }

                    inter_table[key_elem] = coef_elem;
                    elem_count++;
                break;

                case MOL2_LOAD_MODE::bond:
                    //--- check format: 7 parameters, str_list[0] must be integer.
                    if( str_list.size() < 7) continue;
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;

                    coef_bond.r0 = std::stof(str_list[4]);
                    coef_bond.d  = std::stof(str_list[5]);
                    coef_bond.a  = std::stof(str_list[6]);

                    if(       str_list[3] == "anharmonic"){
                        coef_bond.form = IntraFuncForm::anharmonic;
                    } else if(str_list[3] == "harmonic"){
                        coef_bond.form = IntraFuncForm::harmonic;
                    } else {
                        std::cerr << "  file: " << file_name << std::endl;
                        throw std::invalid_argument("undefined form of bond potential.");
                    }

                    key_bond = std::make_tuple( ENUM::which_MolName(model_name),
                                                ENUM::which_AtomName(str_list[1]),
                                                ENUM::which_AtomName(str_list[2]) );
                    //--- check duplication
                    if( bond_table.find(key_bond) != bond_table.cend() ){
                        std::cerr << "WARNIG: 'coef_bond' of " + model_name
                                                       + ": "  + str_list[1]
                                                       + " - " + str_list[2]
                                   + " is overloaded." << std::endl;
                    }

                    bond_table[key_bond] = coef_bond;

                    //--- add reverse shape
                    key_bond = std::make_tuple( ENUM::which_MolName(model_name),
                                                ENUM::which_AtomName(str_list[2]),
                                                ENUM::which_AtomName(str_list[1]) );
                    bond_table[key_bond] = coef_bond;

                    //--- enclement
                    bond_count++;
                break;

                case MOL2_LOAD_MODE::angle:
                    //--- check format: 7 parameters, str_list[0] must be integer.
                    if( str_list.size() < 7) continue;
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;

                    coef_angle.theta0 = std::stof(str_list[5])*Unit::pi/180.0;                           // [degree] -> [rad]
                    coef_angle.k      = std::stof(str_list[6])/std::pow(std::sin(coef_angle.theta0),2);  // make [kcal/mol rad^2]/sin(theta0)^2
                    coef_angle.form   = ENUM::which_IntraFuncForm(str_list[4]);

                    key_angle = std::make_tuple( ENUM::which_MolName(model_name),
                                                 ENUM::which_AtomName(str_list[1]),
                                                 ENUM::which_AtomName(str_list[2]),
                                                 ENUM::which_AtomName(str_list[3]) );
                    //--- check duplication
                    if( angle_table.find(key_angle) != angle_table.cend() ){
                        std::cerr << "WARNIG: 'coef_angle' of " + model_name
                                                        + ": "  + str_list[1]
                                                        + " - " + str_list[2]
                                                        + " - " + str_list[3]
                                   + " is overloaded." << std::endl;
                    }

                    angle_table[key_angle] = coef_angle;
                    angle_count++;
                break;

                case MOL2_LOAD_MODE::torsion:
                    //--- check format: 9 parameters, str_list[0] must be integer.
                    if( str_list.size() < 9) continue;
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;

                    coef_torsion.theta0 = std::stof(str_list[6])*Unit::pi/180.0;   // [degree] -> [rad]
                    coef_torsion.k      = std::stof(str_list[7]);
                    coef_torsion.n_min  = std::stof(str_list[8]);
                    coef_torsion.form   = ENUM::which_IntraFuncForm(str_list[5]);

                    key_torsion = std::make_tuple( ENUM::which_MolName(model_name),
                                                   ENUM::which_AtomName(str_list[1]),
                                                   ENUM::which_AtomName(str_list[2]),
                                                   ENUM::which_AtomName(str_list[3]),
                                                   ENUM::which_AtomName(str_list[4]) );
                    //--- check duplication
                    if( torsion_table.find(key_torsion) != torsion_table.cend() ){
                        std::cerr << "WARNIG: 'coef_angle' of " + model_name
                                                        + ": "  + str_list[1]
                                                        + " - " + str_list[2]
                                                        + " - " + str_list[3]
                                                        + " - " + str_list[4]
                                   + " is overloaded." << std::endl;
                    }

                    torsion_table[key_torsion] = coef_torsion;
                    torsion_count++;
                break;

                default:
                    std::cerr << "  file: " << file_name << std::endl;
                    throw std::invalid_argument("undefined loading mode.");
            }
        }

        //--- check file consistency
        if(n_elem    != elem_count    ||
           n_bond    != bond_count    ||
           n_angle   != angle_count   ||
           n_torsion != torsion_count ){
            std::cerr << "  file: "       << file_name << "\n"
                      << "    n_elem    = " << n_elem    << " / elem_count    = " << elem_count    << "\n"
                      << "    n_bond    = " << n_bond    << " / bond_count    = " << bond_count    << "\n"
                      << "    n_angle   = " << n_angle   << " / angle_count   = " << angle_count   << "\n"
                      << "    n_torsion = " << n_torsion << " / torsion_count = " << torsion_count << std::endl;
            throw std::invalid_argument("number of defined parameters are not match with header.");
        }
    }

    template<class Tptcl, class Tid, class Ttable_inter ,class Ttable_bond, class Ttable_angle, class Ttable_torsion>
    void loading_model_parameter(const std::string                   &model_name,
                                       std::vector<Tptcl>            &atom_list,
                                       std::vector<std::vector<Tid>> &bond_list,
                                       Ttable_inter                  &inter_table,
                                       Ttable_bond                   &bond_table,
                                       Ttable_angle                  &angle_table,
                                       Ttable_torsion                &torsion_table) {

        //--- loading ****.mol2 file
        loading_mol2_file(model_name,
                          atom_list,
                          bond_list );

        //--- loading ****.param file
        loading_param_file(model_name,
                           inter_table,
                           bond_table,
                           angle_table,
                           torsion_table);

        //--- copy VDW coef to atom_list from inter_table
        for(size_t index=0; index<atom_list.size(); index++){
            CoefElement coef_elem;
            try{
                MODEL::KeyElem key_elem = std::make_tuple( atom_list[index].getMolType(),
                                                           atom_list[index].getAtomType() );
                coef_elem = inter_table.at(key_elem);
            }
            catch(std::out_of_range err){
                std::cerr << "  element parameter of " << ENUM::whatis(atom_list[index].getAtomType())
                          << " was not found in '" << model_name + ".param"
                          << "' file." << std::endl;
                throw;
            }
            atom_list[index].setMass(  coef_elem.mass );
            atom_list[index].setVDW_D( coef_elem.vdw_d );
            atom_list[index].setVDW_R( coef_elem.vdw_r );
        }

        //--- copy charge to inter_table from atom_list
        for(auto itr = inter_table.begin(); itr != inter_table.end(); itr++){
            PS::F32 charge_buff;
            bool find_flag = false;
            for(auto atom : atom_list){
                if(std::make_tuple(atom.getMolType(), atom.getAtomType()) == itr->first){
                    #ifdef TEST_MOL_INSTALL
                        if(find_flag){
                            if(charge_buff != atom.getCharge()){
                                std::cerr << "WARNING: charge value of "
                                          << ENUM::whatis(itr->first)
                                          << " is overloaded as "
                                          << std::to_string(atom.getCharge()) << std::endl;
                            }
                        }
                    #endif

                    charge_buff = atom.getCharge();
                    find_flag = true;
                }
            }
            #ifdef TEST_MOL_INSTALL
                if( !find_flag && std::get<0>(itr->first) == ENUM::which_MolName(model_name) ){
                    std::cerr << "WARNING: " << ENUM::whatis(itr->first)
                              << " is not defined in "
                              << model_name + ".param"  << std::endl;
                    std::cerr << "  charge of " << ENUM::whatis(itr->first)
                              << " was set as 0 in inter_table." << std::endl;
                }
            #endif
            itr->second.charge = charge_buff;
        }
    }


    //--- display interface
    template <class Tatom>
    std::string print_atom(const Tatom &atom){
        std::ostringstream oss;

        oss << "    AtomID   : " << atom.getAtomID()   << "\n";
        oss << "    MolID    : " << atom.getMolID()    << "\n";
        oss << "    AtomType : " << atom.getAtomType() << "\n";
        oss << "    MolType  : " << atom.getMolType()  << "\n";
        oss << "    Pos      : ("  << atom.getPos().x
                           << ", " << atom.getPos().y
                           << ", " << atom.getPos().z << ")\n";
        oss << "    Mass     : " << atom.getMass()   << "\n";
        oss << "    charge   : " << atom.getCharge() << "\n";
        oss << "    VDW_R    : " << atom.getVDW_R()  << "\n";
        oss << "    VDW_D    : " << atom.getVDW_D()  << "\n";

        return oss.str();
    }

    template<class Tptcl, class Tid>
    void print_model_template(const std::vector<Tptcl>            &atom_list,
                              const std::vector<std::vector<Tid>> &bond_list){

        std::ostringstream oss;

        for(size_t index=0; index<atom_list.size(); ++index){
            auto& atom = atom_list.at(index);
            auto& bond = bond_list.at(index);

            oss << "  index = " << index << "\n";
            oss << print_atom(atom);
            oss << "    bond     : n=" << bond_list.at(index).size() << ", id=";
            for(auto b : bond){
                oss << " " << b;
            }
            oss << "\n";
            oss << "\n";
        }
        oss << "\n";

        std::cout << oss.str();
    }

    template<class Ttable>
    void print_coef_table(const Ttable &table, const MolName &model){

        std::string str;
        size_t count = 0;
        for(auto coef : table){
            if( std::get<0>(coef.first) != model ) continue;

            count++;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::whatis(coef.first) + "\n";
            str += coef.second.to_string(6);
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }

    template<class Ttable>
    void print_coef_table(const Ttable &table){

        std::string str;
        size_t count = 0;
        for(auto coef : table){
            count++;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::whatis(coef.first) + "\n";
            str += coef.second.to_string(6);
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }
}
