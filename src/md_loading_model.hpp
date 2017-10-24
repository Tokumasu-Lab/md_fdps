//***************************************************************************************
//  This is loading model parameter function.
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

#include <molecular_dynamics_ext.hpp>

#include "atom_class.hpp"
#include "md_coef_table.hpp"


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

    static const std::map<std::string, MOL2_LOAD_MODE> table_str_MOL2_LOAD_MODE{
        {"MOLECULE", MOL2_LOAD_MODE::header  },
        {"ATOM"    , MOL2_LOAD_MODE::atom    },
        {"BOND"    , MOL2_LOAD_MODE::bond    },
        {"ANGLE"   , MOL2_LOAD_MODE::angle   },
        {"TORSION" , MOL2_LOAD_MODE::torsion },
    };

    static const std::map<MOL2_LOAD_MODE, std::string> table_MOL2_LOAD_MODE_str{
        {MOL2_LOAD_MODE::header , "MOLECULE"},
        {MOL2_LOAD_MODE::atom   , "ATOM"    },
        {MOL2_LOAD_MODE::bond   , "BOND"    },
        {MOL2_LOAD_MODE::angle  , "ANGLE"   },
        {MOL2_LOAD_MODE::torsion, "TORSION" },
    };

    MOL2_LOAD_MODE which_MOL2_LOAD_MODE(const std::string &str){
        if(table_str_MOL2_LOAD_MODE.find(str) != table_str_MOL2_LOAD_MODE.end()){
            return table_str_MOL2_LOAD_MODE.at(str);
        } else {
            std::cerr << "  MOL2_LOAD_MODE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in MOL2_LOAD_MODE.");
        }
    }

    std::string what(const MOL2_LOAD_MODE &e){
        if(table_MOL2_LOAD_MODE_str.find(e) != table_MOL2_LOAD_MODE_str.end()){
            return table_MOL2_LOAD_MODE_str.at(e);
        } else {
            using type_base = typename std::underlying_type<MOL2_LOAD_MODE>::type;
            std::cerr << "  MOL2_LOAD_MODE: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in MOL2_LOAD_MODE.");
        }
    }

}

inline std::ostream& operator << (std::ostream& s, const MOL2_LOAD_MODE &e){
    s << ENUM::what(e);
    return s;
}


namespace MODEL {

    const std::string model_dir{"./model/"};

    template<class Tptcl, class Tid>
    void loading_mol2_file(const std::string                   &model_name,
                           const std::string                   &file_name,
                                 std::vector<Tptcl>            &atom_list,
                                 std::vector<std::vector<Tid>> &bond_list){

        std::ifstream file_mol2{file_name};
        if(file_mol2.fail()) throw std::ios_base::failure("file: " + file_name + " was not found.");

        //--- loading mol2 file
        PS::S32 line_count = 0;
        PS::S32 atom_count = 0;
        PS::S32 bond_count = 0;

        PS::S32 line_index = -1;  // illigal value
        PS::S32 n_atom     = -1;
        PS::S32 n_bond     = -1;
        MOL2_LOAD_MODE mode = MOL2_LOAD_MODE::header;

        std::string line;
        while ( getline(file_mol2, line) ){
            STR_TOOL::removeCR(line);

            line_count++;
            std::vector<std::string> str_list = STR_TOOL::split(line, " ");

            //--- skip empty or comment line(start as "!" or "//")
            if(str_list.empty()) continue;
            if(str_list[0].empty()) continue;
            if(str_list[0].substr(0,1) == "!" ||
               str_list[0].substr(0,2) == "//") continue;

            //--- header information
            std::string mark     = "@<TRIPOS>";
            size_t      mark_len = mark.size();
            if(str_list[0].substr(0,mark_len) == mark){
                mode = ENUM::which_MOL2_LOAD_MODE( str_list[0].substr(mark_len) );

                if(mode == MOL2_LOAD_MODE::header){
                    line_index = line_count + 2;
                    if(mode != MOL2_LOAD_MODE::header){
                        throw std::length_error("*.mol2 file must contain single molecule only.");
                    }
                } else if(mode == MOL2_LOAD_MODE::atom){
                    atom_list.clear();
                } else if(mode == MOL2_LOAD_MODE::bond){
                    bond_list.clear();
                    bond_list.resize(atom_count);
                }

                continue;
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
                        std::cerr << "  atom_count = " << atom_count << "\n"
                                  << "  id in file = " << str_list[0] << std::endl;
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
                        std::cerr << "  bond_count = " << bond_count << "\n"
                                  << "  id in file = " << str_list[0] << std::endl;
                        throw std::invalid_argument("bond id in *.mol2 file must be consective noumbers to start with 1.");
                    }

                    id_i = std::stoi(str_list[1]) - 1;
                    id_j = std::stoi(str_list[2]) - 1;

                    bond_list.at(id_i).push_back(id_j);
                    bond_list.at(id_j).push_back(id_i);
                break;

                default:
                    throw std::invalid_argument("undefined loading mode.");
            }
        }

        //--- check file consistency
        if(n_atom != atom_count ||
           n_bond != bond_count ){
            std::cerr << "  n_atom    = " << n_atom << " / atom_count    = " << atom_count << "\n"
                      << "  n_bond    = " << n_bond << " / bond_count    = " << bond_count << std::endl;
            throw std::invalid_argument("number of defined parameters are not match with header.");
        }
    }

    template<class Ttable_inter, class Ttable_res>
    void loading_param_atom(const std::string              &model_name,
                            const std::vector<std::string> &str_list,
                                  Ttable_inter             &atom_table,
                                  Ttable_res               &residue_table){
        MODEL::CoefAtom coef_atom;
        MODEL::KeyAtom  key_atom;

        //--- check format: 6 parameters, str_list[0] must be integer.
        if( str_list.size() < 6 ) return;

        assert(str_list[2].size() <= 3);                  // residue name must be <= 3 characters.

        coef_atom.mass     = 1.e-3*std::stof(str_list[3]);  // convert [g/mol] -> [kg/mol]
        coef_atom.charge   = 0.0;                           // load from *.mol2 file.
        coef_atom.vdw_d    = std::stof(str_list[4]);
        coef_atom.vdw_r    = std::stof(str_list[5]);

        key_atom = std::make_tuple( ENUM::which_MolName(model_name),
                                    ENUM::which_AtomName(str_list[1]) );

        //--- check duplication
        if( atom_table.find(key_atom) != atom_table.cend() ){
            std::cerr << "WARNING: 'coef_atom' of " + ENUM::what(key_atom) + " is overloaded." << std::endl;
        }

        atom_table[key_atom]    = coef_atom;
        residue_table[key_atom] = str_list[2];
    }

    template<class Ttable_bond>
    void loading_param_bond(const std::string              &model_name,
                            const std::vector<std::string> &str_list,
                                   Ttable_bond             &bond_table){
        IntraFuncForm  func_form;
        CoefBond       coef_bond;
        MODEL::KeyBond key_bond;

        //--- check format: total 7 parameters
        //--- set void value in the case of IntraFuncForm == "none".
        if( str_list.size() < 4) return;

        func_form = ENUM::which_IntraFuncForm(str_list[3]);
        if(func_form == IntraFuncForm::none){
            coef_bond.form = func_form;
            coef_bond.r0   = 1.0;
            coef_bond.k    = 0.0;
            coef_bond.a    = 0.0;
        } else if(func_form == IntraFuncForm::harmonic){
            if( str_list.size() < 6) throw std::invalid_argument("invalid format for harmonic bond.");
            coef_bond.form = func_form;
            coef_bond.r0   = std::stof(str_list[4]);
            coef_bond.k    = std::stof(str_list[5]);
            coef_bond.a    = 0.0;
        } else if(func_form == IntraFuncForm::anharmonic){
            if( str_list.size() < 7) throw std::invalid_argument("invalid format for harmonic bond.");
            coef_bond.form = func_form;
            coef_bond.r0   = std::stof(str_list[4]);
            coef_bond.k    = std::stof(str_list[5]);
            coef_bond.a    = std::stof(str_list[6]);
        } else {
            throw std::invalid_argument("undefined form of bond potential.");
        }

        key_bond = std::make_tuple( ENUM::which_MolName(model_name),
                                    ENUM::which_AtomName(str_list[1]),
                                    ENUM::which_AtomName(str_list[2]) );
        //--- check duplication
        if( bond_table.find(key_bond) != bond_table.cend() ){
            std::cerr << "WARNING: 'coef_bond' of " << ENUM::what(key_bond) << " is overloaded." << std::endl;
        }
        bond_table[key_bond] = coef_bond;

        //--- add reverse shape
        key_bond = std::make_tuple( ENUM::which_MolName(model_name),
                                    ENUM::which_AtomName(str_list[2]),
                                    ENUM::which_AtomName(str_list[1]) );
        bond_table[key_bond] = coef_bond;
    }

    template<class Ttable_angle>
    void loading_param_angle(const std::string              &model_name,
                             const std::vector<std::string> &str_list,
                                   Ttable_angle             &angle_table){
        IntraFuncForm   func_form;
        CoefAngle       coef_angle;
        MODEL::KeyAngle key_angle;

        //--- check format: total 7 parameters
        //--- set void value in the case of IntraFuncForm == "none".
        if( str_list.size() < 5) return;

        func_form = ENUM::which_IntraFuncForm(str_list[4]);
        if( func_form == IntraFuncForm::none ){
            coef_angle.form   = ENUM::which_IntraFuncForm(str_list[4]);
            coef_angle.theta0 = 0.0;
            coef_angle.k      = 0.0;
        } else if(func_form == IntraFuncForm::harmonic){
            if( str_list.size() < 7 ) throw std::invalid_argument("invalid format for angle potential parameters.");

            coef_angle.form   = ENUM::which_IntraFuncForm(str_list[4]);
            coef_angle.theta0 = std::stof(str_list[5])*Unit::pi/180.0;                           // [degree] -> [rad]
            coef_angle.k      = std::stof(str_list[6])/std::pow(std::sin(coef_angle.theta0),2);  // make [kcal/mol rad^2]/sin(theta0)^2
        } else {
            throw std::invalid_argument("undefined form of angle potential.");
        }

        key_angle = std::make_tuple( ENUM::which_MolName(model_name),
                                     ENUM::which_AtomName(str_list[1]),
                                     ENUM::which_AtomName(str_list[2]),
                                     ENUM::which_AtomName(str_list[3]) );

        //--- check duplication
        if( angle_table.find(key_angle) != angle_table.cend() ){
            std::cerr << "WARNING: 'coef_angle' of " << ENUM::what(key_angle) << " is overloaded." << std::endl;
        }
        angle_table[key_angle] = coef_angle;

        //--- add reverse shape
        key_angle = std::make_tuple( ENUM::which_MolName(model_name),
                                     ENUM::which_AtomName(str_list[3]),
                                     ENUM::which_AtomName(str_list[2]),
                                     ENUM::which_AtomName(str_list[1]) );
        angle_table[key_angle] = coef_angle;
    }

    template<class Ttable_torsion>
    void loading_param_torsion(const std::string              &model_name,
                               const std::vector<std::string> &str_list,
                                     Ttable_torsion           &torsion_table){
        IntraFuncForm     func_form;
        CoefTorsion       coef_torsion;
        MODEL::KeyTorsion key_torsion;

        //--- check format: total 10 parameters,
        //--- set void value in the case of IntraFuncForm == "none".
        if( str_list.size() < 7) return;

        func_form = ENUM::which_IntraFuncForm(str_list[6]);
        if( func_form == IntraFuncForm::none ){
            coef_torsion.form   = func_form;
            coef_torsion.theta0 = 0.0;
            coef_torsion.k      = 0.0;
            coef_torsion.n_min  = 1;
        } else if( func_form == IntraFuncForm::cos ){
            if( str_list.size() < 10 ){
                throw std::invalid_argument("invalid format for torsion potential parameters.");
            }
            coef_torsion.form   = func_form;
            coef_torsion.theta0 = std::stof(str_list[7])*Unit::pi/180.0;   // [degree] -> [rad]
            coef_torsion.k      = std::stof(str_list[8]);
            coef_torsion.n_min  = std::stoi(str_list[9]);
        } else {
            throw std::invalid_argument("undefined form of torsion potential.");
        }

        key_torsion = std::make_tuple( ENUM::which_MolName(model_name),
                                       ENUM::which_TorsionShape(str_list[1]),
                                       ENUM::which_AtomName(str_list[2]),
                                       ENUM::which_AtomName(str_list[3]),
                                       ENUM::which_AtomName(str_list[4]),
                                       ENUM::which_AtomName(str_list[5])     );

        //--- check duplication
        if( torsion_table.find(key_torsion) != torsion_table.cend() ){
            std::cerr << "WARNING: 'coef_torsion' of " << ENUM::what(key_torsion) << " is overloaded." << std::endl;
        }
        torsion_table[key_torsion] = coef_torsion;

        //--- add reverse shape
        if( ENUM::which_TorsionShape(str_list[1]) == TorsionShape::dihedral ){
            key_torsion = std::make_tuple( ENUM::which_MolName(model_name),
                                           ENUM::which_TorsionShape(str_list[1]),
                                           ENUM::which_AtomName(str_list[5]),
                                           ENUM::which_AtomName(str_list[4]),
                                           ENUM::which_AtomName(str_list[3]),
                                           ENUM::which_AtomName(str_list[2]) );
        } else if( ENUM::which_TorsionShape(str_list[1]) == TorsionShape::improper ){
            key_torsion = std::make_tuple( ENUM::which_MolName(model_name),
                                           ENUM::which_TorsionShape(str_list[1]),
                                           ENUM::which_AtomName(str_list[5]),
                                           ENUM::which_AtomName(str_list[3]),
                                           ENUM::which_AtomName(str_list[4]),
                                           ENUM::which_AtomName(str_list[2]) );
        }
        torsion_table[key_torsion] = coef_torsion;
    }

    template<class Ttable_atom, class Ttable_res,
             class Ttable_bond, class Ttable_angle, class Ttable_torsion>
    void loading_param_file(const std::string    &model_name,
                            const std::string    &file_name,
                                  Ttable_atom    &atom_table,
                                  Ttable_res     &residue_table,
                                  Ttable_bond    &bond_table,
                                  Ttable_angle   &angle_table,
                                  Ttable_torsion &torsion_table){

        std::ifstream file_para{file_name};
        if(file_para.fail()) throw std::ios_base::failure("file: " + file_name + ".para was not found.");

        //--- loading para file
        PS::S32 atom_count    = 0;
        PS::S32 bond_count    = 0;
        PS::S32 angle_count   = 0;
        PS::S32 torsion_count = 0;

        PS::S32 n_atom    = -1;
        PS::S32 n_bond    = -1;
        PS::S32 n_angle   = -1;
        PS::S32 n_torsion = -1;
        MOL2_LOAD_MODE mode = MOL2_LOAD_MODE::header;

        std::string line;
        while (getline(file_para, line)){
            STR_TOOL::removeCR(line);
            std::vector<std::string> str_list = STR_TOOL::split(line, " ");

            //--- skip empty or comment line(start as "!" or "//")
            if(str_list.empty()) continue;
            if(str_list[0].empty()) continue;
            if(str_list[0].substr(0,1) == "!" ||
               str_list[0].substr(0,2) == "//") continue;

            //--- header information
            std::string mark     = "@<PARAM>";
            size_t      mark_len = mark.size();
            if( str_list[0].substr(0,mark_len) == mark ){
                mode = ENUM::which_MOL2_LOAD_MODE( str_list[0].substr(mark_len) );

                if(mode == MOL2_LOAD_MODE::atom   ) n_atom    = std::stoi(str_list.at(1));
                if(mode == MOL2_LOAD_MODE::bond   ) n_bond    = std::stoi(str_list.at(1));
                if(mode == MOL2_LOAD_MODE::angle  ) n_angle   = std::stoi(str_list.at(1));
                if(mode == MOL2_LOAD_MODE::torsion) n_torsion = std::stoi(str_list.at(1));

                continue;
            }

            //--- loading data
            switch (mode) {
                case MOL2_LOAD_MODE::header:
                    continue;
                break;

                case MOL2_LOAD_MODE::atom:
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;
                    loading_param_atom(model_name, str_list, atom_table, residue_table);
                    atom_count++;
                break;

                case MOL2_LOAD_MODE::bond:
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;
                    loading_param_bond(model_name, str_list, bond_table);
                    bond_count++;
                break;

                case MOL2_LOAD_MODE::angle:
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;
                    loading_param_angle(model_name, str_list, angle_table);
                    angle_count++;
                break;

                case MOL2_LOAD_MODE::torsion:
                    if( !STR_TOOL::isInteger(str_list[0]) ) continue;
                    loading_param_torsion(model_name, str_list, torsion_table);
                    torsion_count++;
                break;

                default:
                    throw std::invalid_argument("undefined loading mode.");
            }
        }

        //--- check file consistency
        if(n_atom    != atom_count    ||
           n_bond    != bond_count    ||
           n_angle   != angle_count   ||
           n_torsion != torsion_count ){
            std::cerr << "    n_atom    = " << n_atom    << " / elem_count    = " << atom_count    << "\n"
                      << "    n_bond    = " << n_bond    << " / bond_count    = " << bond_count    << "\n"
                      << "    n_angle   = " << n_angle   << " / angle_count   = " << angle_count   << "\n"
                      << "    n_torsion = " << n_torsion << " / torsion_count = " << torsion_count << std::endl;
            throw std::invalid_argument("number of defined parameters are not match with header.");
        }
    }

    template<class Tptcl, class Tid,
             class Ttable_atom, class Ttable_res,
             class Ttable_bond, class Ttable_angle, class Ttable_torsion>
    void loading_model_parameter(const std::string                   &model_name,
                                       std::vector<Tptcl>            &atom_list,
                                       std::vector<std::vector<Tid>> &bond_list,
                                       Ttable_atom                   &atom_table,
                                       Ttable_res                    &residue_table,
                                       Ttable_bond                   &bond_table,
                                       Ttable_angle                  &angle_table,
                                       Ttable_torsion                &torsion_table){

        std::string file_name;

        //--- loading ****.mol2 file
        if(model_name.find_first_of("/") != std::string::npos){
            file_name = model_name + ".mol2";
        } else {
            file_name = model_dir + model_name + ".mol2";
        }
        try{
            loading_mol2_file(model_name,
                              file_name,
                              atom_list,
                              bond_list );
        }
        catch(...){
            std::cerr << "ERROR at loading file: " << file_name << std::endl;
            throw;
        }

        //--- loading ****.param file
        if(model_name.find_first_of("/") != std::string::npos){
            file_name = model_name + ".param";
        } else {
            file_name = model_dir + model_name + ".param";
        }
        try{
            loading_param_file(model_name,
                               file_name,
                               atom_table,
                               residue_table,
                               bond_table,
                               angle_table,
                               torsion_table);
        }
        catch(...){
            std::cerr << "ERROR at loading file: " << file_name << std::endl;
            throw;
        }

        //--- copy VDW coef to atom_list from inter_table
        for(size_t index=0; index<atom_list.size(); index++){
            CoefAtom coef_atom;
            try{
                MODEL::KeyAtom key_atom = std::make_tuple( atom_list[index].getMolType(),
                                                           atom_list[index].getAtomType() );
                coef_atom = atom_table.at(key_atom);
            }
            catch(std::out_of_range err){
                std::cerr << "  element parameter of " << ENUM::what(atom_list[index].getAtomType())
                          << " was not found in '" << file_name
                          << "' file." << std::endl;
                throw;
            }
            atom_list[index].setMass(  coef_atom.mass );
            atom_list[index].setVDW_D( coef_atom.vdw_d );
            atom_list[index].setVDW_R( coef_atom.vdw_r );
        }

        //--- copy charge to inter_table from atom_list
        for(auto itr = atom_table.begin(); itr != atom_table.end(); itr++){
            PS::F32 charge_buff;
            bool find_flag = false;
            for(auto &atom : atom_list){
                if(std::make_tuple(atom.getMolType(), atom.getAtomType()) == itr->first){
                    #ifdef TEST_MOL_INSTALL
                        if(find_flag){
                            if(charge_buff != atom.getCharge()){
                                std::cerr << "WARNING: charge value of "
                                          << ENUM::what(itr->first)
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
                    std::cerr << "WARNING: " << ENUM::what(itr->first)
                              << " is not defined in "
                              << model_name + ".param"  << "\n"
                              << "  charge of " << ENUM::what(itr->first)
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

    //--- for Coef****::to_str() in md_coef_table.hpp.
    template<class Ttable>
    void print_coef_table(const Ttable &table, const MolName &model){
        std::string str;
        size_t count = 0;
        for(auto &coef : table){
            if( std::get<0>(coef.first) != model ) continue;

            ++count;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::what(coef.first) + "\n";
            str += coef.second.to_str(6);
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
        for(auto &coef : table){
            ++count;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::what(coef.first) + "\n";
            str += coef.second.to_str(6);
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }

    //--- overload for std::unordered_map<key__, std::string>
    template<class Tkey, class Hash, class Pred, class Allocator>
    void print_coef_table( const std::unordered_map<Tkey, std::string, Hash, Pred, Allocator> &table,
                           const MolName &model ){
        std::string str;
        size_t count = 0;
        for(auto &coef : table){
            if( std::get<0>(coef.first) != model ) continue;

            ++count;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::what(coef.first) + "\n";
            str += std::string(6, ' ') + coef.second;
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }

    template<class Tkey, class Hash, class Pred, class Allocator>
    void print_coef_table( const std::unordered_map<Tkey, std::string, Hash, Pred, Allocator> &table ){
        std::string str;
        size_t count = 0;
        for(auto &coef : table){
            ++count;
            str += "  count = " + std::to_string(count) + "\n";
            str += "    key : " + ENUM::what(coef.first) + "\n";
            str += std::string(6, ' ') + coef.second;
            str += "\n";
        }
        str += "  " + std::to_string(count) + " parameters were set.\n";
        str += "\n";

        std::cout << str;
    }
}
