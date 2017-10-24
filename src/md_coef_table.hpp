//***************************************************************************************
//  This is the coefficient data table class for calculate interactions.
//***************************************************************************************
#pragma once

#include <string>
#include <sstream>
#include <unordered_map>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "atom_class_base.hpp"
#include "md_enum.hpp"

//--- grobal object of parameter table
namespace MODEL {

    //--- parameter for intermolecular interactions
    struct CoefAtom {
      public:
        PS::F32 mass;
        PS::F32 charge;
        PS::F32 vdw_d;
        PS::F32 vdw_r;

        inline std::string to_str(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "mass   : " << std::setprecision(8) << this->mass   << "\n";
            oss << std::setw(shift + 9) << "charge : " << std::setprecision(8) << this->charge << "\n";
            oss << std::setw(shift + 9) << "vdw_d  : " << std::setprecision(8) << this->vdw_d  << "\n";
            oss << std::setw(shift + 9) << "vdw_r  : " << std::setprecision(8) << this->vdw_r  << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_str(shift);
        }
    };

    //--- parameter for intramolecular interactions
    struct CoefBond {
      public:
        IntraFuncForm form;
        PS::F32       r0, k, a;

        inline std::string to_str(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 7) << "form : "                         << this->form << "\n";
            oss << std::setw(shift + 7) << "r0   : " << std::setprecision(8) << this->r0   << "\n";
            oss << std::setw(shift + 7) << "k    : " << std::setprecision(8) << this->k    << "\n";
            oss << std::setw(shift + 7) << "a    : " << std::setprecision(8) << this->a    << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_str(shift);
        }
    };

    struct CoefAngle {
      public:
        IntraFuncForm form;
        PS::F32       theta0, k;

        inline std::string to_str(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "form   : "                         << this->form   << "\n";
            oss << std::setw(shift + 9) << "k      : " << std::setprecision(8) << this->k      << "\n";
            oss << std::setw(shift + 9) << "theta0 : " << std::setprecision(8) << this->theta0 << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_str();
        }
    };

    struct CoefTorsion {
      public:
        IntraFuncForm form;
        PS::S32       n_min;
        PS::F32       k, theta0;

        inline std::string to_str(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "form   : "                         << this->form   << "\n";
            oss << std::setw(shift + 9) << "n_min  : " << std::setprecision(8) << this->n_min  << "\n";
            oss << std::setw(shift + 9) << "k      : " << std::setprecision(8) << this->k      << "\n";
            oss << std::setw(shift + 9) << "theta0 : " << std::setprecision(8) << this->theta0 << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_str();
        }
    };

    //--- for intermolecular parameter
    using KeyAtom = std::tuple<MolName,
                               AtomName>;
    std::unordered_map< KeyAtom,
                        CoefAtom,
                        hash_tuple::hash_func<KeyAtom>> coefTable_atom;

    //--- for residue information
    std::unordered_map< KeyAtom,
                        std::string,
                        hash_tuple::hash_func<KeyAtom>> coefTable_residue;

    //--- for intramolecular parameter
    //------ key = (model_name, i_atom, j_atom)
    using KeyBond = std::tuple<MolName,
                               AtomName,
                               AtomName>;
    std::unordered_map< KeyBond,
                        CoefBond,
                        hash_tuple::hash_func<KeyBond>> coefTable_bond;

    //------ key = (model_name, i_atom, j_atom, k_atom)
    using KeyAngle = std::tuple<MolName,
                                AtomName,
                                AtomName,
                                AtomName>;
    std::unordered_map< KeyAngle,
                        CoefAngle,
                        hash_tuple::hash_func<KeyAngle>> coefTable_angle;

    //------ key = (model_name, torsion_shape, i_atom, j_atom, k_atom, l_atom)
    using KeyTorsion = std::tuple<MolName,
                                  TorsionShape,
                                  AtomName,
                                  AtomName,
                                  AtomName,
                                  AtomName    >;
    std::unordered_map< KeyTorsion,
                        CoefTorsion,
                        hash_tuple::hash_func<KeyTorsion>> coefTable_torsion;


    //--- broadcast parameter table from rank=0 to other process
    void broadcast_coefTable(const PS::S32 root = 0){
        COMM_TOOL::broadcast(coefTable_atom,    root);
        COMM_TOOL::broadcast(coefTable_residue, root);
        COMM_TOOL::broadcast(coefTable_bond,    root);
        COMM_TOOL::broadcast(coefTable_angle,   root);
        COMM_TOOL::broadcast(coefTable_torsion, root);
    }

}
