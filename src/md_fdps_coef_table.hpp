//***************************************************************************************
//  This program is the temporary buffer data class for intramolecular force calculation.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <string>
#include <sstream>
#include <unordered_map>

#include <particle_simulator.hpp>

#include "hash_tuple.hpp"
#include "comm_tool.hpp"

#include "md_fdps_atom_class_base.hpp"


//--- enum indicator for intramolecular interaction
enum class IntraFuncForm : int {
    harmonic,
    anharmonic,

    cos,
    none,
};

namespace ENUM {

    std::string whatis(const IntraFuncForm &e){
        switch (e) {
            case IntraFuncForm::harmonic:
                return "harmonic";
            break;

            case IntraFuncForm::anharmonic:
                return "anharmonic";
            break;

            case IntraFuncForm::cos:
                return "cos";
            break;

            case IntraFuncForm::none:
                return "none";
            break;

            default:
                throw std::out_of_range("undefined enum value in IntraFuncForm.");
        }
    }

    IntraFuncForm which_IntraFuncForm(const std::string &str){
        if(        str == "harmonic" ){
            return IntraFuncForm::harmonic;
        } else if( str == "anharmonic" ){
            return IntraFuncForm::anharmonic;
        } else if( str == "cos" ){
            return IntraFuncForm::cos;
        } else if( str == "none" ){
            return IntraFuncForm::none;
        } else {
            std::cerr << "  IntraFuncForm: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in IntraFuncForm.");
        }
    }
}

//--- output function as "std::cout << (enum class::value)"
inline std::ostream& operator << (std::ostream& s, const IntraFuncForm &e){
    s << ENUM::whatis(e);
    return s;
}


//--- grobal object of parameter table
namespace MODEL {

    //--- parameter for intermolecular interactions
    struct CoefElement {
      public:
        PS::F32 mass;
        PS::F32 charge;
        PS::F32 vdw_d;
        PS::F32 vdw_r;

        inline std::string to_string(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "mass   : " << std::setprecision(8) << this->mass   << "\n";
            oss << std::setw(shift + 9) << "charge : " << std::setprecision(8) << this->charge << "\n";
            oss << std::setw(shift + 9) << "vdw_d  : " << std::setprecision(8) << this->vdw_d  << "\n";
            oss << std::setw(shift + 9) << "vdw_r  : " << std::setprecision(8) << this->vdw_r  << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_string(shift);
        }
    };

    //--- parameter for intramolecular interactions
    struct CoefBond {
      public:
        IntraFuncForm form;
        PS::F32 a, d, r0;

        inline std::string to_string(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 7) << "form : "                         << this->form << "\n";
            oss << std::setw(shift + 7) << "d    : " << std::setprecision(8) << this->d    << "\n";
            oss << std::setw(shift + 7) << "a    : " << std::setprecision(8) << this->a    << "\n";
            oss << std::setw(shift + 7) << "r0   : " << std::setprecision(8) << this->r0   << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_string(shift);
        }
    };

    struct CoefAngle {
      public:
        IntraFuncForm form;
        PS::F32 theta0, k;

        inline std::string to_string(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "form   : "                         << this->form   << "\n";
            oss << std::setw(shift + 9) << "k      : " << std::setprecision(8) << this->k      << "\n";
            oss << std::setw(shift + 9) << "theta0 : " << std::setprecision(8) << this->theta0 << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_string();
        }
    };

    struct CoefTorsion {
      public:
        IntraFuncForm form;
        PS::S32 n_min;
        PS::F32 k, theta0;

        inline std::string to_string(const size_t &shift = 0) const {
            std::ostringstream oss;

            oss << std::setw(shift + 9) << "form   : "                         << this->form   << "\n";
            oss << std::setw(shift + 9) << "n_min  : " << std::setprecision(8) << this->n_min  << "\n";
            oss << std::setw(shift + 9) << "k      : " << std::setprecision(8) << this->k      << "\n";
            oss << std::setw(shift + 9) << "theta0 : " << std::setprecision(8) << this->theta0 << "\n";

            return oss.str();
        }

        inline void print(const size_t &shift = 0) const {
            std::cout << this->to_string();
        }
    };

    //--- for intermolecular parameter
    using KeyElem = std::tuple<MolName,
                               AtomName>;
    std::unordered_map< KeyElem,
                        CoefElement,
                        hash_tuple::hash_func<KeyElem>> coefTable_elem;

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

    //------ key = (model_name, i_atom, j_atom, k_atom, l_atom)
    using KeyTorsion = std::tuple<MolName,
                                  AtomName,
                                  AtomName,
                                  AtomName,
                                  AtomName>;
    std::unordered_map< KeyTorsion,
                        CoefTorsion,
                        hash_tuple::hash_func<KeyTorsion>> coefTable_torsion;

    //--- broadcast parameter table from rank=0 to other process
    void broadcast_coefTable(const PS::S32 root = 0){
        COMM_TOOL::broadcast(coefTable_elem,    root);
        COMM_TOOL::broadcast(coefTable_bond,    root);
        COMM_TOOL::broadcast(coefTable_angle,   root);
        COMM_TOOL::broadcast(coefTable_torsion, root);
    }



    //--- intramolecular force parameter manager for PS::ParticleSystem<FP>
    template <typename Tid>
    class IntrarMolecularForcePair_manager {
    private:
        struct Change {
            Tid  i;
            Tid  j;
            bool add;
        };
        std::vector<Change> change_journal;

        void make_intra_mask(  std::vector<Tid> &node_list, const PS::S32 order = 3);
        void make_angle_list(  const Tid &atom_id);
        void make_torsion_list(const Tid &atom_id);

    public:
        using MaskList    = typename std::vector<Tid>;
        using BondList    = typename std::vector<Tid>;
        using AngleList   = typename std::vector<std::tuple<Tid, Tid, Tid>>;
        using TorsionList = typename std::vector<std::tuple<Tid, Tid, Tid, Tid>>;

        std::vector<MaskList>    intra_mask;
        std::vector<BondList>    bond_list;
        std::vector<AngleList>   angle_list;
        std::vector<TorsionList> torsion_list;

        //--- functions
        void setAtomNumber(const Tid &n) {
            this->intra_mask.resize(n);
            this->bond_list.resize(n);
            this->angle_list.resize(n);
            this->torsion_list.resize(n);
        }
        void clear() {
            this->change_journal.clear();

            this->intra_mask.clear();
            this->bond_list.clear();
            this->angle_list.clear();
            this->torsion_list.clear();
        }
        void addBond(const Tid &i, const Tid &j){
            assert(i != j);
            if( std::find(this->bond_list.at(i).begin(), this->bond_list.at(i).end(), j) == this->bond_list.at(i).end()){
                this->bond_list.at(i).emplace_back(j);
            }
            if( std::find(this->bond_list.at(j).begin(), this->bond_list.at(j).end(), i) == this->bond_list.at(j).end()){
                this->bond_list.at(j).emplace_back(i);
            }
        }
        void eraseConnect(const Tid &i, const Tid &j){
            assert(i != j);
            // underdevelop
        }

        void broadcast(const PS::S32 root = 0){
            COMM_TOOL::broadcast(this->bond_list, root);

            this->intra_mask.resize(this->bond_list.size());
            this->angle_list.resize(this->bond_list.size());
            this->torsion_list.resize(this->bond_list.size());

            this->makeIntraList_all();

            this->change_journal.clear();
        }

        void sync(){
            // accumulate this->change_journal
            // then update this->bond_list in each process
            // call this->makeIntraList(id)
            //    at contains this->intra_mask at connection changed atom id.

            // finally this->change_journal.clear();
        }

        void makeIntraList(const Tid &atom_id){
            //--- initialize
            this->intra_mask.at(atom_id).clear();

            this->angle_list.at(atom_id).clear();
            this->torsion_list.at(atom_id).clear();

            //--- construct lists    source: bond_list
            std::vector<Tid> node_list;
            node_list.push_back(atom_id);  // set root node

            this->make_intra_mask(node_list);
            std::sort(this->intra_mask.at(atom_id).begin(), this->intra_mask.at(atom_id).end());

            this->make_angle_list(atom_id);
            this->make_torsion_list(atom_id);
        }

        void makeIntraList_all(){
            for(size_t i=0; i<this->bond_list.size(); ++i){
                this->makeIntraList(i);
            }
        }

        bool checkIntraMask(const Tid &i, const Tid &j){
            //--- general design
            return ( std::find(this->intra_mask.at(i).begin(),
                               this->intra_mask.at(i).end(),
                               j)
                     == this->intra_mask.at(i).end() );

            //--- "intra_mask" must be sorted.
        //    return !std::binary_search(this->intra_mask.at(i).begin(),
        //                               this->intra_mask.at(i).end(),
        //                               j);
        }
    };

    template <typename Tid>
    void IntrarMolecularForcePair_manager<Tid>::make_intra_mask(std::vector<Tid> &node_list,
                                                                const PS::S32     order){
        //--- search in bond_list.at(atom_id)
        Tid root = node_list.front();
        Tid node = node_list.back();
        for(size_t j=0; j<this->bond_list.at(node).size(); ++j){
            Tid next_node = this->bond_list.at(node).at(j);

            //--- search in next node
            if(order > 1){
                node_list.push_back(next_node);
                this->make_intra_mask(node_list,
                                      order-1);
                node_list.pop_back();
            }

            //--- duplication check
            if( next_node == root) continue;
            if( std::find(this->intra_mask.at(root).begin(),
                          this->intra_mask.at(root).end(),
                          next_node)
                != this->intra_mask.at(root).end() ) continue;

            //--- add mask_list
            this->intra_mask.at(root).emplace_back(next_node);
        }
    }

    template <typename Tid>
    void IntrarMolecularForcePair_manager<Tid>::make_angle_list(const Tid &root){
        //--- I-J
        for(size_t j=0; j<this->bond_list.at(root).size(); ++j){
            Tid node_j = this->bond_list.at(root).at(j);

            //std::cout << "node_j = " << node_j << "  , root = " << root << std::endl;
            assert(node_j != root);

            //--- I-J-K
            for(size_t k=0; k<this->bond_list.at(node_j).size(); ++k){
                Tid node_k = this->bond_list.at(node_j).at(k);

                assert(node_k != node_j);
                if(node_k == root) continue;

                this->angle_list.at(root).emplace_back( std::make_tuple(root,
                                                                        node_j,
                                                                        node_k) );
            }

            //--- J-I-K
            for(size_t k=j+1; k<this->bond_list.at(root).size(); ++k){
                Tid node_k = this->bond_list.at(root).at(k);

                assert(node_k != root);
                if(node_k == node_j) continue;

                this->angle_list.at(root).emplace_back( std::make_tuple(node_j,
                                                                        root,
                                                                        node_k) );
            }
        }
    }

    template <typename Tid>
    void IntrarMolecularForcePair_manager<Tid>::make_torsion_list(const Tid &root){
        //--- I-J
        for(size_t j=0; j<this->bond_list.at(root).size(); ++j){
            Tid node_j = this->bond_list.at(root).at(j);

            assert(node_j != root);

            //--- I-J-K
            for(size_t k=0; k<this->bond_list.at(node_j).size(); ++k){
                Tid node_k = this->bond_list.at(node_j).at(k);

                assert(node_k != node_j);
                if(node_k == root) continue;

                //--- I-J-K-L
                for(size_t l=0; l<this->bond_list.at(node_k).size(); ++l){
                    Tid node_l = this->bond_list.at(node_k).at(l);

                    assert(node_l != node_k);
                    if(node_l == node_j ||
                       node_l == root     ) continue;

                    this->torsion_list.at(root).emplace_back( std::make_tuple(root,
                                                                              node_j,
                                                                              node_k,
                                                                              node_l ) );
                }

                //--- I-J<KL
                for(size_t l=k+1; l<this->bond_list.at(node_j).size(); ++l){
                    Tid node_l = this->bond_list.at(node_j).at(l);

                    assert(node_l != node_j);
                    if(node_l == node_k ||
                       node_l == root     ) continue;

                    this->torsion_list.at(root).emplace_back( std::make_tuple(root,
                                                                              node_j,
                                                                              node_k,
                                                                              node_l ) );
                }

                //--- L-I-J-K
                for(size_t l=0; l<this->bond_list.at(root).size(); ++l){
                    Tid node_l = this->bond_list.at(root).at(l);

                    assert(node_l != root);
                    if(node_l == node_j ||
                       node_l == node_k   ) continue;

                    this->torsion_list.at(root).emplace_back( std::make_tuple(node_l,
                                                                              root,
                                                                              node_j,
                                                                              node_k ) );
                }
            }

            //--- improper, I center
            for(size_t k=j+1; k<this->bond_list.at(root).size(); ++k){
                for(size_t l=k+1; l<this->bond_list.at(root).size(); ++l){
                    Tid node_k = this->bond_list.at(root).at(k);
                    Tid node_l = this->bond_list.at(root).at(l);

                    assert(node_k != root);
                    assert(node_l != root);
                    if(node_l == node_j) continue;
                    if(node_l == node_k) continue;

                    this->torsion_list.at(root).emplace_back( std::make_tuple(node_j,
                                                                              root,
                                                                              node_k,
                                                                              node_l ) );
                }
            }
        }
    }

    //--- global object
    IntrarMolecularForcePair_manager<PS::S64> intra_pair_manager;


    //--- display interface
    template <typename Tid>
    void print_connection(const Tid &id){
        std::ostringstream oss;

        oss << " ID = " << id << "\n";

        oss << "  bond =";
        for(auto j : MODEL::intra_pair_manager.bond_list.at(id)){
            oss << " " << j;
        }
        oss << "\n";

        oss << "  mask =";
        for(auto j : MODEL::intra_pair_manager.intra_mask.at(id)){
            oss << " " << j;
        }
        oss << "  / mask_len = " << MODEL::intra_pair_manager.intra_mask.at(id).size() << "\n";

        oss << "  angle =";
        for(auto tuple : MODEL::intra_pair_manager.angle_list.at(id)){
            oss << " " << ENUM::whatis(tuple);
        }
        oss << "\n";

        oss << "  torsion =";
        for(auto tuple : MODEL::intra_pair_manager.torsion_list.at(id)){
            oss << " " << ENUM::whatis(tuple);
        }
        oss << "\n";

        std::cout << oss.str() << std::flush;
    }
}
