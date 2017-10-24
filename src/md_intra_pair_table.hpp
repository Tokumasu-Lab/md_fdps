//***************************************************************************************
//  This is the pair-ID data table class for calculate interactions.
//***************************************************************************************
#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "md_enum.hpp"
#include "atom_class_base.hpp"


//--- grobal object of parameter table
namespace MODEL {

    //--- intramolecular force parameter manager for PS::ParticleSystem<FP>
    template <typename Tid>
    class IntrarMolecularForcePair_manager {
    private:
        struct Change {
            Tid  i, j;
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
        std::vector<TorsionList> dihedral_list;
        std::vector<TorsionList> improper_list;

        //--- functions
        void setAtomNumber(const Tid &n) {
            this->intra_mask.resize(n);
            this->bond_list.resize(n);
            this->angle_list.resize(n);
            this->dihedral_list.resize(n);
            this->improper_list.resize(n);
        }
        void clear() {
            this->change_journal.clear();

            this->intra_mask.clear();
            this->bond_list.clear();
            this->angle_list.clear();
            this->dihedral_list.clear();
            this->improper_list.clear();
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
            throw std::logic_error("underdevelop !");
        }

        void broadcast(const PS::S32 root = 0){
            COMM_TOOL::broadcast(this->bond_list, root);

            this->intra_mask.resize(this->bond_list.size());
            this->angle_list.resize(this->bond_list.size());
            this->dihedral_list.resize(this->bond_list.size());
            this->improper_list.resize(this->bond_list.size());

            this->makeIntraList();

            this->change_journal.clear();
        }

        void sync(){
            throw std::logic_error("underdevelop !");
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
            this->dihedral_list.at(atom_id).clear();
            this->improper_list.at(atom_id).clear();

            //--- construct lists    source: bond_list
            std::vector<Tid> node_list;
            node_list.push_back(atom_id);  // set root node

            this->make_intra_mask(node_list);
            std::sort(this->intra_mask.at(atom_id).begin(), this->intra_mask.at(atom_id).end());

            this->make_angle_list(atom_id);
            this->make_torsion_list(atom_id);
        }

        void makeIntraList(){
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

                //--- dihedral, I-J-K-L
                for(size_t l=0; l<this->bond_list.at(node_k).size(); ++l){
                    Tid node_l = this->bond_list.at(node_k).at(l);

                    assert(node_l != node_k);
                    if(node_l == node_j ||
                       node_l == root     ) continue;

                    this->dihedral_list.at(root).emplace_back( std::make_tuple(root,
                                                                               node_j,
                                                                               node_k,
                                                                               node_l ) );
                }

                //--- improper, I-J<KL
                for(size_t l=k+1; l<this->bond_list.at(node_j).size(); ++l){
                    Tid node_l = this->bond_list.at(node_j).at(l);

                    assert(node_l != node_j);
                    if(node_l == node_k ||
                       node_l == root     ) continue;

                    this->improper_list.at(root).emplace_back( std::make_tuple(root,
                                                                               node_j,
                                                                               node_k,
                                                                               node_l ) );
                }

                //--- dihedral, L-I-J-K
                for(size_t l=0; l<this->bond_list.at(root).size(); ++l){
                    Tid node_l = this->bond_list.at(root).at(l);

                    assert(node_l != root);
                    if(node_l == node_j ||
                       node_l == node_k   ) continue;

                    this->dihedral_list.at(root).emplace_back( std::make_tuple(node_l,
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

                    this->improper_list.at(root).emplace_back( std::make_tuple(node_j,
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
            oss << " " << ENUM::what(tuple);
        }
        oss << "\n";

        oss << "  dohedral torsion =";
        for(auto tuple : MODEL::intra_pair_manager.dihedral_list.at(id)){
            oss << " " << ENUM::what(tuple);
        }
        oss << "\n";

        oss << "  improper torsion =";
        for(auto tuple : MODEL::intra_pair_manager.improper_list.at(id)){
            oss << " " << ENUM::what(tuple);
        }
        oss << "\n";

        std::cout << oss.str() << std::flush;
    }
}
