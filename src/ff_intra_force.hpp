//***************************************************************************************
//  This program is the intramolecular force calculation for "md_fdps_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <tuple>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "md_coef_table.hpp"
#include "md_intra_pair_table.hpp"
#include "ff_intra_force_func.hpp"


template<class Tepi, class Tepj,
         class Ttable, class Tpair_list,
         class Tforce>
void calcForceIntra_IA(const Tepi       &epi,
                       const Tepj       *ep_j,
                       const Ttable     &ID_table,
                       const Tpair_list &intra_pair_list,
                             Tforce     &force_IA        ){

    //--- reference of intra pair table
    const auto& bond_list     = intra_pair_list.bond_list.at(epi.getAtomID());
    const auto& angle_list    = intra_pair_list.angle_list.at(epi.getAtomID());
    const auto& dihedral_list = intra_pair_list.dihedral_list.at(epi.getAtomID());
    const auto& improper_list = intra_pair_list.improper_list.at(epi.getAtomID());

    //--- bond potential
    for(size_t i=0; i<bond_list.size(); ++i){
        PS::S64 id_j    = bond_list.at(i);
        PS::S64 index_j = ID_table.find( id_j, epi.getPos() );

        MODEL::CoefBond bond_prm;
        MODEL::KeyBond  key_bond = std::make_tuple(epi.getMolType(),
                                                   epi.getAtomType(),
                                                   ep_j[index_j].getAtomType());
        try{
            bond_prm = MODEL::coefTable_bond.at(key_bond);
        }
        catch(std::out_of_range err){
            std::cerr << "key_bond = " << ENUM::what(key_bond) << " is not defined." << std::endl;
            throw;
        }

        switch (bond_prm.form) {
            case IntraFuncForm::none:
                //--- free potential
                continue;
            break;

            case IntraFuncForm::anharmonic:
                calcBondForce_anharmonic_IJ(Normalize::realPos(epi.getPos()),
                                            Normalize::realPos(ep_j[index_j].getPos()),
                                            bond_prm,
                                            force_IA);
            break;

            case IntraFuncForm::harmonic:
                calcBondForce_harmonic_IJ(Normalize::realPos(epi.getPos()),
                                          Normalize::realPos(ep_j[index_j].getPos()),
                                          bond_prm,
                                          force_IA);
            break;

            default:
                std::cerr << "  IntraFuncForm = " << ENUM::what(bond_prm.form) << std::endl;
                throw std::invalid_argument("undefined function form: bond force");
        }
    }

    //--- angle potential
    for(size_t i=0; i<angle_list.size(); ++i){
        PS::S64 id_i   = std::get<0>( angle_list.at(i) );
        PS::S64 id_j   = std::get<1>( angle_list.at(i) );
        PS::S64 id_k   = std::get<2>( angle_list.at(i) );
        PS::S64 id_tgt = epi.getAtomID();

        PS::S64 index_i = ID_table.find(id_i, epi.getPos());
        PS::S64 index_j = ID_table.find(id_j, epi.getPos());
        PS::S64 index_k = ID_table.find(id_k, epi.getPos());

        MODEL::CoefAngle angle_prm;
        MODEL::KeyAngle  key_angle = std::make_tuple(epi.getMolType(),
                                                     ep_j[index_i].getAtomType(),
                                                     ep_j[index_j].getAtomType(),
                                                     ep_j[index_k].getAtomType() );
        try{
            angle_prm = MODEL::coefTable_angle.at(key_angle);
        }
        catch(std::out_of_range err){
            std::cerr << "key_angle = " << ENUM::what(key_angle) << " is not defined." << std::endl;
            throw;
        }

        switch (angle_prm.form){
            case IntraFuncForm::none:
                //--- free potential
                continue;
            break;

            case IntraFuncForm::harmonic:
                calcAngleForce_harmonic_IJK(Normalize::realPos( ep_j[index_i].getPos() ),
                                            Normalize::realPos( ep_j[index_j].getPos() ),
                                            Normalize::realPos( ep_j[index_k].getPos() ),
                                            id_i, id_j, id_k,
                                            id_tgt, angle_prm,
                                            force_IA);
            break;

            default:
                std::cerr << "  IntraFuncForm = " << ENUM::what(angle_prm.form) << std::endl;
                throw std::invalid_argument("undefined function form: angle force");
        }
    }

    //--- dihedral torsion potential
    for(size_t i=0; i<dihedral_list.size(); ++i){
        PS::S64 id_i   = std::get<0>( dihedral_list.at(i) );
        PS::S64 id_j   = std::get<1>( dihedral_list.at(i) );
        PS::S64 id_k   = std::get<2>( dihedral_list.at(i) );
        PS::S64 id_l   = std::get<3>( dihedral_list.at(i) );
        PS::S64 id_tgt = epi.getAtomID();

        PS::S64 index_i = ID_table.find(id_i, epi.getPos());
        PS::S64 index_j = ID_table.find(id_j, epi.getPos());
        PS::S64 index_k = ID_table.find(id_k, epi.getPos());
        PS::S64 index_l = ID_table.find(id_l, epi.getPos());

        MODEL::CoefTorsion torsion_prm;
        MODEL::KeyTorsion  key_torsion = std::make_tuple(epi.getMolType(),
                                                        TorsionShape::dihedral,
                                                        ep_j[index_i].getAtomType(),
                                                        ep_j[index_j].getAtomType(),
                                                        ep_j[index_k].getAtomType(),
                                                        ep_j[index_l].getAtomType() );
        try{
            torsion_prm = MODEL::coefTable_torsion.at(key_torsion);
        }
        catch(std::out_of_range err){
            std::cerr << "key_torsion = " << ENUM::what(key_torsion) << " is not defined." << std::endl;
            throw;
        }

        switch (torsion_prm.form){
            case IntraFuncForm::none:
                //--- free potential
                continue;
            break;

            case IntraFuncForm::cos:
                calcTorsionForce_harmonic_IJKL(Normalize::realPos( ep_j[index_i].getPos() ),
                                               Normalize::realPos( ep_j[index_j].getPos() ),
                                               Normalize::realPos( ep_j[index_k].getPos() ),
                                               Normalize::realPos( ep_j[index_l].getPos() ),
                                               id_i, id_j, id_k, id_l,
                                               id_tgt, torsion_prm,
                                               force_IA);
            break;

            default:
                std::cerr << "  IntraFuncForm = " << ENUM::what(torsion_prm.form) << std::endl;
                throw std::invalid_argument("undefined potential type: dihedral torsion force");
        }
    }

    //--- improper torsion potential
    for(size_t i=0; i<improper_list.size(); ++i){
        PS::S64 id_i   = std::get<0>( improper_list.at(i) );
        PS::S64 id_j   = std::get<1>( improper_list.at(i) );
        PS::S64 id_k   = std::get<2>( improper_list.at(i) );
        PS::S64 id_l   = std::get<3>( improper_list.at(i) );
        PS::S64 id_tgt = epi.getAtomID();

        PS::S64 index_i = ID_table.find(id_i, epi.getPos());
        PS::S64 index_j = ID_table.find(id_j, epi.getPos());
        PS::S64 index_k = ID_table.find(id_k, epi.getPos());
        PS::S64 index_l = ID_table.find(id_l, epi.getPos());

        MODEL::CoefTorsion torsion_prm;
        MODEL::KeyTorsion  key_torsion = std::make_tuple(epi.getMolType(),
                                                         TorsionShape::improper,
                                                         ep_j[index_i].getAtomType(),
                                                         ep_j[index_j].getAtomType(),
                                                         ep_j[index_k].getAtomType(),
                                                         ep_j[index_l].getAtomType() );
        try{
            torsion_prm = MODEL::coefTable_torsion.at(key_torsion);
        }
        catch(std::out_of_range){
            std::cerr << "key_torsion = " << ENUM::what(key_torsion) << " is not defined." << std::endl;
            throw;
        }

        switch (torsion_prm.form){
            case IntraFuncForm::none:
                //--- free potential
                continue;
            break;

            case IntraFuncForm::cos:
                calcTorsionForce_harmonic_IJKL(Normalize::realPos( ep_j[index_i].getPos() ),
                                               Normalize::realPos( ep_j[index_j].getPos() ),
                                               Normalize::realPos( ep_j[index_k].getPos() ),
                                               Normalize::realPos( ep_j[index_l].getPos() ),
                                               id_i, id_j, id_k, id_l,
                                               id_tgt, torsion_prm,
                                               force_IA);
            break;

            default:
                std::cerr << "  IntraFuncForm = " << ENUM::what(torsion_prm.form) << std::endl;
                throw std::invalid_argument("undefined potential type: improper torsion force");
        }
    }
}


//--- template function for FDPS interface
template<class Tepj, class Tforce, class Ttree, class Tpsys>
void calcForceIntra(Ttree &tree,
                    Tpsys &psys){

    PS::S64 n_local = psys.getNumberOfParticleLocal();

    //--- pair atom search table
    MD_EXT::ID_Table<PS::S64, PS::S64> ID_table;
    ID_table.reserve(n_local);

    //--- calculate intramolecular force
    for(PS::S64 i=0; i<n_local; ++i){
        ID_table.clear();

        //--- skip no-bonded atoms
        if( MODEL::intra_pair_manager.bond_list.at( psys[i].getAtomID() ).size() == 0 ) continue;

        //--- make ID table for psys[i] neigbor
        Tepj*   ep_j;
        PS::S64 n_neighbor = tree.getNeighborListOneParticle(psys[i], ep_j);
        for(PS::S64 j=0; j<n_neighbor; ++j){
            ID_table.add( ep_j[j].getAtomID(), j, ep_j[j].getPos() );   // ID, array_index, pos
        }

        Tforce force_IA;
               force_IA.clear();
        calcForceIntra_IA(psys[i],
                          ep_j,
                          ID_table,
                          MODEL::intra_pair_manager,
                          force_IA                  );
        psys[i].copyForceIntra(force_IA);
    }
}
