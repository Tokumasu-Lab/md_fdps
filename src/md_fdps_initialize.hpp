//***************************************************************************************
//  This routine is initializer for "md_fdps_main.cpp".
//***************************************************************************************
#pragma once

#include <cmath>
#include <algorithm>
#include <fstream>
#include <random>

#include <particle_simulator.hpp>

#include "boltzmann_dist.hpp"

#include "md_fdps_unit.hpp"
#include "md_fdps_coef_table.hpp"
#include "md_fdps_atom_class.hpp"
#include "md_fdps_loading_condition.hpp"
#include "md_fdps_loading_model.hpp"

namespace Initialize {

    template <typename Tvec>
    Tvec rot_x(const Tvec &v, const PS::F64 &r) {
        Tvec result;
        PS::F64 sin_tmp = std::sin(r);
        PS::F64 cos_tmp = std::cos(r);
        result.x = v.x;
        result.y = cos_tmp*v.y - sin_tmp*v.z;
        result.z = sin_tmp*v.y + cos_tmp*v.z;
        return result;
    }

    template <typename Tvec>
    Tvec rot_y(const Tvec &v, const PS::F64 &r) {
        Tvec result;
        PS::F64 sin_tmp = std::sin(r);
        PS::F64 cos_tmp = std::cos(r);
        result.x = sin_tmp*v.z + cos_tmp*v.x;
        result.y = v.y;
        result.z = cos_tmp*v.z - sin_tmp*v.x;
        return result;
    }

    template <typename Tvec>
    Tvec rot_z(const Tvec &v, const PS::F64 &r) {
        Tvec result;
        PS::F64 sin_tmp = std::sin(r);
        PS::F64 cos_tmp = std::cos(r);
        result.x = cos_tmp*v.x - sin_tmp*v.y;
        result.y = sin_tmp*v.x + cos_tmp*v.y;
        result.z = v.z;
        return result;
    }

    //--- adjast in [0.0,1.0) space
    template <typename Tvec>
    Tvec periodicAdjustNorm_round(const Tvec &pos_norm){
        Tvec pos_new;
        pos_new.x = pos_norm.x - std::round(pos_norm.x) + 0.5;
        pos_new.y = pos_norm.y - std::round(pos_norm.y) + 0.5;
        pos_new.z = pos_norm.z - std::round(pos_norm.z) + 0.5;
        return pos_new;
    }

    template <typename Tvec>
    Tvec periodicAdjustReal_round(const Tvec &pos_real){
        Tvec pos_new;
        pos_new = Normalize::normPos(pos_new);
        pos_new = periodicAdjustNorm_round(pos_new);
        pos_new = Normalize::realPos(pos_new);
        return pos_new;
    }


    template <typename Tpsys, typename Tmask>
    bool check_excluded_vol(const PS::F64vec &pos_new,
                            const Tpsys      &psys,
                            const PS::F64    &r_ex,
                            const PS::S64    &st,
                            const PS::S64    &end,
                            const Tmask      &mask ){
        bool    flag  = true;
        PS::F64 r2_ex = r_ex*r_ex;

        for(PS::S64 i=st; i<end; ++i){
            //--- check mask
            if( std::find(mask.begin(), mask.end(), psys[i].getAtomID() ) != mask.end() ) continue;

            PS::F64vec r_vec = psys[i].getPos() - pos_new;
                       r_vec = Normalize::periodicAdjustNorm(r_vec);
            PS::F64    r2    = r_vec*r_vec;

            if(r2 <= r2_ex){
                flag = false;
                break;
            }
        }

        return flag;
    }

    template <typename Trand, typename Tdist, typename Tblz>
    PS::F64vec calc_vel(const PS::F64 &temperature,
                        const PS::F64 &mass,
                              Trand   &rand,
                              Tdist   &dist,
                              Tblz    &blz         ){

        //--- intensity
        PS::F64    dev = std::sqrt( 2.0*(temperature/Unit::norm_temp)/mass );
        PS::F64vec v   = 0.0;

        v.x = dev*blz.gen( PS::F64(dist(rand)) );

        //--- direction
        v = rot_z(v, Unit::pi*PS::F64(dist(rand)) );  // rotate in z-y-z Euler angle
        v = rot_y(v, Unit::pi*PS::F64(dist(rand)) );
        v = rot_z(v, Unit::pi*PS::F64(dist(rand)) );

        return v;
    }

    template <typename Tpsys, typename Tintralist, typename Tmodel, typename Tbond,
        //      typename Tchecker,
              typename Trand, typename Tdist,      typename Tblz>
    size_t install_molecule(      Tpsys      &psys,
                                  Tintralist &intra_pair_manager,
                            const Tmodel     &model_template,
                            const Tbond      &bond_template,
                            const PS::F64    &ex_r,
                            const size_t     &try_limit,
                            const PS::F64    &temperature,
                                  Trand      &rand,
                                  Tdist      &dist,
                                  Tblz       &blz,
                                  PS::S64    &atom_inclement,
                                  PS::S64    &mol_inclement){

        size_t try_count = 0;

        //--- set molecule properties
        for(size_t i=0; i<model_template.size(); ++i){
            PS::S64 atom_id = atom_inclement + i;
            psys[atom_id].copyFromModelTemplate(atom_id,
                                                mol_inclement,
                                                model_template.at(i) );

            //--- copy and shift bond target ID
            for(auto b : bond_template.at(i)){
                intra_pair_manager.bond_list.at(atom_id).push_back(atom_inclement + b);
            }
        }
        for(size_t i=0; i<model_template.size(); ++i){
            intra_pair_manager.makeIntraList(atom_inclement + i);
        }

        //--- setting position
        for(; try_count<=try_limit; ++try_count){
            PS::F64vec pos_root = { dist(rand), dist(rand), dist(rand) };  // normalized
            PS::F64    rot_1    = Unit::pi*PS::F64(dist(rand));
            PS::F64    rot_2    = Unit::pi*PS::F64(dist(rand));
            PS::F64    rot_3    = Unit::pi*PS::F64(dist(rand));

            //--- set new position of atom (root pos + local_pos)
            bool pos_flag = true;
            for(size_t i=0; i<model_template.size(); ++i){
                PS::F64vec pos_local = Normalize::normPos( model_template.at(i).getPos() );

                pos_local = rot_z(pos_local, rot_1);  // rotate in z-y-z Euler angle
                pos_local = rot_y(pos_local, rot_2);
                pos_local = rot_z(pos_local, rot_3);

                PS::F64vec pos_norm = pos_root + pos_local;
                pos_norm = periodicAdjustNorm_round(pos_norm);

                if( !check_excluded_vol(pos_norm,
                                        psys,
                                        ex_r,
                                        atom_inclement, atom_inclement + i,
                                        intra_pair_manager.intra_mask.at(atom_inclement + i) ) ){
                    pos_flag = false;
                    break;
                }

                //--- add pos list of molecule
                psys[atom_inclement + i].setPos(pos_norm);
            }

            if(try_count == try_limit){ throw std::invalid_argument("box size is too small or molecules are too much."); }

            if( !pos_flag ) {
                continue;  // retry
            } else {
                break;     // position was set
            }
        }

        //--- set velocity
    //    PS::F64 mass = 0.0;
    //    for(size_t i=0; i<model_template.size(); ++i){
    //        mass += model_template.at(i).getMass();
    //    }
    //    PS::F64vec v = calc_vel(temperature, mass, rand, dist, blz);
    //    for(size_t i=0; i<model_template.size(); ++i){
    //        psys[atom_inclement + i].setVel(v);
    //    }
        for(size_t i=0; i<model_template.size(); ++i){
            psys[atom_inclement + i].setVel( calc_vel(temperature, psys[i].getMass(), rand, dist, blz) );
        }

        atom_inclement += model_template.size();
        ++mol_inclement;

        return try_count;
    }


    //--- initialize PS::ParticleSystem<Tptcl>
    template <typename Tpsys>
    void InitParticle(      Tpsys   &psys,
                      const PS::F64 &temperature){

        if( PS::Comm::getRank() != 0 ) return;

        std::cout << " Initialize: temperature = " << temperature << " [K]" << std::endl;
        assert(temperature >= 0.0);

        //--- set array size
        PS::S64 n_total = 0;
        for(size_t index=0; index<System::model_list.size(); ++index){
            PS::S64 n_atom_mol = System::model_template.at(index).size();
            PS::F64 n_insert   = System::model_list.at(index).second;
            n_total += n_atom_mol*n_insert;

            #ifndef TEST_MOL_INSTALL
                if(n_insert == 0) continue;
            #endif

            std::cout << " Initialize: model_name = " << System::model_list.at(index).first << "\n"
                      << "               n_atom_mol = " << std::setw(10) << n_atom_mol
                      <<              " *n_insert = "   << std::setw(10) << n_insert   << std::endl;
        }

        //--- allocate arrays
        psys.setNumberOfParticleLocal(n_total);
        MODEL::intra_pair_manager.setAtomNumber(n_total);

        std::cout << " Initialize: total atoms = " << n_total << std::endl;

        //--- initialize random number & distribution generator
        constexpr int                    seed = std::pow(2, 19) + 1;
        std::mt19937_64                  mt(seed);
        std::uniform_real_distribution<> dist(0.0, 1.0);   // [0.0, 1.0)
        MD_EXT::boltzmann_dist           blz_dist;

        //--- input molecule
        PS::S64 atom_inclement = 0;
        PS::S64 mol_inclement  = 0;
        for(size_t index=0; index<System::model_list.size(); ++index){
            for(PS::S64 num=0; num<System::model_list.at(index).second; ++num){

                #ifdef TEST_MOL_INSTALL
                    std::cout << " installing " << ENUM::whatis(System::model_list.at(index).first) << std::endl;
                #endif

                size_t try_count = install_molecule(psys,
                                                    MODEL::intra_pair_manager,
                                                    System::model_template.at(index),
                                                    System::bond_template.at(index),
                                                    Normalize::normCutOff(System::setting.ex_radius), //  ex_r
                                                    System::setting.try_limit,                        //  try_limit
                                                    temperature,                                      //  temperature
                                                    mt, dist, blz_dist,
                                                    atom_inclement,
                                                    mol_inclement );

                #ifdef TEST_MOL_INSTALL
                    std::cout << "   n_try = "       << try_count << "\n"
                              << "     n_atom_mol: " << System::model_template.at(index).size()
                              << " / MolID: "        << mol_inclement - 1 << std::endl;
                #endif
            }
        }
    }
}
