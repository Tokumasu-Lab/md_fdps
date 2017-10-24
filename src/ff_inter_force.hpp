//***************************************************************************************
//  This program is the intermolecular interactuion of "md_fdps_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <algorithm>

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

#include "md_coef_table.hpp"
#include "ff_inter_force_func.hpp"



//--- Template function calculate Particle-Particle interaction
template <class Tforce, class Tepi, class Tepj>
void calcForceShort_IA(const Tepi    &ep_i,
                       const Tepj    *ep_j,
                       const PS::S32 &n_ep_j,
                       const PS::F64 &r2_cut_LJ,
                       const PS::F64 &r2_cut_coulomb,
                       const PS::F64 &r_cut_coulomb_inv,
                             Tforce  &force_IA){

    for(PS::S32 j=0; j<n_ep_j; ++j){
        //--- skip same atom
        if(ep_i.getAtomID() == ep_j[j].getAtomID() ) continue;

        PS::F32 sign_inter;
        if( MODEL::intra_pair_manager.checkIntraMask( ep_i.getAtomID(), ep_j[j].getAtomID() ) ){
            //--- inter pair
            sign_inter = 1.0;
        } else {
            //--- intra pair
            sign_inter = 0.0;
        }

        //--- intermolecular interaction
        PS::F64vec r_ij = ep_i.getPos() - ep_j[j].getPos();
                   r_ij = Normalize::realPos(r_ij);
        PS::F64    r2   = r_ij*r_ij;

        PS::F32 sign_vdw     = 1.0;
        PS::F32 sign_coulomb = 1.0;
        if( r2 > r2_cut_LJ      ) sign_vdw     = 0.0;
        if( r2 > r2_cut_coulomb ) sign_coulomb = 0.0;

        PS::F64 r2_inv = 1.0/r2;
        PS::F64 r_inv  = sqrt(r2_inv);

        //--- VDW part
        PS::F64 vddm = ep_i.getVDW_D()*ep_j[j].getVDW_D();      // VDW_D values are pre-affected "sqrt"
        PS::F64 vdrm = ep_i.getVDW_R() + ep_j[j].getVDW_R();    // VDW_R values are pre-affected "0.5*"
        PS::F64 sbr6 = vdrm*vdrm*r2_inv;
                sbr6 = sbr6*sbr6*sbr6;                          // coef*(r0/r)^6

                vddm = sign_inter*sign_vdw*vddm;                // affect vdw_sign (intra_mask fillter)

        PS::F64    pot_ij =   0.5*vddm*sbr6*(sbr6-2.0);         // 0.5* for double count
        PS::F64vec f_ij   = (12.0*vddm*sbr6*(sbr6-1.0)*r2_inv)*r_ij;
        force_IA.addPotLJ(    pot_ij );
        force_IA.addForceLJ(  f_ij   );
        force_IA.addVirialLJ( calcVirialEPI(r_ij, f_ij) );

        //--- coulomb part
        PS::F64 r_scale = 2.0*(r2*r_inv)*r_cut_coulomb_inv;
        //--- (1.0 - sign_inter) + 2.0*(sign_inter - 0.5)*S2_cut() = S2_cut       at ( sign == 1.0 )
        //                                                         = 1.0 - S2_cut at ( sign == 0.0 )
        pot_ij =   ( (1.0 - sign_inter) + 2.0*(sign_inter - 0.5)*S2_pcut(r_scale) )*ep_j[j].getCharge()*r_inv;
        f_ij   = ( ( (1.0 - sign_inter) + 2.0*(sign_inter - 0.5)*S2_fcut(r_scale) )*ep_j[j].getCharge()*r2_inv )*r_ij;
        force_IA.addPotCoulomb(   sign_coulomb*pot_ij );
        force_IA.addFieldCoulomb( sign_coulomb*f_ij   );
    }
}

template<class Tforce, class Tepi, class Tepj>
void calcForceShort(      Tepi    *ep_i,
                    const PS::S32  n_ep_i,
                          Tepj    *ep_j,
                    const PS::S32  n_ep_j,
                          Tforce  *force){

    PS::F64 r2_cut_LJ         = Normalize::realCutOff( ep_i->getRcut_LJ() );
            r2_cut_LJ         = r2_cut_LJ*r2_cut_LJ;
    PS::F64 r2_cut_coulomb    = Normalize::realCutOff( ep_i->getRcut_coulomb() );
    PS::F64 r_cut_inv_coulomb = 1.0/r2_cut_coulomb;
            r2_cut_coulomb    = r2_cut_coulomb*r2_cut_coulomb;

    for(PS::S32 i=0; i<n_ep_i; ++i){
        calcForceShort_IA(ep_i[i],
                          ep_j,
                          n_ep_j,
                          r2_cut_LJ,
                          r2_cut_coulomb,
                          r_cut_inv_coulomb,
                          force[i]  );
    }
}
