//***************************************************************************************
//  This program is the intermolecular interactuion of "md_fdps_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <algorithm>

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>

#include "ff_inter_force_func.hpp"


namespace FORCE {

    /*
    *  @breif naive implementation: intramolecular mask is considered in this function.
    */
    struct calcForceShort_naive{
        template <class Tepi, class Tepj, class Tforce>
        void operator () (const Tepi    *ep_i,
                          const PS::S32  n_ep_i,
                          const Tepj    *ep_j,
                          const PS::S32  n_ep_j,
                                Tforce  *force){

            const PS::F64 r_cut_LJ          = Normalize::realCutOff( Tepi::getRcut_LJ() );
            const PS::F64 r2_cut_LJ         = r_cut_LJ*r_cut_LJ;
            const PS::F64 r_cut_coulomb     = Normalize::realCutOff( Tepi::getRcut_coulomb() );
            const PS::F64 r_cut_coulomb_inv = 1.0/r_cut_coulomb;
            const PS::F64 r2_cut_coulomb    = r_cut_coulomb*r_cut_coulomb;

            for(PS::S32 i=0; i<n_ep_i; ++i){
                Tforce force_IA;
                force_IA.clear();
                for(PS::S32 j=0; j<n_ep_j; ++j){
                    Tforce force_ij;
                    force_ij.clear();
                    calcForceShort_IJ_coulombSP_LJ12_6(ep_i[i],
                                                       ep_j[j],
                                                       r2_cut_LJ,
                                                       r2_cut_coulomb,
                                                       r_cut_coulomb_inv,
                                                       force_ij);
                    const auto mask = ep_i[i].find_mask( ep_j[j].getAtomID() );
                    if( mask.is_effective() ){
                        calcForceMask_IJ_coulombSP_LJ12_6(ep_i[i],
                                                          ep_j[j],
                                                          mask,
                                                          force_ij);
                    }
                    force_IA.addPotLJ(       force_ij.getPotLJ()       );
                    force_IA.addForceLJ(     force_ij.getForceLJ()     );
                    force_IA.addVirialLJ(    force_ij.getVirialLJ()    );
                    force_IA.addPotCoulomb(  force_ij.getPotCoulomb()  );
                    force_IA.addFieldCoulomb(force_ij.getFieldCoulomb());
                }
                force[i].copyFromForce(force_IA);

                //--- self consistant term for PM
                force[i].addPotCoulomb( -ep_i[i].getCharge()*(208.0/70.0)*r_cut_coulomb_inv );
            }
        }
    };

    /*
    *  @breif optimized implementation: intramolecular mask is ignored in this function.
    *         when use this function, must consider mask by the function of "calcForceIntraMask()" in below.
    */
    struct calcForceShort{
        template <class Tepi, class Tepj, class Tforce>
        void operator () (const Tepi    *ep_i,
                          const PS::S32  n_ep_i,
                          const Tepj    *ep_j,
                          const PS::S32  n_ep_j,
                                Tforce  *force){

            const PS::F64 r_cut_LJ          = Normalize::realCutOff( Tepi::getRcut_LJ() );
            const PS::F64 r2_cut_LJ         = r_cut_LJ*r_cut_LJ;
            const PS::F64 r_cut_coulomb     = Normalize::realCutOff( Tepi::getRcut_coulomb() );
            const PS::F64 r_cut_coulomb_inv = 1.0/r_cut_coulomb;
            const PS::F64 r2_cut_coulomb    = r_cut_coulomb*r_cut_coulomb;

            for(PS::S32 i=0; i<n_ep_i; ++i){
                Tforce force_IA;
                force_IA.clear();
                for(PS::S32 j=0; j<n_ep_j; ++j){
                    calcForceShort_IJ_coulombSP_LJ12_6(ep_i[i],
                                                       ep_j[j],
                                                       r2_cut_LJ,
                                                       r2_cut_coulomb,
                                                       r_cut_coulomb_inv,
                                                       force_IA);
                }
                force[i].copyFromForce(force_IA);

                //--- self consistant term for PM
                force[i].addPotCoulomb( -ep_i[i].getCharge()*(208.0/70.0)*r_cut_coulomb_inv );
            }
        }
    };


    /*
    *  @breif fuction for intramolecular mask evaluation.
    *         use with the 'calcForceShort()' functor.
    */
    template <class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj,
              class Tpsys>
    void calcForceIntraMask(PS::TreeForForce<TSM,
                                             Tforce,
                                             Tepi,
                                             Tepj,
                                             Tmomloc,
                                             Tmomglb,
                                             Tspj    > &tree,
                            Tpsys                      &atom,
                            std::vector<Tforce>        &pp_force_buff){

        const PS::S64 n_local = atom.getNumberOfParticleLocal();

        assert(static_cast<PS::S64>(pp_force_buff.size()) == n_local );

        //--- calculate intramolecular force
        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            #pragma omp parallel for
        #endif
        for(PS::S64 i=0; i<n_local; ++i){
            const auto& fp_i = atom[i];

            Tforce force_IA;
            force_IA.clear();

            for(const auto& mask : fp_i.mask_list()){
                const auto* ptr_j = tree.getEpjFromId( mask.getId() );
                IntraPair::check_nullptr_ptcl(ptr_j, fp_i.getAtomID(), mask.getId() );

                //--- evaluate mask
                calcForceMask_IJ_coulombSP_LJ12_6(fp_i,
                                                  *ptr_j,
                                                  mask,
                                                  force_IA);
            }
            
            pp_force_buff[i].addPotLJ(        force_IA.getPotLJ()        );
            pp_force_buff[i].addForceLJ(      force_IA.getForceLJ()      );
            pp_force_buff[i].addVirialLJ(     force_IA.getVirialLJ()     );
            pp_force_buff[i].addPotCoulomb(   force_IA.getPotCoulomb()   );
            pp_force_buff[i].addFieldCoulomb( force_IA.getFieldCoulomb() );
        }
    }

}
