/**************************************************************************************************/
/**
* @file  md_force.hpp
* @brief force calculater interface class.
*/
/**************************************************************************************************/
#pragma once

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>
#include <molecular_dynamics_ext.hpp>

#include "atom_class.hpp"
#include "md_coef_table.hpp"
#include "ff_intra_force.hpp"
#include "ff_inter_force.hpp"
#include "ff_pm_wrapper.hpp"
#include "md_setting.hpp"


/**
* @brief force calculater interface class.
*/
class CalcForce {
private:
    //--- FDPS object
    PS::TreeForForceShort<ForceInter<PS::F64>, EP_inter, EP_inter>::Scatter tree_inter;
    PS::TreeForForceShort<ForceIntra<PS::F64>, EP_intra, EP_intra>::Scatter tree_intra;

    //--- ParticleMesh
    #ifdef REUSE_INTERACTION_LIST
        FORCE::PM::CalcForceParticleMesh pm;
    #else
        FORCE::PM::ParticleMesh pm;
    #endif

    //--- intra pair list maker
    struct GetBond {
        MD_EXT::basic_connect<MD_DEFS::ID_type,
                              MD_DEFS::max_bond> operator () (const AtomConnect &atom){
            return atom.bond;
        }
    };
    IntraPair::IntraMaskMaker<  MD_DEFS::ID_type, GetBond> intra_mask_maker;
    IntraPair::AngleListMaker<  MD_DEFS::ID_type, GetBond> angle_list_maker;
    IntraPair::TorsionListMaker<MD_DEFS::ID_type, GetBond> torsion_list_maker;

    //--- result buffer
    std::vector<ForceInter<PS::F64>> inter_force_buff;

public:
    void init(const PS::S64 &n_total){
        this->tree_inter.initialize(n_total,
                                    System::profile.theta,
                                    System::profile.n_leaf_limit,
                                    System::profile.n_group_limit);
        this->tree_intra.initialize(n_total,
                                    System::profile.theta,
                                    System::profile.n_leaf_limit,
                                    System::profile.n_group_limit);

        FORCE::PM::checkMeshSize(n_total);
    }

    /**
    * @brief update cutoff length in normalized space.
    */
    void setRcut(){
        EP_inter::setR_cut_LJ(      Normalize::normCutOff( System::get_cut_off_LJ() ) );
        EP_inter::setR_cut_coulomb( Normalize::normCutOff_PM() );

        EP_intra::setR_cut( Normalize::normCutOff( System::get_cut_off_intra() ) );

        #ifdef REUSE_INTERACTION_LIST
            EP_inter::setR_margin( Normalize::normCutOff( 2.0 ) );
            EP_intra::setR_margin( Normalize::normCutOff( 2.0 ) );
        #endif

        //--- check cut off length
        if(EP_inter::getRcut_LJ() >= 0.5 ||
           EP_inter::getRcut_LJ() <= 0.0 ){
            std::ostringstream oss;
            oss << "RSearch for LJ must be in range of (0.0, 0.5) at normalized space." << "\n"
                << "    EP_inter::getRcut_LJ() = " << EP_inter::getRcut_LJ() << "\n";
            throw std::length_error(oss.str());
        }
        if(EP_intra::getRSearch() >= 0.5 ||
           EP_intra::getRSearch() <= 0.0 ){
            std::ostringstream oss;
            oss << "RSearch for intramolecuar force must be in range of (0.0, 0.5) at normalized space." << "\n"
                << "    EP_intra::getRSearch() = " << EP_intra::getRSearch() << "\n";
            throw std::length_error(oss.str());
        }
    }

    /**
    * @brief update intra pair list at each atom.
    */
    template <class Tptcl, class Tdinfo, class Tmask>
    void update_intra_pair_list(      PS::ParticleSystem<Tptcl> &atom,
                                      Tdinfo                    &dinfo,
                                const Tmask                     &mask_table){

        //--- get neighbor EP_intra information (do not calculate force)
        this->setRcut();
        this->tree_intra.calcForceAll( IntraPair::dummy_func{},
                                       atom,
                                       dinfo                   );

        const PS::S32 n_local = atom.getNumberOfParticleLocal();

        Tptcl::clear_intra_pair_table();

        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            //--- clear & allocate
            for(PS::S32 i=0; i<n_local; ++i){
                atom[i].clear_intra_list();
            }

            //--- search pair
            #pragma omp parallel for
            for(PS::S32 i=0; i<n_local; ++i){
                intra_mask_maker(  atom[i], this->tree_intra, mask_table.at(atom[i].getMolType()),
                                                              atom[i].mask_list());
                angle_list_maker(  atom[i], this->tree_intra, atom[i].angle_list());
                torsion_list_maker(atom[i], this->tree_intra, atom[i].dihedral_list(), atom[i].improper_list());
            }
        #else
            for(PS::S32 i=0; i<n_local; ++i){
                atom[i].clear_intra_list();

                intra_mask_maker(  atom[i], this->tree_intra, mask_table.at(atom[i].getMolType()),
                                                              atom[i].mask_list());
                angle_list_maker(  atom[i], this->tree_intra, atom[i].angle_list());
                torsion_list_maker(atom[i], this->tree_intra, atom[i].dihedral_list(), atom[i].improper_list());
            }
        #endif

    }

    /**
    * @brief update intramolecular force on atom.
    */
    template <class Tpsys, class Tdinfo>
    void update_intra_force(Tpsys  &atom,
                            Tdinfo &dinfo){

        //=================
        // Intra force part
        //=================
        //--- clear intramolecular part
        const PS::S64 n_local = atom.getNumberOfParticleLocal();

        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            #pragma omp parallel for
        #endif
        for(PS::S64 i=0; i<n_local; ++i){
            atom[i].clearForceIntra();
        }

        this->setRcut();

        //--- get neighbor EP_intra information (do not calculate force)
        this->tree_intra.calcForceAll( IntraPair::dummy_func{},
                                       atom,
                                       dinfo                   );  // remake list
        //--- calculate force
        FORCE::calcForceIntra(this->tree_intra, atom);
    }

    /**
    * @brief update intermolecular force on atom (naive version).
    */
    template <class Tpsys, class Tdinfo>
    void update_inter_force_naive(      Tpsys                     &atom,
                                        Tdinfo                    &dinfo,
                                  const PS::INTERACTION_LIST_MODE  reuse_mode = PS::MAKE_LIST){

        //--- clear intermolecular part
        const PS::S64 n_local = atom.getNumberOfParticleLocal();

        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            #pragma omp parallel for
        #endif
        for(PS::S64 i=0; i<n_local; ++i){
            atom[i].clearForceInter();
        }

        this->setRcut();

        //=================
        // PM part
        //=================
        this->pm.setDomainInfoParticleMesh(dinfo);
        this->pm.setParticleParticleMesh(atom, true);   // clear previous charge information
        this->pm.calcMeshForceOnly();
        this->pm.writeBackForce(atom);

        //=================
        // PP part
        //=================
        this->tree_inter.calcForceAll(FORCE::calcForceShort_naive{},
                                      atom,
                                      dinfo,
                                      true,
                                      reuse_mode);
        for(PS::S64 i=0; i<n_local; ++i){
            const auto& result = tree_inter.getForce(i);
            atom[i].addPotLJ(        result.getPotLJ()        );
            atom[i].addForceLJ(      result.getForceLJ()      );
            atom[i].addVirialLJ(     result.getVirialLJ()     );
            atom[i].addPotCoulomb(   result.getPotCoulomb()   );
            atom[i].addFieldCoulomb( result.getFieldCoulomb() );
        }
    }

    /**
    * @brief   update intermolecular force on atom (optimized version).
    * @details delayed evaluation for intramolecular mask. if blanch is removed in P-P calculater kernel.
    */
    template <class Tpsys, class Tdinfo>
    void update_inter_force(      Tpsys                     &atom,
                                  Tdinfo                    &dinfo,
                            const PS::INTERACTION_LIST_MODE  reuse_mode = PS::MAKE_LIST){

        //--- debug use
        #ifdef FORCE_NAIVE_IMPL
            this->update_inter_force_naive(atom, dinfo, reuse_mode);
            return;
        #endif

        //--- clear intermolecular part
        const PS::S64 n_local = atom.getNumberOfParticleLocal();

        this->inter_force_buff.resize(n_local);

        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            #pragma omp parallel for
        #endif
        for(PS::S64 i=0; i<n_local; ++i){
            atom[i].clearForceInter();
        }

        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            #pragma omp parallel for
        #endif
        for(PS::S64 i=0; i<n_local; ++i){
            this->inter_force_buff[i].clear();
        }

        this->setRcut();

        //=================
        // PM part
        //=================
        this->pm.setDomainInfoParticleMesh(dinfo);
        this->pm.setParticleParticleMesh(atom, true);   // clear previous charge information
        this->pm.calcMeshForceOnly();
        this->pm.writeBackForce(atom);

        //=================
        // PP part (without mask)
        //=================
        this->tree_inter.calcForceAll(FORCE::calcForceShort{},
                                      atom,
                                      dinfo,
                                      true,
                                      reuse_mode);
        for(PS::S64 i=0; i<n_local; ++i){
            const auto& result = tree_inter.getForce(i);
                  auto& buf    = this->inter_force_buff.at(i);
            buf.addFieldCoulomb( result.getFieldCoulomb() );
            buf.addPotCoulomb(   result.getPotCoulomb()   );
            buf.addForceLJ(      result.getForceLJ()      );
            buf.addPotLJ(        result.getPotLJ()        );
            buf.addVirialLJ(     result.getVirialLJ()     );
        }
        //=================
        // PP part (evaluate mask)
        //=================
        FORCE::calcForceIntraMask(this->tree_inter,
                                  atom,
                                  this->inter_force_buff);

        //=================
        // PP part (writeback)
        //=================
        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            #pragma omp parallel for
        #endif
        for(PS::S64 i=0; i<n_local; ++i){
            const auto& buf = this->inter_force_buff[i];
            atom[i].addFieldCoulomb( buf.getFieldCoulomb() );
            atom[i].addPotCoulomb(   buf.getPotCoulomb()   );
            atom[i].addForceLJ(      buf.getForceLJ()      );
            atom[i].addPotLJ(        buf.getPotLJ()        );
            atom[i].addVirialLJ(     buf.getVirialLJ()     );
        }
    }

    /**
    * @brief update force on atom (naive version).
    */
    template <class Tpsys, class Tdinfo>
    void update_force_naive(      Tpsys                     &atom,
                                  Tdinfo                    &dinfo,
                            const PS::INTERACTION_LIST_MODE  reuse_mode = PS::MAKE_LIST){

        this->update_inter_force_naive(atom, dinfo, reuse_mode);
        this->update_intra_force(      atom, dinfo);
    }

    /**
    * @brief update force on atom (optimized version).
    */
    template <class Tpsys, class Tdinfo>
    void update_force(      Tpsys                     &atom,
                            Tdinfo                    &dinfo,
                      const PS::INTERACTION_LIST_MODE  reuse_mode = PS::MAKE_LIST){

        this->update_inter_force(atom, dinfo, reuse_mode);
        this->update_intra_force(atom, dinfo);
    }
};
