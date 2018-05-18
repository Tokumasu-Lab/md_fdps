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
    PS::TreeForForceShort<ForceInter<PS::F32>, EP_inter, EP_inter>::Scatter tree_inter;
    PS::TreeForForceShort<ForceIntra<PS::F32>, EP_intra, EP_intra>::Scatter tree_intra;

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
        EP_inter::setRcut_LJ(      Normalize::normCutOff( System::get_cut_off_LJ() ) );
        EP_inter::setRcut_coulomb( Normalize::normCutOff_PM() );

        EP_intra::setRcut( Normalize::normCutOff( System::get_cut_off_intra() ) );

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
    //    this->tree_intra.calcForceAll( IntraPair::dummy_func{},
    //                                   atom,
    //                                   dinfo,
    //                                   false,
    //                                   PS::MAKE_LIST_FOR_REUSE );

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
    * @brief update force on atom.
    */
    template <class Tpsys, class Tdinfo>
    void update_force(      Tpsys                     &atom,
                            Tdinfo                    &dinfo,
                      const PS::INTERACTION_LIST_MODE  reuse_mode = PS::MAKE_LIST){

        //--- clear force
        const PS::S64 n_local = atom.getNumberOfParticleLocal();
        for(PS::S64 i=0; i<n_local; ++i){
            atom[i].clear();
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
        this->tree_inter.calcForceAll(FORCE::calcForceShort<ForceInter<PS::F32>, EP_inter, EP_inter>,
                                      atom,
                                      dinfo,
                                      true,
                                      reuse_mode);
        for(PS::S64 i=0; i<n_local; ++i){
            const auto& result = tree_inter.getForce(i);
            atom[i].addFieldCoulomb( result.getFieldCoulomb() );
            atom[i].addPotCoulomb(   result.getPotCoulomb()   );
            atom[i].addForceLJ(      result.getForceLJ()      );
            atom[i].addPotLJ(        result.getPotLJ()        );
            atom[i].addVirialLJ(     result.getVirialLJ()     );
        }

        //=================
        // Intra force part
        //=================
        //--- get neighbor EP_intra information (do not calculate force)
        #ifdef REUSE_INTERACTION_LIST
            this->tree_intra.calcForceAll( IntraPair::dummy_func{},
                                           atom,
                                           dinfo                   );  // make list in this->update_intra_pair_list()
        //    this->tree_intra.calcForceAll( IntraPair::dummy_func{},
        //                                   atom,
        //                                   dinfo,
        //                                   false,
        //                                   PS::REUSE_LIST          );  // make list in this->update_intra_pair_list()
        #else
            this->tree_intra.calcForceAll( IntraPair::dummy_func{},
                                           atom,
                                           dinfo                   );  // remake list
        #endif
        //--- calculate force
        FORCE::calcForceIntra(this->tree_intra, atom);
    }
};
