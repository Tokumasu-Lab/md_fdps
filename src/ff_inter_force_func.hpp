//***************************************************************************************
//  This program is the intermolecular interactuion of "md_fdps_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

namespace FORCE {

    //--- Cutoff functions  (copy from FDPS-master/sample/c++/p3m/main.cpp)
    inline PS::F64 S2_pcut(const PS::F64 xi) {
       // This is the potential cutoff function where we used Eq.(8.75)
       // in Hockney & Eastwood (1987).

       if (xi <= 1.0) {
          return 1.0 - xi*(208.0
                          +(xi*xi)*(-112.0
                                   +(xi*xi)*(56.0
                                            +xi*(-14.0
                                                +xi*(-8.0
                                                    +3.0*xi)))))/140.0;
       } else if ((1.0 < xi) && (xi < 2.0)) {
          return 1.0 - (12.0
                       +xi*(128.0
                           +xi*(224.0
                               +xi*(-448.0
                                   +xi*(280.0
                                       +xi*(-56.0
                                           +xi*(-14.0
                                               +xi*(8.0
                                                   -xi))))))))/140.0;
       } else {
          return 0.0;
       }
    }

    inline PS::F64 S2_fcut(const PS::F64 xi) {
       // This function returns 1 - R(\xi), where \xi is r/(a/2), a is the
       // scale length of the cutoff function, and R(\xi) is almost the same
       // as the function defined as Eq.(8-72) in Hockney & Eastwood (1987).
       // The only difference is that [1/(r/(a/2))]^2 is factored out
       // in this function from Eq.(8-72).

       if (xi <= 1.0) {
          return 1.0 - (xi*xi*xi)*(224.0
                                  +(xi*xi)*(-224.0
                                           +xi*(70.0
                                               +xi*(48.0-21.0*xi))))/140.0;
       } else if ((1.0 < xi) && (xi < 2.0)) {
          return 1.0 - (12.0
                       +(xi*xi)*(-224.0
                                +xi*(896.0
                                    +xi*(-840.0
                                        +xi*(224.0
                                            +xi*(70.0
                                                +xi*(-48.0+7.0*xi)))))))/140.0;
       } else {
          return 0.0;
       }
    }

    //--- simple functions
    //------ culculate virial value of particle i
    inline PS::F64vec calcVirialEPI(const PS::F64vec &pos, const PS::F64vec &force){
        PS::F64vec tmp = 0.0;
        tmp.x = 0.5 * pos.x * force.x;
        tmp.y = 0.5 * pos.y * force.y;
        tmp.z = 0.5 * pos.z * force.z;
        return tmp;
    }

    //--- basic Particle-Particle function (with cut off)
    template <class Tforce, class Tepi, class Tepj>
    void calcForceShort_IJ_coulombSP_LJ12_6(const Tepi    &ep_i,
                                            const Tepj    &ep_j,
                                            const PS::F64 &r2_cut_LJ,
                                            const PS::F64 &r2_cut_coulomb,
                                            const PS::F64 &r_cut_coulomb_inv,
                                                  Tforce  &force_IJ          ){

        //--- intermolecular interaction
        PS::F64vec r_ij = ep_i.getPos() - ep_j.getPos();
                   r_ij = Normalize::realPos(r_ij);
        PS::F64    r2   = r_ij*r_ij;

        //--- mask for same atom (workaround to zero-devide)
        if( ep_i.getAtomID() == ep_j.getAtomID() ){
            r_ij = PS::F64vec{1e10, 1e10, 1e10};
            r2   = 1e20;  // this value must be > r2_cut.
        }

        PS::F64 r2_inv = 1.0/r2;
        PS::F64 r_inv  = sqrt(r2_inv);

        //--- cut off function for ParticleMesh
        const PS::F64 r_scale        = 2.0*(r2*r_inv)*r_cut_coulomb_inv;
        const PS::F64 factor_S2_pcut = S2_pcut(r_scale);
        const PS::F64 factor_S2_fcut = S2_fcut(r_scale);

        //--- coulomb PP part:
        PS::F64 factor_PM_pot   = factor_S2_pcut;
        PS::F64 factor_PM_force = factor_S2_fcut;

        //--- cut off radius
        PS::F64 factor_LJ = 1.0;
        if( r2 > r2_cut_LJ      ) factor_LJ = 0.0;
        if( r2 > r2_cut_coulomb ){
            factor_PM_pot   = 0.0;
            factor_PM_force = 0.0;
        };

        //--- VDW part
        PS::F64 vddm = ep_i.getVDW_D()*ep_j.getVDW_D();      // VDW_D values are pre-affected "sqrt"
        PS::F64 vdrm = ep_i.getVDW_R() + ep_j.getVDW_R();    // VDW_R values are pre-affected "0.5*"
        PS::F64 sbr6 = vdrm*vdrm*r2_inv;
                sbr6 = sbr6*sbr6*sbr6;                       // (r0/r)^6

                vddm = factor_LJ*vddm;                       // affect scaling mask

        PS::F64    pot_ij =   0.5*vddm*sbr6*(sbr6-2.0);      // 0.5* for double count
        PS::F64vec f_ij   = (12.0*vddm*sbr6*(sbr6-1.0)*r2_inv)*r_ij;
        force_IJ.addPotLJ(    pot_ij );
        force_IJ.addForceLJ(  f_ij   );
        force_IJ.addVirialLJ( calcVirialEPI(r_ij, f_ij) );

        //--- coulomb part
        pot_ij =   factor_PM_pot  *ep_j.getCharge()*r_inv;
        f_ij   = ( factor_PM_force*ep_j.getCharge()*r2_inv )*r_ij;
        force_IJ.addPotCoulomb(   pot_ij );
        force_IJ.addFieldCoulomb( f_ij   );
    }

    //--- basic Particle-Particle mask function
    template <class Tforce, class Tepi, class Tepj, class Tmask>
    void calcForceMask_IJ_coulombSP_LJ12_6(const Tepi    &ep_i,
                                           const Tepj    &ep_j,
                                           const Tmask   &mask_ij,
                                                 Tforce  &force_IJ){

        //--- intermolecular interaction
        PS::F64vec r_ij = ep_i.getPos() - ep_j.getPos();
                   r_ij = Normalize::relativePosAdjustNorm(r_ij);
                   r_ij = Normalize::realPos(r_ij);
        PS::F64    r2   = r_ij*r_ij;

        PS::F64 r2_inv = 1.0/r2;
        PS::F64 r_inv  = sqrt(r2_inv);

        //--- LJ mask ( "-1.0" is cancelation for "calcForceShort_IJ_coulombSP_LJ12_6()")
        PS::F64 factor_LJ = mask_ij.scale_LJ - 1.0;

        //--- coulomb part: ( "-1.0" is cancelation for "calcForceShort_IJ_coulombSP_LJ12_6()" + "ParticleMesh")
        PS::F64 factor_PM_pot   = mask_ij.scale_coulomb - 1.0;
        PS::F64 factor_PM_force = mask_ij.scale_coulomb - 1.0;

        //--- VDW part
        PS::F64 vddm = ep_i.getVDW_D()*ep_j.getVDW_D();      // VDW_D values are pre-affected "sqrt"
        PS::F64 vdrm = ep_i.getVDW_R() + ep_j.getVDW_R();    // VDW_R values are pre-affected "0.5*"
        PS::F64 sbr6 = vdrm*vdrm*r2_inv;
                sbr6 = sbr6*sbr6*sbr6;                       // (r0/r)^6

                vddm = factor_LJ*vddm;                       // affect scaling mask

        PS::F64    pot_ij =   0.5*vddm*sbr6*(sbr6-2.0);      // 0.5* for double count
        PS::F64vec f_ij   = (12.0*vddm*sbr6*(sbr6-1.0)*r2_inv)*r_ij;
        force_IJ.addPotLJ(    pot_ij );
        force_IJ.addForceLJ(  f_ij   );
        force_IJ.addVirialLJ( calcVirialEPI(r_ij, f_ij) );

        //--- coulomb part
        pot_ij =   factor_PM_pot  *ep_j.getCharge()*r_inv;
        f_ij   = ( factor_PM_force*ep_j.getCharge()*r2_inv )*r_ij;
        force_IJ.addPotCoulomb(   pot_ij );
        force_IJ.addFieldCoulomb( f_ij   );
    }

}
