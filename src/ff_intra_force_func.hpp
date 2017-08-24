//***************************************************************************************
//  This program is the intramolecular force calculation for "md_fdps_main.cpp"
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <cmath>
#include <cassert>
#include <stdexcept>

#include <particle_simulator.hpp>

#include "ff_inter_force_func.hpp"


//--- calculate intramolecular potential --------------------------------------------
//------ harmonic bond potential
template <class Tcoef, class Tforce>
void calcBondForce_harmonic_IJ(const PS::F64vec &pos_i,
                               const PS::F64vec &pos_j,
                               const Tcoef      &coef,
                                     Tforce     &force_i){

    PS::F64vec r_ij = pos_i - pos_j;
    PS::F64    r2   = r_ij*r_ij;
    PS::F64    r    = sqrt(r2);

    //--- assert
    #ifndef NDEBUG
        if(r > 4.0){
            std::cerr << "  r = " << r << std::endl;
            throw std::invalid_argument(" bond length is too long.");
        }
    #endif

    //--- potential
    PS::F64 r_diff = r - coef.r0;
    force_i.addPotBond( 0.5*0.5*coef.k*r_diff*r_diff );

    //--- force
    PS::F64vec f_tmp = ( -coef.k*r_diff/r )*r_ij;
    force_i.addForceIntra( f_tmp );

    //--- virial
    force_i.addVirialIntra( calcVirialEPI(r_ij, f_tmp) );
}

//------ anharmonic bond potential
template <class Tcoef, class Tforce>
void calcBondForce_anharmonic_IJ(const PS::F64vec  &pos_i,
                                 const PS::F64vec  &pos_j,
                                 const Tcoef       &coef,
                                       Tforce      &force_i){

    PS::F64vec r_ij = pos_i - pos_j;
    PS::F64    r2   = r_ij*r_ij;
    PS::F64    r    = sqrt(r2);

    //--- assert
    #ifndef NDEBUG
        if(r > 4.0){
            std::cerr << "  r = " << r << std::endl;
            throw std::invalid_argument(" bond length is too long.");
        }
    #endif

    //--- potential
    constexpr PS::F64 factor = 7.0/12.0;  // counteract double count
    PS::F64 ar  = coef.a*(r - coef.r0);
    PS::F64 ar2 = ar*ar;
    force_i.addPotBond(    0.5*coef.k*(  1.0 - ar + factor*ar2)*ar2 );

    //--- force
    PS::F64    ebp   = -coef.a*coef.k*( (2.0 - 3.0*ar) + 4.0*factor*ar2 )*ar/r;
    PS::F64vec f_tmp = ebp*r_ij;
    force_i.addForceIntra( f_tmp );

    //--- virial
    force_i.addVirialIntra( calcVirialEPI(r_ij, f_tmp) );
}

//------ harminic angle potential  (i-j-k form: "j" must be center)
template <class Tid, class Tcoef, class Tforce>
void calcAngleForce_harmonic_IJK(const PS::F64vec &pos_i,
                                 const PS::F64vec &pos_j,
                                 const PS::F64vec &pos_k,
                                 const Tid        &id_i,
                                 const Tid        &id_j,
                                 const Tid        &id_k,
                                 const Tid        &id_tgt,
                                 const Tcoef      &coef,
                                       Tforce     &force_tgt){

    PS::F64vec r_vec_A = pos_i - pos_j;
    PS::F64vec r_vec_B = pos_k - pos_j;

    PS::F64 r_a = std::sqrt(r_vec_A*r_vec_A);
    PS::F64 r_b = std::sqrt(r_vec_B*r_vec_B);

    PS::F64 in_prod  = r_vec_A*r_vec_B;
    PS::F64 r_ab_inv = 1.0/(r_a*r_b);
    PS::F64 cos_tmp  = in_prod*r_ab_inv;
    PS::F64 diff     = cos_tmp - std::cos(coef.theta0);

    //--- potential
    force_tgt.addPotAngle( (1.0/3.0)*0.5*coef.k*diff*diff );

    //--- force
    PS::F64 acoef     = -coef.k*diff;
    PS::F64 r_ab2_inv = r_ab_inv*r_ab_inv;

    PS::F64vec f_a = (acoef*r_ab2_inv)*( (r_b*r_b*cos_tmp)*r_vec_A - (r_a*r_b)*r_vec_B );
    PS::F64vec f_b = (acoef*r_ab2_inv)*( (r_a*r_a*cos_tmp)*r_vec_B - (r_a*r_b)*r_vec_A );

    //--- select output
    if(       id_tgt == id_i){
        force_tgt.addForceIntra(-f_a);
        force_tgt.addVirialIntra( -calcVirialEPI(r_vec_A, f_a) );
    } else if(id_tgt == id_j) {
        force_tgt.addForceIntra(f_a + f_b);
        force_tgt.addVirialIntra(  calcVirialEPI(r_vec_A, f_a)
                                 + calcVirialEPI(r_vec_B, f_b) );
    } else if(id_tgt == id_k) {
        force_tgt.addForceIntra(-f_b);
        force_tgt.addVirialIntra( -calcVirialEPI(r_vec_B, f_b) );
    } else {
        std::cerr << "   id_tgt = " << id_tgt << "\n"
                  << "        i = " << id_i
                  <<       "  j = " << id_j
                  <<       "  k = " << id_k << std::endl;
        throw std::invalid_argument("id_tgt is not match to id_i, id_j, or id_k.");
    }
}

//------ harminic torsion potential  (i-j-k-l form)
template <class Tid, class Tcoef, class Tforce>
void calcTorsionForce_harmonic_IJKL(const PS::F64vec &pos_i,
                                    const PS::F64vec &pos_j,
                                    const PS::F64vec &pos_k,
                                    const PS::F64vec &pos_l,
                                    const Tid        &id_i,
                                    const Tid        &id_j,
                                    const Tid        &id_k,
                                    const Tid        &id_l,
                                    const Tid        &id_tgt,
                                    const Tcoef      &coef,
                                          Tforce     &force_tgt){

    PS::F64vec r_ji = pos_i - pos_j;
    PS::F64vec r_jk = pos_k - pos_j;
    PS::F64vec r_jl = pos_l - pos_j;

    PS::F64vec r_A = pos_j - pos_i;
    PS::F64vec r_B = pos_k - pos_j;
    PS::F64vec r_C = pos_l - pos_k;

    PS::F64vec r_G1 = r_A^r_B;
    PS::F64vec r_G2 = r_B^r_C;
    PS::F64vec r_G3 = r_G1^r_G2;

    PS::F64 r2_g1 = r_G1*r_G1;
    PS::F64 r2_g2 = r_G2*r_G2;
    PS::F64 r_g1  = std::sqrt(r2_g1);
    PS::F64 r_g2  = std::sqrt(r2_g2);
    PS::F64 r_g12 = r_g1*r_g2;
    PS::F64 r_g12_inv = 1.0/r_g12;

    PS::F32 v    = r_G3*r_jk;
    PS::S32 sign = (v > 0.0) - (v < 0.0);   // extract sign, sign = -1 (v<0), 0 (v==0), or 1 (v>0).
    PS::F64 cos_p  = (r_G1*r_G2)*r_g12_inv;
    PS::F64 acos_p = -sign*cos_p;
    PS::F64 sin_p  = std::sin(acos_p);
            cos_p  = std::cos(acos_p);
    //PS::F64 sin_f  = std::sin(PS::F64(coef.n)*acos_p - coef.theta0);
    PS::F64 cos_eq = std::cos(coef.theta0);

    PS::F64 sin_k;
    switch (coef.n_min) {
        case 1:
            sin_k = cos_eq;
        break;

        case 2:
            sin_k = 2.0*cos_p*cos_eq;
        break;

        case 3:
            sin_k = -(4.0*sin_p*sin_p + 3.0)*cos_eq;
        break;

        case 4:
            sin_k = 4.0*cos_p*(2.0*cos_p*cos_p - 1.0)*cos_eq;
        break;

        case 6:
            sin_k = (  ( std::cos(2.0*acos_p)*(4.0*cos_p*cos_p - 2.0) )
                     + ( std::cos(4.0*acos_p) )                         )*2.0*cos_p*cos_eq;
        break;

        default:
            std::cerr << "  n_min = " << coef.n_min << std::endl;
            throw std::invalid_argument("undefined number of equiliblium points.");
    }

    //--- potential
    PS::F64 eng_tmp = 0.0;
    if(        coef.form == IntraFuncForm::cos  ){
        eng_tmp = (1.0/4.0)*0.5*coef.k*( 1.0 - std::cos( PS::F64(coef.n_min)*acos_p - coef.k ) );
    } else if( coef.form == IntraFuncForm::none ){
        eng_tmp = (1.0/4.0)*0.5*coef.k;
    } else {
        std::cerr << "  form = " << ENUM::whatis(coef.form) << std::endl;
        throw std::invalid_argument("undefined potential form for torsion potential.");
    }
    force_tgt.addPotTorsion(eng_tmp);

    PS::F64 pr_aa = r_A*r_A;
    PS::F64 pr_bb = r_B*r_B;
    PS::F64 pr_cc = r_C*r_C;
    PS::F64 pr_ab = r_A*r_B;
    PS::F64 pr_bc = r_B*r_C;
    PS::F64 pr_ca = r_C*r_A;
    PS::F64 t_coef = -0.5*coef.k*PS::F64(coef.n_min)*sin_k*r_g12_inv;

    //--- select output
    PS::F64vec F_tmp;
    if(       id_tgt == id_i){
        F_tmp = (-r2_g2*pr_bb*cos_p)*r_A
              + ( r2_g2*pr_ab*cos_p + r_g12*pr_bc)*r_B
              + (-r_g12*pr_bb)*r_C;
        force_tgt.addForceIntra(t_coef*F_tmp);
        force_tgt.addVirialIntra( 2.0*calcVirialEPI(r_ji, F_tmp) );
    } else if(id_tgt == id_j) {
        F_tmp = (r2_g2*(pr_bb + pr_ab)*cos_p + r_g12*pr_bc)*r_A
              + (   (r2_g2*(pr_aa + pr_ab) + r2_g1*pr_cc)*cos_p
                  + r_g12*(pr_bc + 2.0*pr_ca) )*r_B
              - (r2_g1*pr_bc*cos_p + r_g12*(pr_ab + pr_bb))*r_C;
        force_tgt.addForceIntra(t_coef*F_tmp);
        //force_tgt.addVirialIntra(virial);     virial = 0, considered on other atoms.
    } else if(id_tgt == id_k) {
        F_tmp = (r2_g2*pr_ab*cos_p + r_g12*(pr_bc + pr_bb))*r_A
              + (   (r2_g2*pr_aa + r2_g1*(pr_bc + pr_cc))*cos_p
                  + r_g12*(pr_ab + 2.0*pr_ca) )*r_B
              - (r2_g1*(pr_bc + pr_bb)*cos_p + r_g12*pr_ab)*r_C;
        force_tgt.addForceIntra(t_coef*F_tmp);
        force_tgt.addVirialIntra( 2.0*calcVirialEPI(r_jk, F_tmp) );
    } else if(id_tgt == id_l) {
        F_tmp = (r_g1*r_g2*pr_bb)*r_A
              - (r2_g1*pr_bc*cos_p + r_g12*pr_ab)*r_B
              + (r2_g1*pr_bb*cos_p)*r_C;
        force_tgt.addForceIntra(t_coef*F_tmp);
        force_tgt.addVirialIntra( 2.0*calcVirialEPI(r_jl, F_tmp) );
    } else {
        std::cerr << "   id_tgt = " << id_tgt << "\n"
                  << "        i = " << id_i
                  <<       "  j = " << id_j
                  <<       "  k = " << id_k
                  <<       "  l = " << id_l << std::endl;
        throw std::invalid_argument("id_tgt is not match to id_i, id_j, id_k, or id_l.");
    }
}
