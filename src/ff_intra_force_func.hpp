//***************************************************************************************
//  This is the intramolecular force calculation for "md_fdps_main.cpp"
//***************************************************************************************
#pragma once

#include <cmath>
#include <cassert>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "ff_inter_force_func.hpp"


//--- calculate intramolecular potential --------------------------------------------
//------ harmonic bond potential
template <class Tcoef, class Tforce>
void calcBondForce_harmonic_IJ(const PS::F64vec &pos_i,
                               const PS::F64vec &pos_j,
                               const Tcoef      &coef,
                                     Tforce     &force_i){

    PS::F64vec R_ij = pos_i - pos_j;
    PS::F64    r2   = R_ij*R_ij;
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
    PS::F64vec F_ij = ( -coef.k*r_diff/r )*R_ij;
    force_i.addForceIntra( F_ij );

    //--- virial
    force_i.addVirialIntra( calcVirialEPI(F_ij, R_ij) );
}

//------ anharmonic bond potential
template <class Tcoef, class Tforce>
void calcBondForce_anharmonic_IJ(const PS::F64vec  &pos_i,
                                 const PS::F64vec  &pos_j,
                                 const Tcoef       &coef,
                                       Tforce      &force_i){

    PS::F64vec R_ij = pos_i - pos_j;
    PS::F64    r2   = R_ij*R_ij;
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
    PS::F64vec F_ij = ebp*R_ij;
    force_i.addForceIntra( F_ij );

    //--- virial
    force_i.addVirialIntra( calcVirialEPI(F_ij, R_ij) );
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

    PS::F64vec R_ij = pos_i - pos_j;
    PS::F64vec R_kj = pos_k - pos_j;

    PS::F64 r_a = VEC_EXT::norm(R_ij);
    PS::F64 r_b = VEC_EXT::norm(R_kj);

    PS::F64 in_prod  = R_ij*R_kj;
    PS::F64 r_ab_inv = 1.0/(r_a*r_b);
    PS::F64 cos_tmp  = in_prod*r_ab_inv;
    PS::F64 diff     = cos_tmp - std::cos(coef.theta0);

    //--- potential
    force_tgt.addPotAngle( (1.0/3.0)*0.5*coef.k*diff*diff );

    //--- force
    PS::F64 acoef     = -coef.k*diff;
    PS::F64 r_ab2_inv = r_ab_inv*r_ab_inv;

    PS::F64vec F_ij = (acoef*r_ab2_inv)*( (r_b*r_b*cos_tmp)*R_ij - (r_a*r_b)*R_kj );
    PS::F64vec F_kj = (acoef*r_ab2_inv)*( (r_a*r_a*cos_tmp)*R_kj - (r_a*r_b)*R_ij );

    //--- select output
    if(       id_tgt == id_i){
        force_tgt.addForceIntra(-F_ij);
        force_tgt.addVirialIntra( -calcVirialEPI(F_ij, R_ij) );
    } else if(id_tgt == id_j) {
        force_tgt.addForceIntra(F_ij + F_kj);
        force_tgt.addVirialIntra(  calcVirialEPI(F_ij, R_ij)
                                 + calcVirialEPI(F_kj, R_kj) );
    } else if(id_tgt == id_k) {
        force_tgt.addForceIntra(-F_kj);
        force_tgt.addVirialIntra( -calcVirialEPI(F_kj, R_kj) );
    } else {
        std::cerr << "   id_tgt = " << id_tgt << "\n"
                  << "        i = " << id_i
                  <<       "  j = " << id_j
                  <<       "  k = " << id_k << std::endl;
        throw std::invalid_argument("id_tgt is not match to id_i, id_j, or id_k.");
    }
}

//------ harminic torsion potential  (i-j-k-l form)
//  shape:      i
//              |
//              j---k     ---- axis ----
//              |   |
//   (improper) l   l (dihedral)
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

    if( std::abs(coef.theta0)            <= 1.e-5 ||
        std::abs(coef.theta0 - Unit::pi) <= 1.e-5 ){
        //--- pass checking
    } else {
        throw std::invalid_argument("theta0 = " + std::to_string(coef.theta0)
                                    + ", must be 0.0 or pi (0.0 or 180.0 in degree).");
    }

    PS::F64vec R_ij = pos_i - pos_j;
    PS::F64vec R_ik = pos_i - pos_k;
    PS::F64vec R_jl = pos_j - pos_l;
    PS::F64vec R_kj = pos_k - pos_j;
    PS::F64vec R_kl = pos_k - pos_l;

    PS::F64vec R_a     = VEC_EXT::cross(R_ij, R_kj);
    PS::F64    r_a     = VEC_EXT::norm(R_a);
    PS::F64    r_a_inv = 1.0/r_a;

    PS::F64vec R_b     = VEC_EXT::cross(R_kj, R_kl);
    PS::F64    r_b     = VEC_EXT::norm(R_b);
    PS::F64    r_b_inv = 1.0/r_b;

    PS::F64 cos_tmp = VEC_EXT::dot(R_a, R_b)*(r_a_inv*r_b_inv);
    if(cos_tmp < -1.0) cos_tmp = -1.0;
    if(cos_tmp >  1.0) cos_tmp =  1.0;

    PS::F32 sign = VEC_EXT::dot(R_kj, VEC_EXT::cross(R_a, R_b));
    if(sign < 0.0){
        sign = -1.0;
    } else {
        sign = 1.0;
    }

    PS::F64 phi = -sign*std::abs(std::acos(cos_tmp));
    PS::F64 sin_p = std::sin(phi);
    PS::F64 cos_p = std::cos(phi);

    PS::F64 cos_eq = std::cos(coef.theta0);
    PS::F64 sin_tmp, sin_1, sin_2;
    switch (coef.n_min) {
        case 1:
            sin_tmp = cos_eq;
        break;

        case 2:
            sin_tmp = 2.0*cos_p*cos_eq;
        break;

        case 3:
            sin_tmp = (-4.0*sin_p*sin_p + 3.0)*cos_eq;
        break;

        case 4:
            sin_tmp = 4.0*cos_p*(2.0*cos_p*cos_p - 1.0)*cos_eq;
        break;

        case 6:
            sin_1   = 4.0*cos_p*(2.0*cos_p*cos_p - 1.0)*std::cos(2.0*phi);
            sin_2   = 2.0*cos_p*std::cos(4.0*phi);
            sin_tmp = (sin_1 + sin_2)*cos_eq;
        break;

        default:
            sin_tmp = std::sin( PS::F64(coef.n_min)*phi - coef.theta0 )/sin_p;
            //sin_tmp = std::sin( PS::F64(coef.n_min)*(phi - coef.theta0) )/sin_p;
    }

    PS::F64 eng_tmp   = (1.0/4.0)*0.5*coef.k*(1.0 - std::cos( PS::F64(coef.n_min)*phi - coef.theta0 ));
    //PS::F64 eng_tmp   = (1.0/4.0)*0.5*coef.k*(1.0 - std::cos( PS::F64(coef.n_min)*(phi - coef.theta0) ));
    PS::F64 force_tmp = -0.5*coef.k*PS::F64(coef.n_min)*sin_tmp;

    PS::F64vec F_a = (R_b*r_b_inv - R_a*r_a_inv*cos_p)*r_a_inv;
    PS::F64vec F_b = (R_a*r_a_inv - R_b*r_b_inv*cos_p)*r_b_inv;

    force_tgt.addPotTorsion(eng_tmp);

    PS::F64vec F_tmp;
    if(       id_tgt == id_i){
        F_tmp = force_tmp*VEC_EXT::cross(F_a, R_kj);
        force_tgt.addForceIntra(F_tmp);
        force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, R_ij) );
    } else if(id_tgt == id_j){
        F_tmp = force_tmp*(VEC_EXT::cross(F_a, R_ik) + VEC_EXT::cross(-F_b, R_kl));
        force_tgt.addForceIntra(F_tmp);
        //--- virial = 0.
    } else if(id_tgt == id_k){
        F_tmp = force_tmp*(VEC_EXT::cross(-F_a, R_ij) + VEC_EXT::cross(F_b, R_jl));
        force_tgt.addForceIntra(F_tmp);
        force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, R_kj) );
    } else if(id_tgt == id_l){
        F_tmp = force_tmp*VEC_EXT::cross(F_b, R_kj);
        force_tgt.addForceIntra(F_tmp);
        force_tgt.addVirialIntra( 2.0*calcVirialEPI(F_tmp, -R_jl) );
    } else {
        std::cerr << "   id_tgt = " << id_tgt << "\n"
                  << "        i = " << id_i
                  <<       ", j = " << id_j
                  <<       ", k = " << id_k
                  <<       ", l = " << id_l << std::endl;
        throw std::invalid_argument("id_tgt is not match to id_i, id_j, id_k, or id_l.");
    }
}
