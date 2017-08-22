//***************************************************************************************
//  This program is the normalize interface of position data.
//    convert pos{0~boxdh__x, 0~boxdh__y, 0~boxdh__z} and pos{0~1, 0~1, 0~1}.
//    This code is used by "md_fdps_main.cpp"
//***************************************************************************************
#pragma once

#include <cmath>
#include <algorithm>

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>


namespace Normalize {

    //--- system box size
    PS::F64vec boxdh_{1.0, 1.0, 1.0};

    //--- inverse value of system box size
    PS::F64vec boxdh_inv_{1.0, 1.0, 1.0};

    //--- parameter for Particle Mesh
    //------ mesh size in normalized spase for Particle Mesh
    //------ "SIZE_OF_MESH" is defined in "$(PS_DIR)/src/particle_mesh/param_fdps.h".
    //------ recommended value of "SIZE_OF_MESH" is N^(1/3)/2. (N is the total number of particles)
    PS::F64vec coefPMForce{1.0, 1.0, 1.0};
    PS::F64vec mesh_dSize{ 1.0/SIZE_OF_MESH,
                           1.0/SIZE_OF_MESH,
                           1.0/SIZE_OF_MESH };

    //--- internal functions
    void setNormalizePosition_() {
        boxdh_inv_.x = 1.0/boxdh_.x;
        boxdh_inv_.y = 1.0/boxdh_.y;
        boxdh_inv_.z = 1.0/boxdh_.z;
    }
    void setCoefPM_(){
        coefPMForce.x = 1.0*(boxdh_inv_.x*boxdh_inv_.x);
        coefPMForce.y = 1.0*(boxdh_inv_.y*boxdh_inv_.y);
        coefPMForce.z = 1.0*(boxdh_inv_.z*boxdh_inv_.z); // squared
    }

    //--- access functions (must be call only below functions)
    //------ initialize or setting the real domain size
    void setBoxSize(const PS::F64vec & box){
        boxdh_ = box;
        setNormalizePosition_();
        setCoefPM_();
    }

    PS::F64vec getBoxSize(){
        return boxdh_;
    }
    PS::F64vec getBoxSizeInv(){
        return boxdh_inv_;
    }
    PS::F64 getVol(){
        return (boxdh_.x*boxdh_.y*boxdh_.z);
    }
    PS::F64 getVolInv(){
        return (boxdh_inv_.x*boxdh_inv_.y*boxdh_inv_.z);
    }

    void broadcast_boxSize(const PS::S32 &ref){
        PS::Comm::broadcast(&boxdh_, 1, ref);
        setNormalizePosition_();
        setCoefPM_();
    }

    //------ convert real pos to norm pos
    inline PS::F64vec normPos(const PS::F64vec &pos){
        PS::F64vec tmp = 0.0;
                   tmp.x = pos.x * boxdh_inv_.x;
                   tmp.y = pos.y * boxdh_inv_.y;
                   tmp.z = pos.z * boxdh_inv_.z;
        return tmp;
    }
    //------ convert norm pos to real pos
    inline PS::F64vec realPos(const PS::F64vec &pos){
        PS::F64vec tmp = 0.0;
                   tmp.x = pos.x * boxdh_.x;
                   tmp.y = pos.y * boxdh_.y;
                   tmp.z = pos.z * boxdh_.z;
        return tmp;
    }
    //------ convert norm PM force to real PM force
    inline PS::F64vec realPMForce(const PS::F64vec &Fpm){
        PS::F64vec tmp   = 0.0;
                   tmp.x = Fpm.x * coefPMForce.x;
                   tmp.y = Fpm.y * coefPMForce.y;
                   tmp.z = Fpm.z * coefPMForce.z;
        return tmp;
    }
    inline PS::F64 realPMPotential(const PS::F64 &Fpm){
        //--- treat cubic system only
        assert(boxdh_.x == boxdh_.y &&
               boxdh_.x == boxdh_.z);
        return boxdh_inv_.x*Fpm;
    }
    //------ convert real cutoff length to norm cutoff length
    inline PS::F64 normCutOff(const PS::F64 &realCutOff){
        PS::F64 len_inv = std::max( {boxdh_inv_.x,
                                     boxdh_inv_.y,
                                     boxdh_inv_.z},
                                     std::greater<PS::F64>() );
        return realCutOff*len_inv;
    }
    //------ convert norm cutoff length to real cutoff length
    inline PS::F64 realCutOff(const PS::F64 &normCutOff){
        PS::F64 len = std::min( {boxdh_.x,
                                 boxdh_.y,
                                 boxdh_.z},
                                 std::less<PS::F64>() );
        return normCutOff*len;
    }
    inline PS::F64 normCutOff_PM(){
        PS::F64 r = std::max( {mesh_dSize.x,
                               mesh_dSize.y,
                               mesh_dSize.z},
                              std::greater<PS::F64>() );
        return r*3.0;
    }

    //------ adjustment position into normalized periodic boundary system
    inline PS::F64vec periodicAdjustNorm(const PS::F64vec &pos_norm){
        PS::F64vec pos_new = pos_norm;
        if(pos_new.x < 0.0)  pos_new.x += 1.0;
        if(pos_new.y < 0.0)  pos_new.y += 1.0;
        if(pos_new.z < 0.0)  pos_new.z += 1.0;
        if(pos_new.x >= 1.0) pos_new.x -= 1.0;
        if(pos_new.y >= 1.0) pos_new.y -= 1.0;
        if(pos_new.z >= 1.0) pos_new.z -= 1.0;
      //  pos_new.x -= round(pos_new.x - 0.5000000000001);
      //  pos_new.y -= round(pos_new.y - 0.5000000000001);
      //  pos_new.z -= round(pos_new.z - 0.5000000000001);
        return pos_new;
    }
    inline PS::F64vec periodicAdjustReal(const PS::F64vec &pos_real){
        PS::F64vec pos_new = pos_real;
        if(pos_new.x < 0.0)       pos_new.x += boxdh_.x;
        if(pos_new.y < 0.0)       pos_new.y += boxdh_.y;
        if(pos_new.z < 0.0)       pos_new.z += boxdh_.z;
        if(pos_new.x >= boxdh_.x) pos_new.x -= boxdh_.x;
        if(pos_new.y >= boxdh_.y) pos_new.y -= boxdh_.y;
        if(pos_new.z >= boxdh_.z) pos_new.z -= boxdh_.z;
        return pos_new;
    }

    //------ check the position is out of space or not
    inline bool checkPosInSpace(const PS::F64vec &pos_norm){
        if(pos_norm.x <  0.0 ||
           pos_norm.y <  0.0 ||
           pos_norm.z <  0.0 ||
           pos_norm.x >= 1.0 ||
           pos_norm.y >= 1.0 ||
           pos_norm.z >= 1.0   ){
            return false;
        } else {
            return true;
        }
    }

    //------ drift in normalized space
    inline PS::F64vec normDrift(const PS::F64vec &move){
        return normPos(move);
    }

    //------ length converter
    inline PS::F64 normXLen(const PS::F64 &x_real){
        PS::F64 x_norm = x_real * boxdh_inv_.x;
        return x_norm;
    }
    inline PS::F64 realXLen(const PS::F64 &x_norm){
        PS::F64 x_real = x_norm * boxdh_.x;
        return x_real;
    }
    inline PS::F64 normYLen(const PS::F64 &y_real){
        PS::F64 y_norm = y_real * boxdh_inv_.y;
        return y_norm;
    }
    inline PS::F64 realYLen(const PS::F64 &y_norm){
        PS::F64 y_real = y_norm * boxdh_.y;
        return y_real;
    }
    inline PS::F64 normZLen(const PS::F64 &z_real){
        PS::F64 z_norm = z_real * boxdh_inv_.z;
        return z_norm;
    }
    inline PS::F64 realZLen(const PS::F64 &z_norm){
        PS::F64 z_real = z_norm * boxdh_.z;
        return z_real;
    }
}
