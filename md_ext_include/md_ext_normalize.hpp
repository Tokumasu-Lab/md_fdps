/**************************************************************************************************/
/**
* @file  md_ext_normalize.hpp
* @brief normalize interface for position, potential, and force.
*/
/**************************************************************************************************/
#pragma once

#include <cmath>
#include <algorithm>

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>


/**
* @brief normalize interface for position, potential, and force.
*/
namespace Normalize {

    //! @brief system box size.
    PS::F64vec boxdh_{1.0, 1.0, 1.0};

    //! @brief inverse value of system box size.
    PS::F64vec boxdh_inv_{1.0, 1.0, 1.0};

    //! @brief coefficient for force from PS::PM::ParticleMesh.
    PS::F64vec coefPMForce{1.0, 1.0, 1.0};

    /**
    * @brief volume of a mesh.
    * @details "SIZE_OF_MESH" is defined in "$(PS_DIR)/src/particle_mesh/param_fdps.h".
    * @details The recommended value of "SIZE_OF_MESH" is N^(1/3)/2. (N is the total number of particles)
    */
    PS::F64vec mesh_dSize{ 1.0/SIZE_OF_MESH,
                           1.0/SIZE_OF_MESH,
                           1.0/SIZE_OF_MESH };

   /**
   * @brief set boxdh_inv_. internal use only.
   */
    void setNormalizePosition_() {
        boxdh_inv_.x = 1.0/boxdh_.x;
        boxdh_inv_.y = 1.0/boxdh_.y;
        boxdh_inv_.z = 1.0/boxdh_.z;
    }

    /**
    * @brief set coefPMForce. internal use only.
    */
    void setCoefPM_(){
        coefPMForce.x = 1.0*(boxdh_inv_.x*boxdh_inv_.x);
        coefPMForce.y = 1.0*(boxdh_inv_.y*boxdh_inv_.y);
        coefPMForce.z = 1.0*(boxdh_inv_.z*boxdh_inv_.z); // squared
    }


    //--- access functions (must be call only below functions)
    /**
    * @brief setting the real domain size.
    * @param[in] box size of domain.
    */
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

    /**
    * @brief broadcast setting of domain size.
    * @param[in] ref ID of source process.
    */
    void broadcast_boxSize(const PS::S32 &ref){
        PS::Comm::broadcast(&boxdh_, 1, ref);
        setNormalizePosition_();
        setCoefPM_();
    }

    //! @brief convert real pos to normalized pos.
    inline PS::F64vec normPos(const PS::F64vec &pos){
        PS::F64vec tmp = 0.0;
                   tmp.x = pos.x * boxdh_inv_.x;
                   tmp.y = pos.y * boxdh_inv_.y;
                   tmp.z = pos.z * boxdh_inv_.z;
        return tmp;
    }
    //! @brief convert normalized pos to real pos.
    inline PS::F64vec realPos(const PS::F64vec &pos){
        PS::F64vec tmp = 0.0;
                   tmp.x = pos.x * boxdh_.x;
                   tmp.y = pos.y * boxdh_.y;
                   tmp.z = pos.z * boxdh_.z;
        return tmp;
    }
    //! @brief convert normalized PM force to real PM force.
    inline PS::F64vec realPMForce(const PS::F64vec &Fpm){
        PS::F64vec tmp   = 0.0;
                   tmp.x = Fpm.x * coefPMForce.x;
                   tmp.y = Fpm.y * coefPMForce.y;
                   tmp.z = Fpm.z * coefPMForce.z;
        return tmp;
    }
    //! @brief convert normalized PM potential to real PM potential.
    //! @details the domain must be cubic (X=Y=Z).
    inline PS::F64 realPMPotential(const PS::F64 &Fpm){
        //--- treat cubic system only
        assert(boxdh_.x == boxdh_.y &&
               boxdh_.x == boxdh_.z);
        return boxdh_inv_.x*Fpm;
    }
    //! @brief convert real cutoff length to normalized cutoff length.
    inline PS::F64 normCutOff(const PS::F64 &realCutOff){
        PS::F64 len_inv = std::max( {boxdh_inv_.x,
                                     boxdh_inv_.y,
                                     boxdh_inv_.z},
                                     std::greater<PS::F64>() );
        return realCutOff*len_inv;
    }
    //! @brief convert normalized cutoff length to real cutoff length.
    inline PS::F64 realCutOff(const PS::F64 &normCutOff){
        PS::F64 len = std::min( {boxdh_.x,
                                 boxdh_.y,
                                 boxdh_.z},
                                 std::less<PS::F64>() );
        return normCutOff*len;
    }
    //! @brief get cutoff length for short-range interaction with PS::PM::ParticleMesh.
    //! @details cutoff length is 3.0/SIZE_OF_MESH.
    inline PS::F64 normCutOff_PM(){
        PS::F64 r = std::max( {mesh_dSize.x,
                               mesh_dSize.y,
                               mesh_dSize.z},
                              std::greater<PS::F64>() );
        return r*3.0;
    }

    //! @brief adjustment position into normalized periodic boundary system.
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
    //! @brief adjustment position into real periodic boundary system
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

    //! @brief check the position is out of domain or not
    //! @param[in] pos_norm normalized position.
    //! @return bool True: out of domain. False: exists in domain.
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

    //! @brief convert real drift dispracement to normalized drift.
    inline PS::F64vec normDrift(const PS::F64vec &move){
        return normPos(move);
    }

    //------ length converter
    //! @brief convert real X length to normalized X length.
    inline PS::F64 normXLen(const PS::F64 &x_real){
        PS::F64 x_norm = x_real * boxdh_inv_.x;
        return x_norm;
    }
    //! @brief convert normalized X length to real X length.
    inline PS::F64 realXLen(const PS::F64 &x_norm){
        PS::F64 x_real = x_norm * boxdh_.x;
        return x_real;
    }
    //! @brief convert real Y length to normalized Y length.
    inline PS::F64 normYLen(const PS::F64 &y_real){
        PS::F64 y_norm = y_real * boxdh_inv_.y;
        return y_norm;
    }
    //! @brief convert normalized Y length to real Y length.
    inline PS::F64 realYLen(const PS::F64 &y_norm){
        PS::F64 y_real = y_norm * boxdh_.y;
        return y_real;
    }
    //! @brief convert real Z length to normalized Z length.
    inline PS::F64 normZLen(const PS::F64 &z_real){
        PS::F64 z_norm = z_real * boxdh_inv_.z;
        return z_norm;
    }
    //! @brief convert normalized Z length to real Z length.
    inline PS::F64 realZLen(const PS::F64 &z_norm){
        PS::F64 z_real = z_norm * boxdh_.z;
        return z_real;
    }
}
