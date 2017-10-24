//***************************************************************************************
//  This code is base class for FullParticle and EssentialParticle.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <vector>
#include <string>

#include <particle_simulator.hpp>

#include "md_enum.hpp"


//--- base classes for particle data ----------------------------------------------------
//------ identifier
class AtomID{
protected:
    PS::S64 atom_id = -1;  // atom id
    PS::S64 mol_id  = -1;  // molecular id

public:

    //--- member functions
    void setAtomID(const PS::S64 &id){ this->atom_id = id; }
    void setMolID( const PS::S64 &id){ this->mol_id  = id; }
    inline PS::S64 getAtomID() const { return this->atom_id; }
    inline PS::S64 getMolID()  const { return this->mol_id; }

    //--- copy function
    template<class Tptcl>
    void copyAtomID(const Tptcl &fp){
        this->atom_id  = fp.getAtomID();
        this->mol_id   = fp.getMolID();
    }
};

//------ sub identifier
class AtomType{
protected:
    AtomName atom_type = static_cast<AtomName>(-1);  // illigal value
    MolName  mol_type  = static_cast<MolName>(-1);

public:
    void setAtomType(const AtomName &type){ this->atom_type = type; }
    void setAtomType(const std::string &str){
        this->setAtomType(ENUM::which_AtomName(str));
    }
    void setMolType(const MolName &type){ this->mol_type = type; }
    void setMolType(const std::string &str){
        this->setMolType(ENUM::which_MolName(str));
    }

    inline AtomName getAtomType() const { return this->atom_type; }
    inline MolName  getMolType()  const { return this->mol_type; }

    template<class Tptcl>
    void copyAtomType(const Tptcl &fp){
        this->atom_type = fp.getAtomType();
        this->mol_type  = fp.getMolType();
    }
};


//------ position
class AtomPos{
protected:
    PS::F64vec pos = 0.0;

public:

    //--- member functions
    void setPos(const PS::F64vec &pos_new){ this->pos = pos_new; }
    inline PS::F64vec getPos() const { return this->pos; }

    //--- copy function
    template<class Tptcl>
    void copyAtomPos(const Tptcl &fp){
        this->pos = fp.getPos();
    }
};

//------ mass & velocity
class AtomVel{
protected:
    PS::F64    mass = 0.0;  // mass of atom
    PS::F64vec vel  = 0.0;  // velocity of atom

public:

    void setMass(const PS::F64    &mass_new){ this->mass = mass_new; }
    void setVel( const PS::F64vec &vel_new){  this->vel  = vel_new; }
    inline PS::F64    getMass() const { return this->mass; }
    inline PS::F64vec getVel()  const { return this->vel; }

    inline PS::F64 getKinetic() const {
        PS::F64 v2 = (this->vel*this->vel);
        return 0.5*v2/this->mass;
    }
};

//------ permanent charge
class AtomCharge{
protected:
    PS::F64 charge = 0.0;

public:

    void setCharge(const PS::F64 &q){ this->charge = q; }
    inline PS::F64 getCharge()             const { return this->charge; }
    inline PS::F64 getChargeParticleMesh() const { return this->charge; }

    template<class Tptcl>
    void copyAtomCharge(const Tptcl &fp){
        this->charge = fp.getCharge();
    }
};

//------ VDW parameters
class AtomVDW{
protected:
    PS::F64 vdw_r = 0.0;
    PS::F64 vdw_d = 0.0;

public:

    void setVDW_R(const PS::F64 &vdw_r){ this->vdw_r = vdw_r; }
    void setVDW_D(const PS::F64 &vdw_d){ this->vdw_d = vdw_d; }
    inline PS::F64 getVDW_R() const { return this->vdw_r; }
    inline PS::F64 getVDW_D() const { return this->vdw_d; }

    template<class Tptcl>
    void copyAtomVDW(const Tptcl &fp){
        this->vdw_r = fp.getVDW_R();
        this->vdw_d = fp.getVDW_D();
    }
};


//--- Force class -----------------------------------------------------------------------
//------ This class has the result of force. (used as temporary data)
//------ This class is the subset of Full Particle class.

//------ LJ interaction with cutoff length (simple cutoff)
class ForceLJ {
protected:
    PS::F64vec force_LJ  = 0.0;
    PS::F64vec virial_LJ = 0.0;
    PS::F64    pot_LJ    = 0.0;

public:
    void clearForceLJ(){
        this->force_LJ  = 0.0;
        this->virial_LJ = 0.0;
        this->pot_LJ    = 0.0;
    }
    void clear(){ this->clearForceLJ(); }

    inline PS::F64vec getForceLJ()  const { return this->force_LJ;  }
    inline PS::F64vec getVirialLJ() const { return this->virial_LJ; }
    inline PS::F64    getPotLJ()    const { return this->pot_LJ;    }
    inline void addForceLJ( const PS::F64vec &force){  this->force_LJ  += force;  }
    inline void addVirialLJ(const PS::F64vec &virial){ this->virial_LJ += virial; }
    inline void addPotLJ(   const PS::F64    &pot){    this->pot_LJ    += pot;    }

    template <class T>
    void copyForceLJ(const T &f){
        this->force_LJ  = f.getForceLJ();
        this->virial_LJ = f.getVirialLJ();
        this->pot_LJ    = f.getPotLJ();
    }
    void copyFromForce(const ForceLJ &f){
        this->copyForceLJ(f);
    }
};

//------ coulomb interaction
class ForceCoulomb {
protected:
    PS::F64vec field_coulomb = 0.0;
    PS::F64      pot_coulomb = 0.0;

public:
    void clearForceCoulomb(){
        field_coulomb = 0.0;
          pot_coulomb = 0.0;
    }
    void clear(){ this->clearForceCoulomb(); }

    inline PS::F64vec getFieldCoulomb() const { return this->field_coulomb; }
    inline PS::F64    getPotCoulomb()   const { return this->pot_coulomb;   }
    inline void addFieldCoulomb(const PS::F64vec &f){ this->field_coulomb += f; }
    inline void addPotCoulomb(  const PS::F64    &p){ this->pot_coulomb   += p; }

    template <class T>
    void copyForceCoulomb(const T &f){
        this->field_coulomb  = f.getFieldCoulomb();
        this->pot_coulomb    = f.getPotCoulomb();
    }
    void copyFromForce(const ForceCoulomb &f){
        this->copyForceCoulomb(f);
    }
};

//--- Intramoleecular interaction
class ForceIntra {
protected:
    PS::F64vec force_intra  = 0.0;
    PS::F64vec virial_intra = 0.0;
    PS::F64    pot_bond     = 0.0;
    PS::F64    pot_angle    = 0.0;
    PS::F64    pot_torsion  = 0.0;

public:
    void clearForceIntra(){
        this->force_intra  = 0.0;
        this->virial_intra = 0.0;
        this->pot_bond     = 0.0;
        this->pot_angle    = 0.0;
        this->pot_torsion  = 0.0;
    }
    void clear(){ this->clearForceIntra(); }

    inline PS::F64vec getForceIntra()  const { return this->force_intra;  }
    inline PS::F64vec getVirialIntra() const { return this->virial_intra; }
    inline PS::F64    getPotBond()     const { return this->pot_bond;     }
    inline PS::F64    getPotAngle()    const { return this->pot_angle;    }
    inline PS::F64    getPotTorsion()  const { return this->pot_torsion;  }
    inline void addForceIntra( const PS::F64vec &f){ this->force_intra  += f; }
    inline void addVirialIntra(const PS::F64vec &v){ this->virial_intra += v; }
    inline void addPotBond(   const PS::F64 &p){ this->pot_bond    += p; }
    inline void addPotAngle(  const PS::F64 &p){ this->pot_angle   += p; }
    inline void addPotTorsion(const PS::F64 &p){ this->pot_torsion += p; }

    template <class T>
    void copyForceIntra(const T &f){
        this->force_intra  = f.getForceIntra();
        this->virial_intra = f.getVirialIntra();
        this->pot_bond     = f.getPotBond();
        this->pot_angle    = f.getPotAngle();
        this->pot_torsion  = f.getPotTorsion();
    }
    void copyFromForce(const ForceIntra &f){
        this->copyForceIntra(f);
    }
};
