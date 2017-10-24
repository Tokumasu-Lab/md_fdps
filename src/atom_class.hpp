//***************************************************************************************
//  This code is the FullParticle class and EssentialParticle class.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "unit.hpp"
#include "md_enum.hpp"
#include "atom_class_base.hpp"


//--- derived force class
class Force_FP :
  public ForceLJ,
  public ForceCoulomb,
  public ForceIntra {
  public:

    //--- copy results of calculated force
    template <class Tforce>
    void copyFromForce(const Tforce &force){
        this->copyForceLJ(force);
        this->copyForceCoulomb(force);
        this->copyForceIntra(force);
    }

    //--- clear dinamic data
    void clear() {
        this->clearForceIntra();
        this->clearForceLJ();
        this->clearForceCoulomb();
    }
};

//--- Full Particle class ---------------------------------------------------------------
//------ This class should contain all properties of the particle.
//------ This class is based on base classes.
class Atom_FP :
  public AtomID,
  public AtomType,
  public AtomPos,
  public AtomVel,
  public AtomCharge,
  public AtomVDW,
  public Force_FP {
  public:

    //--- output interaction result
    inline PS::F64vec getForce() const {
        return   this->getForceIntra()
               + this->getForceLJ()
               + this->getFieldCoulomb()*this->getCharge();
    }
    inline PS::F64vec getVirial() const {
        return   this->getVirialIntra()
               + this->getVirialLJ()
               + this->getPotCoulomb()*this->getCharge();
    }

    //--- copy model property from molecular model template (using FP class)
    void copyFromModelTemplate(const PS::S64 &a_id,
                               const PS::S64 &m_id,
                               const Atom_FP &t){
        assert(a_id >= 0 && m_id >= 0);
        *this = t;

        //--- id
        this->setAtomID(a_id);
        this->setMolID(m_id);

        //--- preprocess of interaction parameters
        this->setMass(  t.getMass()/Unit::mass_C );
        this->setVDW_R( 0.5*t.getVDW_R() );
        this->setVDW_D( std::sqrt( t.getVDW_D() ) );
        this->setCharge( Unit::coef_coulomb*t.getCharge() );

    //    std::cout << "  ID= "   << this->getAtomID()
    //              << "  type= " << this->getMolType()
    //              << ", "       << this->getAtomType() << std::endl;
    }

    //--- prototype for file I/O
    void writeAscii(FILE *fp) const;
    void readAscii(FILE *fp);
};


//--- Essential Particle class ----------------------------------------------------------
//------ These classes have the meta-data for force calculation.
//------ They are the subset of Full Particle class.
//------ They are based on base classes.

class Atom_EP :
  public AtomID,
  public AtomType,
  public AtomPos,
  public AtomCharge,
  public AtomVDW {
  private:
    static PS::F64 Rcut;
    static PS::F64 Rcut_LJ;
    static PS::F64 Rcut_coulomb;

  public:
    static void setRcut_LJ(const PS::F64 &r){
        Atom_EP::Rcut_LJ = r;
        Atom_EP::Rcut    = std::max(r, Atom_EP::Rcut_coulomb);
    }
    static void setRcut_coulomb(const PS::F64 &r){
        Atom_EP::Rcut_coulomb = r;
        Atom_EP::Rcut         = std::max(r, Atom_EP::Rcut_LJ);
    }

    static PS::F64 getRSearch()      { return Atom_EP::Rcut; }
    static PS::F64 getRcut_LJ()      { return Atom_EP::Rcut_LJ; }
    static PS::F64 getRcut_coulomb() { return Atom_EP::Rcut_coulomb; }

    template <class T>
    void copyFromFP(const T &fp){
        this->copyAtomID(fp);
        this->copyAtomType(fp);
        this->copyAtomPos(fp);
        this->copyAtomCharge(fp);
        this->copyAtomVDW(fp);
    }
};
PS::F64 Atom_EP::Rcut         = 0.0;
PS::F64 Atom_EP::Rcut_LJ      = 0.0;
PS::F64 Atom_EP::Rcut_coulomb = 0.0;
