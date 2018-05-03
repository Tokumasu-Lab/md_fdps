//***************************************************************************************
//  This is basic kick and drift routine for atom.
//***************************************************************************************
#pragma once

#include <pair>
#include <tuple>
#include <string>
#include <sstream>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "atom_class_base.hpp"
#include "file_IO_pos.hpp"

namespace analysis {

    class Atom_FP :
      public AtomType,
      public AtomID,
      public AtomPos<PS::F32> {
      private:
        PS::F32vec trj;

      public:

        inline void       clearTrj()                           { this->trj = 0.0;   }
        inline void       addTrj(const PS::F32vec &move)       { this->trj += move; }
        inline PS::F32vec getTrj()                       const { return this->trj;  }

        std::string str() const {
            std::ostringstream oss;

            oss << "    AtomID   : " << this->getAtomID()   << "\n";
            oss << "    MolID    : " << this->getMolID()    << "\n";
            oss << "    AtomType : " << this->getAtomType() << "\n";
            oss << "    MolType  : " << this->getMolType()  << "\n";
            oss << "    Pos      : ("  << this->getPos().x
                               << ", " << this->getPos().y
                               << ", " << this->getPos().z << ")\n";

            return oss.str();
        }

        //--- interface for file I/O (read only)
        //------ for pos
        void read_pos_ascii( FILE *fp) { FILE_IO::Pos::read_atom(fp, *this); }
    };

}
