//***************************************************************************************
//  This is basic kick and drift routine for atom.
//***************************************************************************************
#pragma once

#include <utility>
#include <tuple>
#include <string>
#include <sstream>

#include <particle_simulator.hpp>

#include "enum_model.hpp"

namespace Analysis {

    //--- value filter
    template <class T>
    class TargetValue {
    private:
        bool effective = false;
        bool inverse   = false;
        T    data;

    public:
        using value_type = T;

        void assign(const T &value, const bool inverse = false){
            this->effective = true;
            this->inverse   = inverse;
            this->data      = value;
        }
        void clear(){ this->effective = false; }

        bool isEffective() const { return this->effective; }
        bool isInverse()   const { return this->inverse;   }
        T    get()         const { return this->data;      }

        TargetValue() = default;
        TargetValue(const T &value, const bool inverse = false){
            this->assign(value, inverse);
        }
        TargetValue(const std::pair<T,bool> &tgt){
            this->assign(tgt.first, tgt.second);
        }
        TargetValue(const TargetValue &tgt){
            this->assign(tgt.get(), tgt.isInverse());
        }

        TargetValue& operator = (const T &value){
            this->assign(value);
            return *this;
        }
        TargetValue& operator = (const std::pair<T,bool> &tgt){
            this->assign(tgt.first, tgt.second);
            return *this;
        }
        TargetValue& operator = (const TargetValue &tgt){
            this->assign(tgt.get(), tgt.isInverse());
            return *this;
        }

        bool isMatch(const T &value) const {
            bool match = (this->data == value);
            if(this->inverse) match = not match;
            return (!(this->effective) || match);
        }

        //--- for duplication check
        bool operator == (const TargetValue &tgt) const {
            return (this->effective == tgt.isEffective() &&
                    this->inverse   == tgt.isInverse()   &&
                    this->data      == tgt.get()           );
        }
        bool operator != (const TargetValue &tgt) const {
            return !(*this == tgt);
        }

        std::string str() const {
            std::ostringstream oss;
            if(this->isEffective()){
                if(this->isInverse()){
                    oss << "not ";
                }
                oss << this->get();
            } else {
                oss << "[wildcard]";
            }
            return oss.str();
        }

    };

    template <class T>
    inline std::ostream& operator << (std::ostream& s, const TargetValue<T> &v){
        s << v.str();
        return s;
    }


    //--- generic atom target
    class TargetAtom {
    public:
        TargetValue<AtomName> atom_type;
        TargetValue<MolName>  mol_type;

        template <class Tptcl>
        bool isMatch(const Tptcl &atom) const {
            return (this->atom_type.isMatch( atom.getAtomType() ) &&
                    this->mol_type.isMatch(  atom.getMolType()  )   );
        }

        std::string str() const {
            std::ostringstream oss;
            oss << "atom: "
                << "AtomName = " << this->atom_type
                << ", "
                << "MolName = "  << this->mol_type;
            return oss.str();
        }
    };

    inline std::ostream& operator << (std::ostream& s, const TargetAtom &atom){
        s << atom.str();
        return s;
    }

}
