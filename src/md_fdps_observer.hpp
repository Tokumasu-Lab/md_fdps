//***************************************************************************************
//  This is observer functions.
//***************************************************************************************
#pragma once

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cassert>

#include <particle_simulator.hpp>

#include "md_ext_normalize.hpp"

#include "md_fdps_unit.hpp"


namespace Observer {

    class FilePrinter {
    private:
        PS::S32       rank = 0;
        std::string   name;
        std::ofstream ofs;

    public:
        void file_init(const std::string &name, const PS::S32 &rank){
            if( this->ofs.is_open() ){ this->ofs.close(); }
            this->rank = rank;
            this->name = name;

            if(PS::Comm::getRank() != this->rank) return;
            this->ofs.open(name, std::ios::trunc);
        }

        void print(const std::string &str){
            if(PS::Comm::getRank() != this->rank) return;
            this->ofs << str;
        }

        std::string file_name() const { return this->name; }
        PS::S32     getRank()   const { return this->rank; }
    };

    class Energy {
    private:
        FilePrinter printer;
        PS::S64     start;

    public:
        PS::F64 bond;
        PS::F64 angle;
        PS::F64 torsion;

        PS::F64 vdw;
        PS::F64 coulomb;

        PS::F64    kin;
        PS::F64    ext_sys;
        PS::F64vec virial;

        PS::F64 density;
        PS::F64 n_atom;    // buffer for property calculation

        //--- manipulator
        void clear(){
            this->bond    = 0.0;
            this->angle   = 0.0;
            this->torsion = 0.0;

            this->vdw     = 0.0;
            this->coulomb = 0.0;

            this->kin     = 0.0;
            this->ext_sys = 0.0;
            this->virial  = 0.0;

            this->density = 0.0;
            this->n_atom  = 0.0;
        }
        void plus(const Energy &rv, const PS::F64 sign = 1.0){
            this->bond    += sign*rv.bond;
            this->angle   += sign*rv.angle;
            this->torsion += sign*rv.torsion;

            this->vdw     += sign*rv.vdw;
            this->coulomb += sign*rv.coulomb;

            this->kin     += sign*rv.kin;
            this->ext_sys += sign*rv.ext_sys;
            this->virial  += sign*rv.virial;

            this->density += sign*rv.density;
            this->n_atom  += sign*rv.n_atom;
        }
        void assign(const Energy &rv){
            this->clear();
            this->plus(rv);
        }
        void multiple(const PS::F64 &r){
            Energy buf;
            buf.plus(*this, r);
            this->assign(buf);
        }

        //--- constructor
        Energy(){
            this->clear();
        }
        Energy(const Energy &rv){
            this->assign(rv);
        }

        //--- operator
        Energy& operator = (const Energy &rv){
            this->assign(rv);
            return *this;
        }
        Energy operator + (const Energy &rv)        { this->plus(rv,  1.0); return *this; }
        Energy operator - (const Energy &rv)        { this->plus(rv, -1.0); return *this; }
        const Energy& operator += (const Energy &rv){ this->plus(rv,  1.0); return *this; }
        const Energy& operator -= (const Energy &rv){ this->plus(rv, -1.0); return *this; }

        PS::F64 inter() const {
            return   this->vdw
                   + this->coulomb;
        }
        PS::F64 intra() const {
            return   this->bond
                   + this->angle
                   + this->torsion;
        }
        PS::F64 kinetic() const {
            return this->kin;
        }
        PS::F64 total() const {
            return   this->inter()
                   + this->intra()
                   + this->kinetic()
                   + this->ext_sys;
        }

        //--- sampling from PS::ParticleSystem<FP>
        template <class Tpsys>
        void getEnergy(const Tpsys &psys){
            Energy  buf;
                    buf.clear();
            PS::S64 n_local  = psys.getNumberOfParticleLocal();

            for(PS::S64 i=0; i<n_local; ++i){
                buf.bond    += psys[i].getPotBond();
                buf.angle   += psys[i].getPotAngle();
                buf.torsion += psys[i].getPotTorsion();

                buf.vdw     += psys[i].getPotLJ();
                buf.coulomb += psys[i].getPotCoulomb()*psys[i].getCharge();

                PS::F64vec v  = psys[i].getVel();
                buf.kin      += 0.5*psys[i].getMass()*(v*v);

            //    buf.ext_sys  = 0.0;

                buf.virial  += psys[i].getVirial();

                buf.density += psys[i].getMass();   // total mass
            //    buf.n_atom   = 0;
            }
            buf.n_atom = PS::F64(n_local);

            //--- sumation
            buf.bond    = PS::Comm::getSum(buf.bond);
            buf.angle   = PS::Comm::getSum(buf.angle);
            buf.torsion = PS::Comm::getSum(buf.torsion);
            buf.vdw     = PS::Comm::getSum(buf.vdw);
            buf.coulomb = PS::Comm::getSum(buf.coulomb);
            buf.kin     = PS::Comm::getSum(buf.kin);
        //    buf.ext_sys = PS::Comm::getSum(buf.ext_sys);
            buf.virial  = PS::Comm::getSum(buf.virial);
            buf.density = PS::Comm::getSum(buf.density);
            buf.n_atom  = PS::Comm::getSum(buf.n_atom);

            buf.ext_sys = this->ext_sys;

            buf.density = buf.density*Normalize::getVolInv();   // mass -> density
            this->assign(buf);
        }

        std::string header() const {
            std::ostringstream oss;
            const size_t word_length = 16;

            oss << std::left << std::setw(12) << "  Time[fs]";
            oss << std::left << std::setw(12) << "  step";

            oss << std::left << std::setw(word_length) << "  Total";
            oss << std::left << std::setw(word_length) << "  kinetic";
            oss << std::left << std::setw(word_length) << "  bond";
            oss << std::left << std::setw(word_length) << "  angle";
            oss << std::left << std::setw(word_length) << "  torsion";
            oss << std::left << std::setw(word_length) << "  vdw";
            oss << std::left << std::setw(word_length) << "  coulomb";
            oss << std::left << std::setw(word_length) << "  virial_x";
            oss << std::left << std::setw(word_length) << "  virial_y";
            oss << std::left << std::setw(word_length) << "  virial_z";
            oss << std::left << std::setw(word_length) << "  ext_sys";

            oss << "\n";

            return oss.str();
        }
        std::string to_string(const PS::F64 &dt,
                              const PS::S64 &i_step) const {
            std::ostringstream oss;
            const size_t word_length = 16;

            oss << std::setw(12) << std::setprecision(8) << dt*PS::F64(i_step)*Unit::norm_time/Unit::femto_second;
            oss << std::setw(12) << i_step;

            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->total();
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->kinetic();
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->bond;
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->angle;
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->torsion;
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->vdw;
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->coulomb;
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->virial.x;
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->virial.y;
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->virial.z;
            oss << std::setw(word_length) << std::setprecision(8) << std::scientific << this->ext_sys;

            oss << "\n";

            return oss.str();
        }
        void display(const PS::S64 &i_step) const {
            if(PS::Comm::getRank() != 0) return;

            std::cout << "step " << i_step;

            std::cout << " | energy total= " << this->total();
            std::cout << " kin= "            << this->kinetic();
            std::cout << " intra= "          << this->intra();
            std::cout << " inter= "          << this->inter();

            std::cout << std::endl;
        }

        void file_init(const std::string &name,
                       const PS::S64      start = 0,
                       const PS::S32      proc  = 0){
            this->printer.file_init(name, proc);
            this->printer.print( this->header() );
            this->start = start;
        }
        void record(const PS::F64 &dt,
                    const PS::S64 &i_step) {
            if(i_step < this->start) return;
            this->printer.print( this->to_string(dt, i_step) );
        }
    };

    class Property {
    private:
        FilePrinter printer;
        PS::S64     start;

    public:
        PS::F64 temperature;
        PS::F64 density;
        PS::F64 pressure;

        PS::F64 press_temperature;
        PS::F64 press_virial;

        //--- manipulator
        void clear(){
            this->temperature = 0.0;
            this->density     = 0.0;
            this->pressure    = 0.0;

            this->press_temperature = 0.0;
            this->press_virial      = 0.0;
        }
        void plus(const Property &rv, const PS::F64 sign = 1.0){
            this->temperature += sign*rv.temperature;
            this->density     += sign*rv.density;
            this->pressure    += sign*rv.pressure;

            this->press_temperature += sign*rv.press_temperature;
            this->press_virial      += sign*rv.press_virial;
        }
        void assign(const Property &rv){
            this->clear();
            this->plus(rv);
        }
        void multiple(const PS::F64 &r){
            Property buf;
            buf.plus(*this, r);
            this->assign(buf);
        }

        //--- constructor
        Property() = default;
        Property(const Property& rv){
            this->assign(rv);
        }

        template <class Teng>
        void getProperty(const Teng &eng){
            //--- temperature [K]
            this->temperature = eng.kinetic()*(2.0/(3.0*eng.n_atom - 1.0))*Unit::norm_temp;

            //--- density [g/cm^3] = 10^-3 [kg/m^3]
            this->density     = eng.density*Unit::norm_dens*1.e-3;

            //--- pressure [Pa]
            this->press_temperature = Normalize::getVolInv()*Unit::norm_press*eng.kinetic()*(2.0/3.0);
            this->press_virial      = Normalize::getVolInv()*Unit::norm_press*(  eng.virial.x
                                                                               + eng.virial.y
                                                                               + eng.virial.z)/3.0;

            this->pressure = this->press_temperature
                           + this->press_virial;
        }

        std::string header(){
            std::ostringstream oss;

            oss << std::left << std::setw(12) << "  Time[fs]";
            oss << std::left << std::setw(12) << "  step";

            oss << std::left << std::setw(17) << "  Temperature[K]";
            oss << std::left << std::setw(17) << "  density[g/cm^3]";
            oss << std::left << std::setw(17) << "  pressure[Pa]";
            oss << std::left << std::setw(17) << "  press_temp[Pa]";
            oss << std::left << std::setw(17) << "  press_virial[Pa]";

            oss << "\n";

            return oss.str();
        }
        std::string to_string(const PS::F64 &dt,
                              const PS::S64 &i_step){
            std::ostringstream oss;

            oss << std::setw(12) << std::setprecision(8) << dt*PS::F64(i_step)*Unit::norm_time/Unit::femto_second;
            oss << std::setw(12) << i_step;

            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->temperature;
            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->density;
            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->pressure;
            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->press_temperature;
            oss << std::setw(17) << std::setprecision(7) << std::scientific << this->press_virial;

            oss << "\n";

            return oss.str();
        }
        void display(const PS::S64 &i_step){
            return;  // under developing
        }

        void file_init(const std::string &name,
                       const PS::S64      start = 0,
                       const PS::S32      proc  = 0){
            this->printer.file_init(name, proc);
            this->printer.print( this->header() );
            this->start = start;
        }
        void record(const PS::F64 &dt,
                    const PS::S64 &i_step) {
            if(i_step < this->start) return;
            this->printer.print( this->to_string(dt, i_step) );
        }
    };


    template <class Tprop>
    class MovingAve {
    private:
        Tprop ref;
        Tprop sum;

        PS::S64 start;
        PS::S64 cycle;
        PS::S64 count = 0;

    public:
        MovingAve<Tprop>(){
            this->ref.clear();
            this->sum.clear();
        }
        void clear(){
            this->sum.clear();
            this->count = 0;
        }
        void clearReference(){
            this->ref.clear();
        }
        void file_init(const std::string &name,
                       const PS::S64      start = 0,
                       const PS::S64      cycle = 1,
                       const PS::S32      rank  = 0 ){
            assert(start >= 0);
            assert(cycle >  0);
            assert(rank  >= 0);
            this->sum.file_init(name, 0, rank);
            this->start = start;
            this->cycle = cycle;

            this->clear();
            this->clearReference();
        }

        void setReference(const Tprop &prop){
            this->ref.assign(prop);
        }
        void record(const PS::F64 &dt,
                    const PS::S64 &i_step,
                    const Tprop   &prop   ){

            if(i_step < this->start) return;

            //--- add data
            this->sum.plus(prop, 1.0);
            ++(this->count);

            if(this->count < this->cycle) return;

            //--- output avarage value
            this->sum.multiple(1.0/this->count);
            this->sum.plus(this->ref, -1.0);
            this->sum.record(dt, i_step);

            this->sum.clear();
            this->count = 0;
        }
    };
}
