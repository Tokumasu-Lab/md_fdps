//***************************************************************************************
//  This is basic kick and drift routine for atom.
//***************************************************************************************
#pragma once

#include <sstream>
#include <vector>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

namespace ATOM_MOVE {

    template <class Tpsys>
    void kick(const PS::F64 &dt,
                    Tpsys   &psys){

        //--- kick atom
        const PS::S64 n_local = psys.getNumberOfParticleLocal();

        PS::F64vec v_barycentric = 0.0;
        PS::F64    mass_total    = 0.0;
        for(PS::S64 i=0; i<n_local; ++i){
            PS::F64vec v_new = psys[i].getVel()
                             + psys[i].getForce()*( dt/psys[i].getMass() );
            psys[i].setVel(v_new);

            v_barycentric += psys[i].getMass()*v_new;
            mass_total    += psys[i].getMass();
        }

        //--- cancel barycentric velocity
        mass_total    = PS::Comm::getSum(mass_total);
        v_barycentric = PS::Comm::getSum(v_barycentric);
        v_barycentric = v_barycentric*(1.0/mass_total);
        for(PS::S64 i=0; i<n_local; ++i){
            psys[i].setVel( psys[i].getVel() - v_barycentric );
        }
    }

    template <class Tptcl>
    void check_move_runaway(const PS::F64                   &dt,
                                  PS::ParticleSystem<Tptcl> &psys){
        const PS::S64 n_local    = psys.getNumberOfParticleLocal();
        const PS::F64 move_limit = 0.5;

        std::vector<Tptcl> report_vec;
        report_vec.reserve(10);
        report_vec.clear();

        PS::F64 max_move = 0.0;
        for(PS::S64 i=0; i<n_local; ++i){
            const PS::F64vec move   = Normalize::normDrift(psys[i].getVel()*dt);
            const PS::F64    move_r = std::sqrt(move*move);
            max_move = std::max(max_move, move_r);

            if(move_r > move_limit){
                report_vec.push_back(psys[i]);
            }
        }

        //--- show dettailed report
        for(PS::S32 i_proc=0; i_proc<PS::Comm::getNumberOfProc(); ++i_proc){
            if(PS::Comm::getRank() == i_proc    &&
               max_move            >  move_limit  ){
                std::ostringstream oss;
                oss << "  proc = " << i_proc << ", max_move = " << max_move << "\n";
                for(const auto& atom : report_vec){
                    oss << "    id = "  << atom.getId()
                        << ", pos = "   << atom.getPos()
                        << ", force = " << atom.getForce() << "\n";
                }
                std::cerr << oss.str() << std::flush;
            }
            COMM_TOOL::barrier();
        }

        //--- stopper
        if(max_move > move_limit){
            throw std::logic_error("atoms speed runaway");
        }
    }

    //--- after calling drift(), must call psys.adjustPositionIntoRootDomain(dinfo);
    template <class Tpsys>
    PS::F64 drift(const PS::F64 &dt,
                        Tpsys   &psys){

        const PS::S64 n_local    = psys.getNumberOfParticleLocal();
        const PS::F64 move_limit = 0.5;

        #ifndef NDEBUG
            //--- detailed error check
            check_move_runaway(dt, psys);
        #endif

        PS::F64 max_move = 0.0;

        for(PS::S64 i=0; i<n_local; ++i){
            PS::F64vec move    = psys[i].getVel()*dt;
            PS::F64vec pos_new = psys[i].getPos() + Normalize::normDrift(move);
            psys[i].addTrj( move );
            psys[i].setPos( pos_new );

            //--- check largest move at step
            move     = Normalize::normDrift(move);
            max_move = std::max(max_move, move*move);
        }
        max_move = std::sqrt(max_move);

        //--- simple error check
        if(max_move > move_limit){
            std::ostringstream oss;
            oss << "  proc = " << PS::Comm::getRank() << ", max_move = " << max_move << "\n";
            std::cerr << oss.str() << std::flush;
            throw std::logic_error("atoms speed runaway");
        }

        return max_move;
    }

}
