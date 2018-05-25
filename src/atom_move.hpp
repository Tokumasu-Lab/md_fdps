//***************************************************************************************
//  This is basic kick and drift routine for atom.
//***************************************************************************************
#pragma once

#include <sstream>
#include <vector>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>


enum class RESPA_MODE {
    all,
    intra,
    inter,
};

namespace ENUM {
    static const std::map<std::string, RESPA_MODE> table_str_RESPA_MODE{
        {"all"  , RESPA_MODE::all  },
        {"intra", RESPA_MODE::intra},
        {"inter", RESPA_MODE::inter},
    };
    static const std::map<RESPA_MODE, std::string> table_RESPA_MODE_str{
        {RESPA_MODE::all  , "all"  },
        {RESPA_MODE::intra, "intra"},
        {RESPA_MODE::inter, "inter"},
    };

    RESPA_MODE which_RESPA_MODE(const std::string &str){
        if(table_str_RESPA_MODE.find(str) != table_str_RESPA_MODE.end()){
            return table_str_RESPA_MODE.at(str);
        } else {
            std::cerr << "  RESPA_MODE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in RESPA_MODE.");
        }
    }
    std::string what(const RESPA_MODE &e){
        if(table_RESPA_MODE_str.find(e) != table_RESPA_MODE_str.end()){
            return table_RESPA_MODE_str.at(e);
        } else {
            using type_base = typename std::underlying_type<RESPA_MODE>::type;
            std::cerr << "  RESPA_MODE: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in RESPA_MODE.");
        }
    }
}

namespace std {
    inline string to_string(const RESPA_MODE &e){ return ENUM::what(e); }
}
inline std::ostream& operator << (std::ostream& s, const RESPA_MODE &e){
    s << ENUM::what(e);
    return s;
}


namespace ATOM_MOVE {

    namespace _Impl{

        template <class FGetForce, class Tpsys>
        void kick_atom(const PS::F64    &dt,
                             Tpsys      &psys,
                             PS::F64vec &v_barycentric){

            const PS::S64 n_local = psys.getNumberOfParticleLocal();

                    v_barycentric = 0.0;
            PS::F64 mass_total    = 0.0;

            #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                #pragma omp parallel for reduction(+: v_baricentric, mass_total)
            #endif
            for(PS::S64 i=0; i<n_local; ++i){
                const PS::F64vec acc   = FGetForce()(psys[i])*( dt/psys[i].getMass() );
                const PS::F64vec v_new = acc + psys[i].getVel();
                psys[i].setVel(v_new);

                v_barycentric += psys[i].getMass()*v_new;
                mass_total    += psys[i].getMass();
            }

            mass_total    = PS::Comm::getSum(mass_total);
            v_barycentric = PS::Comm::getSum(v_barycentric);
            v_barycentric = v_barycentric*(1.0/mass_total);
        }
    }

    struct GetForceInter {
        template <class Tptcl>
        decltype(declval<Tptcl>().getForceInter()) operator () (const Tptcl &ptcl) const {
            return ptcl.getForceInter();
        }
    };

    struct GetForceIntra {
        template <class Tptcl>
        decltype(declval<Tptcl>().getForceIntra()) operator () (const Tptcl &ptcl) const {
            return ptcl.getForceIntra();
        }
    };

    struct GetForceTotal {
        template <class Tptcl>
        decltype(declval<Tptcl>().getForce()) operator () (const Tptcl &ptcl) const {
            return ptcl.getForce();
        }
    };

    template <class Tpsys>
    void kick(const PS::F64    &dt,
                    Tpsys      &psys,
              const RESPA_MODE  respa_mode = RESPA_MODE::all){

        //--- kick atom
        const PS::S64 n_local = psys.getNumberOfParticleLocal();

        PS::F64vec v_barycentric = 0.0;

        switch(respa_mode){
            case RESPA_MODE::all:
                _Impl::kick_atom<GetForceTotal>(dt, psys, v_barycentric);
            break;

            case RESPA_MODE::intra:
                _Impl::kick_atom<GetForceIntra>(dt, psys, v_barycentric);
            break;

            case RESPA_MODE::inter:
                _Impl::kick_atom<GetForceInter>(dt, psys, v_barycentric);
            break;

            default:
                throw std::invalid_argument("undefined RESPA_MODE: " + ENUM::what(respa_mode));
        }

        //--- cancel barycentric velocity
        #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
            #pragma omp parallel for
        #endif
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
