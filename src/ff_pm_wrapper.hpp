/**************************************************************************************************/
/**
* @file  ff_pm_wrapper.hpp
* @brief wrapper for PS::ParticleMesh. To work around the limitation that a particle must be in local domain.
*/
/**************************************************************************************************/
#pragma once

#include <particle_simulator.hpp>
#include <particle_mesh.hpp>
#include <molecular_dynamics_ext.hpp>


namespace FORCE {
    namespace PM {

    /*
    * @brief mesh size calculator.
    */
    PS::S32 recommendMeshSize(const PS::S64 &n_total){
        assert(n_total > 0);
        const PS::S32 min_mesh_size = 16; // default of FDPS

        const PS::F32 ideal_mesh_size = std::pow(n_total, 1.0/3.0)*0.5;  // N^(1/3)/2 is recommended by FDPS spec.
        const PS::S32 p               = static_cast<PS::S32>( std::ceil( std::log2(ideal_mesh_size) ) );

        const PS::S32 recommend_mesh_size = std::pow(2, p);
        return std::max(min_mesh_size, recommend_mesh_size);
    }

    /*
    * @brief mesh size checker.
    * @details check "SIZE_OF_MESH" value is appropriate or not. when it not good, sugest recommended value.
    */
    void checkMeshSize(const PS::S64 &n_total){
        const PS::S32 recommend_mesh_size = recommendMeshSize(n_total);

        if(SIZE_OF_MESH   < recommend_mesh_size ||
           SIZE_OF_MESH/2 > recommend_mesh_size){
            if(PS::Comm::getRank() == 0){
                std::ostringstream oss;
                oss << "\n"
                    << "WARNING: the value of 'SIZE_OF_MESH' is not appropriate." << "\n"
                    << "    SIZE_OF_MESH = " << SIZE_OF_MESH << "\n"
                    << "    recommended  = " << recommend_mesh_size << "\n"
                    << "  edit value in the file of '(FDPS_DIR)/src/particle_mesh/param_fdps.h'." << "\n"
                    << "\n";
                std::cerr << oss.str() << std::flush;
            }
        }
    }

    /*
    * @brief temporary data class for CalcForceParticleMesh (internal use)
    */
    class EP_ParticleMesh {
    private:
        PS::F32vec pos;
        PS::F32    charge;

    public:
        inline PS::F32vec getPos()                const { return this->pos;    }
        inline PS::F32    getChargeParticleMesh() const { return this->charge; }

        inline void setPos(const PS::F32vec &pos_new) { this->pos = pos_new; }

        template <class Tptcl>
        void copyFromFP(const Tptcl &ptcl){
            this->pos    = ptcl.getPos();
            this->charge = ptcl.getChargeParticleMesh();
        }
    };

    /*
    * @brief temporary data class for CalcForceParticleMesh (internal use)
    */
    class Result_ParticleMesh {
    private:
        PS::S32    proc_id = -1;
        PS::S32    index   = -1;

        PS::F32vec pos;

        PS::F32    pot;
        PS::F32vec field;

    public:
        void clear(){
            this->pot   = 0.0;
            this->field = PS::F32vec{0.0, 0.0, 0.0};
        }

        template <class Tptcl>
        void copyFromFP(const Tptcl &fp){
            this->proc_id = PS::Comm::getRank();
            this->pos     = fp.getPos();
            this->clear();
        }

        inline void setIndex(const PS::S64 &index){ this->index = index; }

        inline PS::S32    getProc()  const { return this->proc_id; }
        inline PS::S32    getIndex() const { return this->index;   }

        inline PS::F32vec getPos()  const { return this->pos;     }

        inline PS::F32    getPot()   const { return this->pot;   }
        inline PS::F32vec getField() const { return this->field; }

        inline void setPos(  const PS::F32vec &pos)  { this->pos    = pos;   }

        inline void addPotParticleMesh(  const PS::F32    &pot)  { this->pot   += pot;   }
        inline void addFieldParticleMesh(const PS::F32vec &field){ this->field += field; }
    };


    /*
    * @brief wrapper for PS::ParticleMesh.
    * @details workaround the limitation which the position of point charge must be inside the local domain.
    * @details [tradeoff] = increasing CPU & MPI load for convert & communicate intermediate particles.
    */
    class CalcForceParticleMesh {
    public:
        //--- data type
        using EP_type     = EP_ParticleMesh;
        using Result_type = Result_ParticleMesh;

    private:
        //--- FDPS object
        PS::PM::ParticleMesh pm;

        //--- buffer object for PM
        PS::DomainInfo                  pm_dinfo;
        PS::ParticleSystem<EP_type>     pm_ep_buff;
        PS::ParticleSystem<Result_type> pm_result_buff;

        //--- buffer for MPI communication
        std::vector<std::vector<Result_type>> result_send_buff;
        std::vector<std::vector<Result_type>> result_recv_buff;

        //--- buffer for serialize
        std::vector<Result_type> result_seq;

        bool dinfo_flag = true;

        template <class Tptcl, class Tep>
        void _impl_copy_into_buffer(PS::ParticleSystem<Tptcl> &psys_fp,
                                    PS::ParticleSystem<Tep>   &psys_ep){

            const PS::S64 n_local = psys_fp.getNumberOfParticleLocal();
            psys_ep.setNumberOfParticleLocal(n_local);

            #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                #pragma omp parallel for
            #endif
            for(PS::S64 i=0; i<n_local; ++i){
                psys_ep[i].copyFromFP(psys_fp[i]);
            }
        }

        void _impl_check_n_local(const PS::S64 &n) const {
            if(n > std::numeric_limits<PS::S32>::max()){
                std::ostringstream oss;
                oss << "n_local in process is too large." << "\n"
                    << "    n_local = " << n << " > max value of PS::S32." << "\n";
                throw std::length_error(oss.str());
            }
        }

    public:
        CalcForceParticleMesh(){
            this->pm_ep_buff.initialize();
            this->pm_result_buff.initialize();

            this->pm_dinfo.initialize( 1.0f );    // exchange particle by newest position only
            this->pm_dinfo.setBoundaryCondition(PS::BOUNDARY_CONDITION_PERIODIC_XYZ);
            this->pm_dinfo.setPosRootDomain( PS::F32vec{0.0, 0.0, 0.0},
                                             PS::F32vec{1.0, 1.0, 1.0} );  // fixed size for using PS::ParticleMesh

            this->result_send_buff.resize(PS::Comm::getNumberOfProc());
            this->result_recv_buff.resize(PS::Comm::getNumberOfProc());
        }
        CalcForceParticleMesh(const CalcForceParticleMesh&) = delete;
        CalcForceParticleMesh& operator = (const CalcForceParticleMesh&) = delete;
        ~CalcForceParticleMesh() = default;

        /*
        * @breif dummy function. (for compativility with PS::PM::ParticleMesh)
        * @details this class use internal domaininfo.
        */
        template <class Tdinfo>
        void setDomainInfoParticleMesh(Tdinfo &dinfo){ return; }

        template <class Tptcl>
        void setParticleParticleMesh(PS::ParticleSystem<Tptcl> &psys,
                                     const bool                 clear_flag = true){

            //--- copy into buffer
            this->_impl_copy_into_buffer(psys, this->pm_ep_buff);

            //--- update domaininfo
            if(clear_flag || this->dinfo_flag){
                this->pm_dinfo.decomposeDomainAll(this->pm_ep_buff);
                this->pm.setDomainInfoParticleMesh(this->pm_dinfo);

                this->dinfo_flag = false;
            }

            //--- update buffer
            this->pm_ep_buff.adjustPositionIntoRootDomain(this->pm_dinfo);
            this->pm_ep_buff.exchangeParticle(this->pm_dinfo);

            //--- set charge into PM
            this->pm.setParticleParticleMesh(this->pm_ep_buff, clear_flag);
        }

        void calcMeshForceOnly(){
            this->pm.calcMeshForceOnly();
        }

        template <class Tptcl>
        void writeBackForce(PS::ParticleSystem<Tptcl> &psys){
            const PS::S64 n_local = psys.getNumberOfParticleLocal();

            this->_impl_check_n_local(n_local);

            //--- copy into buffer
            this->pm_result_buff.setNumberOfParticleLocal(n_local);

            #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                #pragma omp parallel for
            #endif
            for(PS::S32 i=0; i<n_local; ++i){
                this->pm_result_buff[i].copyFromFP(psys[i]);
                this->pm_result_buff[i].setIndex(i);
            }

            //--- update buffer
            this->pm_result_buff.adjustPositionIntoRootDomain(this->pm_dinfo);
            this->pm_result_buff.exchangeParticle(this->pm_dinfo);

            //--- get fleld from PM
            const PS::S64 n_local_buff = this->pm_result_buff.getNumberOfParticleLocal();
            this->_impl_check_n_local(n_local_buff);

            #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                #pragma omp parallel for
            #endif
            for(PS::S32 i=0; i<n_local_buff; ++i){
                const auto& pos_i = this->pm_result_buff[i].getPos();
                this->pm_result_buff[i].addPotParticleMesh(   Normalize::realPMPotential( -this->pm.getPotential( pos_i ) ) );
                this->pm_result_buff[i].addFieldParticleMesh( Normalize::realPMForce(     -this->pm.getForce(     pos_i ) ) );
            }

            //--- send result
            for(auto& vec : this->result_send_buff){ vec.clear(); }
            for(PS::S32 i=0; i<n_local_buff; ++i){
                const auto& i_proc = this->pm_result_buff[i].getProc();
                this->result_send_buff[i_proc].push_back( this->pm_result_buff[i] );
            }
            COMM_TOOL::allToAll(this->result_send_buff, this->result_recv_buff);

            //--- check sum
            PS::S32 recv_total = 0;
            for(const auto& v : this->result_recv_buff){
                recv_total += v.size();
            }
            if(n_local != recv_total){
                std::ostringstream oss;
                oss << "error in virtual particle management." << "\n"
                    << "  n_local = " << n_local    << "\n"
                    << "  n_recv  = " << recv_total << ", must be same." << "\n";
                throw std::logic_error(oss.str());
            }

            #ifdef CHECK_FORCE_STRENGTH
                std::vector< std::tuple<decltype(declval<Tptcl>().getId()),
                                        PS::F32vec,
                                        PS::F32vec,
                                        PS::F32vec                         > > report_force_err;
            #endif

            //--- writeback result (with serialize)
            this->result_seq.clear();
            this->result_seq.reserve(n_local);
            for(const auto& recv_vec : this->result_recv_buff){
                for(const auto& result : recv_vec){
                    this->result_seq.emplace_back(result);
                }
            }
            #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                #pragma omp parallel for
            #endif
            for(const auto& result : this->result_seq){
                const auto& index = result.getIndex();

                psys[index].addPotParticleMesh(   result.getPot()   );
                psys[index].addFieldParticleMesh( result.getField() );

                #ifdef CHECK_FORCE_STRENGTH
                    const PS::F64 dt_val = 0.00295207;  // normalized dt at dt = 0.5 [fs]

                    const auto field_pm    = result.getField();
                    const auto force_pm    = field_pm*psys[index].getCharge();
                    const auto acc_pm      = force_pm*(dt_val/psys[index].getMass());
                    const auto vec_pm_norm = Normalize::normDrift(acc_pm*dt_val);
                    const auto move_r_norm = std::sqrt(vec_pm_norm*vec_pm_norm);

                    if(move_r_norm > 0.5){
                        report_force_err.push_back( std::make_tuple(psys[index].getId(),
                                                                    psys[index].getPos(),
                                                                    field_pm,
                                                                    force_pm)             );
                    }
                #endif
            }

            #ifdef CHECK_FORCE_STRENGTH
                for(PS::S32 i_proc=0; i_proc<PS::Comm::getNumberOfProc(); ++i_proc){
                    if(PS::Comm::getRank() == i_proc &&
                       report_force_err.size() >  0    ){
                        std::ostringstream oss;
                        oss << "  proc = " << i_proc << ", PM field is too strong." << "\n";
                        for(const auto& rep : report_force_err){
                            oss << "    id = "     << std::get<0>(rep)
                                << ", pos = "      << std::get<1>(rep)
                                << ", field_pm = " << std::get<2>(rep)
                                << ", force_pm = " << std::get<3>(rep) << "\n";
                        }
                        std::cerr << oss.str() << std::flush;
                    }
                    COMM_TOOL::barrier();
                }
            #endif


            //--- writeback result (without serialize & OpenMP)
            /*
            for(const auto& recv_vec : this->result_recv_buff){
                for(const auto& result : recv_vec){
                    const auto& index = result.getIndex();

                    psys[index].addPotCoulomb(   result.getPot()   );
                    psys[index].addFieldCoulomb( result.getField() );
                }
            }
            */
        }
    };


    /*
    * @brief wrapper for PS::ParticleMesh.
    * @details provide same interface with CalcForceParticleMesh class.
    * @details this implementation is raw PS::ParticleMesh. the 'coef_ema' must be 1.0f at PS::DomainInfo::initialize().
    */
    class ParticleMesh {
    private:
        PS::PM::ParticleMesh pm;

    public:
        template <class Tdinfo>
        void setDomainInfoParticleMesh(Tdinfo &dinfo){
            this->pm.setDomainInfoParticleMesh(dinfo);
        }

        template <class Tptcl>
        void setParticleParticleMesh(PS::ParticleSystem<Tptcl> &psys,
                                     const bool                 clear_flag = true){
            this->pm.setParticleParticleMesh(psys, clear_flag);
        }

        void calcMeshForceOnly(){
            this->pm.calcMeshForceOnly();
        }

        template <class Tptcl>
        void writeBackForce(PS::ParticleSystem<Tptcl> &psys){
            const PS::S64 n_local = psys.getNumberOfParticleLocal();

            #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                #pragma omp parallel for
            #endif
            for(PS::S64 i=0; i<n_local; ++i){
                const auto& pos_i = psys[i].getPos();
                psys[i].addFieldParticleMesh( Normalize::realPMForce(     -this->pm.getForce(     pos_i ) ) );
                psys[i].addPotParticleMesh(   Normalize::realPMPotential( -this->pm.getPotential( pos_i ) ) );
            }
        }
    };


    }
}
