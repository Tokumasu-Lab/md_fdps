//***************************************************************************************
//  This is generic Radial Distribution Function (RDF) analyzer.
//***************************************************************************************
#pragma once

#include <string>
#include <sstream>
#include <cassert>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>


namespace Analysis {

    namespace RDF_DEFS {
        constexpr PS::F64 default_range_min  =  0.0;
        constexpr PS::F64 default_range_max  = 10.0;
        constexpr size_t  default_resolution = 100;
    }

    template <class Target>
    class RDF_Sampler {
    private:
        FS_TOOL::FilePrinter printer;
        std::string          file_name;
        PS::S32              rank         = -1;
        bool                 setting_flag = false;
        bool                 initialized  = false;

        std::vector<Target> sight_list;
        std::vector<Target> tgt_list;

        PS::F64 range_min;
        PS::F64 range_max;
        size_t  resolution;

        size_t  array_range;

        PS::F64 dr     = 0.0;
        PS::F64 dr_inv = 1.0;

        std::vector<std::vector<size_t>> rdf_count_buff;

        size_t sight_count;
        size_t tgt_count;

        size_t sample_count;

        std::vector<PS::F64> rdf_raw;

        size_t total_sight_count;
        size_t total_tgt_count;

        void _check_setting() const {
            if(!this->setting_flag){
                std::ostringstream oss;
                oss << "the setting is invalid."
                    << " set root as the process that call 'input_setting()'." << "\n";
                throw std::logic_error(oss.str());
            }
        }

        void _check_init() const {
            if(!this->initialized){
                std::ostringstream oss;
                oss << "the RDF_Sampler is not initialized."
                    << " call 'broadcast_setting()'." << "\n";
                throw std::logic_error(oss.str());
            }
        }

        PS::S32 _get_thread_id() const {
            #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                const PS::S32 thread_id = omp_get_thread_num();
            #else
                const PS::S32 thread_id = 0;
            #endif
            return thread_id;
        }

    public:
        //--- accessor for OpenMP buffer
        std::vector<size_t>& rdf_count()  { return this->rdf_count_buff.at(  this->_get_thread_id()); }

        //--- initializer
        void input_setting(const std::string         &file,
                           const std::vector<Target> &sight_list,
                           const std::vector<Target> &tgt_list,
                           const PS::F64              range_min  = RDF_DEFS::default_range_min,
                           const PS::F64              range_max  = RDF_DEFS::default_range_max,
                           const size_t               resolution = RDF_DEFS::default_resolution,
                           const PS::S32              rank       = 0                            ){

            //--- arguments check
            if(file.size() <= 0){
                throw std::invalid_argument("output file name is invalid.");
            }
            if(sight_list.size() <= 0 ||
               tgt_list.size()   <= 0   ){
                std::ostringstream oss;
                oss << "sight_list & tgt_list must have more than 1 target." << "\n"
                    << "   sight_list: " << sight_list.size() << " targets." << "\n"
                    << "   tgt_list  : " << tgt_list.size()   << " targets." << "\n";
                throw std::invalid_argument(oss.str());
            }
            if(range_min <  0.0 ||
               range_max <= 0.0 ||
               range_min >= range_max){
                std::ostringstream oss;
                oss << "range setting must be '0 <= range_min < range_max'." << "\n"
                    << "   range_min = " << range_min << "\n"
                    << "   range_max = " << range_max << "\n";
                throw std::invalid_argument(oss.str());
            }
            if(resolution <= 1){
                std::ostringstream oss;
                oss << "resolution must be > 1." << "\n"
                    << "   resolution = " << resolution << "\n";
                throw std::invalid_argument(oss.str());
            }
            COMM_TOOL::check_proc_rank(rank);

            this->rank = rank;

            this->file_name  = file;
            this->sight_list = sight_list;
            this->tgt_list   = tgt_list;
            this->range_min  = range_min;
            this->range_max  = range_max;
            this->resolution = resolution;

            this->setting_flag = true;
            this->initialized  = false;
        }

        void broadcast_setting(const PS::S32 rank = -1){
            PS::S32 source_proc = 0;
            if(rank < 0){
                //--- auto detect
                source_proc = PS::Comm::getMaxValue(this->rank);
            } else {
                //--- use presented rank
                source_proc = rank;
            }

            COMM_TOOL::broadcast(this->setting_flag, source_proc);
            this->_check_setting();

            this->rank = source_proc;
            COMM_TOOL::broadcast(this->file_name , source_proc);
            COMM_TOOL::broadcast(this->sight_list, source_proc);
            COMM_TOOL::broadcast(this->tgt_list  , source_proc);
            COMM_TOOL::broadcast(this->range_min , source_proc);
            COMM_TOOL::broadcast(this->range_max , source_proc);
            COMM_TOOL::broadcast(this->resolution, source_proc);

            this->printer.file_init(file_name, this->rank);
            this->array_range = this->resolution + 1;

            this->dr     = (this->range_max - this->range_min)/static_cast<PS::F64>(this->resolution);
            this->dr_inv = 1.0/this->dr;

            this->rdf_count_buff.clear();

            //--- alloc buffer for OpenMP
            const PS::S32 n_thread = PS::Comm::getNumberOfThread();
            this->rdf_count_buff.resize(n_thread);
            for(auto& v : this->rdf_count_buff){
                v.resize(this->array_range, 0);  // .at(0) = total value in range < this->range_min
            }

            this->rdf_raw.clear();
            this->rdf_raw.resize(this->array_range, 0.0);

            this->sample_count = 0;

            this->sight_count = 0;
            this->tgt_count   = 0;

            this->total_sight_count = 0;
            this->total_tgt_count   = 0;

            this->initialized = true;
        }

        RDF_Sampler() = default;
        RDF_Sampler(const std::string         &file,
                    const std::vector<Target> &sight_list,
                    const std::vector<Target> &tgt_list,
                    const PS::F64              range_min  = RDF_DEFS::default_range_min,
                    const PS::F64              range_max  = RDF_DEFS::default_range_max,
                    const size_t               resolution = RDF_DEFS::default_resolution,
                    const PS::S32              rank       = 0                            ){
            this->input_setting(file,
                                sight_list,
                                tgt_list,
                                range_min,
                                range_max,
                                resolution,
                                rank       );
        }

        std::string setting() const {
            this->_check_setting();

            std::ostringstream oss;
            oss << "  RDF_Sampler:" << "\n"
                << "    output file: "     << this->printer.file_name() << "\n"
                << "    process    : at "  << this->rank << "\n"
                << "    range      : r = " << this->range_min << " ~ " << this->range_max << "\n"
                << "    resolution : n = " << this->resolution << "\n";

            oss << "    sight filter:" << "\n";
            for(const auto& s : this->sight_list){
                oss << "      " << s << "\n";
            }

            oss << "    target filter:" << "\n";
            for(const auto& t : this->tgt_list){
                oss << "      " << t << "\n";
            }

            return oss.str();
        }

        template <class Tptcl>
        bool isSight(const Tptcl &atom_i) const {
            for(const auto s : this->sight_list){
                if(s.isMatch(atom_i)) return true;
            }
            return false;
        }
        template <class Tptcl>
        bool isTarget(const Tptcl &atom_j) const {
            for(const auto t : this->tgt_list){
                if(t.isMatch(atom_j)) return true;
            }
            return false;
        }

        size_t get_sight_count() const { return this->sight_count; }
        size_t get_tgt_count()   const { return this->tgt_count;   }
        void set_sight_count(const size_t n){ this->sight_count = n; }
        void set_tgt_count  (const size_t n){ this->tgt_count   = n; }

        template <class Tptcl>
        void add_ptcl_j(const Tptcl &atom_j, const PS::F64 r){
            #ifndef NDEBUG
            this->_check_init();
            #endif

            //--- check window
            PS::S64 int_r     = 0;
            PS::S64 inclement = 1;
            if(r > this->range_max){
                int_r     = 0;
                inclement = 0;
            }

            //--- count up
            int_r = std::floor( (r - this->range_min)*this->dr_inv ) + 1;

            int_r = std::max(int_r, static_cast<PS::S64>(0));
            int_r = std::min(int_r, static_cast<PS::S64>(this->resolution));
            this->rdf_count().at(int_r) += inclement;
        }

        template <class Tptcl,
                  class TSM, class Tforce, class Tepi, class Tepj, class Tmomloc, class Tmomglb, class Tspj>
        void sampling_IA(Tptcl                      &ptcl,
                         PS::TreeForForce<TSM,
                                          Tforce,
                                          Tepi,
                                          Tepj,
                                          Tmomloc,
                                          Tmomglb,
                                          Tspj    > &tree){
            if( !this->isSight(ptcl) ) return;

            Tepj *ptr;
            const PS::S32 n_list = tree.getNeighborListOneParticle(ptcl, ptr);
            for(PS::S32 j=0; j<n_list; ++j){
                const PS::F32vec r_vec = Normalize::realPos(ptr[j].getPos() - ptcl.getPos());
                const PS::F32    r     = std::sqrt(r_vec*r_vec);

                if(ptcl.getId() == ptr[j].getId()) continue;

                if( this->isTarget(ptr[j]) ){
                    this->add_ptcl_j(ptr[j], r);
                }
            }
        }

        template <class Tpsys, class Ttree, class Tdinfo>
        void sampling(Tpsys  &psys,
                      Ttree  &tree,
                      Tdinfo &dinfo){

            const PS::S64 n_local  = psys.getNumberOfParticleLocal();

            tree.calcForceAll(IntraPair::dummy_func{}, psys, dinfo);

            #ifdef PARTICLE_SIMULATOR_THREAD_PARALLEL
                #pragma omp parallel for
            #endif
            for(PS::S64 i=0; i<n_local; ++i){
                this->sampling_IA(psys[i], tree);
            }

            //--- count up target particle
            size_t n_sight = 0;
            size_t n_tgt   = 0;
            for(PS::S64 i=0; i<n_local; ++i){
                if( this->isSight( psys[i]) ) ++(n_sight);
                if( this->isTarget(psys[i]) ) ++(n_tgt  );
            }
            this->set_sight_count(n_sight);
            this->set_tgt_count(  n_tgt  );
        }

        void convert_count_to_raw(const PS::F64 vol){
            this->_check_init();
            ++(this->sample_count);

            if(vol <= 0.0){
                std::ostringstream oss;
                oss << "vol must be > 0.0." << "\n"
                    << "   vol = " << vol << "\n";
                throw std::invalid_argument(oss.str());
            }

            //--- accumulate # of sight and target
            size_t n_sight = PS::Comm::getSum(this->get_sight_count());
            size_t n_tgt   = PS::Comm::getSum(this->get_tgt_count()  );

            if(n_sight == 0 ||
               n_tgt   == 0   ){
                //--- sight or tgt was not found.
            } else {
                //--- non-zero rdf count
                const PS::F64 n_sight_inv = 1.0/n_sight;
                const PS::F64 n_tgt_inv   = 1.0/n_tgt;
                const PS::F64 coef        = n_sight_inv*n_tgt_inv*vol;
                for(const auto& v : this->rdf_count_buff){
                    for(size_t i=0; i<=this->resolution; ++i){
                        this->rdf_raw[i] += static_cast<PS::F64>(v[i])*coef;
                    }
                }
            }

            //--- clear rdf count
            for(auto& v : this->rdf_count_buff){
                std::fill(v.begin(), v.end(), 0);
            }

            //--- clear # of sight and target
            this->total_sight_count += n_sight;
            this->total_tgt_count   += n_tgt;

            this->sight_count = 0;
            this->tgt_count   = 0;
        }

        void output(){
            this->_check_init();

            //--- collect data
            const auto sample_count_gather = COMM_TOOL::gather(this->sample_count, this->rank);
            const auto rdf_raw_gather      = COMM_TOOL::gather(this->rdf_raw     , this->rank);

            if(this->rank != PS::Comm::getRank()) return;

            //--- check data consistency
            for(const auto& c : sample_count_gather){
                if(c != this->sample_count){
                    std::ostringstream oss;
                    oss << "invalid sample count. it must be same between all process." << "\n"
                        << "   note: sample count is count up in RDF_Sampler.convert_count_to_raw()." << "\n"
                        << "         call this function in all process at same timing." << "\n";
                    size_t i_rank=0;
                    for(const auto& count : sample_count_gather){
                        oss << "      in rank = " << i_rank << ", sample_count = " << count << "\n";
                        ++i_rank;
                    }
                    throw std::logic_error(oss.str());
                }
            }
            for(const auto& v : rdf_raw_gather){
                if(v.size() != this->rdf_raw.size()){
                    std::ostringstream oss;
                    oss << "invalid resolution. it must be same between all process." << "\n";
                    size_t i_rank=0;
                    for(const auto& vec : rdf_raw_gather){
                        oss << "   in rank = " << i_rank << ", resolution = " << vec.size() << "\n";
                        ++i_rank;
                    }
                    throw std::logic_error(oss.str());
                }
            }

            //--- check data count
            if(this->sample_count <= 0){
                std::ostringstream oss;
                oss << "WARNING: RDF data is empty." << "\n"
                    << "         the file: '" << this->printer.file_name() << "' is not generated." << "\n";
                std::cout << oss.str() << std::flush;
                return;
            }

            //--- accumulate result
            const PS::F64        n_sample_inv = 1.0/static_cast<PS::F64>(this->sample_count);
            std::vector<PS::F64> rdf_final;
            rdf_final.resize(this->array_range, 0.0);
            for(const auto& vec : rdf_raw_gather){
                for(size_t i=0; i<=this->resolution; ++i){
                    rdf_final[i] += vec[i]*n_sample_inv;
                }
            }

            //--- output result
            std::ostringstream oss;
            //------ header
            oss.str("");
            oss << "r" << "\t" << "RDF" << "\t" << "coordinate" << "\n";
            this->printer.print(oss.str());

            //------ body
            PS::F64 rdf        = 0.0;
            PS::F64 coordinate = rdf_final[0];
            for(size_t int_r=1; int_r<=this->resolution; ++int_r){
                PS::F64 r         = this->range_min + this->dr*static_cast<PS::F64>(int_r);
                PS::F64 vol_shell = 4.0*Unit::pi*r*r*this->dr;

                if(vol_shell != 0.0){
                    rdf = rdf_final[int_r]/vol_shell;
                } else {
                    rdf = 0.0;
                }
                coordinate += rdf_final[int_r];

                oss.str("");
                oss << std::setprecision(6) << r          << "\t"
                    << std::setprecision(6) << rdf        << "\t"
                    << std::setprecision(6) << coordinate << "\n";
                this->printer.print(oss.str());
            }

            //--- finish
            oss.str("");
            oss << "   RDF_Sampler: the result file '" << this->printer.file_name() << "' was generated." << "\n"
                << "      average target count: " << static_cast<PS::F64>(this->total_tgt_count  )/this->sample_count << "\n"
                << "      average sight  count: " << static_cast<PS::F64>(this->total_sight_count)/this->sample_count << "\n";
            std::cout << oss.str() << std::flush;
        }

    };

}
