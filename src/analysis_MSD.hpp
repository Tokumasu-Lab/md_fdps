//***************************************************************************************
//  This is generic Mean Square Displacement (MSD) analyzer.
//***************************************************************************************
#pragma once

#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <tuple>
#include <unordered_map>
#include <stdexcept>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

namespace Analysis {
    enum class MSD_SAMPLING_MODE {
        raw,
        resampling,
    };
}

namespace ENUM {

    //====================================
    //  enum interface for Analysis::MSD_SAMPLING_MODE
    //====================================
    static const std::unordered_map<std::string, Analysis::MSD_SAMPLING_MODE> table_str_msd_sampling_mode{
        {"raw"       , Analysis::MSD_SAMPLING_MODE::raw       },
        {"resampling", Analysis::MSD_SAMPLING_MODE::resampling},
    };

    static const std::unordered_map<Analysis::MSD_SAMPLING_MODE, std::string> table_msd_sampling_mode_str{
        {Analysis::MSD_SAMPLING_MODE::raw       , "raw"        },
        {Analysis::MSD_SAMPLING_MODE::resampling, "resampling"},
    };

    Analysis::MSD_SAMPLING_MODE which_MSD_SAMPLING_MODE(const std::string &str){
        if(table_str_msd_sampling_mode.find(str) != table_str_msd_sampling_mode.end()){
            return table_str_msd_sampling_mode.at(str);
        } else {
            std::cerr << "  Analysis::MSD_SAMPLING_MODE: input = " << str << std::endl;
            throw std::out_of_range("undefined enum value in Analysis::MSD_SAMPLING_MODE.");
        }
    }

    std::string what(const Analysis::MSD_SAMPLING_MODE &e){
        if(table_msd_sampling_mode_str.find(e) != table_msd_sampling_mode_str.end()){
            return table_msd_sampling_mode_str.at(e);
        } else {
            using type_base = typename std::underlying_type<Analysis::MSD_SAMPLING_MODE>::type;
            std::cerr << "  Analysis::MSD_SAMPLING_MODE: input = " << static_cast<type_base>(e) << std::endl;
            throw std::out_of_range("undefined enum value in Analysis::MSD_SAMPLING_MODE.");
        }
    }
}

//--- specialize for std::to_string() & std::ostream
namespace std {
    inline string to_string(const Analysis::MSD_SAMPLING_MODE &e){ return ENUM::what(e); }

    inline ostream& operator << (ostream& s, const Analysis::MSD_SAMPLING_MODE &e){
        s << ENUM::what(e);
        return s;
    }
}

namespace Analysis {

    template <class Tid, class Tmove = PS::F32vec>
    class MSD_Tracer {
    public:
        struct ID_range {
            Tid begin;
            Tid end;

            bool isInRange(const Tid &id) const {
                return (this->begin <= id && id < this->end);
            }
        };

    private:
        Tid                   n_piece;
        std::vector<ID_range> id_range_table;

        bool initialized = false;

        std::vector<std::vector<std::pair<Tid, Tmove>>> send_vec;
        std::vector<std::vector<std::pair<Tid, Tmove>>> recv_vec;

        std::unordered_map<Tid, std::vector<Tmove>> move_trace;
        std::vector<PS::F64>                        time_history;


        void _check_init() const {
            if(!this->initialized){
                throw std::logic_error("the MSD_Tracer is not initialized. call 'init()' at first.");
            }
        }

    public:
        void init(const Tid &n_total){
            const auto n_list = COMM_TOOL::allGather(n_total);
            for(const auto& n : n_list){
                if(n != n_total){
                    std::ostringstream oss;
                    oss << "n_total must be same in all MPI process." << "\n";
                    for(size_t i=0; i<n_list.size(); ++i){
                        oss << "  Proc = " << i << ", n_total = " << n_list[i] << "\n";
                    }
                    throw std::invalid_argument(oss.str());
                }
            }

            const PS::S32 n_proc = PS::Comm::getNumberOfProc();
            this->n_piece = (n_total/n_proc) + 1;

            this->id_range_table.resize(n_proc);
            for(PS::S32 i=0; i<n_proc; ++i){
                this->id_range_table[i].begin = static_cast<Tid>(n_piece*i);
                this->id_range_table[i].end   = std::min(n_piece*(i+1), n_total);
            }

            this->send_vec.resize(n_proc);
            this->recv_vec.resize(n_proc);

            this->move_trace.clear();
            this->move_trace.max_load_factor(0.7);
            this->move_trace.reserve(n_piece);

            this->time_history.clear();

            this->initialized = true;
        }

        //--- O(N) function. safe & slow.
        PS::S32 which_proc_naive(const Tid &id) const {
            for(PS::S32 i_proc=0; i_proc<static_cast<PS::S32>(this->id_range_table.size()); ++i_proc){
                if( this->id_range_table[i_proc].isInRange(id) ){
                    return i_proc;
                }
            }
            //--- illigal value
            if(id < 0 || id >= this->id_range_table.back().end){
                std::ostringstream oss;
                oss << "id must be in range: 0 <= id < n_total." << "\n"
                    << "   id      = " << id << "\n"
                    << "   n_total = " <<  this->id_range_table.back().end - 1 << "\n"
                    << " note: it is recommended the id (return value of ptcl.getId()) is a sequencial positive integer value." << "\n";
                throw std::out_of_range(oss.str());
            }
            return -1;
        }

        //--- O(1) function. be careful for the data structure of this->id_range_table
        PS::S32 which_proc(const Tid &id) const {
            //--- predict
            PS::S32 i_proc = static_cast<PS::S32>(id/this->n_piece);
            i_proc = std::max(i_proc, 0);
            i_proc = std::min(i_proc, PS::Comm::getNumberOfProc()-1);
            //--- correct search
            if(this->id_range_table[i_proc].isInRange(id)) return i_proc;

            if(this->id_range_table[i_proc].end <= id){
                //--- forword iterate
                for(PS::S32 i=i_proc+1; i<static_cast<PS::S32>(this->id_range_table.size()); ++i){
                    if(this->id_range_table[i].isInRange(id)) return i;
                }
            } else {
                //--- reverse iterate
                for(PS::S32 i=i_proc-1; i>=0; --i){
                    if(this->id_range_table[i].isInRange(id)) return i;
                }
            }

            //--- illigal value
            if(id < 0 || id >= this->id_range_table.back().end){
                std::ostringstream oss;
                oss << "id must be in range: 0 <= id < n_total." << "\n"
                    << "   id      = " << id << "\n"
                    << "   n_total = " <<  this->id_range_table.back().end - 1 << "\n"
                    << " note: it is recommended the id (return value of ptcl.getId()) is a sequencial positive integer value." << "\n";
                throw std::out_of_range(oss.str());
            }
            return -1;
        }

        template <class Tprof, class Tpsys>
        void add_trace(const Tprof &profile,
                             Tpsys &psys    ){

            this->_check_init();

            //--- clear communicate buffers
            for(auto& v : this->send_vec){
                v.clear();
            }

            //--- input move into local buffer
            const PS::S64 n_local = psys.getNumberOfParticleLocal();
            for(PS::S64 i=0; i<n_local; ++i){
                const PS::S32 i_proc = this->which_proc(psys[i].getId());
                this->send_vec[i_proc].push_back( std::make_pair(psys[i].getId(),
                                                                 psys[i].getTrj()) );
            }

            //--- collect move for trace
            COMM_TOOL::allToAll(this->send_vec, this->recv_vec);
            for(const auto& v : this->recv_vec){
                for(const auto& atom : v){
                    auto& trace = this->move_trace[atom.first];
                    if( trace.empty() ){
                        //--- first step
                        trace.push_back(atom.second);
                    } else {
                        //--- after second step (calc total move)
                        Tmove move_total = trace.back() + atom.second;
                        trace.push_back(move_total);
                    }
                }
            }

            //--- time record
            PS::F64 time_now = profile.get_trj_time();
            if( !this->time_history.empty() ){
                time_now += this->time_history.back();
            }
            this->time_history.push_back(time_now);

            //--- check sum
            const size_t n_ref = this->time_history.size();
            for(const auto& trace_IA : this->move_trace){
                if(trace_IA.second.size() != n_ref){
                    std::ostringstream oss;
                    oss << "missing the trajectory data." << "\n"
                        << "   n_sample  = " << n_ref          << "\n"
                        << "   id        = " << trace_IA.first << "\n";
                    throw std::invalid_argument(oss.str());
                }
            }
        }

        //--- interface for sampler
        const std::vector<PS::F64>&                        get_time_history() const { return this->time_history; }
        const std::unordered_map<Tid, std::vector<Tmove>>& get_move_trace()   const { return this->move_trace;   }

        bool isRegularInterval(const PS::F64 eps = 1.e-6) const {
            if(this->time_history.empty()) return false;

            const auto dt = this->time_history.front();
            if(dt == 0.0) return false;

            const auto dt_inv = 1.0/dt;
            for(size_t i=0; i<this->time_history.size()-1; ++i){
                const PS::F64 dt_i = this->time_history[i+1] - this->time_history[i];
                if( (dt_i - dt)*dt_inv > eps ) return false;
            }

            return true;
        }
    };

    template <class Target, class Tid = PS::S64>
    class MSD_Sampler {
    private:
        std::string         file_name;
        std::vector<Target> tgt_list;
        PS::S32             rank = -1;

        bool setting_flag = false;
        bool sync_flag    = false;
        bool target_flag  = false;

        size_t           n_tgt_total;
        std::vector<Tid> id_list;

        MSD_SAMPLING_MODE sampling_mode = MSD_SAMPLING_MODE::raw;

        struct R2_DATA {
            PS::F64 r2_3d = 0.0;
            PS::F64 r2_x  = 0.0;
            PS::F64 r2_y  = 0.0;
            PS::F64 r2_z  = 0.0;

            void clear(){
                this->r2_3d = 0.0;
                this->r2_x  = 0.0;
                this->r2_y  = 0.0;
                this->r2_z  = 0.0;
            }
            template <class Tf>
            void add_sample(const PS::Vector3<Tf> &vec){
                const Tf x2 = vec.x*vec.x;
                const Tf y2 = vec.y*vec.y;
                const Tf z2 = vec.z*vec.z;
                this->r2_3d += (x2 + y2 + z2);
                this->r2_x  += x2;
                this->r2_y  += y2;
                this->r2_z  += z2;
            }
            void apply_coef(const PS::F64 &r){
                this->r2_3d = r*this->r2_3d;
                this->r2_x  = r*this->r2_x;
                this->r2_y  = r*this->r2_y;
                this->r2_z  = r*this->r2_z;
            }

            R2_DATA& operator += (const R2_DATA &rhs){
                this->r2_3d += rhs.r2_3d;
                this->r2_x  += rhs.r2_x;
                this->r2_y  += rhs.r2_y;
                this->r2_z  += rhs.r2_z;

                return *this;
            }
        };

        template <class Tptcl>
        bool isTarget(const Tptcl               &ptcl,
                      const std::vector<Target> &tgt_list){

            for(const auto& tgt : tgt_list){
                if( tgt.isMatch(ptcl) ) return true;
            }
            return false;
        }

    public:
        void clear(){
            this->setting_flag = false;
            this->sync_flag    = false;
            this->target_flag  = false;

            this->file_name.clear();
            this->tgt_list.clear();
        }
        void input_setting(const std::string          &file_name,
                           const std::vector<Target>  &tgt_list,
                           const MSD_SAMPLING_MODE     sampling_mode = MSD_SAMPLING_MODE::resampling,
                           const PS::S32               rank = 0 ){

            //--- sequence check
            if(this->setting_flag){
                std::ostringstream oss;
                oss << " setting is already input." << "\n"
                    << "    file_name    : " << this->file_name << "\n"
                    << "    target filter: " << this->tgt_list.size() << " factor." << "\n"
                    << "    sampling mode: " << this->sampling_mode << "\n"
                    << "    rank         : " << this->rank << "\n";
                throw std::logic_error(oss.str());
            }

            //--- arguments check
            if(file_name.size() <= 0){
                throw std::invalid_argument("output file name is invalid.");
            }
            if(tgt_list.size() <= 0){
                std::ostringstream oss;
                oss << "tgt_list must have more than 1 target." << "\n";
                throw std::invalid_argument(oss.str());
            }
            COMM_TOOL::check_proc_rank(rank);

            this->file_name     = file_name;
            this->tgt_list      = tgt_list;
            this->sampling_mode = sampling_mode;
            this->rank          = rank;

            this->setting_flag = true;
        }
        MSD_Sampler() = default;
        MSD_Sampler(const std::string          &file_name,
                    const std::vector<Target>  &tgt_list,
                    const MSD_SAMPLING_MODE     sampling_mode = MSD_SAMPLING_MODE::resampling,
                    const PS::S32               rank = 0){
            this->input_setting(file_name,
                                tgt_list,
                                sampling_mode,
                                rank          );
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

            this->rank = source_proc;
            COMM_TOOL::broadcast(this->file_name    , source_proc);
            COMM_TOOL::broadcast(this->tgt_list     , source_proc);
            COMM_TOOL::broadcast(this->sampling_mode, source_proc);

            if(this->file_name.size() < 1){
                throw std::logic_error(" output file name is empty.");
            }
            if(this->tgt_list.size() < 1){
                throw std::logic_error(" target list is empty.");
            }

            this->sync_flag = true;
        }

        template <class Tpsys, class Tracer>
        void make_id_table(      Tpsys  &psys,
                           const Tracer &tracer){

            const PS::S32 n_proc  = PS::Comm::getNumberOfProc();
            const PS::S64 n_local = psys.getNumberOfParticleLocal();

            std::vector<std::vector<Tid>> send_id_list;

            send_id_list.resize(n_proc);

            for(auto& v : send_id_list){
                v.clear();
            }

            for(PS::S64 i=0; i<n_local; ++i){
                if( !this->isTarget(psys[i], this->tgt_list) ) continue;

                const Tid     id_i   = psys[i].getId();
                const PS::S32 i_proc = tracer.which_proc(id_i);
                send_id_list[i_proc].push_back(id_i);
            }

            const auto recv_id_list = COMM_TOOL::allToAll(send_id_list);

            size_t n_tgt_recv = 0;
            this->id_list.clear();
            for(const auto& recv_v : recv_id_list){
                n_tgt_recv += recv_v.size();
                for(const auto& tgt_id : recv_v){
                    this->id_list.push_back(tgt_id);
                }
            }

            this->n_tgt_total = PS::Comm::getSum(n_tgt_recv);
            this->target_flag = true;
        }

    private:
        void _impl_check_data_effectivity(const PS::S32 n_sample) const {
            if( !this->sync_flag ){
                throw std::logic_error(" the MSD_Sampler must call 'broadcast_setting()' before calling 'output()'.");
            }
            if( !this->target_flag ){
                throw std::logic_error(" the MSD_Sampler must call 'make_id_table()' before calling 'output()'.");
            }
            if(n_sample < 2){
                throw std::logic_error(" sample data is not enough to calculate MSD (more than 2 samples required).");
            }
        }

        template <class Tdata_trace>
        void _impl_calc_msd_local_resampling(const Tdata_trace          &move_trace,
                                             const PS::S32               n_sample,
                                                   std::vector<R2_DATA> &data_local) const {
            data_local.clear();

            R2_DATA tmp;
            tmp.clear();

            for(PS::S32 i_span=1; i_span<n_sample; ++i_span){
                tmp.clear();

                for(const auto& id_atom : this->id_list){
                    const auto& trace_atom = move_trace.at(id_atom);

                    //--- interpolation
                    for(PS::S32 i_begin=0; i_begin<n_sample-i_span; ++i_begin){
                        const PS::F64vec r_vec = trace_atom.at(i_begin + i_span) - trace_atom.at(i_begin);
                        tmp.add_sample(r_vec);
                    }
                }
                const PS::F64 sample_coef = 1.0/static_cast<PS::F64>(this->n_tgt_total*(n_sample-i_span));
                tmp.apply_coef(sample_coef);

                data_local.push_back(tmp);
            }
        }

        template <class Tdata_trace>
        void _impl_calc_msd_local(const Tdata_trace          &move_trace,
                                  const PS::S32               n_sample,
                                        std::vector<R2_DATA> &data_local) const {
            data_local.clear();

            R2_DATA tmp;
            tmp.clear();

            for(PS::S32 i_span=1; i_span<n_sample; ++i_span){
                tmp.clear();

                for(const auto& id_atom : this->id_list){
                    const auto&      trace_atom = move_trace.at(id_atom);
                    const PS::F64vec r_vec      = trace_atom.at(i_span) - trace_atom.front();
                    tmp.add_sample(r_vec);
                }
                const PS::F64 sample_coef = 1.0/static_cast<PS::F64>(this->n_tgt_total);
                tmp.apply_coef(sample_coef);

                data_local.push_back(tmp);
            }
        }

        template <class Thistory>
        void _impl_write_result(const Thistory             &time_history,
                                const std::vector<R2_DATA> &data_local   ) const {

            //--- accumulate data
            const auto data_all = COMM_TOOL::gather(data_local, this->rank);

            if(PS::Comm::getRank() != this->rank) return;

            //------ alloc & clear buffer
            std::vector<R2_DATA> data_final;
            data_final.resize( data_all[0].size() );
            for(auto& value : data_final){
                value.clear();
            }
            //------ accumulate
            for(const auto& v_recv : data_all){
                for(size_t i=0; i<v_recv.size(); ++i){
                    data_final.at(i) += v_recv.at(i);
                }
            }

            //--- output file
            FS_TOOL::FilePrinter printer{this->file_name, this->rank};

            //------ header
            std::ostringstream oss;
            oss << "Time[fs]" << "\t"
                << "MSD_3D"   << "\t"
                << "MSD_x"    << "\t"
                << "MSD_y"    << "\t"
                << "MSD_z"    << "\n";
            printer.print(oss.str());

            oss << std::scientific;

            oss.str("");
            oss << std::setprecision(6) << 0.0 << "\t"
                << std::setprecision(6) << 0.0 << "\t"
                << std::setprecision(6) << 0.0 << "\t"
                << std::setprecision(6) << 0.0 << "\t"
                << std::setprecision(6) << 0.0 << "\n";
            printer.print(oss.str());

            //------ body
            const PS::S32 n_sample = time_history.size();
            for(PS::S32 i_sample=0; i_sample<n_sample-1; ++i_sample){
                oss.str("");

                oss << std::setprecision(6) << time_history[i_sample]     << "\t"
                    << std::setprecision(6) << data_final[i_sample].r2_3d << "\t"
                    << std::setprecision(6) << data_final[i_sample].r2_x  << "\t"
                    << std::setprecision(6) << data_final[i_sample].r2_y  << "\t"
                    << std::setprecision(6) << data_final[i_sample].r2_z  << "\n";

                printer.print(oss.str());
            }

            //--- output report
            oss.str("");
            oss << "   MSD_Sampler: the result file '" << printer.file_name() << "' was generated." << "\n"
                << "                number of trace target: " << this->n_tgt_total << "\n";
            std::cout << oss.str() << std::flush;
        }

    public:
        template <class Tracer>
        void output(const Tracer &tracer) const {
            const auto&   move_trace   = tracer.get_move_trace();
            const auto&   time_history = tracer.get_time_history();
            const PS::S32 n_sample     = time_history.size();

            this->_impl_check_data_effectivity(n_sample);

            std::vector<R2_DATA> data_local;
            data_local.reserve(n_sample);
            data_local.clear();

            switch (this->sampling_mode){
                case MSD_SAMPLING_MODE::resampling:
                    if( tracer.isRegularInterval() ){
                        this->_impl_calc_msd_local_resampling(move_trace,
                                                              n_sample,
                                                              data_local);
                    } else {
                        if(PS::Comm::getRank() == this->rank){
                            std::ostringstream oss;
                            oss << " warning: the trajectory data is not regular intervals." << "\n"
                                << "          MSD result of '" << this->file_name << "' without resampling." << "\n";
                            std::cout << oss.str() << std::flush;
                        }
                        this->_impl_calc_msd_local(move_trace,
                                                   n_sample,
                                                   data_local);
                    }
                break;

                case MSD_SAMPLING_MODE::raw:
                    this->_impl_calc_msd_local(move_trace,
                                               n_sample,
                                               data_local);
                break;
            }

            this->_impl_write_result(time_history,
                                     data_local   );
        }

        //--- high level API
        template <class Tpsys, class Tracer>
        void output(      Tpsys               &psys,
                    const Tracer              &tracer,
                    const std::vector<Target> &tgt_list,
                    const std::string         &file_name,
                    const MSD_SAMPLING_MODE   &sampling_mode = MSD_SAMPLING_MODE::resampling,
                    const PS::S32              rank          = 0 ){
            this->input_setting(file_name, rank, sampling_mode);
            this->broadcast_setting(rank);
            this->make_id_table(psys, tracer);
            this->output(tracer);
        }

    };

}
