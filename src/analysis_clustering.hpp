//***************************************************************************************
//  This is generic Clustering analyzer.
//***************************************************************************************
#pragma once

#include <string>
#include <sstream>
#include <unordered_map>
#include <limits>

#include <particle_simulator.hpp>
#include <molecular_dynamics_ext.hpp>

#include "file_Out_VMD.hpp"


namespace Analysis {

    namespace CLUSTERING_DEFS {

        template <class ID_type>
        class TreeResult {
        private:
            ID_type cluster_id;

        public:
            ID_type getClusterId() const { return this->cluster_id; }

            void clearId(){
                this->cluster_id = std::numeric_limits<ID_type>::max();
            }
            void clear(){ this->clearId(); }
            void update_id(const ID_type id){
                this->cluster_id = std::min(this->cluster_id, id);
            }

            template <class Tptcl>
            void copyFromForce(const Tptcl &ptcl){
                this->cluster_id = ptcl.getClusterId();
            }
        };

        template <class Tptcl>
        class TreePtcl {
        private:

            using ID_type  = decltype(declval<Tptcl>().getId());
            using PosType  = decltype(declval<Tptcl>().getPos());
            using AtomType = decltype(declval<Tptcl>().getAtomType());

            Tptcl   atom;
            ID_type cluster_id;
            bool    is_target          = false;
            bool    is_largest_cluster = false;

            static PS::F64 r_search;

        public:
            ID_type getId()        const { return this->atom.getId(); }
            ID_type getClusterId() const { return this->cluster_id;   }

            PosType getPos() const { return this->atom.getPos(); }

            void setClusterId(const ID_type id){
                this->cluster_id = id;
            }

            bool isTarget() const { return this->is_target; }
            void setTargetFlag(const bool flag) { this->is_target = flag; }

            bool isLargestCluster() const { return this->is_largest_cluster; }
            void setLargestClusterFlag(const bool flag) { this->is_largest_cluster = flag; }

            //--- static func
            static void    setRSearch(const PS::F64 r){ r_search = r; }
            static PS::F64 getRSearch()    { return r_search;          }
            static PS::F64 getRSearch_sq() { return r_search*r_search; }

            //--- wrapper for Tptcl
            void setPos(const PosType pos_new){
                this->atom.setPos(pos_new);
            }

            void clearId(){
                this->cluster_id         = atom.getId();
                this->is_largest_cluster = false;
            }
            void copyFromFP(const Tptcl    &ptcl){
                this->atom.copyFromFP(ptcl);
                this->clearId();
            }
            void copyFromFP(const TreePtcl &ptcl){
                this->atom.copyFromFP(ptcl.atom);
                this->cluster_id         = ptcl.getClusterId();
                this->is_target          = ptcl.isTarget();
                this->is_largest_cluster = ptcl.isLargestCluster();
            }
            void copyFromForce(const TreeResult<ID_type> &result){
                this->cluster_id = result.getClusterId();
            }

            //--- interface for FILE_IO::VMD::write_atom()
            ID_type     getAtomID()   const { return this->atom.getAtomID();   }
            AtomType    getAtomType() const { return this->atom.getAtomType(); }
            ID_type     getMolID()    const {
                if(this->is_target){
                    return this->getClusterId();
                } else {
                    return -1;
                }
            }
            std::string getResidue()  const {
                std::string s;
                if( this->isLargestCluster() ){
                    s = "LRG";
                } else {
                    if( this->isTarget() ){
                        s = "ELS";
                    } else {
                        s = "nt ";
                    }
                }
                return s;
            }
            void write_VMD_ascii(FILE *fp) const { FILE_IO::VMD::write_atom(fp, *this); }
        };
        template <class Tptcl>
        PS::F64 TreePtcl<Tptcl>::r_search;

        template <class Tresult, class Tepi, class Tepj>
        void JudgeCluster(const Tepi    *ep_i,
                          const PS::S32  n_ep_i,
                          const Tepj    *ep_j,
                          const PS::S32  n_ep_j,
                                Tresult *result ){

            const PS::F64 r2_search = ep_j->getRSearch_sq();

            Tresult tmp;
            for(PS::S32 i=0; i<n_ep_i; ++i){
                tmp.clear();

                //--- non target
                if(!ep_i[i].isTarget()){
                    result[i].copyFromForce(tmp);
                    continue;
                }

                //--- target
                tmp.update_id( ep_i[i].getClusterId() );
                for(PS::S32 j=0; j<n_ep_j; ++j){
                    if( !ep_j[j].isTarget() ) continue;

                    const PS::F64vec r_vec = ep_j[j].getPos() - ep_i[i].getPos();
                    const PS::F64    r2    = r_vec*r_vec;

                    if(r2 <= r2_search){
                        tmp.update_id( ep_j[j].getClusterId() );
                    }
                }
                result[i].copyFromForce(tmp);
            }
        }

    }

    template <class Target, class Tptcl>
    class Clustering {
    private:

        using ID_type = decltype(declval<Tptcl>().getId());

        using TreeResult = CLUSTERING_DEFS::TreeResult<ID_type>;
        using TreePtcl   = CLUSTERING_DEFS::TreePtcl<Tptcl>;


        std::string          statistics_file;
        std::string          histogram_file;
        FS_TOOL::FilePrinter printer;

        PS::S32              rank = -1;

        PS::F64             range;
        PS::S64             cycle_limit;
        std::vector<Target> tgt_list;

        bool setting_flag = false;
        bool sync_flag    = false;

        bool cycle_flag = false;

        bool ps_initialized = false;
        PS::ParticleSystem<TreePtcl>                          psys;
        typename PS::TreeForForceShort< TreeResult,
                                        TreePtcl,
                                        TreePtcl   >::Scatter tree;


        std::vector<ID_type>                     forword_table;  // index      -> cluster_id
        std::unordered_map<ID_type,
                           std::vector<ID_type>> inverse_map;    // cluster_id -> std::vector<index>

        std::unordered_map<ID_type, PS::U64>     size_map;       // cluster_id -> size
        std::pair<std::vector<ID_type>, PS::U64> largest_cluster;

        std::vector<            std::pair<ID_type, PS::U64>>  cls_size_send_buff;
        std::vector<std::vector<std::pair<ID_type, PS::U64>>> cls_size_recv_buff;


        size_t                         n_sample = 0;
        size_t                         n_pin_10x;
        MD_EXT::logspace_array<size_t> histogram_array;

        void _check_sync_flag() const {
            if(!this->sync_flag){
                std::ostringstream oss;
                oss << " MPI sync error. call 'broadcast_setting()' before execute analyze." << "\n";
                throw std::logic_error(oss.str());
            }
        }
        void _check_output_flag() const {
            if(this->n_sample <= 0){
                std::ostringstream oss;
                oss << " result data is empty. call 'add_sample()' before output result." << "\n";
                throw std::logic_error(oss.str());
            }
        }

    public:
        template <class Tprof>
        void init_local_buffer(const Tprof   &profile,
                               const PS::S64 &n_total){

            if(this->ps_initialized){
                throw std::logic_error("ERROR: internal buffer of Clustering was already initialized. call 'init_local_buffer()' at only once. ");
            }

            this->psys.initialize();
            this->tree.initialize(n_total,
                                  profile.theta,
                                  profile.n_leaf_limit,
                                  profile.n_group_limit);
            this->ps_initialized = true;

            this->inverse_map.max_load_factor(0.7);
            this->size_map.max_load_factor(0.7);

            const PS::S64 n_piece = static_cast<PS::S64>( (n_total)*1.5/PS::Comm::getNumberOfProc() );
            this->inverse_map.reserve(n_piece);
            this->size_map.reserve(n_piece);

            this->histogram_array.init( PS::F64(1),
                                        PS::F64(10),
                                        this->n_pin_10x );
            this->histogram_array.resize( static_cast<PS::F64>(n_total) );
            this->histogram_array.fill(0);
        }

        void input_setting(const std::string         &statistics_file,
                           const std::string         &histogram_file,
                           const std::vector<Target> &tgt_list,
                           const PS::F64              range,
                           const size_t               n_pin_10x = 6,
                           const PS::S64              cycle_limit = 1000,
                           const PS::S32              rank = 0    ){

            //--- argument check
            if(statistics_file.size() <= 0){
                throw std::invalid_argument("statistics_file name is invalid.");
            }
            if(histogram_file.size() <= 0){
                throw std::invalid_argument("histogram_file name is invalid.");
            }
            if(tgt_list.size() <= 0){
                std::ostringstream oss;
                oss << "tgt_list must have more than 1 target." << "\n";
                throw std::invalid_argument(oss.str());
            }
            if(range <= 0.0){
                std::ostringstream oss;
                oss << "range must be > 0.0" << "\n"
                    << "   range = " << range << "\n";
                throw std::invalid_argument(oss.str());
            }
            if(n_pin_10x < 2){
                std::ostringstream oss;
                oss << "n_pin_10x must be > 1." << "\n"
                    << "   n_pin_10x = " << n_pin_10x << "\n";
                throw std::invalid_argument(oss.str());
            }
            if(cycle_limit < 1){
                std::ostringstream oss;
                oss << "cycle_limit must be > 1." << "\n"
                    << "   cycle_limit = " << cycle_limit << "\n";
                throw std::invalid_argument(oss.str());
            }
            COMM_TOOL::check_proc_rank(rank);

            //--- set value
            this->statistics_file = statistics_file;
            this->histogram_file  = histogram_file;
            this->tgt_list        = tgt_list;
            this->range           = range;
            this->n_pin_10x       = n_pin_10x;
            this->cycle_limit     = cycle_limit;
            this->rank            = rank;

            this->setting_flag = true;

            this->sync_flag   = false;
        }
        Clustering() = default;
        Clustering(const std::string         &statistics_file,
                   const std::string         &histogram_file,
                   const std::vector<Target> &tgt_list,
                   const PS::F64              range,
                   const size_t               n_pin_10x = 6,
                   const PS::S64              cycle_limit = 1000,
                   const PS::S32              rank = 0    ){
            this->input_setting(statistics_file,
                                histogram_file,
                                tgt_list,
                                range,
                                n_pin_10x,
                                cycle_limit,
                                rank           );
        }
        void broadcast_setting(const PS::S32 rank = -1){
            if(this->n_sample > 0){
                throw std::logic_error("the Clustering analyzer has result data. do not reuse for other settings.");
            }

            PS::S32 source_proc = 0;
            if(rank < 0){
                //--- auto detect
                source_proc = PS::Comm::getMaxValue(this->rank);
            } else {
                //--- use presented rank
                source_proc = rank;
            }

            this->rank = source_proc;
            COMM_TOOL::broadcast(this->statistics_file, source_proc);
            COMM_TOOL::broadcast(this->histogram_file , source_proc);
            COMM_TOOL::broadcast(this->tgt_list       , source_proc);
            COMM_TOOL::broadcast(this->range          , source_proc);
            COMM_TOOL::broadcast(this->n_pin_10x      , source_proc);
            COMM_TOOL::broadcast(this->cycle_limit    , source_proc);
            COMM_TOOL::broadcast(this->setting_flag   , source_proc);

            if(!this->setting_flag){
                std::ostringstream oss;
                oss << "the setting is invalid."
                    << " call 'input_setting()' before calling 'broadcast_setting()' at source rank." << "\n"
                    << "   source_proc = " << source_proc << "\n";
                throw std::logic_error(oss.str());
            }

            //--- init printer
            this->printer.file_init(this->statistics_file, source_proc);
            std::ostringstream oss;
            oss << "time"      << "\t"
                << "n_cluster" << "\t"
                << "max_size"  << "\t"
                << "ave_size"  << "\t"
                << "ave_size_withoutMax" << "\n";
            this->printer.print(oss.str());

            this->n_sample = 0;

            this->sync_flag = true;
        }

        bool isTarget(const Tptcl &ptcl) const {
            bool result = false;
            for(const auto& tgt : this->tgt_list){
                if(tgt.isMatch(ptcl)){
                    result = true;
                    break;
                }
            }
            return result;
        }

    private:
        template <class Tpsys>
        void _impl_copy_atom_into_buffer(Tpsys &atom){
            const PS::S64 n_local = atom.getNumberOfParticleLocal();
            this->psys.setNumberOfParticleLocal(n_local);

            for(PS::S64 i=0; i<n_local; ++i){
                const bool tgt_flag = this->isTarget(atom[i]);

                this->psys[i].copyFromFP(atom[i]);
                this->psys[i].clearId();
                this->psys[i].setTargetFlag( tgt_flag );
                this->psys[i].setLargestClusterFlag( false );
            }
        }

        template <class Tdinfo>
        void _impl_cluster_search(Tdinfo &dinfo){
            const PS::S64 n_local = this->psys.getNumberOfParticleLocal();

            this->tree.calcForceAll(CLUSTERING_DEFS::JudgeCluster<TreeResult, TreePtcl, TreePtcl>,
                                    this->psys,
                                    dinfo);

            for(PS::S64 i=0; i<n_local; ++i){
                const auto result = this->tree.getForce(i);
                if( this->psys[i].getClusterId() != result.getClusterId() ){
                    this->cycle_flag = true;
                    this->psys[i].copyFromForce( result );
                }
            }
        }

        void _impl_deformation_forword_map(      std::unordered_map<ID_type, ID_type> &f_map,
                                           const ID_type                               id_prev,
                                                 ID_type                              &id_root ){

            //--- at root (terminate in other proc, root id is not exist in table)
            if( f_map.count(id_prev) == 0 ){
                id_root = id_prev;
                return;
            }

            ID_type& id_next = f_map.at(id_prev);
            assert(id_next <= id_prev);

            //--- at root (terminate in local)
            if(id_next == id_prev){
                id_root = id_next;
                return;
            }

            //--- not at root (recursive search)
            this->_impl_deformation_forword_map(f_map, id_next, id_root);

            //--- update pointer
            if(id_next != id_root){
                id_next = id_root;
                this->cycle_flag = true;
            }
        }
        void _impl_update_local_table(){
            const PS::S64 n_local = this->psys.getNumberOfParticleLocal();

            //--- make forword_map
            std::unordered_map<ID_type, ID_type> forword_map;
            forword_map.max_load_factor(0.7);
            forword_map.reserve(n_local);

            for(PS::S64 i=0; i<n_local; ++i){
                if( !this->psys[i].isTarget() ) continue;

                const ID_type atom_id    = this->psys[i].getId();
                const ID_type cluster_id = this->psys[i].getClusterId();
                forword_map[atom_id] = cluster_id;
            }

            //--- deformation forword_map as 2-level table
            for(auto& elem : forword_map){
                ID_type tmp;
                this->_impl_deformation_forword_map(forword_map, elem.first, tmp);
            }

            //--- update forword_table
            this->forword_table.clear();
            this->forword_table.resize(n_local, std::numeric_limits<ID_type>::max());
            for(PS::S64 i=0; i<n_local; ++i){
                if( this->psys[i].isTarget() ){
                    this->forword_table[i] = forword_map.at( this->psys[i].getId() );
                }
            }

            //--- update psys & inverse_map
            assert(static_cast<PS::S64>(this->forword_table.size()) == n_local);
            this->inverse_map.clear();
            for(PS::S64 i=0; i<n_local; ++i){
                const PS::S64& index      = i;
                const ID_type& cluster_id = this->forword_table[i];

                if( !this->psys[i].isTarget() ) continue;

                this->psys[i].setClusterId( cluster_id );

                if( this->inverse_map.count(cluster_id) == 0 ){
                    auto& list = this->inverse_map[cluster_id];
                    list.clear();
                    list.push_back(index);
                } else {
                    this->inverse_map.at(cluster_id).push_back(index);
                }
            }
        }

        template <class Tdinfo>
        void _impl_cluster_search_iterate(Tdinfo &dinfo){
            const PS::S64 n_local = this->psys.getNumberOfParticleLocal();

            assert( static_cast<PS::S64>(this->forword_table.size()) == n_local );
            assert( !this->inverse_map.empty() );

            this->tree.calcForceAll(CLUSTERING_DEFS::JudgeCluster<TreeResult, TreePtcl, TreePtcl>,
                                    this->psys,
                                    dinfo);

            //--- extract diff
            std::unordered_map<ID_type, ID_type> update_map;
            update_map.max_load_factor(0.7);
            update_map.reserve(n_local);

            for(PS::S64 i=0; i<n_local; ++i){
                const ID_type old_cls_id = this->psys[i].getClusterId();
                const ID_type new_cls_id = this->tree.getForce(i).getClusterId();

                if( new_cls_id == old_cls_id ) continue;

                this->cycle_flag = true;

                if( update_map.count(old_cls_id) == 0 ){
                    update_map[old_cls_id] = new_cls_id;
                } else {
                    auto& itr = update_map.at(old_cls_id);
                    itr = std::min(itr, new_cls_id);
                }
            }

            //--- update internal table
            for(const auto& diff : update_map){
                const ID_type old_cls_id = diff.first;
                const ID_type new_cls_id = diff.second;

                //--- alloc new list
                auto& old_list = this->inverse_map.at(old_cls_id);
                auto& new_list = this->inverse_map[new_cls_id];
                new_list.reserve(new_list.size() + old_list.size());

                //--- move pointer
                for(const auto& index : old_list){
                    //--- forword_table
                    this->forword_table.at(index) = new_cls_id;

                    //--- inverse_map
                    new_list.push_back(index);
                }
                old_list.clear();
                this->inverse_map.erase(old_cls_id);
            }

            //--- update psys
            for(PS::S64 i=0; i<n_local; ++i){
                if( !this->psys[i].isTarget() ) continue;

                this->psys[i].setClusterId( this->forword_table.at(i) );
            }
        }

        void _impl_collect_cluster_size(){
            //------ sum up in local
            this->size_map.clear();
            for(const auto& cls : this->inverse_map){
                const ID_type cls_id   = cls.first;
                const PS::U64 cls_size = cls.second.size();

                if( cls_size == 0 ) continue;

                if( this->size_map.count(cls_id) == 0 ){
                    this->size_map[cls_id]     = cls_size;
                } else {
                    this->size_map.at(cls_id) += cls_size;
                }
            }

            //------ collect inter proc
            this->cls_size_send_buff.clear();
            this->cls_size_send_buff.reserve( this->size_map.size() );
            for(const auto& elem : this->size_map){
                if(elem.second > 0){
                    this->cls_size_send_buff.push_back(elem);
                }
            }

            COMM_TOOL::allGather(this->cls_size_send_buff, this->cls_size_recv_buff);

            this->size_map.clear();
            for(const auto& v_proc : this->cls_size_recv_buff){
                for(const auto& cls : v_proc){
                    const auto cls_id   = cls.first;
                    const auto cls_size = cls.second;

                    if(this->size_map.count(cls_id) == 0){
                        this->size_map[cls_id]     = cls_size;
                    } else {
                        this->size_map.at(cls_id) += cls_size;
                    }
                }
            }

            //--- get max cluster id & size
            auto& max_cls_id_list = this->largest_cluster.first;
            auto& max_cls_size    = this->largest_cluster.second;
            max_cls_id_list = { std::numeric_limits<ID_type>::max() };
            max_cls_size    = 0;
            for(const auto& cls : this->size_map){
                const auto cls_id   = cls.first;
                const auto cls_size = cls.second;

                if(       max_cls_size == cls_size ){
                    max_cls_id_list.push_back(cls_id);
                } else if(max_cls_size <  cls_size){
                    max_cls_id_list = { cls_id };
                    max_cls_size    = cls_size;
                }
            }

            //--- set max cluster flag
            for(const auto& cls_id : this->largest_cluster.first){
                if(this->inverse_map.count(cls_id) == 0) continue;  // the cls is not in this proc.

                const auto& index_list = this->inverse_map.at(cls_id);
                for(const auto& index : index_list){
                    this->psys[index].setLargestClusterFlag( true );
                }
            }
        }

        template <class Tprof>
        void _impl_calc_statistics(const Tprof &profile){
            PS::S64 n_cls          = 0;
            PS::F64 size_ave_all   = 0.0;
            PS::F64 size_ave_woMax = 0.0;

            for(const auto& cls : this->size_map){
                if(cls.second > 0){
                    ++n_cls;
                    size_ave_all += cls.second;

                    //--- count up histogram
                    this->histogram_array.at(cls.second) += 1;
                }
            }

            const auto& max_cls_list = this->largest_cluster.first;
            const auto& max_cls_size = this->largest_cluster.second;

            size_ave_woMax = size_ave_all - static_cast<PS::F64>( max_cls_list.size()*max_cls_size );

            size_ave_all   = size_ave_all/static_cast<PS::F64>(n_cls);

            if( n_cls == static_cast<PS::S64>(max_cls_list.size()) ){
                // all clusters are same size (= max size)
                size_ave_woMax = 0.0;
            } else {
                size_ave_woMax = size_ave_woMax/static_cast<PS::F64>(n_cls - max_cls_list.size());
            }

            //--- write file
             std::ostringstream oss;
             oss << profile.get_time() << "\t"
                 << n_cls              << "\t"
                 << max_cls_size       << "\t"
                 << size_ave_all       << "\t"
                 << size_ave_woMax     << "\n";
             this->printer.print(oss.str());
        }

    public:
        template <class Tprof, class Tpsys, class Tdinfo>
        void add_sample(const Tprof  &profile,
                              Tpsys  &atom,
                              Tdinfo &dinfo   ){

            this->_check_sync_flag();

            if(!this->ps_initialized){
                const PS::U64 n_total = atom.getNumberOfParticleGlobal();
                this->init_local_buffer(profile, n_total);
            }

            //--- copy ptcl into buffer
            this->_impl_copy_atom_into_buffer(atom);

            //--- set RSearch
            TreePtcl::setRSearch( Normalize::normCutOff(this->range) );

            //--- clustring first step
            this->forword_table.clear();
            this->inverse_map.clear();
            this->cycle_flag = false;
            this->_impl_cluster_search(dinfo);
            this->_impl_update_local_table();

            //--- clustering iterate
            PS::S64 n_cycle = 0;
            while (this->cycle_flag){
                ++n_cycle;

                this->cycle_flag = false;
                this->_impl_cluster_search_iterate(dinfo);
                this->cycle_flag = PS::Comm::synchronizeConditionalBranchOR(this->cycle_flag);

                if(n_cycle > this->cycle_limit){
                    std::ostringstream oss;
                    oss << "did not completed in cycle_limit." << "\n"
                        << "   cycle_limit = " << this->cycle_limit << "\n";
                    throw std::logic_error(oss.str());
                }
            }

            this->_impl_collect_cluster_size();
            this->_impl_calc_statistics(profile);

            ++(this->n_sample);

            if(PS::Comm::getRank() == this->rank){
                std::ostringstream oss;
                oss << "  Clustering: iteration = " << n_cycle << "\n";
                std::cout << oss.str() << std::flush;
            }
        }

        //--- output interface
        const ID_type              getClusterId(const PS::S64 i) const { return this->psys[i].getClusterId(); }
        const std::vector<ID_type> getLargestClusterId()         const { return this->largest_cluster.first;      }
        const PS::U64              getLargestClusterSize()       const { return this->largest_cluster.second;     }

        void output_particle(const std::string &ptcl_file,
                             const bool         extract_flag = false){
            this->_check_sync_flag();
            this->_check_output_flag();

            if(extract_flag){
                //--- output target atoms only
                const PS::S64 n_local = this->psys.getNumberOfParticleLocal();
                PS::ParticleSystem<TreePtcl> psys_tmp;

                psys_tmp.setNumberOfParticleLocal(n_local);
                PS::S64 n_tmp = 0;
                for(PS::S64 i=0; i<n_local; ++i){
                    if( this->psys[i].isTarget() ){
                        psys_tmp[n_tmp].copyFromFP( this->psys[i] );
                        ++n_tmp;
                    }
                }
                //--- shrink buffer
                psys_tmp.setNumberOfParticleLocal(n_tmp);

                //--- output
                psys_tmp.writeParticleAscii( ptcl_file.c_str(), &TreePtcl::write_VMD_ascii );

            } else {
                //--- output all atoms
                this->psys.writeParticleAscii( ptcl_file.c_str(), &TreePtcl::write_VMD_ascii );
            }
        }

        void output_histogram() const {
            this->_check_sync_flag();
            this->_check_output_flag();

            if(PS::Comm::getRank() != this->rank) return;

            FS_TOOL::FilePrinter histogram_printer;
            histogram_printer.file_init(this->histogram_file, this->rank);

            std::ostringstream oss;

            //--- header
            oss << "cls_size"    << "\t"
                << "prpbability" << "\n";
            histogram_printer.print(oss.str());

            //--- body
            const PS::F64 n_sample_inv = 1.0/static_cast<PS::F64>(this->n_sample);
            oss << std::scientific;
            for(const auto& cls_count : this->histogram_array){
                //const PS::U64 cls_size    = static_cast<PS::U64>(cls_count.first) + 1;
                const auto    cls_size    = cls_count.first;
                const PS::F64 probability = static_cast<PS::F64>(cls_count.second)*n_sample_inv;

                oss.str("");
                oss << cls_size    << "\t"
                    << probability << "\n";
                histogram_printer.print(oss.str());
            }

            //--- report
            oss.str("");
            oss << "  Clustering: logspace histogram" << "\n"
                << "      range_min = " << this->histogram_array.range().first  << "\n"
                << "      range_max = " << this->histogram_array.range().second << "\n"
                << "      pin_size  = " << this->histogram_array.ratio()
                                        << " | devided by " << this->n_pin_10x << " at 10x space." "\n"
                << "      n_pin     = " << this->histogram_array.size()  << "\n";

            oss << "    the result file '" << this->histogram_file << "' was generated." << "\n";
            std::cout << oss.str() << std::flush;
        }



        //--- naive implementation
        template <class Tprof, class Tpsys, class Tdinfo>
        void add_sample_naive(const Tprof  &profile,
                                    Tpsys  &atom,
                                    Tdinfo &dinfo   ){

            this->_check_sync_flag();

            if(!this->ps_initialized){
                const PS::U64 n_total = atom.getNumberOfParticleGlobal();
                this->init_local_buffer(profile, n_total);
            }

            //--- copy ptcl into buffer
            this->_impl_copy_atom_into_buffer(atom);

            //--- set RSearch
            TreePtcl::setRSearch( Normalize::normCutOff(this->range) );

            //--- clustring first step
            this->forword_table.clear();
            this->inverse_map.clear();
            this->cycle_flag = false;
            this->_impl_cluster_search(dinfo);

            //--- clustering iterate
            PS::S64 n_cycle = 0;
            while (this->cycle_flag){
                ++n_cycle;

                this->cycle_flag = false;
                this->_impl_cluster_search(dinfo);
                this->cycle_flag = PS::Comm::synchronizeConditionalBranchOR(this->cycle_flag);

                if(n_cycle > this->cycle_limit){
                    std::ostringstream oss;
                    oss << "did not completed in cycle_limit." << "\n"
                        << "   cycle_limit = " << this->cycle_limit << "\n";
                    throw std::logic_error(oss.str());
                }
            }

            //--- make inverse_map for count clusters
            const PS::S64 n_local = this->psys.getNumberOfParticleLocal();
            this->inverse_map.clear();
            for(PS::S64 i=0; i<n_local; ++i){
                const PS::S64& index      = i;
                const ID_type& cluster_id = this->psys[i].getClusterId();

                if( !this->psys[i].isTarget() ) continue;

                if( this->inverse_map.count(cluster_id) == 0 ){
                    auto& list = this->inverse_map[cluster_id];
                    list.clear();
                    list.push_back(index);
                } else {
                    this->inverse_map.at(cluster_id).push_back(index);
                }
            }

            this->_impl_collect_cluster_size();
            this->_impl_calc_statistics(profile);

            ++(this->n_sample);

            if(PS::Comm::getRank() == this->rank){
                std::ostringstream oss;
                oss << "  Clustering: iteration = " << n_cycle << " (naive)" << "\n";
                std::cout << oss.str() << std::flush;
            }
        }

    };

}
