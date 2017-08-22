//***************************************************************************************
//  This program is the intra pair table for "md_fdps_main.cpp"
//    search nearest intra pair particle in "CalcForce***" from EPJ particle array.
//    This code is using the Framework for Developing Particle Simulator (FDPS).
//    https://github.com/FDPS
//***************************************************************************************
#pragma once

#include<cassert>
#include<climits>
#include<unordered_map>
#include<unordered_set>

namespace MD_EXT {

    //--- data buffer for searching image paricle
    template<class Tindex>
    class ImageEPJ_ {
      private:
        std::pair<Tindex, PS::F64vec> pair;

      public:
        //--- constructor
        ImageEPJ_(const Tindex i, const PS::F64vec & r){
            this->pair = std::make_pair(i, r);
        }

        inline void setIndex(const Tindex i)     { this->pair.first  = i; }
        inline void setPos(const PS::F64vec & r) { this->pair.second = r; }
        inline Tindex     getIndex() const { return this->pair.first;  }
        inline PS::F64vec getPos()   const { return this->pair.second; }
    };

    //--- image particle list for same atomID
    template<class Tindex, std::size_t list_size>
    class ImageEPJ_List_ {
      private:
        std::vector<ImageEPJ_<Tindex>> list;

      public:
        //--- constructor
        ImageEPJ_List_(){
            this->list.reserve(list_size);
            //for(PS::S32 i=0; i<size; i++){
            //    this->list.emplace_back(-1, PS::F64vec(-1.0, -1.0, -1.0));
            //}
            //this->list.clear();
        }

       // //--- move interface
       // ImageEPJ_List_(ImageEPJ_List_&& value){
       //     this->list = std::move(value);
       // }
       // ImageEPJ_List_& operator=(ImageEPJ_List_&& value){
       //     return std::move(this->list);
       // }

        //--- access functions
        inline void add(const Tindex i, const PS::F64vec & pos){
            this->list.emplace_back(i, pos);
        }
        inline Tindex find(const PS::F64vec & pos_ref) const {
            PS::F64 r2_ref = std::numeric_limits<PS::F64>::max();  // largest F64 value
            Tindex  tgt;
            //--- search nearest image particle
            for(auto & buf : this->list){
                PS::F64vec r_ij   = buf.getPos() - pos_ref;
                PS::F64    r2_now = r_ij*r_ij;
                if(r2_now <= r2_ref){
                    tgt    = buf.getIndex();
                    r2_ref = r2_now;
                }
            }
            //--- return nearest EPJ index
            return tgt;
        }
        PS::S32 size() const { return this->list.size(); }

        //--- iterator interface
        using const_iterator = typename std::vector<ImageEPJ_<Tindex>>::const_iterator;
        const_iterator cbegin() const {
            return this->list.cbegin();
        }
        const_iterator cend() const {
            return this->list.cend();
        }
    };


    //--- interface class
    template<class Tindex, class Tid>
    class IntraPairTable {
      private:
        std::unordered_multimap<Tid, ImageEPJ_<Tindex>> table;

      public:
        //--- clear all data
        void clear(){
            this->table.clear();
        }
        //--- reserve table size
        void reserve(std::size_t size){
            this->table.reserve(size);
        }

        //--- constructor
        IntraPairTable(){
            this->table.clear();
            this->table.max_load_factor(0.7);
        }
        IntraPairTable(std::size_t size){
            this->table.clear();
            this->table.max_load_factor(0.7);
            this->reserve(size);
        }

        //--- access functions
        inline void add(const Tid id, const Tindex j, const PS::F64vec & pos){
            this->table.insert( std::pair<Tid, ImageEPJ_<Tindex>>(id, ImageEPJ_<Tindex>(j, pos)) );
        }
        inline Tindex find(const Tid id, const PS::F64vec & pos_ref) const {
            PS::F64 r2_ref = std::numeric_limits<PS::F64>::max();  // largest F64 value
            Tindex  tgt    = std::numeric_limits<Tindex>::max();   // illigal value

            assert( this->table.count(id) != 0 );

            auto range = this->table.equal_range(id);

            //--- search nearest pair
            for(auto itr = range.first; itr != range.second; itr++){
                PS::F64vec r_ij   = itr->second.getPos() - pos_ref;
                PS::F64    r2_now = r_ij*r_ij;
                if(r2_now < r2_ref){
                    tgt    = itr->second.getIndex();
                    r2_ref = r2_now;
                }
            }
            return tgt;
        }
        PS::S32 count(const Tid id) const {
            return this->table.count(id);
        }
        std::unordered_set<Tid> getKeyList() const {
            std::unordered_set<Tid> keys;
            std::for_each(this->table.cbegin(), this->table.cend(),
                          [&keys](std::pair<Tid, ImageEPJ_<Tindex>> elem) { keys.insert(elem.first); } );
            return std::move(keys);
        }

        //--- iterator interface
        using const_itr = typename std::unordered_multimap<Tid, ImageEPJ_<Tindex>>::const_iterator;
        const_itr cbegin() const {
            return this->table.cbegin();
        }
        const_itr cend() const {
            return this->table.cend();
        }
        std::pair<const_itr, const_itr> equal_range(const Tid id) const {
            return this->table.equal_range(id);
        }
    };

}
