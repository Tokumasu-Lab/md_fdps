/**************************************************************************************************/
/**
* @file  logspace_array.hpp
* @brief the array with logarism space index.
*/
/**************************************************************************************************/
#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>


namespace MD_EXT {

    template <class T, class Texp = double>
    class logspace_array {
    private:

        struct Info {
            Texp begin = 0.0;
            Texp end   = 0.0;

            Texp log_begin = 1.0;
            Texp log_end   = 1.0;

            size_t n = 0;

            Texp r_log     = 1.0;
            Texp r_log_inv = 1.0;

            Texp r_real = 1.0;

            bool init_flag = false;
        };

        Info                            _info;
        std::vector<std::pair<Texp, T>> _data;

        void _check_init() const {
            #ifndef NDEBUG
            if( !this->_info.init_flag ){
                throw std::logic_error("the logspace_array is not initialized. call 'init()' at first.");
            }
            #endif
        }

        void _M_range_check(const Texp f) const {
            this->_check_init();

            if( f < this->_info.begin || this->_info.end <= f ){
                std::ostringstream oss;
                oss << "f must be 'begin <= f < end'." << "\n"
                    << "   begin = " << this->_info.begin << "\n"
                    << "   end   = " << this->_info.end   << "\n"
                    << "   f     = " << f << "\n";
                throw std::out_of_range(oss.str());
            }
        }

        int_fast64_t _get_index(const Texp f) const {
            const Texp log_f = std::log(f);
            return size_t( std::floor( (log_f - this->_info.log_begin)*this->_info.r_log_inv ) );
        }

    public:
        using reference              = typename decltype(_data)::reference;
        using const_reference        = typename decltype(_data)::const_reference;
        using iterator               = typename decltype(_data)::iterator;
        using const_iterator         = typename decltype(_data)::const_iterator;
        using size_type              = typename decltype(_data)::size_type;
        using difference_type        = typename decltype(_data)::difference_type;
        using value_type             = typename decltype(_data)::value_type;
        using pointer                = typename decltype(_data)::pointer;
        using const_pointer          = typename decltype(_data)::const_pointer;
        using reverse_iterator       = typename decltype(_data)::reverse_iterator;
        using const_reverse_iterator = typename decltype(_data)::const_reverse_iterator;

        //--- iterators
        iterator         begin()  noexcept { return this->_data.begin();  }
        iterator         end()    noexcept { return this->_data.end();    }
        reverse_iterator rbegin() noexcept { return this->_data.rbegin(); }
        reverse_iterator rend()   noexcept { return this->_data.rend();   }

        const_iterator   begin()   const noexcept { return this->_data.begin();  }
        const_iterator   end()     const noexcept { return this->_data.end();    }
        const_iterator   cbegin()  const noexcept { return this->_data.cbegin(); }
        const_iterator   cend()    const noexcept { return this->_data.cend();   }
        const_reverse_iterator crbegin() const noexcept { return this->crbegin(); }
        const_reverse_iterator crend()   const noexcept { return this->crend();   }


        T& at(const Texp f){
            this->_M_range_check(f);
            return this->_data[ this->_get_index(f) ].second;
        }
        const T& at(const Texp f) const {
            this->_M_range_check(f);
            return this->_data[ this->_get_index(f) ].second;
        }

        T&       operator [] (const Texp f)       { this->_data[ this->_get_index(f) ].second; }
        const T& operator [] (const Texp f) const { this->_data[ this->_get_index(f) ].second; }

        T&       front()       { return this->_data.front().second; }
        const T& front() const { return this->_data.front().second; }
        T&       back()        { return this->_data.back().second;  }
        const T& back()  const { return this->_data.back().second;  }

        size_type size()  const noexcept { return this->_data.size(); }
        bool      empty() const noexcept { return (this->_data.size() == 0); }

        Texp ratio()     const noexcept { return this->_info.r_real; }
        Texp min_range() const noexcept { return this->_info.begin;  }
        Texp max_range() const noexcept { return this->_info.end;    }
        std::pair<Texp, Texp> range() const noexcept {
            return std::make_pair(this->min_range(),
                                  this->max_range() );
        }
        std::pair<Texp, Texp> range(const Texp f) const {
            this->_M_range_check(f);
            const Texp begin = this->_data[ this->_get_index(f) ].first;
            const Texp end   = begin*this->ratio();
            return std::make_pair(begin, end);
        }

        void fill(const T &u){
            this->_check_init();
            for(auto& p : this->_data){
                p.second = u;
            }
        }
        void swap(logspace_array &rhs){
            const auto tmp = this->_info;
            this->_info = rhs._info;
            rhs._info   = tmp;

            this->_data.swap(rhs._data);
        }

        void init(const Texp   begin,
                  const Texp   end,
                  const size_t n     ){

            if( !(begin > 0.0) ||
                !(end   > 0.0) ||
                !(begin < end)   ){
                std::ostringstream oss;
                oss << "range is must be '0 <= begin < end'." << "\n"
                    << "   begin = " << begin << "\n"
                    << "   end   = " << end   << "\n";
                throw std::invalid_argument(oss.str());
            }
            if( !(n > 0) ){
                std::ostringstream oss;
                oss << "n must be > 0." "\n"
                    << "   n = " << n << "\n";
                throw std::invalid_argument(oss.str());
            }

            this->_info.begin = begin;
            this->_info.end   = end;
            this->_info.n     = n;

            this->_info.log_begin = std::log(begin);
            this->_info.log_end   = std::log(end);
            this->_info.r_log     = (this->_info.log_end - this->_info.log_begin)/static_cast<Texp>(n);
            this->_info.r_log_inv = 1.0/this->_info.r_log;

            this->_info.r_real = std::exp(this->_info.r_log);

            this->_data.clear();
            this->_data.resize(n);

            for(size_t i=0; i<this->_data.size(); ++i){
                this->_data[i].first = std::exp(  this->_info.log_begin
                                                + this->_info.r_log*static_cast<Texp>(i) );
            }

            this->_info.init_flag = true;
        }

        void resize_array(const size_t n_new){
            this->_check_init();
            if( !(n_new > 1) ){
                std::ostringstream oss;
                oss << "n_new must be > 1." "\n"
                    << "   n_new = " << n_new << "\n";
                throw std::invalid_argument(oss.str());
            }
            const size_t n_old = this->_data.size();

            this->_info.log_end = this->_info.log_begin + this->_info.r_log*static_cast<Texp>(n_new);
            this->_info.end     = std::exp(this->_info.log_end);

            this->_info.n = n_new;
            this->_data.resize(n_new);

            for(size_t i=n_old; i<this->_data.size(); ++i){
                this->_data[i].first = std::exp(  this->_info.log_begin
                                                + this->_info.r_log*static_cast<Texp>(i) );
            }
        }
        template <class Tin>
        void resize(const Tin x_max_in){
            const Texp x_max = static_cast<Texp>(x_max_in);
            if( !(x_max > this->min_range()) ){
                std::ostringstream oss;
                oss << "x_max must be > min_range()." "\n"
                    << "   x_max     = " << x_max << "\n"
                    << "   min_range = " << this->min_range() << "\n";
                throw std::length_error(oss.str());
            }

            const auto index_nead = this->_get_index(x_max);
            this->resize_array( static_cast<size_t>(index_nead + 1) );
        }

        logspace_array() = default;
        logspace_array(const Texp   begin,
                       const Texp   end,
                       const size_t n     ){
            this->init(begin, end, n);
        }
        logspace_array& operator = (const logspace_array& rhs) noexcept {
            this->_info = rhs._info;
            this->_data = rhs._data;
            return *this;
        }
        logspace_array(const logspace_array &rhs){
            *this = rhs;
        }

        ~logspace_array() = default;
    };

}
