//***************************************************************************************
//  This program is loading model parameter function.
//***************************************************************************************
#pragma once

#include <vector>
#include <cmath>

#include <particle_simulator.hpp>

namespace MD_EXT {

    class boltzmann_dist {
    private:
        PS::F64 range_min = 0.0;
        PS::F64 range_max = 4.0;

        std::vector<PS::F64> table;
        PS::F64              dv;

        PS::F64 f_bar(const PS::F64 &v){
            constexpr PS::F64 pi_sqrtinv = 1.0/std::sqrt(3.141592653589793);
            return 4.0*pi_sqrtinv*v*v*std::exp(-v*v);
        }

    public:
        //--- initialize function
        void init(const PS::F64 range_min,
                  const PS::F64 range_max,
                  const size_t  resolution,
                  const size_t  n_integral ){

            assert(range_min >= 0.0);
            assert(range_max > range_min);
            assert(resolution >= 10);
            assert(n_integral >= 10);

            //--- save parameters
            this->range_min = range_min;
            this->range_max = range_max;
            this->dv        = (range_max - range_min)/resolution;

            //--- make distribution table
            for(size_t i=0; i<resolution-1; ++i){
                PS::F64 v = range_min + this->dv*PS::F64(i);

                //--- integrate by simpson equation
                PS::F64 tmp = 0.0;
                PS::F64 dj = (v - range_min)/PS::F64(2*n_integral);
                for(size_t j=1; j<=n_integral; ++j){
                    PS::F64 x = range_min + dj*(2.0*j - 1);
                    tmp += ( f_bar(x - dj) + 4.0*f_bar(x) + f_bar(x + dj) )/3.0;
                }
                tmp = tmp*dj;

                this->table.emplace_back(tmp);
            }
            this->table.emplace_back(1.0);
        }

        //--- constructor
        boltzmann_dist(){
            this->init(0.0, 4.0, 1023, 512);
        }
        boltzmann_dist(const PS::F64 range_min,
                       const PS::F64 range_max){
            this->init(range_min, range_max, 1023, 512);
        }
        boltzmann_dist(const PS::F64 range_min,
                       const PS::F64 range_max,
                       const PS::S32 resolution,
                       const PS::S32 n_integral){
            this->init(range_min, range_max, resolution, n_integral);
        }


        PS::F64 gen(const PS::F64 r){
            //--- input range = [0.0, 1.0)
            assert(r >= 0.0 && r < 1.0);

            PS::F64 result;
            for(size_t i=0; i<this->table.size(); ++i){
                if(this->table[i] <= r && r < this->table[i+1]){
                    PS::F64 x = PS::F64(i);
                            x += (r - this->table[i])/(this->table[i+1] - this->table[i]);  // linear interpolation
                    result = this->range_min + this->dv*x;
                    break;
                }
            }
            return result;
        }
        PS::F64 operator () (const PS::F64 r){ return this->gen(r); }

        void get_cumulative_dist(std::vector<PS::F64> &dist,
                                 PS::F64              &range_min,
                                 PS::F64              &dv       ){
            dist.clear();

            dist      = this->table;
            range_min = this->range_min;
            dv        = this->dv;
        }
    };

}
