//***************************************************************************************
//  This is unit test of loading model parameter function.
//***************************************************************************************

#include <fstream>
#include <iostream>
#include <vector>
#include <random>

#include <particle_simulator.hpp>

#include "boltzmann_dist.hpp"

int main(int argc, char *argv[]) {
    PS::Initialize(argc, argv);

    constexpr size_t n         = 400000;
    constexpr size_t seed      = 54321;
    std::string      file_name = "blz_dist_result.txt";

    std::mt19937_64                  mt(seed);
    std::uniform_real_distribution<> dist(0.0, 1.0);   // [0.0, 1.0)

    //--- instance of boltzmann_dist
    boltzmann_dist  blz_dist;


    std::vector<PS::F64> cumulative_dist;
    PS::F64              range_min, dv;

    std::vector<PS::S32> blz_count;

    //--- get table in boltzmann_dist
    blz_dist.get_cumulative_dist(cumulative_dist, range_min, dv);


    //--- count blz_dist result
    blz_count.resize(cumulative_dist.size());
    for(size_t i=0; i<n; ++i){
        PS::F64 blz_value = blz_dist( PS::F64(dist(mt)) );

        PS::S32 index    =  PS::S32( (blz_value - range_min)/dv );
        blz_count[index] += 1;
    }


    //--- output result
    std::ofstream file(file_name);

    file << "x,            "
         << "blz_count     "
         << "blz_cumulative" << std::endl;

    for(size_t i=0; i<cumulative_dist.size(); ++i){
        file << range_min + PS::F64(i)*dv << " "
             << blz_count[i]              << " "
             << cumulative_dist[i]        << std::endl;
    }

    file.close();


    //--- show instruction
    std::cout << std::endl
              << "  --- type the command in below to show test result. --- " << std::endl
              << "  $ gnuplot" << std::endl
              << "  > plot '" << file_name << "' u 1:2 w p" << std::endl;

    PS::Finalize();
    return 0;
}
