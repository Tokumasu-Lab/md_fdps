/**************************************************************************************************/
/**
* @file  comm_tool.hpp
* @brief STL container adapter for PS::Comm
*/
/**************************************************************************************************/
#pragma once

#include "comm_tool_SerDes.hpp"

#include "comm_tool_broadcast.hpp"
#include "comm_tool_gather.hpp"
#include "comm_tool_scatter.hpp"
#include "comm_tool_allGather.hpp"
#include "comm_tool_allToAll.hpp"


/**
* @brief wrappers for PS::Comm
* @details STL container adapter for MPI communication.
* @details If you need more higher peformance or more generalized function,
* @details The boost.MPI library and boost.Serialization library are recommended.
*/
namespace COMM_TOOL {

    /**
    * @brief process rank checker.
    */
    void check_proc_rank(const PS::S32 rank){
        if(rank < 0 || rank >= PS::Comm::getNumberOfProc()){
            std::ostringstream oss;
            oss << " rank = " << rank << " is invalid." << "\n"
                << "   exist process: 0 - " << PS::Comm::getNumberOfProc()-1 << "\n";
            throw std::out_of_range(oss.str());
        }
    }

}
