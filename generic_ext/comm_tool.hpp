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
* @details The boost.MPI library and boost.Serialization library are recommended.
*/
namespace COMM_TOOL {

    //--- implemented in "comm_tool_***.hpp"

}
