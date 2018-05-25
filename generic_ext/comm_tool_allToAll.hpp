/**************************************************************************************************/
/**
* @file  comm_tool_allToAll.hpp
* @brief STL container wrapper for PS::Comm::allToAll()
*/
/**************************************************************************************************/
#pragma once

#include <cassert>
#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <unordered_map>

#include <particle_simulator.hpp>

#include "comm_tool_SerDes.hpp"


namespace COMM_TOOL {

    //=====================
    //  wrapper functions
    //=====================

    /**
    * @brief wrapper for static size class.
    * @param[in]  send_vec send target.
    *                      MUST NOT contain pointer member.
    * @param[out] recv_vec recieve data.
    * @details    send/recv_vec[i] is send/recieve data to/from process i.
    */
    template <class T>
    void allToAll(const std::vector<T> &send_vec,
                        std::vector<T> &recv_vec,
                        MPI_Comm        comm = MPI_COMM_WORLD){

        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** alltoall_<std::vector<T>>" << std::endl;
        #endif

        const int n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec.size() >= static_cast<size_t>(n_proc));

        recv_vec.resize(n_proc);

        #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Alltoall(&send_vec[0], 1, PS::GetDataType<T>(),
                         &recv_vec[0], 1, PS::GetDataType<T>(),
                          comm);
        #else
            recv_vec[0] = send_vec[0];
        #endif
    }

    /**
    * @brief specialization for std::vector<T>.
    * @param[in]  send_vec_vec send target.
    *                         MUST NOT contain pointer member.
    * @param[out] recv_vec_vec recieve data.
    * @details    send/recv_vec[i] is send/recieve data to/from process i.
    */
    template <class T>
    void allToAll(const std::vector<std::vector<T>> &send_vec_vec,
                        std::vector<std::vector<T>> &recv_vec_vec,
                        MPI_Comm                     comm = MPI_COMM_WORLD){

        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** allToAll_<std::vector<std::vector<T>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_vec.size() >= static_cast<size_t>(n_proc));

        #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            std::vector<PS::S32> n_send;
            std::vector<PS::S32> n_send_disp;
            std::vector<PS::S32> n_recv;
            std::vector<PS::S32> n_recv_disp;

            n_send.resize(n_proc);
            n_recv.resize(n_proc);
            n_send_disp.resize(n_proc+1);
            n_recv_disp.resize(n_proc+1);

            n_send_disp[0] = 0;
            for(PS::S32 i=0; i<n_proc; ++i){
                n_send[i] = send_vec_vec[i].size();
            }
            allToAll(n_send, n_recv, comm);

            n_recv_disp[0] = 0;
            for(PS::S32 i=0; i<n_proc; ++i){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }

            std::vector<T> send_data;
            std::vector<T> recv_data;

            serialize_vector_vector(send_vec_vec, send_data, n_send_disp);
            recv_data.resize(n_recv_disp.back());

            MPI_Alltoallv(&send_data[0], &n_send[0], &n_send_disp[0], PS::GetDataType<T>(),
                          &recv_data[0], &n_recv[0], &n_recv_disp[0], PS::GetDataType<T>(),
                           comm);

            deserialize_vector_vector(recv_data, n_recv_disp, recv_vec_vec);
        #else
            recv_vec_vec.resize(1);
            recv_vec_vec[0] = send_vec_vec[0];
        #endif
    }

    /**
    * @brief specialization for std::string.
    * @param[in]  send_vec_str send target.
    *                         MUST NOT contain pointer member.
    * @param[out] recv_vec_str recieve data.
    * @details    send/recv_vec[i] is send/recieve data to/from process i.
    */
    void allToAll(const std::vector<std::string> &send_vec_str,
                        std::vector<std::string> &recv_vec_str,
                        MPI_Comm                  comm = MPI_COMM_WORLD){

        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** allToAll_<std::vector<std::string>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_str.size() >= static_cast<size_t>(n_proc));

        std::vector<std::vector<char>> send_vec_vec_char;
        std::vector<std::vector<char>> recv_vec_vec_char;

        send_vec_vec_char.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            serialize_string(send_vec_str[i], send_vec_vec_char[i]);
        }

        allToAll(send_vec_vec_char, recv_vec_vec_char, comm);

        recv_vec_str.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            deserialize_string(recv_vec_vec_char.at(i),
                               recv_vec_str.at(i)      );
        }
    }

    /**
    * @brief specialization for std::vector<std::string>.
    * @param[in]  send_vec_str send target.
    *                         MUST NOT contain pointer member.
    * @param[out] recv_vec_str recieve data.
    * @details    send/recv_vec[i] is send/recieve data to/from process i.
    */
    void allToAll(const std::vector<std::vector<std::string>> &send_vec_vec_str,
                        std::vector<std::vector<std::string>> &recv_vec_vec_str,
                        MPI_Comm                               comm = MPI_COMM_WORLD){

        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** allToAll_<std::vector<std::string>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_vec_str.size() >= static_cast<size_t>(n_proc));

        std::vector<std::vector<char>> send_vec_vec_char;
        std::vector<std::vector<char>> recv_vec_vec_char;

        send_vec_vec_char.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            serialize_vector_string(send_vec_vec_str[i], send_vec_vec_char[i]);
        }

        allToAll(send_vec_vec_char, recv_vec_vec_char, comm);

        recv_vec_vec_str.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            deserialize_vector_string(recv_vec_vec_char[i], recv_vec_vec_str[i]);
        }
    }

    /**
    * @brief specialization for std::vector<std::vector<T>>.
    * @param[in] send_vec_vv send target.
    *                         MUST NOT contain pointer member.
    * @param[out] recv_vec_vv recieve data.
    * @details    send/recv_vec[i] is send/recieve data to/from process i.
    * @details    devide std::vector<std::vector<T>> into std::vector<T> and std::vector<index>, then call alltoallV() for std::vector<T>.
    * @details    class T accepts std::vector<>, std::string, or user-defined class WITHOUT pointer member.
    * @details    If 4 or more nested std::vector<...> is passed, this function will call itself recurcively.
    */
    template <class T>
    void allToAll(const std::vector<std::vector<std::vector<T>>> &send_vec_vv,
                        std::vector<std::vector<std::vector<T>>> &recv_vec_vv,
                        MPI_Comm                                  comm = MPI_COMM_WORLD){

        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** allToAll_<std::vector<std::vector<std::vector<T>>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_vv.size() >= static_cast<size_t>(n_proc));

        std::vector<std::vector<T>>      send_data;
        std::vector<std::vector<size_t>> send_index;

        send_data.resize(n_proc);
        send_index.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            serialize_vector_vector(send_vec_vv[i],
                                    send_data[i], send_index[i]);
        }

        std::vector<std::vector<T>>      recv_data;
        std::vector<std::vector<size_t>> recv_index;

        allToAll(send_data,  recv_data , comm);
        allToAll(send_index, recv_index, comm);

        recv_vec_vv.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            deserialize_vector_vector(recv_data[i], recv_index[i],
                                      recv_vec_vv[i]              );
        }
    }

    /**
    * @brief specialization for std::vector<std::pair<Ta, Tb>>.
    * @param[in] send_vec_vp send target.
    *                         MUST NOT contain pointer member.
    * @param[out] recv_vec_vp recieve data.
    * @details    send/recv_vec[i] is send/recieve data to/from process i.
    * @details    devide std::vector<std::pair<Ta, Tb>> into std::vector<Ta> and std::vector<Tb>, then call alltoallV() for std::vector<T>.
    * @details    class Ta and Tb accept std::vector<>, std::string, or user-defined class WITHOUT pointer member.
    */
    template <class Ta, class Tb>
    void allToAll(const std::vector<std::vector<std::pair<Ta, Tb>>> &send_vec_vp,
                        std::vector<std::vector<std::pair<Ta, Tb>>> &recv_vec_vp,
                        MPI_Comm                                     comm = MPI_COMM_WORLD){

        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** allToAll_<std::vector<std::pair<Ta, Tb>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_vp.size() >= static_cast<size_t>(n_proc));

        std::vector<std::vector<Ta>> send_vec_1st;
        std::vector<std::vector<Tb>> send_vec_2nd;

        send_vec_1st.resize(n_proc);
        send_vec_2nd.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            split_vector_pair(send_vec_vp[i],
                              send_vec_1st[i], send_vec_2nd[i]);
        }

        std::vector<std::vector<Ta>> recv_vec_1st;
        std::vector<std::vector<Tb>> recv_vec_2nd;
        allToAll(send_vec_1st, recv_vec_1st, comm);
        allToAll(send_vec_2nd, recv_vec_2nd, comm);

        recv_vec_vp.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            combine_vector_pair(recv_vec_1st[i], recv_vec_2nd[i],
                                recv_vec_vp[i]                   );
        }
    }


    /**
    * @brief specialization for returning type interface.
    * @param[in] send_vec send target.
    * @return    recv_vec recieved data.
    * @details   wrapper for the use case: "auto recv_vec = COMM_TOOL::allToAll(vec);"
    * @details   recv_vec.size() = PS::Comm::getNumberOfProc();
    * @details   send/recv_vec[i] is send/recieve data to/from process i.
    */
    template <class T>
    std::vector<T> allToAll(const std::vector<T> &send_vec,
                                  MPI_Comm        comm = MPI_COMM_WORLD){

        std::vector<T> recv_vec;
        allToAll(send_vec, recv_vec, comm);
        return recv_vec;
    }

}
