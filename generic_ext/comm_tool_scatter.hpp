/**************************************************************************************************/
/**
* @file  comm_tool_scatter.hpp
* @brief STL container wrapper for PS::Comm::scatter()
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
    * @param[out] recv_data recieve data.
    * @details    send_vec[i] is send data to process rank i.
    */
    template <class T>
    void scatter(const std::vector<T> &send_vec,
                       T              &recv_data,
                 const PS::S32         root     ,
                       MPI_Comm        comm = MPI_COMM_WORLD){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** scatter_<std::vector<T>>" << std::endl;
        #endif

        const int n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec.size() >= static_cast<size_t>(n_proc));
        assert(0 <= root && root < n_proc);

        #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            MPI_Scatter(&send_vec[0], 1, PS::GetDataType<T>(),
                        &recv_data  , 1, PS::GetDataType<T>(),
                         root,
                         comm);
        #else
            recv_data = send_vec[0];
        #endif
    }

    /**
    * @brief specialization for std::vector<T>.
    * @param[in]  send_vec_vec send target.
    *                          MUST NOT contain pointer member.
    * @param[out] recv_vec recieve data.
    * @details    send_vec_vec[i] is send data to process rank i.
    */
    template <class T>
    void scatter(const std::vector<std::vector<T>> &send_vec_vec,
                       std::vector<T>              &recv_vec,
                 const PS::S32                      root,
                       MPI_Comm                     comm = MPI_COMM_WORLD){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** scatter_<std::vector<std::vector<T>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_vec.size() >= static_cast<size_t>(n_proc));
        assert(0 <= root && root < n_proc);

        #ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            std::vector<PS::S32> n_send;

            n_send.resize(n_proc);
            for(PS::S32 i=0; i<n_proc; ++i){
                n_send[i] = send_vec_vec[i].size();
            }

            PS::S32 n_recv;
            scatter(n_send, n_recv, root);
            recv_vec.resize(n_recv);

            std::vector<T>       send_data;
            std::vector<PS::S32> n_send_disp;

            serialize_vector_vector(send_vec_vec, send_data, n_send_disp);

            MPI_Scatterv(&send_data[0], &n_send[0], &n_send_disp[0], PS::GetDataType<T>(),
                         &recv_vec[0] ,  n_recv   ,                  PS::GetDataType<T>(),
                          root,
                          comm);
        #else
            recv_vec.resize(1);
            recv_vec[0] = send_vec_vec[0];
        #endif
    }

    /**
    * @brief specialization for std::string.
    * @param[in]  send_vec_str send target.
    * @param[out] recv_str recieve data.
    * @details    send_vec_str[i] is send data to process rank i.
    */
    void scatter(const std::vector<std::string> &send_vec_str,
                       std::string              &recv_str,
                 const PS::S32                   root,
                       MPI_Comm                  comm = MPI_COMM_WORLD){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** scatter_<std::vector<std::string>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_str.size() >= static_cast<size_t>(n_proc));
        assert(0 <= root && root < n_proc);

        std::vector<std::vector<char>> send_vec_vec_char;
        std::vector<char>              recv_vec_char;

        send_vec_vec_char.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            serialize_string(send_vec_str[i], send_vec_vec_char[i]);
        }

        scatter(send_vec_vec_char, recv_vec_char, root, comm);

        deserialize_string(recv_vec_char, recv_str);
    }

    /**
    * @brief specialization for std::vector<std::string>.
    * @param[in]  send_vec_vec_str send target.
    * @param[out] recv_vec_str recieve data.
    * @details    send_vec_vec_str[i] is send data to process rank i.
    */
    void scatter(const std::vector<std::vector<std::string>> &send_vec_vec_str,
                       std::vector<std::string>              &recv_vec_str,
                 const PS::S32                                root,
                       MPI_Comm                               comm = MPI_COMM_WORLD){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** scatter_<std::std::vector<vector<std::string>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_vec_str.size() >= static_cast<size_t>(n_proc));
        assert(0 <= root && root < n_proc);

        std::vector<std::vector<char>> send_vec_vec_char;
        std::vector<char>              recv_vec_char;

        send_vec_vec_char.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            serialize_vector_string(send_vec_vec_str[i], send_vec_vec_char[i]);
        }

        scatter(send_vec_vec_char, recv_vec_char, root, comm);

        deserialize_vector_string(recv_vec_char, recv_vec_str);
    }

    /**
    * @brief specialization for std::vector<std::vector<T>>.
    * @param[in] send_vec_vv send target.
    *                         MUST NOT contain pointer member.
    * @param[out] recv_vec_v recieve data.
    * @details    send_vec_vec_str[i] is send data to process rank i.
    * @details    devide std::vector<std::vector<T>> into std::vector<T> and std::vector<index>, then call scatterV() for std::vector<>.
    * @details    class T accepts std::vector<>, std::string, or user-defined class WITHOUT pointer member.
    * @details    If 4 or more nested std::vector<...> is passed, this function will call itself recurcively.
    */
    template <class T>
    void scatter(const std::vector<std::vector<std::vector<T>>> &send_vec_vv,
                       std::vector<std::vector<T>>              &recv_vec_v,
                 const PS::S32                                   root,
                       MPI_Comm                                  comm = MPI_COMM_WORLD){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** scatter_<std::vector<std::vector<std::vector<T>>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_vv.size() >= static_cast<size_t>(n_proc));
        assert(0 <= root && root < n_proc);

        std::vector<std::vector<T>>      send_data;
        std::vector<std::vector<size_t>> send_index;

        send_data.resize(n_proc);
        send_index.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            serialize_vector_vector(send_vec_vv[i],
                                    send_data[i], send_index[i]);
        }

        std::vector<T>      recv_data;
        std::vector<size_t> recv_index;

        scatter(send_data,  recv_data , root, comm);
        scatter(send_index, recv_index, root, comm);

        deserialize_vector_vector(recv_data, recv_index,
                                  recv_vec_v            );
    }

    /**
    * @brief specialization for std::vector<std::pair<Ta, Tb>>.
    * @param[in] send_vec_vp send target.
    *                         MUST NOT contain pointer member.
    * @param[out] recv_vec_pair recieve data.
    * @details    send_vec_vec_str[i] is send data to process rank i.
    * @details    devide std::vector<std::pair<Ta, Tb>> into std::vector<Ta> and std::vector<Tb>, then call scatterV() for std::vector<T>.
    * @details    class Ta and Tb accept std::vector<>, std::string, or user-defined class WITHOUT pointer member.
    */
    template <class Ta, class Tb>
    void scatter(const std::vector<std::vector<std::pair<Ta, Tb>>> &send_vec_vp,
                       std::vector<std::pair<Ta, Tb>>              &recv_vec_pair,
                 const PS::S32                                      root,
                       MPI_Comm                                     comm = MPI_COMM_WORLD){
        #ifdef DEBUG_COMM_TOOL
            if(PS::Comm::getRank() == 0 )  std::cout << " *** scatter_<std::vector<std::vector<std::pair<Ta, Tb>>>>" << std::endl;
        #endif

        const PS::S32 n_proc = PS::Comm::getNumberOfProc();
        assert(send_vec_vp.size() >= static_cast<size_t>(n_proc));
        assert(0 <= root && root < n_proc);

        std::vector<std::vector<Ta>> send_vec_1st;
        std::vector<std::vector<Tb>> send_vec_2nd;

        send_vec_1st.resize(n_proc);
        send_vec_2nd.resize(n_proc);
        for(PS::S32 i=0; i<n_proc; ++i){
            split_vector_pair(send_vec_vp[i],
                              send_vec_1st[i], send_vec_2nd[i]);
        }

        std::vector<Ta> recv_vec_1st;
        std::vector<Tb> recv_vec_2nd;
        scatter(send_vec_1st, recv_vec_1st, root, comm);
        scatter(send_vec_2nd, recv_vec_2nd, root, comm);

        combine_vector_pair(recv_vec_1st, recv_vec_2nd,
                            recv_vec_pair              );
    }


    /**
    * @brief specialization for returning type interface.
    * @param[in] send_vec send target.
    * @return    recv_vec recieved data.
    * @detail    wrapper for the use case: "auto recv_vec = COMM_TOOL::allToAll(vec);"
    * @detail    recv_vec.size() = PS::Comm::getNumberOfProc();
    * @detail    send/recv_vec[i] is send/recieve data to/from process i.
    */
    template <class T>
    T scatter(const std::vector<T> &send_vec,
              const PS::S32         root,
                    MPI_Comm        comm = MPI_COMM_WORLD){
        T recv_data;
        scatter(send_vec, recv_data, root, comm);
        return recv_data;
    }

}
