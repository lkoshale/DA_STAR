
#ifdef __NVCC__

#ifndef A_STAR_CUH
#define A_STAR_CUH

#include "dynamic_graph.cuh"
// #include "a_star_kernels.cuh"

#include "a_star.h"

template < class T, class U >
class GPU_A_Star
{
    private:
        GPU_Dynamic_Graph<T> *graph;
        GPU_Dynamic_Graph<T> rev_graph;
        U* Hx;
        U* Cx;

        unsigned int start_node;
        unsigned int end_node;

        unsigned int* PQ;
        unsigned int* PQ_size;
        unsigned int num_pq;
        int* parent;
        int* open_list;

        bool is_set_hx;

        //device pointers
        U* d_Hx;
        U* d_Cx;
        unsigned int* d_PQ;
        unsigned int* d_PQ_size;
        int* d_parent;
        int* d_open_list;

        void __alloc_gpu();

    public:
        GPU_A_Star(GPU_Dynamic_Graph<T> *graph, unsigned int start,unsigned int end, unsigned int K );

        void set_huiristics(U* hx);

        std::vector<int> get_path();

        void free_gpu();
                           

};

#include "a_star.cu"


#endif

#endif