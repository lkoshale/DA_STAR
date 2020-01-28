
#ifdef __NVCC__

#ifndef DIFF_GRAPH_CUH
#define DIFF_GRAPH_CUH

#include "graph.cuh"

template < class T = unsigned int>
class GPU_Diff_Graph : GPU_Graph<T>
{
    private:
        int id;
    
    public:
        GPU_Diff_Graph(int id);
};


/*************************
*       Diff Graph
*
**************************/
template< class T >
GPU_Diff_Graph<T> :: GPU_Diff_Graph(int id){
    this->id = id;
    this->N = 0;
    this->E = 0;
}


#endif

#endif