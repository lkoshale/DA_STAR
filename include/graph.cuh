#ifdef __NVCC__

#ifndef GRAPH_CUH
#define GRAPH_CUH

#include <cuda.h>

#include "utils.cuh"
#include "graph.h"


template < class T >
class GPU_Graph : public Graph<T>
{
    private:
        int* d_offsets;
        int* d_edges;
        T* d_weights;
    
    public:

        GPU_Graph();
        GPU_Graph(Graph<T>* graph_cpu);
         
        void __gpu_alloc();

        void __gpu_free();

        int* get_offsets();
        int* get_edges();
        T* get_weight();

};


/*
* Implemenattion
* Here
*/


template< class T >
GPU_Graph<T>:: GPU_Graph()
{
    this->N = 0;
    this->E = 0;
}

template< class T >
GPU_Graph<T>:: GPU_Graph(Graph<T>* graph_cpu)
{
    this->N = graph_cpu->get_num_nodes();
    this->E = graph_cpu->get_num_edges();
    this->offsets = graph_cpu->get_offsets();
    this->edges = graph_cpu->get_edges();
    this->weights = graph_cpu->get_weights();

    __gpu_alloc();
}

template< class T >
void GPU_Graph<T>:: __gpu_alloc()
{
    gpuErrchk ( cudaMalloc(&d_offsets,sizeof(int)*(this->N) ) );
    gpuErrchk ( cudaMalloc(&d_edges,sizeof(int)*(this->E) ) );
    gpuErrchk ( cudaMalloc(&d_weights,sizeof(T)*(this->E) ) );

    gpuErrchk ( cudaMemcpy(d_offsets,this->offsets,sizeof(int)*(this->N),cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(d_edges,this->edges,sizeof(int)*(this->E),cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(d_weights,this->weights,sizeof(T)*(this->E),cudaMemcpyHostToDevice) );

}

template< class T >
int* GPU_Graph<T>:: get_offsets(){ return this->d_offsets; }

template< class T >
int* GPU_Graph<T>:: get_edges(){ return this->d_edges; }

template< class T >
T* GPU_Graph<T>:: get_weight(){ return this->d_weights; }


template< class T >
void  GPU_Graph<T>:: __gpu_free(){
    gpuErrchk ( cudaFree(d_offsets) );
    gpuErrchk ( cudaFree(d_edges) );
    gpuErrchk ( cudaFree(d_weights) );
}



#endif

#endif
