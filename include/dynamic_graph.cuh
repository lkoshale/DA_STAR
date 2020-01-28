#ifdef __NVCC__

#ifndef DYNAMIC_GRAPH_CUH
#define DYNAMIC_GRAPH_CUH

#include "diff_graph.cuh"

template < class T >
class GPU_Dynamic_Graph
{
    private:
        GPU_Graph<T> *graph;
        GPU_Diff_Graph<T> *diff;
        bool is_diff = false;

    public:
        GPU_Dynamic_Graph();

        GPU_Dynamic_Graph(GPU_Graph<T> * graph, GPU_Diff_Graph<T> * diff);

        GPU_Dynamic_Graph(GPU_Graph<T> * graph);

        GPU_Graph<T> get_graph();
        GPU_Diff_Graph<T> get_diff_graph();

};



template < class T >
GPU_Dynamic_Graph<T> :: GPU_Dynamic_Graph()
{
    this->graph = new GPU_Graph<T>();
    this->diff = new GPU_Diff_Graph<T>(0);
    is_diff = false;
}


template < class T >
GPU_Dynamic_Graph<T> :: GPU_Dynamic_Graph(GPU_Graph<T> * graph, GPU_Diff_Graph<T> * diff)
{
    this->graph = graph;
    this->diff = diff;
    is_diff = false;
}

template < class T >
GPU_Dynamic_Graph<T> :: GPU_Dynamic_Graph(GPU_Graph<T> * graph)
{
    this->graph = graph;
    is_diff = false;
}

template < class T >
GPU_Graph<T> GPU_Dynamic_Graph<T> ::  get_graph() { return *(this->graph); }

template < class T >
GPU_Diff_Graph<T> GPU_Dynamic_Graph<T> :: get_diff_graph(){ return *(this->diff); }


#endif

#endif