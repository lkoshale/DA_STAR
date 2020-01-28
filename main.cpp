#include <iostream>
#include <cstdlib>
#include <string>

#include "include/graph.h"
#include "include/diff_graph.h"
#include "include/dynamic_graph.h"
#include "include/a_star.h"

#ifdef __NVCC__
    #include <cuda.h>
    #include "include/graph.cuh"
    #include "include/a_star.cuh"
#endif


int main(){
    Graph<unsigned int > graph;
    std::string filename = "graph.txt";
    graph.read_graph(filename);

    Dynamic_Graph<unsigned int> dgraph(&graph);
    
   
    int* edge = dgraph.get_graph().get_edges();
    int* off = dgraph.get_graph().get_offsets();
    unsigned int* w = dgraph.get_graph().get_weights();
 
    
    for(int i=0;i<graph.get_num_edges();i++){
        std::cout<<edge[i]<<" ";
    }
    std::cout<<"\n";

    for(int i=0;i<graph.get_num_nodes();i++){
        std::cout<<off[i]<<" ";
    }
    std::cout<<"\n";

    for(int i=0;i<graph.get_num_edges();i++){
        std::cout<<w[i]<<" ";
    }
    std::cout<<"\n";

    int start= 0;
    int end = 7;

    int* hx = (int*)malloc(sizeof(int)*graph.get_num_nodes());
    for(int i=0;i<graph.get_num_nodes();i++)
        hx[i] = 0;

    A_Star <unsigned int,int>* algo = new A_Star<unsigned int,int>(&dgraph,start,end,1);

    algo->set_huiristics(hx);

    std::vector<int> path = algo->get_path();

    for(int i=0;i<path.size();i++){
        printf("%d ",path[i]);
    }
    printf("\n");
    
    path.clear();

    #ifdef __NVCC__
        GPU_Graph<unsigned int> g(&graph);
        printf("allocated\n");
        GPU_Dynamic_Graph<unsigned int> dg(&g);
        GPU_A_Star < unsigned int, int > g_algo(&dg,start,end,1);
        g_algo.set_huiristics(hx);
        
        path = g_algo.get_path();
        
        for(int i=0;i<path.size();i++){
             printf("%d ",path[i]);
        }
        printf("\n");
    #endif



    return 0;
}