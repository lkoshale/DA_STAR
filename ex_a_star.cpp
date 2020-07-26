#include <iostream>
#include <cstdlib>
#include <string>

#include "include/graph.h"
#include "include/diff_graph.h"
#include "include/dynamic_graph.h"
#include "include/a_star.h"

#include <cuda.h>
#include "include/graph.cuh"
#include "include/a_star.cuh"

/*****************
 *  graph.txt Format
 *  N E
 *  edge array 
 *  offset array
 *  weight array
 * 
 * ex: adjacency list, all edge wieghts= 1
 * 0 -> 1, 2, 3
 * 1 -> 0, 2, 3
 * 2 -> 0, 1, 3
 * 3 -> 0, 1, 2
 * 
 * graph.txt
 * 4, 12
 * 1 2 3 0 2 3 0 1 3 0 1 2
 * 0 3 6 9
 * 1 1 1 1 1 1 1 1 1 1 1 1
 * **************/

int main(){

    // create CPU graph from file
    Graph<unsigned int > graph;
    std::string filename = "graph.txt";
    graph.read_graph(filename);

    // read or set heuristics value for each node
    int* hx = (int*)malloc(sizeof(int)*graph.get_num_nodes());
    for(int i=0;i<graph.get_num_nodes();i++)
        hx[i] = 0;

    // GPU graph from cpu graph
    GPU_Graph<unsigned int> gpu_graph(&graph);

    // dynamic GPU graph object from graph
    GPU_Dynamic_Graph<unsigned int> dgpu_graph(&gpu_graph);

    // sorce , destination, and number of parallel priority queues
    int start= 0;
    int end = 7;
    int parallel_pqs = 2;

    // create algorithm object
    GPU_A_Star < unsigned int, int > gpu_algo(&dgpu_graph,start,end,parallel_pqs);

    // set heruristics value
    gpu_algo.set_heuristics(hx);
    
    // find path on GPU
    std::vector<int> path = gpu_algo.get_path();
    
    printf("Path:");
    for(int i=0;i<path.size();i++){
            printf("%d ",path[i]);
    }
    printf("\n");

    return 0;
}