# DA_STAR

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![build:passing](https://img.shields.io/badge/build-passing-brightgreen.svg)](https://github.com/lkoshale/DDP)
<!-- [![forthebadge](https://forthebadge.com/images/badges/made-with-c-plus-plus.svg)](https://forthebadge.com) -->

This repository contains the generic implementation of the A*(star) algorithm for dynamic graphs.

## About
A\* is one of the widely used path planning and shortest path approximation algorithms. It is applied in a diverse set of problems from path-finding in video games and robotics to codon optimization in genetics. In this work, we focus on A\* for graphs that are subjected to update operation, where an update operation refers to the insertion or
deletion of an edge. Instead of performing A\* again from start each time graph is subject to update, our algorithm process the sub-graphs which are affected by the update. For temporal graph available at [SNAP](http://snap.stanford.edu/data.)
for 100 updates we got 25x-40x of performance improvement over repeated A* search. [More details](#dynamic-graphs)


## Features
- Supports both directed and undirected positively weighted graphs 
- Supports all types of heuristics function 
- Supports insertions/deletions of edges onto the graph
- Uses Diff CSR to store dynamic graphs
- improved performance from repeated A* search

## System Requirements

To use this library ensure you meet the following requirements:

**Linux**

* [gcc and g++](https://gcc.gnu.org/)  with std c++11

* [CUDA](https://developer.nvidia.com/cuda-toolkit) >= 10.0

## How to Build
This library is self-contained in the header files at the [include](https://github.com/lkoshale/DA_STAR/tree/master/include) directory of this repository with all required classes and functions defined.

- You can copy the headers of this library to your projects include directory :
    - ```cp -r DA_STAR/include/   <your-project>/include/```

- Or You can copy the headers of this library to your root/CUDA include folder

    - ```sudo cp -r DA_STAR/include/   /usr/include/```
    - ```cp -r DA_STAR/include/     <CUDA_PATH>/include/```

## Examples
```c
#include "include/dynamic_graph.cuh"
#include "include/d_a_star.cuh"
#include "include/utils.cuh"

int main(){

    std::string input_graph = "graph.txt";
    std::string input_updates = "updates.txt";
    
    //reuires weight type
    GPU_Dynamic_Graph<unsigned int> graph;

    graph.read(input_graph);

    int start_node = 0;             //start node
    int end_node = 7;               //end node
    int num_pq = 10;                //number of parallel priority queue
    
    //requires weight type and heuristics/Cost type
    GPU_D_A_Star<unsigned int, int> algo(&graph,start_node,end_node,num_pq);

    graph.set_update_file(input_update);

    while(graph.update()){

        std::vector<int> path = algo.get_path();
        print_path(path);
        
    }

    return 0;
}

```


## Dynamic Graphs


## Credits
