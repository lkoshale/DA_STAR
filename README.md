# DA_STAR

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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


## Dynamic Graphs


## Credits