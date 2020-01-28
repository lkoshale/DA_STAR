
#ifndef GRAPH_H
#define GRAPH_H


#include <unordered_map>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>


template <class DataType > 
class Graph 
{
    public:
        unsigned int N;         //nodes from 0 to N-1
        unsigned int E;
        int* offsets;
        int* edges;
        DataType* weights;

        std::unordered_map<long,int> vertex_map;     //if vertex doesn't follow 0-1 naming
 
        Graph();

        void __alloc__(int n, int e);

        //void build_graph( FuncType function );

        int get_num_nodes();
        int get_num_edges();

        int* get_offsets();
        int* get_edges();
        DataType* get_weights();

        void read_graph(std::string filename);

        void free();

        
};



/*
* Implementation Here as 
* templates cant have diff files for implementation
* (or move to another file .tpp nand include it here)
*/

template <class DataType >
Graph<DataType >::Graph()
{
    N = 0;
    E = 0;
}

template <class DataType >
void Graph<DataType > :: __alloc__(int n, int e){
    N = n;
    E = e;
    this->offsets = (int*)malloc(sizeof(int)*N);
    this->edges = (int*)malloc(sizeof(int)*E);
    this->weights = (DataType*)malloc(sizeof(DataType)*E);
}


// template <class DataType >
// void Graph<DataType > :: build_graph(FuncType builder )
// {
//     builder(this->offsets,this->edges,this->weights);
// }

template <class DataType >
int Graph<DataType > :: get_num_nodes(){ return N; }

template <class DataType >
int Graph<DataType > :: get_num_edges(){ return E; }

template <class DataType >
int* Graph<DataType > :: get_offsets(){ return this->offsets; }

template <class DataType >
int*  Graph<DataType > :: get_edges(){ return this->edges; }

template <class DataType >
DataType*  Graph<DataType > :: get_weights(){ return this->weights; }


template <class DataType >
void Graph<DataType > :: read_graph(std::string filename){
    std:: ifstream infile;
    infile.open(filename);
    if (!infile){
        std::cout<<"[ERROR] Couldn't open graph file\n";
        exit(1);
    } 

    int n,e;
    infile >> n >> e;

    this->__alloc__(n,e);
    for(int i=0;i<e;i++){
       infile >> this->edges[i];
    }

    for(int i=0;i<n;i++){
       infile >> this->offsets[i];
    }

    for(int i=0;i<e;i++){
       infile >> this->weights[i];
    }


}


template <class DataType >
void  Graph<DataType > :: free()
{
    free(this->offsets);
    free(this->edges);
    free(this->weights);
}


#endif