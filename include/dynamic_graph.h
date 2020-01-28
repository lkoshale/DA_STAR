
#ifndef DYNAMIC_GRAPH_H
#define DYNAMIC_GRAPH_H

#include "diff_graph.h"

template <class T = unsigned int>
class Dynamic_Graph
{
    private:
        Graph<T> *graph;
        Diff_Graph<T> *diff;
        bool is_diff = false;
    
    public:
        Dynamic_Graph();

        Dynamic_Graph(Graph<T> *graph, Diff_Graph<T> * diff);

        Dynamic_Graph(Graph<T> *graph);

        Graph<T> get_graph();
        Diff_Graph<T> get_diff_graph();

        void merge();
        void make_diff();
        void update_Del();

};


/*
* Implemetation
*/

template < class T > 
Dynamic_Graph<T> :: Dynamic_Graph()
{
    this->graph = new Graph<T>();
    this->diff = new Diff_Graph<T>(0);
    is_diff = true;
}

template < class T >
Dynamic_Graph<T> :: Dynamic_Graph(Graph<T> *graph)
{
    this->graph = graph;
    is_diff = false;
}

template < class T >
Dynamic_Graph<T> :: Dynamic_Graph(Graph<T> *graph, Diff_Graph<T> * diff)
{
    this->graph = graph;
    this->diff = diff;
    is_diff = true;
}

template < class T >
Graph<T> Dynamic_Graph<T> :: get_graph(){return *(this->graph); }

template < class T >
Diff_Graph<T> Dynamic_Graph<T> :: get_diff_graph(){ return *(this->diff); }



#endif