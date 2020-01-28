#ifndef DIFF_GRAPH_H
#define DIFF_GRAPH_H

#include "graph.h"

template <class T = unsigned int>
class Diff_Graph : public Graph<T>
{
    private:
        int id;             //diff graph id
    
    public:
        Diff_Graph();
        Diff_Graph(int i);
        
        ~Diff_Graph();
};


/*
* Implementation
*/

template < class T >
Diff_Graph<T> :: Diff_Graph(){}


template <class T >
Diff_Graph<T> :: Diff_Graph(int id)
{
    this->id = id;
}

#endif
