#ifndef UTILS_H
#define UTILS_H

#include "graph_builder.h"
#include <utility>
#include <cstdio>
#include <cstdlib>






/*
void build_fn( std::unordered_map<unsigned int, Node<unsigned int>* > *graph,unsigned int *N,unsigned int* E){
    int n =0;
    int e = 0;

    std::unordered_map<unsigned int, Node<unsigned int>* > Graph = *graph;

    FILE* fptr = fopen("pair.txt","r");

    unsigned int a,b,c;
    while( fscanf(fptr,"%u %u %u\n",&a,&b,&c)!=EOF){

        std::unordered_map< unsigned int,Node<unsigned int>* >:: iterator itr;
        itr = Graph.find(a);
        if(itr!=Graph.end()){
            Node<unsigned int>* n = itr->second;
            std::unordered_map< unsigned int,Node<unsigned int>* >:: iterator it;
            it = Graph.find(b);
            if(it!=Graph.end()){
                Node<unsigned int>* v = it->second;
                n->addEdge(v,c);
            }
            else{
                Node<unsigned int>* v = new Node<unsigned int>(b);
                n->addEdge(v,c);
                Graph.insert(std::pair<unsigned int,Node<unsigned int>*>(b,v));
            }

        }
        else{
            Node <unsigned int>* n =new Node<unsigned int>(a);
            
            Graph.insert(std::pair<unsigned int,Node<unsigned int>*>(a,n));

            std::unordered_map<unsigned int,Node<unsigned int>*>:: iterator it;

            it = Graph.find(b);
            if(it!=Graph.end()){
                Node<unsigned int>* v = it->second;
                n->addEdge(v,c);
            }
            else{
                Node<unsigned int>* v = new Node<unsigned int>(b);
                n->addEdge(v,c);
                Graph.insert(std::pair<unsigned int,Node<unsigned int>*>(b,v));
            }
        }

        if( a > n ) n=a;
        if( b > n ) n=b;

        e++;
    
    }
    n++;
    
    *N = n;
    *E = e;

    return;
}

*/



#endif