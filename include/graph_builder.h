#ifndef GRAPH_BUILDER_H
#define GRAPH_BUILDER_H

#include <vector>
#include <unordered_map>
#include <iostream>

template < class T >
class Node
{
    public:
        unsigned int val;
        std::vector< Node<T>* > Edges;
        std::vector<T>weights;
    
         Node(int val){
            this->val = val;
         }

        void addEdge(Node<T> *v,T w){
            this->Edges.push_back(v);
            this->weights.push_back(w);
        }
};

template < class T >
class Builder
{
    private:
        std::unordered_map<unsigned int, Node<T>* > Graph;
        unsigned int N;
        unsigned int E;
    public:
        Builder(){
            N = 0;
            E = 0;
        }

        void build_graph( void (*build_fn)(std::unordered_map<unsigned int, Node<T>* > *Graph,unsigned int *N,unsigned int* E)){
            build_fn(&(this->Graph),&N,&E);
        }

        void operator()(int* offset,int* edges, T* weights){
            std:: vector<unsigned int> off;
            off.push_back(0);
            int k =0;
            std::cout<<Graph.size()<<" \n";
            for(int i=0;i<N;i++){
                typename std:: unordered_map<unsigned int,Node<T>* >:: iterator itr;
                itr = Graph.find(i); 
                if(itr!=Graph.end()){
                    Node<T> *n = itr->second;
                    for(int j=0;j<n->Edges.size();j++){
                        if(n->Edges[j]!=NULL){
                            edges[k] = n->Edges[j]->val ;
                            weights[k] = n->weights[j];
                            k++;
                        }
                    }
                    off.push_back(k);
                }
                else{
                    off.push_back(k);
                    std::cout<<k<<" "<<i<<"\n";
                }
            }

            
            for(int j=0;j<off.size()-1;j++){
                offset[j] = off[j];
            }
        }

        unsigned int get_num_nodes(){return N;}
        unsigned int get_num_edges(){return E;}

};




#endif