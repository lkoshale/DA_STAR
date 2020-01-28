#ifndef A_STAR_H
#define A_STAR_H

#include <climits>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <algorithm>

#include "dynamic_graph.h"


template <class T, class U >
class A_Star
{
    private:
        Dynamic_Graph<T> *graph;
        Dynamic_Graph<T> rev_graph;
        U* Hx;
        U* Cx;

        unsigned int start_node;
        unsigned int end_node;

        unsigned int* PQ;
        unsigned int* PQ_size;
        unsigned int num_pq;
        int* parent;
        int* open_list;

        bool is_set_hx;
    
    public:
        A_Star(Dynamic_Graph<T> *graph, unsigned int start,unsigned int end, unsigned int K );

        void set_huiristics(U* hx);


        void insert_pq(int node);

        int get_min_pq();

        void update_pq(int node);

        std::vector<int> get_path();

        std::vector<int> get_path_gpu();

        
    
};


template <class T, class U >
A_Star< T, U> :: A_Star(Dynamic_Graph<T> *graph, unsigned int start_node,unsigned int end_node, unsigned int K )
{
    this->graph = graph;
    this->num_pq = K;
    this->start_node = start_node;
    this->end_node = end_node;
    
    int N = this->graph->get_graph().get_num_nodes();

    this->PQ = (unsigned int*)malloc(sizeof(unsigned int)*N );
    this->PQ_size = (unsigned int*)malloc(sizeof(unsigned int)*K);

    this->Cx = (U*)malloc(sizeof(U)*N);
    this->Hx = (U*)malloc(sizeof(U)*N);
    this->parent = (int*)malloc(sizeof(int)*N);
    this->open_list = (int*)malloc(sizeof(int)*N);

    memset(this->parent,-1,sizeof(int)*N);
    memset(this->open_list,-1,sizeof(int)*N);

    is_set_hx = false;
    *(this->PQ_size)=0;
    
    //todo make it memset
    for(int i=0;i<N;i++){
        this->Cx[i] = INT_MAX;
    }

}

template <class T, class U >
void A_Star< T, U> :: set_huiristics(U* hx)
{
    this->Hx = hx;
    is_set_hx = true;
}


template <class T, class U >
void A_Star< T, U>  :: insert_pq(int node)
{

    if(*PQ_size == 0){
        PQ[*PQ_size]=node;
        (*PQ_size)++;
    }
    else{
        PQ[*PQ_size]=node;
        (*PQ_size)++;

        int i = *PQ_size - 1;
        while(i>0){
            if(Cx[PQ[(i-1)/2]] > Cx[PQ[i]] ){
                int swap = PQ[i];
                PQ[i] = PQ[(i-1)/2];
                PQ[(i-1)/2] = swap;
                i = (i-1)/2;
            }
            else
                break;
        }
    }
}

template <class T, class U >
int A_Star< T, U> :: get_min_pq()
{
    if(*PQ_size==0)
        return -1;
    
    int min = PQ[0];

    //swap with last
    PQ[0]= PQ[*PQ_size -1];
    (*PQ_size)--;

    int i=0;
    while(2*i+1< *PQ_size){
        if(2*i+2 >= *PQ_size){
            if(Cx[PQ[2*i+1]] < Cx[PQ[i]]){
                int swap = PQ[2*i+1];
                PQ[2*i+1]=PQ[i];
                PQ[i] = swap;
                i = 2*i +1;
            }
            else 
                break;
        }
        else{
            if( Cx[PQ[i]] > Cx[PQ[2*i+1]] && Cx[PQ[2*i+1]]<= Cx[PQ[2*i+2]]){
                int swap = PQ[2*i+1];
                PQ[2*i+1]=PQ[i];
                PQ[i] = swap;
                i = 2*i +1;
            }
            else if(Cx[PQ[i]] > Cx[PQ[2*i+2]] && Cx[PQ[2*i+2]]<= Cx[PQ[2*i+1]]){
                int swap = PQ[2*i+2];
                PQ[2*i+2]=PQ[i];
                PQ[i] = swap;
                i = 2*i +2;
            }
            else
                break;
        }

    }

    return min;
}


template <class T, class U >
void A_Star< T,  U>:: update_pq(int node)             //update happens when Cx is decreased
{
    
    int index=-1;
    for(int i=0;i<*PQ_size;i++){
        if(PQ[i]==node)
            index = i;
    }

    if(index>0){
        int i = index;
        while(i>0){
            if(Cx[PQ[(i-1)/2]] > Cx[PQ[i]] ){
                int swap = PQ[i];
                PQ[i] = PQ[(i-1)/2];
                PQ[(i-1)/2] = swap;
                i = (i-1)/2;
            }
            else
                break;
        }
    }
}

template <class T, class U >
std::vector<int> A_Star< T,  U>::get_path()
{

    std::vector<int> path;

   // Hx not read
    if(!is_set_hx){
        path.push_back(-1);
        return path;
    }

    int N = this->graph->get_graph().get_num_nodes();
    int E = this->graph->get_graph().get_num_edges();

    int* offsets = this->graph->get_graph().get_offsets();
    int* edges = this->graph->get_graph().get_edges(); 
    T* weights = this->graph->get_graph().get_weights();

    

    Cx[this->start_node]=Hx[this->start_node];

    insert_pq(this->start_node);

    open_list[this->start_node]=1;
    
    bool flag_found = false;
    int curr_node ;

    while( *PQ_size >0 && !flag_found){

        curr_node = get_min_pq();
        //remove from openList
        open_list[curr_node]=-1;

        //got dest
        if(curr_node==end_node){
            flag_found = true;
            break;
        }

        //iterate 
        int begin = offsets[curr_node];
        int end = E-1;
        if(curr_node!=N-1)
            end = offsets[curr_node+1];
        
        for(int j=begin;j<end;j++){
            int child = edges[j];
            if( Cx[child] > Cx[curr_node]- Hx[curr_node]+ weights[j]+ Hx[child]){
                Cx[child] = Cx[curr_node] - Hx[curr_node] + weights[j]+Hx[child];
                if(open_list[child]==1){
                    update_pq(child);
                }
                else{
                    insert_pq(child);
                    open_list[child]=1;
                }
                parent[child]=curr_node;
            }
        }
    }

    

    if(flag_found){
        // printf("cost %d\n",H_cx[endNode]);

        int c_node = end_node;
        path.push_back(c_node);

        //printf("%d ",c_node);
        while(parent[c_node]!=-1){
           // printf("%d ",H_parent[c_node]);
            c_node = parent[c_node];
            path.push_back(c_node);
        }
        // printf("\n");
    }

    std::reverse(path.begin(),path.end());

    return path;
}



#endif