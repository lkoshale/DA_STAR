/******************
* Geneartes CSR fromat graph from edge list
* creates directed graph
* currently reads from stdin
*
* usage: ./executable < input_edge_file.txt
* edge u->v of weight w1 
* edge v->y of weight w2
* input_edge_file.txt format:
* u v w1
* v y w2
*******************/

#include <iostream>
#include <vector>
#include <string>
#include <unordered_map>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stack>

using namespace std;

class Node{
public:
    unsigned int val;
    vector<unsigned int> weights;
    vector<Node*> Edges;
    vector<Node*> backEdge;
    vector<unsigned int> backWeight;

    Node(int val){
        this->val = val;
    }

    void addEdge(Node* v,unsigned int w){
        this->Edges.push_back(v);
        this->weights.push_back(w);
    }

};

void printGraphCSR(int N,int E,unordered_map<unsigned int,Node*>& Graph, string filename);

int main(){

    srand(42);

    unordered_map<unsigned int,Node*> Graph;

    unsigned int N = 0;
    unsigned int a,b,c;

    int E = 0;

    while( scanf("%u %u %u\n",&a,&b,&c)!=EOF){
        // c = rand()%100;
        
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(a);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            unordered_map<unsigned int,Node*>:: iterator it;
            it = Graph.find(b);
            if(it!=Graph.end()){
                Node* v = it->second;
                n->addEdge(v,c);
            }
            else{
                Node* v = new Node(b);
                n->addEdge(v,c);
                Graph.insert(pair<unsigned int,Node*>(b,v));
            }

        }
        else{
            Node* n =new Node(a);
            Graph.insert(pair<unsigned int,Node*>(a,n));

            unordered_map<unsigned int,Node*>:: iterator it;
            it = Graph.find(b);
            if(it!=Graph.end()){
                Node* v = it->second;
                n->addEdge(v,c);
            }
            else{
                Node* v = new Node(b);
                n->addEdge(v,c);
                Graph.insert(pair<unsigned int,Node*>(b,v));
            }

        }

        if( a > N ) N=a;
        if( b > N ) N=b;

        E++;
    }
    N++;

    printGraphCSR(N,E,Graph,"graph.txt");

    return 0;
}


void printGraphCSR(int N,int E,unordered_map<unsigned int,Node*>& Graph, string filename){
    FILE* fptr = fopen(filename.c_str(),"w");

    fprintf(fptr,"%d %d\n",N,E);
    
    vector<unsigned int> off;
    off.push_back(0);
    int k =0;
    for(int i=0;i<N;i++){
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(i);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            for(int j=0;j<n->Edges.size();j++){
                if(n->Edges[j]!=NULL){
                    fprintf(fptr,"%d ",n->Edges[j]->val);
                    k++;
                }
            }
            off.push_back(k);
        }
        else{
            off.push_back(k);
        }
    }

    fprintf(fptr,"\n");

    for(int j=0;j<off.size()-1;j++){
        fprintf(fptr,"%d ",off[j]);
    }

    fprintf(fptr,"\n");
    

    for(int i=0;i<N;i++){
        unordered_map<unsigned int,Node*>:: iterator itr;
        itr = Graph.find(i);
        if(itr!=Graph.end()){
            Node* n = itr->second;
            for(int j=0;j<n->weights.size();j++){
                if(n->Edges[j]!=NULL){
                    fprintf(fptr,"%d ",n->weights[j]);
                }
            }
        }
    }
    
    fprintf(fptr,"\n");
    fclose(fptr);

}