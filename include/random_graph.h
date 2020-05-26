#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <vector>

void gen_random_graph(){

    int N,E;
    printf("Enter Number of Nodes: ");
    scanf("%d",&N);
    printf("Enter number of Edges: ");
    scanf("%d",&E);

    FILE* fptr = fopen("random_graph.txt","w");

    for(int i=0;i<E;i++){
        int u = rand()%N;
        int v = rand()%N;
        int w = rand()%90+10;

        fprintf(fptr,"%d %d %d\n",u,v,w);
    }

    fclose(fptr);

}

void gen_csr(std::string filename){
    
}

