#ifdef __NVCC__


template < class U >
__global__ void extractMin(unsigned int* PQ, unsigned int* PQ_size, int* expandNodes,int* expandNodes_size,U* Cx,int* openList,int N,int K);


template < class T, class U >
__global__ void A_star_expand(int* off,int* edge,T* W, U* Hx,int* parent,volatile U* Cx,
    int* expandNodes,int* expandNodes_size, int* lock ,int* flagfound,int* openList,
    int N,int E, int K,int dest,int* nVFlag );


template < class U >
__global__ void keepHeapPQ(unsigned int* PQ, unsigned int* PQ_size,U* Cx,int N,int K);


__global__ void setNV(int* nextFlag,int* nextV,int* nvSize,int N);

template <class U >
__global__ void insertPQ(unsigned int* PQ,unsigned int* PQS,int* nextV,int* nVsize,U* Cx,int K,int N,int* openList);


template < class U >
__global__ void checkMIN(unsigned int* PQ, unsigned int* PQ_size,int* flagEnd,U* Cx,int dest,int N,int K);

template <class U>
__global__ void getCx(U* Cx,int dest,U* val);


#include "kernels/a_star_kernels.cu"

#include "a_star.cuh"

#ifdef DEBUG
    #include <cstdio>
#endif


template <class T, class U >
GPU_A_Star< T, U> :: GPU_A_Star(GPU_Dynamic_Graph<T> *graph, unsigned int start_node,unsigned int end_node, unsigned int K )
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
    memset(this->PQ_size,0,sizeof(int)*K);

    is_set_hx = false;
    
    //todo make it memset
    for(int i=0;i<N;i++){
        this->Cx[i] = INT_MAX;
    }

}

template <class T, class U >
void GPU_A_Star< T, U> :: __alloc_gpu()
{
    int N = this->graph->get_graph().get_num_nodes();

    gpuErrchk ( cudaMalloc(&d_Cx,sizeof(U)*N ) );

    gpuErrchk ( cudaMalloc(&d_parent,sizeof(int)*N ) );
    gpuErrchk ( cudaMalloc(&d_open_list,sizeof(int)*N ) );

    gpuErrchk ( cudaMalloc(&d_PQ,sizeof(unsigned int)*N ) );
    gpuErrchk ( cudaMalloc(&d_PQ_size,sizeof(unsigned int)*num_pq ) );

    gpuErrchk ( cudaMemcpy(d_Cx,Cx,sizeof(U)*N,cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpy(d_PQ_size,PQ_size,sizeof(unsigned int)*num_pq,cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpy(d_parent,parent,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(d_open_list,open_list,sizeof(int)*N,cudaMemcpyHostToDevice) );

}



template <class T, class U >
void GPU_A_Star< T, U> :: set_huiristics(U* hx)
{

    this->Hx = hx;
    is_set_hx = true;
    int N = this->graph->get_graph().get_num_nodes();
    gpuErrchk ( cudaMalloc(&d_Hx,sizeof(U)*N ) );
    gpuErrchk ( cudaMemcpy(d_Hx,Hx,sizeof(U)*N,cudaMemcpyHostToDevice) );

}


template <class T, class U >
std::vector<int>  GPU_A_Star< T, U>:: get_path()
{

    int N = this->graph->get_graph().get_num_nodes();
    int E = this->graph->get_graph().get_num_edges();
    int K = this->num_pq;


    //init Host var
    int* flag_end = (int*)malloc(sizeof(int));
    int* flag_found = (int*)malloc(sizeof(int));
    int* __a0 = (int*)malloc(sizeof(int));
    *__a0 = 0;



    //required coz if many tries to add same in diff threads high low lower
    int* next_vertices_flag = (int*)malloc(sizeof(int)*N);
    memset(next_vertices_flag,-1,sizeof(int)*N);    

    *flag_end = 0;
    *flag_found = 0;

    //insert startNode in PQ[0]
    Cx[this->start_node] = Hx[this->start_node];
    PQ[0] = this->start_node;
    PQ_size[0]=1;
    open_list[this->start_node]=0;  

    //alloc
    __alloc_gpu();

    //next nodes flag
    int* d_next_vertices_flag;
    //next nodes array to insert PQ
    int* d_next_vertices;
    int* d_next_vertices_size;
    
    //nodes to be expanded ( extracted from PQ )
    int* d_expand_nodes;
    int* d_expand_nodes_size;
    
    //flag to end while loop and found the destination
    int* d_flag_end;
    int* d_flag_found;

    //cost of endNode
    U* d_dest_cost;

    //lock for nodes
    int* d_lock;


    gpuErrchk ( cudaMalloc(&d_lock,sizeof(int)*N) );

    gpuErrchk ( cudaMalloc(&d_dest_cost,sizeof(U)) );

    
    //for next set of vertices to add in PQ
    gpuErrchk ( cudaMalloc(&d_next_vertices,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&d_next_vertices_size,sizeof(int)) );
    gpuErrchk ( cudaMalloc(&d_next_vertices_flag,sizeof(int)*N) );

    //next nodes to expand
    gpuErrchk ( cudaMalloc(&d_expand_nodes,sizeof(int)*K) );  //changed to K
    gpuErrchk ( cudaMalloc(&d_expand_nodes_size,sizeof(int)) );

    //flag to end search
    gpuErrchk( cudaMalloc(&d_flag_end,sizeof(int)) );
    gpuErrchk( cudaMalloc(&d_flag_found,sizeof(int)) );
    

    gpuErrchk ( cudaMemcpy(d_flag_end,flag_end,sizeof(int),cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(d_flag_found,flag_found,sizeof(int),cudaMemcpyHostToDevice) );
   
    gpuErrchk ( cudaMemcpy(d_next_vertices_flag,next_vertices_flag,sizeof(int)*N,cudaMemcpyHostToDevice) );

    // gpuErrchk ( cudaMemcpy(d_next_vertices_size,__a0,sizeof(int),cudaMemcpyHostToDevice) );
    // gpuErrchk ( cudaMemcpy(d_expand_nodes_size,__a0,sizeof(int),cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemset(d_next_vertices_size,0,sizeof(int)) );
    gpuErrchk ( cudaMemset(d_expand_nodes_size,0,sizeof(int)) );
    
    gpuErrchk ( cudaMemset(d_lock,0,sizeof(int)*N) );

    int flag_PQ_not_empty = 0;
    for(int i=0;i<K;i++){
        if(PQ_size[i]>0)
            flag_PQ_not_empty=1;
    }


    int numThreads = 512;
    int numBlocks = (K+numThreads-1)/numThreads;
    int N_numBlocks = (N+numThreads-1)/numThreads;

    #ifdef DEBUG
        printf("[INFO] A* started\n");
    #endif
   
    //DO A* initailly on whole graph
    while(*flag_end==0 && flag_PQ_not_empty==1){
        
        //extract min
        extractMin <U> <<<numBlocks,numThreads>>>( d_PQ, d_PQ_size, d_expand_nodes, d_expand_nodes_size, d_Cx, d_open_list, N, K );
        
        gpuErrchk(cudaPeekAtLastError() );

        cudaDeviceSynchronize();

        
        A_star_expand < T, U > <<<numBlocks,numThreads>>> (
            this->graph->get_graph().get_offsets(),this->graph->get_graph().get_edges(),this->graph->get_graph().get_weight(),
            d_Hx,d_parent,d_Cx,
            d_expand_nodes,d_expand_nodes_size, d_lock ,d_flag_found,d_open_list,
            N,E,K,this->end_node,d_next_vertices_flag 
        );
        
        gpuErrchk(cudaPeekAtLastError() );
        
        cudaDeviceSynchronize();


        keepHeapPQ < U > <<<numBlocks,numThreads>>>( d_PQ, d_PQ_size, d_Cx, N, K );
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        
        //gen from flag d_next_vertices
        //for N in parallel
        setNV<<<N_numBlocks,numThreads>>>(d_next_vertices_flag, d_next_vertices, d_next_vertices_size, N );
        
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        

        insertPQ < U > <<<numBlocks,numThreads>>>( d_PQ, d_PQ_size, d_next_vertices, d_next_vertices_size, d_Cx, K, N, d_open_list );
        
        gpuErrchk(cudaPeekAtLastError() );
        cudaDeviceSynchronize();
       
        //cpy flagend and flagEmpty
        gpuErrchk( cudaMemcpy(flag_found,d_flag_found, sizeof(int),cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(PQ_size,d_PQ_size, sizeof(int)*K,cudaMemcpyDeviceToHost) );
        
        //reset nVFlag
        gpuErrchk( cudaMemcpy(d_next_vertices_flag,next_vertices_flag,sizeof(int)*N,cudaMemcpyHostToDevice) );

        //reset next insert array
        gpuErrchk( cudaMemcpy(d_next_vertices_size, __a0,sizeof(int),cudaMemcpyHostToDevice) );
        gpuErrchk( cudaMemcpy(d_expand_nodes_size, __a0,sizeof(int),cudaMemcpyHostToDevice) );
        

        flag_PQ_not_empty = 0;
        for(int i=0;i<K;i++){
            if(PQ_size[i]>0)
                flag_PQ_not_empty=1;
        }

        //check for mins
        if( *flag_found==1 && flag_PQ_not_empty==1){
            //end 
            gpuErrchk( cudaMemcpy(d_flag_end,flag_found,sizeof(int),cudaMemcpyHostToDevice) );

            checkMIN < U > <<< numBlocks,numThreads >>>(d_PQ, d_PQ_size, d_flag_end, d_Cx, this->end_node, N, K );
            
            gpuErrchk(cudaPeekAtLastError() );
            cudaDeviceSynchronize();
            gpuErrchk( cudaMemcpy(flag_end,d_flag_end, sizeof(int),cudaMemcpyDeviceToHost) );
        }
   
    }

    getCx < U > <<<1,1>>>( d_Cx, this->end_node,d_dest_cost);

    U dest_cost;
    gpuErrchk( cudaMemcpy(&dest_cost,d_dest_cost, sizeof(U),cudaMemcpyDeviceToHost) );
 
    gpuErrchk( cudaMemcpy(parent,d_parent, sizeof(int)*N,cudaMemcpyDeviceToHost) );


    std::vector<int> Path;
    if(dest_cost != INT_MAX){
        int p = this->end_node;
        while(parent[p]!=-1){
            Path.push_back(p); 
            p = parent[p];
        }
        Path.push_back(p); 
    }

    std::reverse(Path.begin(),Path.end());

    return Path;

}

#endif