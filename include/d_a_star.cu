#ifdef __NVCC__

template <class U>
__global__ void extractMin(unsigned int* PQ, unsigned int* PQ_size, int* expandNodes,int* expandNodes_size,U* Cx,int* openList,int N,int K);


template <class T,class U>
__global__ void A_star_expand(int* off, int* edge, T* W, U* Hx, int* parent, volatile U* Cx,
    int* expandNodes, int* expandNodes_size, int* lock , int* flagfound, int* openList, int* nVFlag,
    int N, int E, int K, int dest,
    int flagDiff, int dE,
    int* diff_off, int* diff_edge, unsigned int* diff_weight );


template <class U>
__global__ void keepHeapPQ(unsigned int* PQ,unsigned int* PQ_size,U* Cx,int N,int K);

__global__ void setNV(int* nextFlag,int* nextV,int* nvSize,int N);

template <class U>
__global__ void insertPQ(unsigned int* PQ,unsigned int* PQS,int* nextV,int* nVsize,U* Cx,int K,int N,int* openList);

template <class U>
__global__ void checkMIN(unsigned int* PQ, unsigned int* PQ_size,int* flagEnd,U* Cx,int dest,int N,int K);

template <class T, class U>
__global__ void propogateDel(int* delEdgesV,int delEdge, volatile U* Cx,
                int* rev_offset,int* rev_edges,T* rev_weight,int N,int E,
                U* Hx,volatile int* parent,int* parent_old,int* addFlag,
                int* rev_diff_offset,int* rev_diff_edges,T* rev_diff_weight,int dE);


template <class T, class U>
__global__ void propogateAdd(int* diff_off, int* diff_edges,T* diff_W,U* Hx,int* addFlag,
            volatile U* Cx,int* lock, int* parent, int* parent_old, int N, int dE);


template <class T,class U>
__global__ void insert_propagate(int* nodes, int* size, int* off, int* edge,T* W,U* Hx,
            int N,int E,volatile U* Cx,int* lock, int* parent,int* addFlag,
            int* diff_off,int* diff_edge,T* diff_W,int dE);



template <class T,class U>
__global__ void delete_propagate(int* nodes, int* size, int* off, int* edge,T* W,U* Hx,
                    int N,int E,volatile U* Cx,int* lock, int* parent,int* parent_old,int* addFlag,
                    int* diff_off,int* diff_edge,T* diff_W,int dE,
                    int* rev_offset,int* rev_edges,T* rev_weight,
                    int* rev_diff_offset,int* rev_diff_edges,T* rev_diff_weight);

template <class U>
__global__ void insertDest(int* PQ_size,U* Cx,int dest,int* openList);


template <class U>
__global__ void getCx(U* Cx,int dest,U* val);


///////////////////////////////////////////

#include "kernels/d_a_star_kernels.cu"

// #include "d_a_star.cuh"

#ifdef DEBUG
    #include <cstdio>
#endif

template <class T,class U>
GPU_D_A_Star<T,U>:: GPU_D_A_Star(GPU_Dynamic_Graph<T> *graph, unsigned int start,unsigned int end, unsigned int K )
{
    this->graph = graph;
    this->start_node = start;
    this->end_node = end;
    this->num_pq = K;

    this->flag_end = 0;
    this->flag_found = 0;

    this->is_set_hx = false;
    this->next_vertices_size = 0;

    this->num_threads = 512;

    this->num_updates = 0;

    this->update_file = NULL;

}

template <class T,class U>
void GPU_D_A_Star<T,U>:: set_heuristics(U* hx)
{
    this->Hx = hx;
    is_set_hx = true;
    int N = this->graph->get_graph().get_num_nodes();
    gpuErrchk ( cudaMalloc(&d_Hx,sizeof(U)*N ) );
    gpuErrchk ( cudaMemcpy(d_Hx,Hx,sizeof(U)*N,cudaMemcpyHostToDevice) );

    __alloc_cpu();

}


template <class T,class U>
void GPU_D_A_Star<T,U>:: __alloc_cpu()
{
    int N = this->graph->get_graph().get_num_nodes();
    int K = this->num_pq;

    this->PQ = (unsigned int*)malloc(sizeof(unsigned int)*N );
    this->PQ_size = (unsigned int*)malloc(sizeof(unsigned int)*K);

    this->Cx = (U*)malloc(sizeof(U)*N);
   
    
    this->open_list = (int*)malloc(sizeof(int)*N);

    this->parent = (int*)malloc(sizeof(int)*N);
    this->parent_old = (int*)malloc(sizeof(int)*N);

    this->next_vertices_flag = (int*)malloc(sizeof(int)*N);
    this->next_vertices = (int*)malloc(sizeof(int)*N);

    memset(this->parent,-1,sizeof(int)*N);
    memset(this->parent_old,-1,sizeof(int)*N);

    memset(this->open_list,-1,sizeof(int)*N);
    memset(this->PQ_size,0,sizeof(int)*K);

    memset(this->next_vertices_flag,-1,sizeof(int)*N);
    
    //todo make it memset
    for(int i=0;i<N;i++){
        this->Cx[i] = INT_MAX;
    }

    //append start node
    //insert startNode in PQ[0]
    Cx[start_node]=Hx[start_node];
    PQ[0]=start_node;
    PQ_size[0]=1;
    open_list[start_node]=0;


    //allocate on GPU

    __alloc_gpu();

}

template <class T,class U>
void GPU_D_A_Star<T,U>:: __alloc_gpu()
{
    int N = this->graph->get_graph().get_num_nodes();

    gpuErrchk ( cudaMalloc(&d_Cx,sizeof(U)*N ) );

    gpuErrchk ( cudaMalloc(&d_parent,sizeof(int)*N ) );
    gpuErrchk ( cudaMalloc(&d_parent_old,sizeof(int)*N ) );
    
    gpuErrchk ( cudaMalloc(&d_open_list,sizeof(int)*N ) );

    gpuErrchk ( cudaMalloc(&d_PQ,sizeof(unsigned int)*N ) );
    gpuErrchk ( cudaMalloc(&d_PQ_size,sizeof(unsigned int)*num_pq ) );


    gpuErrchk ( cudaMalloc(&d_lock,sizeof(int)*N) );

    //for next set of vertices to add in PQ
    gpuErrchk ( cudaMalloc(&d_next_vertices,sizeof(int)*N) );
    gpuErrchk ( cudaMalloc(&d_next_vertices_size,sizeof(int)) );
    gpuErrchk ( cudaMalloc(&d_next_vertices_flag,sizeof(int)*N) );

    //next nodes to expand
    gpuErrchk ( cudaMalloc(&d_expand_nodes,sizeof(int)*num_pq) );  //changed to K
    gpuErrchk ( cudaMalloc(&d_expand_nodes_size,sizeof(int)) );

    //flag to end search
    gpuErrchk( cudaMalloc(&d_flag_end,sizeof(int)) );
    gpuErrchk( cudaMalloc(&d_flag_found,sizeof(int)) );
    
    gpuErrchk ( cudaMemset(d_next_vertices_size,0,sizeof(int)) );
    gpuErrchk ( cudaMemset(d_expand_nodes_size,0,sizeof(int)) );
    gpuErrchk ( cudaMemset(d_lock,0,sizeof(int)*N) );

    //copy from cpu to gpu
    gpuErrchk ( cudaMemcpy(d_Cx,Cx,sizeof(U)*N,cudaMemcpyHostToDevice) );
    

    gpuErrchk ( cudaMemcpy(d_PQ_size,PQ_size,sizeof(unsigned int)*num_pq,cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpy(d_parent,parent,sizeof(int)*N,cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(d_open_list,open_list,sizeof(int)*N,cudaMemcpyHostToDevice) );

    gpuErrchk ( cudaMemcpy(d_flag_end,&flag_end,sizeof(int),cudaMemcpyHostToDevice) );
    gpuErrchk ( cudaMemcpy(d_flag_found,&flag_found,sizeof(int),cudaMemcpyHostToDevice) );
   
    gpuErrchk ( cudaMemcpy(d_next_vertices_flag,next_vertices_flag,sizeof(int)*N,cudaMemcpyHostToDevice) );

}


template <class T,class U>
void GPU_D_A_Star<T,U>:: extract_min()
{
    //K parallel
    int N = this->graph->get_graph().get_num_nodes();
    int num_blocks = (num_pq+num_threads-1)/num_threads;

    extractMin < U > <<< num_blocks,num_threads >>>( d_PQ, d_PQ_size, d_expand_nodes, d_expand_nodes_size, d_Cx, d_open_list, N, num_pq);
        
    gpuErrchk(cudaPeekAtLastError() );

    cudaDeviceSynchronize();

}


template <class T,class U>
void GPU_D_A_Star<T,U>:: expand()
{
    int N = this->graph->get_graph().get_num_nodes();
    int E = this->graph->get_graph().get_num_edges();
    int num_blocks = (num_pq+num_threads-1)/num_threads;

    A_star_expand < T, U > <<<num_blocks,num_threads>>>(
        this->graph->get_graph().get_offsets(), this->graph->get_graph().get_edges(), this->graph->get_graph().get_weights(),
        d_Hx, d_parent, d_Cx,
        d_expand_nodes,d_expand_nodes_size,
        d_lock ,d_flag_found,d_open_list,d_next_vertices_flag,
        N, E, num_pq, end_node,
        false,0,
        this->graph->get_diff_graph().get_offsets(),this->graph->get_diff_graph().get_edges(),this->graph->get_diff_graph().get_weights()
    );

    gpuErrchk(cudaPeekAtLastError() );

    cudaDeviceSynchronize();
}


template <class T,class U>
void GPU_D_A_Star<T,U>::  maintain_heap()
{
    
    int N = this->graph->get_graph().get_num_nodes();

    int num_blocks = (num_pq+num_threads-1)/num_threads;

    keepHeapPQ < U > <<< num_blocks, num_threads>>>(d_PQ, d_PQ_size, d_Cx, N, num_pq);

    gpuErrchk(cudaPeekAtLastError() );

    cudaDeviceSynchronize();

}


template <class T,class U>
void GPU_D_A_Star<T,U>:: set_flags()
{
    int N = this->graph->get_graph().get_num_nodes();
    int num_blocks = (N+num_threads-1)/num_threads;           //N_num_blocks

    setNV <<<num_blocks,num_threads>>>(d_next_vertices_flag, d_next_vertices, d_next_vertices_size, N);
    gpuErrchk(cudaPeekAtLastError() );

    cudaDeviceSynchronize();
}

template <class T,class U>
void GPU_D_A_Star<T,U>:: insert()
{
    int N = this->graph->get_graph().get_num_nodes();

    int num_blocks = (num_pq+num_threads-1)/num_threads;

    insertPQ< U > <<<num_blocks, num_threads >>>(d_PQ, d_PQ_size, d_next_vertices, d_next_vertices_size, d_Cx, num_pq, N, d_open_list);
    gpuErrchk(cudaPeekAtLastError() );

    cudaDeviceSynchronize();

}

template <class T,class U>
void GPU_D_A_Star<T,U>:: check_all_min_pq()
{
    int N = this->graph->get_graph().get_num_nodes();

    int num_blocks = (num_pq+num_threads-1)/num_threads;

    checkMIN < U > <<< num_blocks, num_threads >>>(d_PQ, d_PQ_size, d_flag_end, d_Cx, end_node, N, num_pq);

    gpuErrchk(cudaPeekAtLastError() );

    cudaDeviceSynchronize();

}

template <class T, class U>
bool GPU_D_A_Star<T,U>:: is_empty_pq_cpu()
{
    bool is_not_empty = false;
    for(int i=0;i<num_pq;i++){
        if(PQ_size[i]>0)
            is_not_empty=true;
    }

    return !is_not_empty;
}

template <class T, class U>
std::vector<int>  GPU_D_A_Star<T,U>:: initial_path()
{

    if(!is_set_hx){
        printf("[INFO] heuristics value not set\n");
        printf("[DO] set heuristics to init algorithm\n");
        exit(0);
    }

    int N = this->graph->get_graph().get_num_nodes();

    while(flag_end == 0 && !is_empty_pq_cpu())
    {
        extract_min();

        expand();

        maintain_heap();

        set_flags();

        insert();

        //copy
        gpuErrchk( cudaMemcpy(&flag_found, d_flag_found, sizeof(int), cudaMemcpyDeviceToHost) );
        gpuErrchk( cudaMemcpy(PQ_size, d_PQ_size, sizeof(int)*num_pq, cudaMemcpyDeviceToHost) );

        //reset
        gpuErrchk( cudaMemcpy(d_next_vertices_flag, next_vertices_flag, sizeof(int)*N,cudaMemcpyHostToDevice) );

        //reset next insert array
        gpuErrchk ( cudaMemset(d_next_vertices_size,0,sizeof(int)) );
        gpuErrchk ( cudaMemset(d_expand_nodes_size,0,sizeof(int)) );

        if( !is_empty_pq_cpu() && flag_found==1)
        {
            gpuErrchk( cudaMemcpy(d_flag_end, &flag_found,sizeof(int),cudaMemcpyHostToDevice) );

            check_all_min_pq();

            gpuErrchk( cudaMemcpy(&flag_end,d_flag_end, sizeof(int),cudaMemcpyDeviceToHost) );

        }

    }

    gpuErrchk( cudaMemcpy(parent, d_parent, sizeof(int)*N,cudaMemcpyDeviceToHost) );

    U dest_cost;
    gpuErrchk( cudaMemcpy(&dest_cost,d_Cx+end_node, sizeof(U),cudaMemcpyDeviceToHost) );

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



template <class T, class U>
std::vector<int>  GPU_D_A_Star<T,U>:: get_path()
{
    std::vector<int> path;
    if(num_updates==0){
        path = initial_path();
        num_updates++;
    }
    else
    {
        //check fo update file
        //call dynamic prop

    }
    
    return path;
}

template <class T, class U>
std::vector<int>  GPU_D_A_Star<T,U>:: updated_path()
{
    std::vector<int> path;

    return path;
}


template <class T, class U>
void  GPU_D_A_Star<T,U>:: set_update_file(FILE* fptr)
{
    this->update_file = fptr;
}


#endif