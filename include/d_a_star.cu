#ifdef __NVCC__

template <class U>
__global__ void extractMin(int* PQ_size, int* expandNodes,int* expandNodes_size,U* Cx,int* openList,int N,int K);


template <class T,class U>
__global__ void A_star_expand(int* off,int* edge,unsigned T* W,U* Hx,int* parent,volatile U* Cx,
    int* expandNodes,int* expandNodes_size, int* lock ,int* flagfound,int* openList,
    int N,int E, int K,int dest,int* nVFlag,int* PQ_size,
    int flagDiff,int* diff_off,int* diff_edge,unsigned int* diff_weight,int dE );


template <class U>
__global__ void keepHeapPQ(int* PQ_size,U* Cx,int N,int K);

__global__ void setNV(int* nextFlag,int* nextV,int* nvSize,int N);

template <class U>
__global__ void insertPQ(int* PQS,int* nextV,int* nVsize,U* Cx,int K,int N,int* openList);

template <class U>
__global__ void checkMIN(int* PQ_size,int* flagEnd,U* Cx,int dest,int N,int K);

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

#include "d_a_star.cuh"

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

    this->num_updated_paths = 0;

    __alloc_cpu();

}

template <class T,class U>
void GPU_D_A_Star<T,U>:: set_huiristics(U* hx)
{
    this->Hx = hx;
    is_set_hx = true;
    int N = this->graph->get_graph().get_num_nodes();
    gpuErrchk ( cudaMalloc(&d_Hx,sizeof(U)*N ) );
    gpuErrchk ( cudaMemcpy(d_Hx,Hx,sizeof(U)*N,cudaMemcpyHostToDevice) );

}


template <class T,class U>
void GPU_D_A_Star<T,U>:: __alloc_cpu()
{
    int N = this->graph->get_graph().get_num_nodes();
    int K = this->num_pq;

    this->PQ = (unsigned int*)malloc(sizeof(unsigned int)*N );
    this->PQ_size = (unsigned int*)malloc(sizeof(unsigned int)*K);

    this->Cx = (U*)malloc(sizeof(U)*N);
    this->Hx = (U*)malloc(sizeof(U)*N);
    
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
    gpuErrchk ( cudaMalloc(&d_expand_nodes,sizeof(int)*K) );  //changed to K
    gpuErrchk ( cudaMalloc(&d_expand_nodes_size,sizeof(int)) );

    //flag to end search
    gpuErrchk( cudaMalloc(&d_flag_end,sizeof(int)) );
    gpuErrchk( cudaMalloc(&d_flag_found,sizeof(int)) );
    
    gpuErrchk ( cudaMemset(d_next_vertices_size,0,sizeof(int)) );
    gpuErrchk ( cudaMemset(d_expand_nodes_size,0,sizeof(int)) );
    gpuErrchk ( cudaMemset(d_lock,0,sizeof(int)*N) );


    // gpuErrchk ( cudaMemcpy(d_Cx,Cx,sizeof(U)*N,cudaMemcpyHostToDevice) );

    // gpuErrchk ( cudaMemcpy(d_PQ_size,PQ_size,sizeof(unsigned int)*num_pq,cudaMemcpyHostToDevice) );

    // gpuErrchk ( cudaMemcpy(d_parent,parent,sizeof(int)*N,cudaMemcpyHostToDevice) );
    // gpuErrchk ( cudaMemcpy(d_open_list,open_list,sizeof(int)*N,cudaMemcpyHostToDevice) );

}




#endif