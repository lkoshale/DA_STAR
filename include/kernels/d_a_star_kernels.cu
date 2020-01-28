#ifdef __NVCC__

__device__ volatile int PQ[MAX_NODE];


//K in parallel
template <class U>
__global__ void extractMin(int* PQ_size, int* expandNodes,int* expandNodes_size,U* Cx,int* openList,int N,int K){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id<K && PQ_size[id]>0){

        //extract min from PQ
        int front = id* ( (N+K-1)/K );
        int node = PQ[front];

        // restructure the heap
        PQ[front]=PQ[front+PQ_size[id]-1];
        PQ_size[id]-=1;
        int pqIndex = 0;

        while(2*pqIndex+1 < PQ_size[id]){
            if(2*pqIndex+2 >= PQ_size[id]){
                if( Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+1]]){
                    int swap = PQ[front + 2*pqIndex+1];
                    PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+1;
                }
                else
                    break;
            }
            else{
                if( Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+1]] && Cx[PQ[front+2*pqIndex+1]] <= Cx[PQ[front+2*pqIndex+2]] ){
                    int swap = PQ[front + 2*pqIndex+1];
                    PQ[front + 2*pqIndex+1] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+1;
                }
                else if(Cx[PQ[front+pqIndex]] > Cx[PQ[front+2*pqIndex+2]] && Cx[PQ[front+2*pqIndex+2]] <= Cx[PQ[front+2*pqIndex+1]] ){
                    int swap = PQ[front + 2*pqIndex+2];
                    PQ[front + 2*pqIndex+2] = PQ[front +pqIndex];
                    PQ[front + pqIndex] = swap;
                    pqIndex = 2*pqIndex+2;                    
                } 
                else{
                    break;
                }
            }
            
        }

        //removed from openList
        openList[node] = -1;

        //added to expand next
        int len = atomicAdd(expandNodes_size,1);
        expandNodes[len]=node;
    }

} 


//for K in parallel
template <class T,class U>
__global__ void A_star_expand(int* off,int* edge,unsigned T* W,U* Hx,int* parent,volatile U* Cx,
    int* expandNodes,int* expandNodes_size, int* lock ,int* flagfound,int* openList,
    int N,int E, int K,int dest,int* nVFlag,int* PQ_size,
    int flagDiff,int* diff_off,int* diff_edge,unsigned int* diff_weight,int dE ){
       
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id< *expandNodes_size ){

        int node = expandNodes[id];
        
        //reach dest
        if(node == dest){
            atomicOr(flagfound,1);
        }

        // expand
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
        while(start < end){ 
            int child = edge[start];
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }

            //array L initilaized with 0
            //get the lock for child to update C(x)
            //loop till acquire the lock
            bool leaveLoop = false;

            while(leaveLoop==false){

                if(atomicCAS(&lock[child],0,1)==0){
                    //critical section
                    if( Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        if(openList[child]==-1){
                            nVFlag[child]=1;
                            //add only once
                        }
                    }

                    //end critical section
                    leaveLoop = true;

                    atomicCAS(&lock[child],1,0);

                }

                __syncthreads();

            }

            start++;
        }

        //diff expand
        if(flagDiff){

            start = diff_off[node];
            end = dE;
            if(node!=N-1)
                end = diff_off[node+1];

            while(start<end){ 
                int child = diff_edge[start];
                
                //deleted edges
                if(child<0){
                    start++;
                    continue;
                }
    
                //array L initilaized with 0
                //get the lock for child to update C(x)
                bool leaveLoop = false;

                while(!leaveLoop){

                    if(atomicCAS(&lock[child],0,1)==0){
                        //critical section
                        if( Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                            Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                            __threadfence();
                            parent[child] = node;
            
                            if(openList[child]==-1){
                                nVFlag[child]=1;
                                //add only once
                            }
                        }

                        //end critical section
                        leaveLoop = true;

                        atomicCAS(&lock[child],1,0);

                    }

                    __syncthreads();
                    
                }
                
                start++;
            }
            
        }
        //end diff
        
    }//end 

}


//K in parallel -- O(N)
template <class U>
__global__ void keepHeapPQ(int* PQ_size,U* Cx,int N,int K){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < K && PQ_size[id] > 0){
        int front  = id*( (N+K-1)/K );
        int size = PQ_size[id];
        
        for(int i=front;i<front+size;i++){
            if(2*i+2 < front+size){
                int cost = Cx[PQ[i]];
                int costLeft = Cx[PQ[2*i+1]];
                int costRight = Cx[PQ[2*i+2]]; 
                if( cost > costLeft  ||  cost > costRight  ){
                    int index ;
                    if(costLeft <= costRight)
                        index = 2*i+1;
                    else
                        index = 2*i+2;
                    
                    while(index > front){
                        if( Cx[PQ[(index-1)/2]] > Cx[PQ[index]] ){
                            int swap = PQ[index];
                            PQ[index] = PQ[(index-1)/2];
                            PQ[(index-1)/2] = swap;
                            index = (index-1)/2;
                        }
                        else
                            break;
                    }
                }
            }
            else if(2*i+1 < front+size){
                if(Cx[PQ[i]] > Cx[PQ[2*i+1]]){
                    int index = 2*i+1;
                    while(index > front){
                        if( Cx[PQ[(index-1)/2]] > Cx[PQ[index]] ){
                            int swap = PQ[index];
                            PQ[index] = PQ[(index-1)/2];
                            PQ[(index-1)/2] = swap;
                            index = (index-1)/2;
                        }
                        else
                            break;
                    }
                }
            }
        }
    }
}

//N threads
__global__ void setNV(int* nextFlag,int* nextV,int* nvSize,int N){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < N){
        if(nextFlag[id]==1){
            int index = atomicAdd(nvSize,1);
            nextV[index]=id;
        }
    }
}


//for K in parallel
template <class U>
__global__ void insertPQ(int* PQS,int* nextV,int* nVsize,U* Cx,int K,int N,int* openList){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < K){

        int front = id*( (N+K-1)/K );
        int i = id;
        
        while(i<*nVsize){            
            //if not already present
            if(openList[nextV[i]]!=-1){
                i+=K;
                continue;
            }

            PQ[front+PQS[id]]= nextV[i];
            PQS[id]+=1;

            //add in openList
            openList[nextV[i]] = id;

            if(PQS[id]>1){
                int index = PQS[id]-1;
                while(index>0){
                    if(Cx[PQ[front+ (index-1)/2]] > Cx[PQ[front+index]]){
                        int swap = PQ[front+index];
                        PQ[front+index]=PQ[front+ (index-1)/2];
                        PQ[front+ (index-1)/2] = swap;
                        index = (index-1)/2;
                    }
                    else
                        break;
                }
            }
            i += K;
        }
    }
}


//for K in parallel
template <class U>
__global__ void checkMIN(int* PQ_size,int* flagEnd,U* Cx,int dest,int N,int K){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id < K && PQ_size[id] > 0 ){
        int front = id* ( (N+K-1)/K );
        int node = PQ[front];
        //check if atleast one min, dont end the a*
        if( Cx[node] < Cx[dest] ){
            atomicAnd(flagEnd,0);
        }
    }
}

template <class T, class U>
__global__ void propogateDel(int* delEdgesV,int delEdge, volatile U* Cx,
                int* rev_offset,int* rev_edges,T* rev_weight,int N,int E,
                U* Hx,volatile int* parent,int* parent_old,int* addFlag,
                int* rev_diff_offset,int* rev_diff_edges,T* rev_diff_weight,int dE){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;

    if(id<delEdge){
        int node = delEdgesV[id];
        //check for the parent and add to nextflag and update the cost
        

        int start = rev_offset[node];
        int end = E;
        if(node!=N-1)
            end = rev_offset[node+1];
        
        //no parent
        // write in parent read always from old_parent
        parent[node] = -1;
        Cx[node]=INT_MAX;
        addFlag[node]=1;

        int cost = INT_MAX;
        int opt_parent = -1;

        //if any parent can change the cost 
        while(start< end){
            int p = rev_edges[start];
            
            //del edges
            if(p<0 || p==node){
                start++;
                continue;
            }
            
            int weight = rev_weight[start];
            int flag_cycle = false;
            
            //check parent doesn't contain node
            int ancestor = parent_old[p];
    
            while(ancestor>0){
                if(ancestor==node){
                    flag_cycle = true;
                    break;
                }
                ancestor = parent_old[ancestor];
                
            }
            
            
            //no need to lock only single parent so only one node in array so one node per thread
            if(!flag_cycle && Cx[p]!=INT_MAX && cost > (Cx[p]-Hx[p])+weight+Hx[node] ){
                cost = (Cx[p]-Hx[p] )+weight+Hx[node];
                opt_parent = p;
            }

            start++;
        }

        start = rev_diff_offset[node];
        end = dE;
        if(node!=N-1)
            end = rev_diff_offset[node+1];
        
        while(start< end){
            int p = rev_diff_edges[start];
            
            //del edges
            if(p<0 || p==node){
                start++;
                continue;
            }
            
            int weight = rev_diff_weight[start];
            int flag_cycle = false;
            
            //check parent doesn't contain node
            int ancestor = parent_old[p];
            
            while(ancestor!=-1){
                if(ancestor==node){
                    flag_cycle = true;
                    break;
                }
                ancestor = parent_old[ancestor];
                
            }

            //no need to lock only single parent so only one node in array so one node per thread
            if(!flag_cycle && Cx[p]!=INT_MAX && cost > (Cx[p]-Hx[p])+weight+Hx[node] ){
                cost = (Cx[p]-Hx[p] )+weight+Hx[node];
                opt_parent = p;
            }

            start++;
        }

        //write here
        if(cost!=INT_MAX){
            Cx[node]=cost;
            parent[node]=opt_parent;
        }

    }

}

//add inserted edges to propogate
template <class T, class U>
__global__ void propogateAdd(int* diff_off, int* diff_edges,T* diff_W,U* Hx,int* addFlag,
            volatile U* Cx,int* lock, int* parent, int* parent_old, int N, int dE){
    
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    
    if(id < N){
        int node = id;
        
        int start = diff_off[node];
        int end = dE;
        if(node!=N-1)
            end = diff_off[node+1];
        
        while(start < end ){
            int child = diff_edges[start];
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }

            //array L initilaized with 0
            //get the lock for child to update C(x)
            //loop till acquire the lock
            bool leaveLoop = false;

            while(!leaveLoop){

                if(atomicCAS(&lock[child],0,1)==0){
                    //critical section
                    bool flag_cycle = false;

                    int ancestor = node;
                    while(ancestor > 0){
                        if(ancestor==child){
                            flag_cycle = true;
                            break;
                        }
                        ancestor = parent_old[ancestor];
                        
                    }
                   
                    if(!flag_cycle && Cx[node] != INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child] ){
                        
                        Cx[child] =  (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child];
                
                        parent[child] = node;
                        __threadfence();

                        addFlag[child]=1;

                    }

                    //end critical section
                    leaveLoop = true;

                    atomicCAS(&lock[child],1,0);
                }

                __syncthreads();
            }

            start++;
        }
    }    

}


template <class T,class U>
__global__ void insert_propagate(int* nodes, int* size, int* off, int* edge,T* W,U* Hx,
            int N,int E,volatile U* Cx,int* lock, int* parent,int* addFlag,
            int* diff_off,int* diff_edge,T* diff_W,int dE){

    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < *size){
        int node = nodes[id];
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        while(start < end ){
            int child = edge[start];
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }
            bool leaveLoop = false;
                
            while(!leaveLoop){

                if(atomicExch(&lock[child],1)==0){
                    
                    if(Cx[node]!=INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }

                    leaveLoop = true;    
                    atomicExch(&lock[child],0);
                }
                __syncthreads();
            }
            
            start++;
        }

        start = diff_off[node];
        end = dE;
        if(node!=N-1)
            end = diff_off[node+1];
        
        while(start < end ){
            int child = diff_edge[start];
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }
            bool leaveLoop = false;

            while(!leaveLoop){

                if(atomicCAS(&lock[child],0,1)==0){
                    //critical section

                    if(Cx[node]!=INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ diff_W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                         
                    }

                    //end critical section
                    leaveLoop = true;
    
                    atomicCAS(&lock[child],1,0);

                }

                __syncthreads();

            }

            start++;
        }
    
    }

}

template <class T,class U>
__global__ void delete_propagate(int* nodes, int* size, int* off, int* edge,T* W,U* Hx,
                    int N,int E,volatile U* Cx,int* lock, int* parent,int* parent_old,int* addFlag,
                    int* diff_off,int* diff_edge,T* diff_W,int dE,
                    int* rev_offset,int* rev_edges,T* rev_weight,
                    int* rev_diff_offset,int* rev_diff_edges,T* rev_diff_weight){

    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id < *size){
        int node = nodes[id];
        int start = off[node];
        int end = E;
        if(node!=N-1)
            end = off[node+1];
        
        while(start < end ){
            int child = edge[start];
            if(child<0){
                start++;
                continue;
            }

            bool leaveLoop = false;
                
            while(!leaveLoop){

                if(atomicExch(&lock[child],1)==0){
                    if(Cx[node]!=INT_MAX && Cx[child]!=INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }
                    else 
                    if( (Cx[node]==INT_MAX && parent[child]==node ) || ( parent[child]==node && (Cx[child] < Cx[node] - Hx[node]+ W[start]+ Hx[child]) )  ){
                        //use back edges
                        int rstart = rev_offset[child];
                        int rend = E;
                        if(child!=N-1)
                            rend = rev_offset[child+1];
                        
                        //there is always one parent that is node.
                        Cx[child] = INT_MAX;
                        parent[child]=-1;

                        while(rstart < rend){
                            int p = rev_edges[rstart]; 
                            if(p<0 || p == child){
                                rstart++;
                                continue;
                            }

                            int weight = rev_weight[rstart];
                            bool flag_cycle = false;
                            
                            //check parent doesn't contain child
                        
                            int ancestor = parent_old[p];
                            
                            while(ancestor > 0){
                                if(ancestor==child){
                                    flag_cycle = true;
                                    break;
                                }
                                ancestor = parent_old[ancestor];
                            }
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }
                        
                        rstart =  rev_diff_offset[child];
                        rend = dE;
                        if(child!=N-1)
                            rend = rev_diff_offset[child+1];
                    
                        while(rstart < rend){
                            int p = rev_diff_edges[rstart]; 
                            
                            if(p<0 || p==child){
                                rstart++;
                                continue;
                            }

                            int weight = rev_diff_weight[rstart];
                            int flag_cycle = false;
                            
                            //check parent doesn't contain child
                            int ancestor = parent_old[p];
                            while(ancestor!=-1){
                                if(ancestor==child){
                                    flag_cycle = true;
                                    break;
                                }
                                    
                                ancestor = parent_old[ancestor];
                            }
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }

                        addFlag[child]=1;
                    }
                        

                    leaveLoop = true;
    
                    atomicExch(&lock[child],0);

                }

                __syncthreads();
            }

            start++;

        }


        start = diff_off[node];
        end = dE;
        if(node!=N-1)
            end = diff_off[node+1];
        
        while(start < end ){
            int child = diff_edge[start];
            
            //deleted edges
            if(child<0){
                start++;
                continue;
            }

            //array L initilaized with 0
            //get the lock for child to update C(x)
            //loop till acquire the lock
            bool leaveLoop = false;

            while(!leaveLoop){

                if(atomicCAS(&lock[child],0,1)==0){
                    if(Cx[node]!=INT_MAX && Cx[child]!=INT_MAX && Cx[child] > (Cx[node] - Hx[node])+ W[start]+ Hx[child] ){
                        Cx[child]  = (Cx[node] - Hx[node])+ W[start]+ Hx[child];
                        __threadfence();
                        parent[child] = node;
        
                        addFlag[child]=1;

                    }
                    else 
                    if((Cx[node]==INT_MAX && parent[child]==node )|| ( parent[child]==node && (Cx[child] < Cx[node] - Hx[node]+ diff_W[start]+ Hx[child]) )  ){
                        //use back edges
                        int rstart = rev_offset[child];
                        int rend = E;
                        if(child!=N-1)
                            rend = rev_offset[child+1];
                        
                        //there is always one parent that is node.
                        Cx[child] = INT_MAX;
                        parent[child]=-1;
 
                        while(rstart < rend){
                            int p = rev_edges[rstart]; 
                            
                            if(p<0 || p ==child){
                                rstart++;
                                continue;
                            }
 
                            int weight = rev_weight[rstart];
                            int flag_cycle = false;
                            
                            //check parent doesn't contain child
                            int ancestor = parent_old[p];
                            while(ancestor!=-1){
                                if(ancestor==child)
                                    flag_cycle = true;
                                ancestor = parent_old[ancestor];
                            }
                            
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }
 
                        rstart =  rev_diff_offset[child];
                        rend = dE;
                        if(child!=N-1)
                            rend = rev_diff_offset[child+1];
                    
                        while(rstart < rend){
                            int p = rev_diff_edges[rstart]; 
                            
                            if(p<0 || p==child){
                                rstart++;
                                continue;
                            }
 
                            int weight = rev_diff_weight[rstart];
                            int flag_cycle = false;
                            
                            //check parent doesn't contain child
                            int ancestor = parent_old[p];
                            while(ancestor!=-1){
                                if(ancestor==child){
                                 flag_cycle = true;
                                 break;
                                }
                                    
                                ancestor = parent_old[ancestor];
                            }
                            
                            if(!flag_cycle && Cx[p]!=INT_MAX && Cx[child] > (Cx[p]-Hx[p])+weight+Hx[child] ){
                                Cx[child] = (Cx[p]-Hx[p] )+weight+Hx[child];
                                parent[child] = p;
                            }
                            
                            rstart++;
                        }
 
                        addFlag[child]=1;
                    }
                    
                    
                    //end critical section
                    leaveLoop = true;
    
                    atomicCAS(&lock[child],1,0);
                }

                __syncthreads();

            }

            start++;
        }
        
    }

}

//do in 1 thread
template <class U>
__global__ void insertDest(int* PQ_size,U* Cx,int dest,int* openList){
    int id = 0;
    int front = 0;
    if(openList[dest]==-1){
        PQ[front+PQ_size[id]]= dest;
        PQ_size[id]+=1;

        //add in openList
        openList[dest] = id;

        if(PQ_size[id]>1){
            int index = PQ_size[id]-1;
            while(index>0){
                if(Cx[PQ[front+ (index-1)/2]] > Cx[PQ[front+index]]){
                    int swap = PQ[front+index];
                    PQ[front+index]=PQ[front+ (index-1)/2];
                    PQ[front+ (index-1)/2] = swap;
                    index = (index-1)/2;
                }
                else
                    break;
            }
        }
    }
    
}

template <class U>
__global__ void getCx(U* Cx,int dest,U* val){
    int id = blockIdx.x*blockDim.x+threadIdx.x;
    if(id==0){
        *val = Cx[dest];
    }
}



#endif