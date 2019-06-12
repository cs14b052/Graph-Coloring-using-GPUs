#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../include/cuda_launch_config.hpp"
#include "../include/sirgl32.h"
#include <thrust/scan.h>
#include <thrust/execution_policy.h>
#include <thrust/sequence.h>
#include <sys/time.h>
#include <curand.h>
#include <curand_kernel.h>
#define MAX_COLOR 1024
#define BLOCKSIZE 1024
#define ITERATIONS 1
#define FULL_MASK 0xffffffff
using namespace std;


#define CUDA_SAFE_CALL(ans) { cudaSafeCheck((ans), __FILE__, __LINE__);}
inline void cudaSafeCheck(cudaError_t call, const char *file, int line, bool abort=true){
  if (call != cudaSuccess){
    printf("Error: %s in file: %s at line: %d\n", cudaGetErrorString(call), file, line);
    if (abort)
      exit(call);
  }
}

__global__ void firstFit(int* d_numVertices, const int* __restrict__ d_nonZeroIndices, const int* __restrict__ d_numColIndices, int* d_globalColoring, int* maxColor, int* adjColors){
  
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  int totalThreads = blockDim.x * gridDim.x;
  int VerticesPerThread = *d_numVertices/totalThreads;
  if (*d_numVertices%totalThreads != 0)
      VerticesPerThread++;
  int maxColorVal = (int)ceilf((*maxColor)/32.0);
  for (int id = threadID; VerticesPerThread > 0 && id < *d_numVertices; id += totalThreads, VerticesPerThread--){
      // get colors of neighbours
      int start = __ldg(d_numColIndices + id);
      int end = __ldg(d_numColIndices + id + 1);
      for (int i = start; i < end; i++){
        int color = d_globalColoring[__ldg(d_nonZeroIndices + i)];
        int adjColorIndex = blockIdx.x*(maxColorVal)*blockDim.x + threadIdx.x + blockDim.x*(color/32);
        if (color == *maxColor + 1)
          continue;
        adjColors[adjColorIndex] = adjColors[adjColorIndex] & ~(1 << (color%32));
      }
      for (int i = 0; i < maxColorVal; i++){
        int adjColorIndex = blockIdx.x*maxColorVal*blockDim.x + threadIdx.x + blockDim.x*i;
        int color = __ffs(adjColors[adjColorIndex]);
        if (color){
          d_globalColoring[id] = i * 32 + color - 1;
          break;
        }
      }

      for (int j = 0; j < maxColorVal; j++){
          int adjColorIndex = blockIdx.x*maxColorVal*blockDim.x + threadIdx.x + blockDim.x*j;
          adjColors[adjColorIndex] = 0xffffffff;
      }
  }
}

extern __shared__ int array[];

__global__ void conflictResolve(const int* __restrict__ d_nonZeroIndices,const int* __restrict__ d_numColIndices, int* d_globalColoring,const int* __restrict__ d_nextVertexIndices, int* inlist, int* inNumV, int* outlist, int* globalIndex, int* maxColor, int* adjColors){
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  int totalThreads = blockDim.x * gridDim.x;
  int VerticesPerThread = *inNumV/totalThreads;
  if (*inNumV%totalThreads != 0)
      VerticesPerThread++;
  int *threadBlockPrefixSum = array;
  int threadOutNumV = 0;
  int maxColorVal = (int)ceilf((*maxColor)/32.0);
  for (int vertex = threadID; VerticesPerThread > 0; vertex += totalThreads, VerticesPerThread--){
      // get colors of neighbours
      int recolorVertex;
      threadOutNumV = 0;
      if (vertex < *inNumV){
        int id = inlist[vertex];
        int nextVertex  = __ldg(d_nextVertexIndices + id);
        int end = __ldg(d_numColIndices + id + 1);
        if (nextVertex == -1)
          nextVertex = end;
        int index;
        for (index = nextVertex; index < end; index++){
          if (d_globalColoring[id] == d_globalColoring[__ldg(d_nonZeroIndices + index)]){
            for (int j = 0; j < maxColorVal; j++){
                int adjColorIndex = blockIdx.x*maxColorVal*blockDim.x + threadIdx.x + blockDim.x*j;
                adjColors[adjColorIndex] = 0xffffffff;
            }

            int start = __ldg(d_numColIndices + id);
            for (int i = start; i < end; i++){
              int color = d_globalColoring[__ldg(d_nonZeroIndices + i)];
              int adjColorIndex = blockIdx.x*maxColorVal*blockDim.x + threadIdx.x + blockDim.x*(color/32);
              adjColors[adjColorIndex] = adjColors[adjColorIndex] & ~(1 << (color%32));
            }
            
            for (int i = 0; i < maxColorVal; i++){
              int adjColorIndex = blockIdx.x*maxColorVal*blockDim.x + threadIdx.x + blockDim.x*i;
              int color = __ffs(adjColors[adjColorIndex]);
              if (color){
                d_globalColoring[id] = i * 32 + color - 1;
                recolorVertex = id;
                threadOutNumV++;
                break;
              }
            }
            break;
          }
        } 
      }
      threadBlockPrefixSum[threadIdx.x] = threadOutNumV;
      int PrefixSumIndex = 0;
      int reduceOffset = 0;
      int lastThreadNum = (blockDim.x > 32) ? 31 : (blockDim.x-1);
      if (threadIdx.x%32 == lastThreadNum){
          thrust::exclusive_scan(thrust::device, threadBlockPrefixSum + threadIdx.x - lastThreadNum, threadBlockPrefixSum + threadIdx.x + 1, threadBlockPrefixSum + threadIdx.x - lastThreadNum);
          PrefixSumIndex = atomicAdd(globalIndex, threadBlockPrefixSum[threadIdx.x] + threadOutNumV);
      }
      __syncwarp(FULL_MASK);
      reduceOffset = threadBlockPrefixSum[threadIdx.x];
      threadBlockPrefixSum[threadIdx.x] = PrefixSumIndex;
      PrefixSumIndex = thrust::reduce(thrust::device, threadBlockPrefixSum + threadIdx.x - threadIdx.x%32, threadBlockPrefixSum + threadIdx.x - threadIdx.x%32 + 32, 0, thrust::maximum<int>());
      int arrayindex = PrefixSumIndex + reduceOffset;
      
      if (threadOutNumV == 1)
        outlist[arrayindex] = recolorVertex;
  }
}


__global__ void initializeColorsL32(int* d_globalColoring, int* d_numVertices, int* maxDeg){
  int threadID = blockIdx.x * blockDim.x + threadIdx.x;
  int totalThreads = blockDim.x * gridDim.x;
  int VerticesPerThread = *d_numVertices/totalThreads;
  if (*d_numVertices%totalThreads != 0)
      VerticesPerThread++;
  for (int id = threadID; VerticesPerThread > 0 && id < *d_numVertices; id += blockDim.x, VerticesPerThread--){
    int color = (*maxDeg + 1);
    d_globalColoring[id] = color;
  } 
}

__global__ void swapPointersL32(int* inNumV, int* globalIndex){
  *inNumV = *globalIndex;
  *globalIndex = 0;
}


int* colorGraphSIRGL32(int numVertices, int* nonZeroIndices, int* numColIndices, int numElems, int* nextVertexIndices, int maxDeg){

  double kernel_runtime[ITERATIONS];
  int numColors[ITERATIONS];
  int* globalColors;
  globalColors = (int *)malloc(sizeof(int)*numVertices);
  int maxColor;
  int blocks;
  for(int i = 0; i < ITERATIONS; i++){
    kernel_runtime[i] = 0;
    numColors[i] = 0;
    int *d_numVertices;
    int *d_nonZeroIndices, *d_numColIndices,*d_nextVertexIndices;
    int *d_globalColoring;
    cudaEvent_t cudaStart, cudaEnd;
    cudaEventCreate(&cudaStart);
    cudaEventCreate(&cudaEnd);
    cudaMalloc((void **)&d_numVertices, sizeof(int));
    cudaMemcpy((void *)d_numVertices, (const void *) &numVertices, sizeof(int), cudaMemcpyHostToDevice);
    cudaMalloc((void **)&d_nonZeroIndices, sizeof(int)*numElems);
    cudaMemcpy((void *)d_nonZeroIndices, (const void *) nonZeroIndices, sizeof(int)*numElems, cudaMemcpyHostToDevice);
    cudaMalloc((void **)&d_numColIndices, sizeof(int)*(numVertices + 1));
    cudaMemcpy((void *)d_numColIndices, (const void *) numColIndices, sizeof(int)*(numVertices + 1), cudaMemcpyHostToDevice);
    cudaMalloc((void **)&d_nextVertexIndices, sizeof(int)*numVertices);
    cudaMemcpy((void *)d_nextVertexIndices, (const void *) nextVertexIndices, sizeof(int)*numVertices, cudaMemcpyHostToDevice);
    cudaMalloc((void **)&d_globalColoring, sizeof(int)*numVertices);

    cudaDeviceProp deviceProperty;
    cudaGetDeviceProperties(&deviceProperty, 0);
    int numSMs = deviceProperty.multiProcessorCount;

    blocks = maximum_residency(conflictResolve, BLOCKSIZE, 0) * numSMs;
    blocks = numSMs;
    int* inlist, *outlist;
    cudaMalloc((void **)&inlist, sizeof(int)*numVertices);
    thrust::sequence(thrust::device, inlist, inlist + numVertices);
    cudaMalloc((void **)&outlist, sizeof(int)*numVertices);

    int *inNumV;
    cudaMalloc((void **)&inNumV, sizeof(int));
    cudaMemcpy((void *)inNumV, (const void *) &numVertices, sizeof(int), cudaMemcpyHostToDevice);

    int *globalIndex;
    cudaMalloc((void **)&globalIndex, sizeof(int));


    maxColor = pow(2, ceil(log(maxDeg + 1)/log(2)));
    int* d_maxColor;
    cudaMalloc((void **)&d_maxColor, sizeof(int));
    cudaMemcpy((void *)d_maxColor, (const void *)&maxColor, sizeof(int), cudaMemcpyHostToDevice);
    
    int *d_adjColors;
    cudaMalloc((void **)&d_adjColors, sizeof(int)*BLOCKSIZE*blocks*(ceil(maxColor/32.0)));
    cudaMemset((void *)d_adjColors, 0xffffffff, sizeof(int)*BLOCKSIZE*blocks*(ceil(maxColor/32.0)));
    cudaEventRecord(cudaStart, 0);

    initializeColorsL32<<<blocks, BLOCKSIZE>>>(d_globalColoring, d_numVertices, d_maxColor);
    firstFit<<<blocks, BLOCKSIZE>>> (d_numVertices, d_nonZeroIndices, d_numColIndices, d_globalColoring, d_maxColor, d_adjColors);
    int hInNumV = numVertices;
    int iterNum = 0;
    while(hInNumV > 0){
      if ((iterNum & 1) == 0){
        conflictResolve<<<blocks, BLOCKSIZE, BLOCKSIZE*4>>> (d_nonZeroIndices, d_numColIndices, d_globalColoring, d_nextVertexIndices, inlist, inNumV, outlist, globalIndex, d_maxColor, d_adjColors);
      }
      else{
        conflictResolve<<<blocks, BLOCKSIZE, BLOCKSIZE*4>>> (d_nonZeroIndices, d_numColIndices, d_globalColoring, d_nextVertexIndices, outlist, inNumV, inlist, globalIndex, d_maxColor, d_adjColors);
      }
       swapPointersL32<<<1,1>>>(inNumV, globalIndex);
       cudaMemcpy((void *)&hInNumV, (const void *)inNumV, sizeof(int), cudaMemcpyDeviceToHost);
      iterNum++;
    }
    cudaDeviceSynchronize();
    cudaEventRecord(cudaEnd, 0);
    cudaFree(inlist);
    cudaFree(outlist);
    cudaFree(inNumV);
    cudaFree(globalIndex);
    cudaFree(d_maxColor);
    cudaFree(d_adjColors);
    cudaMemcpy((void *)globalColors,(const void *) d_globalColoring, sizeof(int)*numVertices, cudaMemcpyDeviceToHost);
    cudaFree((void *)d_globalColoring);
    cudaFree(d_numVertices);
    cudaFree(d_nonZeroIndices);
    cudaFree(d_numColIndices);
    cudaFree(d_nextVertexIndices);
    cudaDeviceSynchronize();
    float time;
    cudaEventElapsedTime(&time, cudaStart, cudaEnd);
    kernel_runtime[i] = time;
    numColors[i] = countColors(globalColors, numVertices);
    cudaEventDestroy(cudaStart);
    cudaEventDestroy(cudaEnd);
  }
    
  double averageKernelRuntime = 1.0;
  double averageColors = 1.0;

  for(int i = 0; i < ITERATIONS; i++){
    averageKernelRuntime *= kernel_runtime[i];
    averageColors *= numColors[i];
//    cout << "[ " << kernel_runtime[i] << " " << numColors[i] << " ] ";
  }
//  cout << endl;
  averageKernelRuntime = pow(averageKernelRuntime, 1.0/((double)ITERATIONS));
  averageColors = pow(averageColors, 1.0/((double)ITERATIONS));
  cout << fixed;
  cout << averageKernelRuntime << " " << averageColors << " " << sizeof(int)*BLOCKSIZE*blocks*(ceil(maxColor/32.0)) << endl;
//  cout << "Average Runtime: " << averageKernelRuntime << endl;
//  cout << "Average Colors: " << averageColors << endl;
//  cout << "Memory Used for adjColors array: " << sizeof(int)*BLOCKSIZE*blocks*(ceil(maxColor/32.0)) << " bytes" << endl;
  return globalColors;
}
