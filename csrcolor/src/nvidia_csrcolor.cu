#include <stdio.h>
#include <iostream>
#include <math.h>
#include "../include/nvidia_csrcolor.h"
#include <sys/time.h>
#include <cusparse.h>
#define BLOCKSIZE 1024
using namespace std;

int* colorGraph(int numVertices, int* nonZeroIndices, int* numColIndices, int numElems, int* nextVertexIndices, int maxDeg){

  int* globalColors;
  globalColors = (int *)malloc(sizeof(int)*numVertices);

  int *d_nonZeroIndices, *d_numColIndices;
  int *d_globalColoring;
  float* d_csrVal;
  int numColors = 0;
  float fractionToColor = 1.0;
  cudaEvent_t cudaStart, cudaEnd;
  cudaEventCreate(&cudaStart);
  cudaEventCreate(&cudaEnd);
  cudaMalloc((void **)&d_nonZeroIndices, sizeof(int)*numElems);
  cudaMemcpy((void *)d_nonZeroIndices, (const void *) nonZeroIndices, sizeof(int)*numElems, cudaMemcpyHostToDevice);
  cudaMalloc((void **)&d_numColIndices, sizeof(int)*(numVertices + 1));
  cudaMemcpy((void *)d_numColIndices, (const void *) numColIndices, sizeof(int)*(numVertices + 1), cudaMemcpyHostToDevice);
  cudaMalloc((void **)&d_globalColoring, sizeof(int)*numVertices);
  cudaMemset((void *)d_globalColoring, -1, sizeof(int)*numVertices);
  cudaMalloc((void **)&d_csrVal, sizeof(float)*numElems);
  
  cusparseStatus_t status;

  cusparseHandle_t handle;
  status = cusparseCreate(&handle);
  cusparseMatDescr_t descriptor;
  status = cusparseCreateMatDescr(&descriptor);
  cusparseColorInfo_t info;
  status = cusparseCreateColorInfo(&info);
  
  cudaEventRecord(cudaStart, 0);
  status = cusparseScsrcolor(handle, numVertices, numElems, descriptor, d_csrVal, d_numColIndices, d_nonZeroIndices, &fractionToColor, &numColors, d_globalColoring, NULL, info);
  cudaDeviceSynchronize();
  cudaEventRecord(cudaEnd, 0);

  cudaMemcpy((void *)globalColors,(const void *) d_globalColoring, sizeof(int)*numVertices, cudaMemcpyDeviceToHost);
  
  cudaFree((void *)d_globalColoring);
  cudaFree(d_nonZeroIndices);
  cudaFree(d_numColIndices);
  cudaFree(d_csrVal);
  cudaDeviceSynchronize();
  
  float time;
  cudaEventElapsedTime(&time, cudaStart, cudaEnd);
  /*
  cout << "Cuda Kernel Runtime: " << time << endl;
  cout << "Number of Colors: " << numColors << endl;
  */
  cout << time << " " << numColors << endl;
  return globalColors;

} 
