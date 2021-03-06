// Copyright 2016, National University of Defense Technology
// Authors: Xuhao Chen <cxh@illinois.edu> and Pingfan Li <lipingfan@163.com>
#include <cub/cub.cuh>
#include <iostream>
#include "gbar.cuh"
#include "cuda_launch_config.hpp"
#include "cutil_subset.h"
#include "common.h"
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/count.h>
#include <thrust/reduce.h>
#include <thrust/functional.h>
#include <thrust/execution_policy.h>
#include "worklistc.h"
#define	SCRATCHSIZE BLKSIZE
#define	MAXCOLOR 1024
#define FUSION
typedef cub::BlockScan<int, BLKSIZE> BlockScan;
using namespace std;

__global__ void initialize(int *coloring, int m) {
	unsigned id = blockIdx.x * blockDim.x + threadIdx.x;
	if (id < m) {
		coloring[id] = MAXCOLOR;
	}   
}

__device__ __forceinline__ void assignColor(unsigned *forbiddenColors, int *coloring, int vertex) {
	int vertex_color;
	for (vertex_color = 0; vertex_color < MAXCOLOR/32; vertex_color++) {
		int pos = __ffs(forbiddenColors[vertex_color]);
		if(pos) {
			coloring[vertex] = vertex_color * 32 + pos - 1;
			break;
		}
	}
	assert(vertex_color < MAXCOLOR);
}
#if 0
__device__ __forceinline__ void maskByCta(int tid, int m, int vertex, int size, int *row_offsets, int *column_indices, int *degree, Worklist2 &inwl, int *coloring, /*bool forbiddenColors[BLKSIZE][MAXCOLOR+1]*/ unsigned forbiddenColors[BLKSIZE][MAXCOLOR/32+1]) {
	__shared__ int owner;
	__shared__ int sh_vertex;
	owner = -1;
	//int total_inputs = (*inwl.dindex + gridDim.x * blockDim.x - 1)/(gridDim.x * blockDim.x);
	//for (int id = tid; total_inputs > 0; id += blockDim.x * gridDim.x, total_inputs--) {
	int id = tid;
	int owner_tid;
	int my_size = size;
	//int vertex;
	//if(inwl.pop_id(id, vertex)) {
		//row_begin = row_offsets[vertex];
		//row_begin = __ldg(row_offsets + vertex);
		//row_end = row_offsets[vertex + 1];
		//row_end = __ldg(row_offsets + vertex + 1);
		//size = row_end - row_begin;
		//size = degree[vertex];
	//}
	while(true) {
		if(my_size >= BLKSIZE)
			owner = threadIdx.x;
		__syncthreads();
		owner_tid = owner;
		if(owner_tid == -1)
			break;
		__syncthreads();
		if(owner == threadIdx.x) {
			sh_vertex = vertex;
			// mark this vertex as processed already
			//inwl.dwl[id] = -1;
			cub::ThreadStore<cub::STORE_CG>(inwl.dwl + id, -1);
			owner = -1;
			my_size = 0;
		}
		__syncthreads();
		int row_begin = row_offsets[sh_vertex];
		//row_begin = __ldg(row_offsets + sh_vertex);
		//row_end = row_offsets[sh_vertex + 1];
		//row_end = __ldg(row_offsets + sh_vertex + 1);
		//int neighbor_size = row_end - row_begin;
		int neighbor_size = degree[sh_vertex];
		int num = ((neighbor_size - 1) / BLKSIZE + 1) * BLKSIZE;
		for(int i = threadIdx.x; i < num; i += BLKSIZE) {
			if(i < neighbor_size) {
				int edge = row_begin + i;
				int dst = column_indices[edge];
				//int dst = __ldg(column_indices + edge);
				int color = coloring[dst];
				//forbiddenColors[owner_tid][color] = true;
				atomicAnd(&forbiddenColors[owner_tid][color/32], ~(1<<(color%32)));
			}
		}
	}
	//}
}
#endif
/*
__global__ void firstFit(int m, int nnz, int *csrRowPtr, int *csrColInd, Worklist2 inwl, int *coloring)
{
	int global_id = blockIdx.x * blockDim.x + threadIdx.x;
	int local_id = threadIdx.x;
	typedef cub::BlockScan<int, BLKSIZE> BlockScan;
	__shared__ BlockScan::TempStorage temp_storage;
	__shared__ int gather_offsets[SCRATCHSIZE];
	__shared__ int srcIndex[BLKSIZE];	
	__shared__ bool forbiddenColors[BLKSIZE][MAXCOLOR+1];
	//__shared__ unsigned int forbiddenColors[BLKSIZE][MAXCOLOR/32];
	gather_offsets[local_id] = 0;
	int vertex;
	int neighbor_size = 0;
	int row_begin = 0;
	int row_end = 0;
	int scratch_offset = 0;
	int total_edges = 0;
	if (inwl.pop_id(global_id, vertex)) {
		if (vertex != -1) {
			row_begin = csrRowPtr[vertex];
			//row_begin = __ldg(csrRowPtr + vertex);
			row_end = csrRowPtr[vertex + 1];
			//row_end = __ldg(csrRowPtr + vertex + 1);
			neighbor_size = row_end - row_begin;
		}
	}
	BlockScan(temp_storage).ExclusiveSum(neighbor_size, scratch_offset, total_edges);
	//for (int i = 0; i < MAXCOLOR/32; i++)
		//forbiddenColors[local_id][i] = 0xffffffff;
	for (int i = 0; i < MAXCOLOR; i++)
		forbiddenColors[local_id][i] = false;
	int done = 0;
	int neighborsdone = 0;
	while (total_edges > 0) {
		__syncthreads();
		int i;
		for (i = 0; (neighborsdone + i) < neighbor_size && (scratch_offset + i - done) < SCRATCHSIZE; i++) {
			gather_offsets[scratch_offset + i - done] = row_begin + neighborsdone + i;
			srcIndex[scratch_offset + i - done] = local_id;
		}
		neighborsdone += i;
		scratch_offset += i;
		__syncthreads();

		if (threadIdx.x < total_edges) {
			int edge = gather_offsets[local_id];
			int dst = csrColInd[edge];
			//int dst = __ldg(csrColInd + edge);
			int index = srcIndex[local_id];
			int color = coloring[dst];
			//int color = cub::ThreadLoad<cub::LOAD_CG>(coloring + dst);
			//atomicAnd(&forbiddenColors[index][color/32], ~(1<<(color%32)));
			forbiddenColors[index][color] = true;
		}
		total_edges -= BLKSIZE;
		done += BLKSIZE;
	}
	__syncthreads();
	if (neighbor_size) {
		//assignColor(forbiddenColors[local_id], coloring, vertex);
		int vertex_color;
		for (vertex_color = 0; vertex_color < MAXCOLOR; vertex_color++) {
			if (!forbiddenColors[local_id][vertex_color]) {
				coloring[vertex] = vertex_color;
				break;
			}
		}
		assert(vertex_color < MAXCOLOR);
	}
}
//*/
__device__ void firstFit(int m, const int* __restrict__ csrRowPtr, const int* __restrict__ csrColInd, Worklist2 &inwl, int *coloring) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned forbiddenColors[MAXCOLOR/32+1];
	int id = tid;
	int total_inputs = (*inwl.dindex + gridDim.x * blockDim.x - 1)/(gridDim.x * blockDim.x);
	for (int id = tid; total_inputs > 0; id += blockDim.x * gridDim.x, total_inputs--) {
		int vertex;
		if (inwl.pop_id(id, vertex)) {
			int row_begin = __ldg(csrRowPtr + vertex);
			int row_end= __ldg(csrRowPtr + vertex + 1);
			for (int j = 0; j < MAXCOLOR/32; j++)
				forbiddenColors[j] = 0xffffffff;
			for (int offset = row_begin; offset < row_end; offset ++) {
				int neighbor = __ldg(csrColInd + offset);
				int color = coloring[neighbor];
				forbiddenColors[color / 32] &= ~(1 << (color % 32));
			}
			assignColor(forbiddenColors, coloring, vertex);
		}
	}
}

/*
__global__ void firstFit(int m, int *row_offsets, int *column_indices, int *degree, Worklist2 inwl, int *coloring) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tx = threadIdx.x;
	//__shared__ bool forbiddenColors[BLKSIZE][MAXCOLOR+1];
	__shared__ unsigned forbiddenColors[BLKSIZE][MAXCOLOR/32+1];
	//for (int j = 0; j < MAXCOLOR; j++)
	//	forbiddenColors[tx][j] = false;
	for (int j = 0; j < MAXCOLOR/32; j++)
		forbiddenColors[tx][j] = 0xffffffff;
	__syncthreads();
	int old_vertex;
	int my_degree = 0;
	if (inwl.pop_id(tid, old_vertex)) {
		my_degree = degree[old_vertex];
	}
	maskByCta(tid, m, old_vertex, my_degree, row_offsets, column_indices, degree, inwl, coloring, forbiddenColors);

	//int total_inputs = (*inwl.dindex + gridDim.x * blockDim.x - 1)/(gridDim.x * blockDim.x);
	//for (int id = tid; total_inputs > 0; id += blockDim.x * gridDim.x, total_inputs--) {
		int id = tid;
		int vertex;
		if (inwl.pop_id(id, vertex)) {
			if (vertex != -1) {
				int row_begin = row_offsets[vertex];
				int row_end = row_offsets[vertex + 1];
				for (int offset = row_begin; offset < row_end; offset ++) {
					int neighbor = column_indices[offset];
					int color = coloring[neighbor];
					forbiddenColors[tx][color / 32] &= ~(1 << (color % 32));
					//forbiddenColors[tx][color] = true;
				}
			}
			assignColor(forbiddenColors[tx], coloring, vertex);
			#if 0
			int vertex_color;
			for (vertex_color = 0; vertex_color < MAXCOLOR; vertex_color ++) {
				if (!forbiddenColors[tx][vertex_color]) {
					coloring[old_vertex] = vertex_color;
					break;
				}
			}
			assert(vertex_color < MAXCOLOR);
			#endif
		}
	//}
}
//*/

__device__ __forceinline__ void conflictDetect1(int src, int dst, int *coloring, bool &is_conflict) {
	if (coloring[src] == coloring[dst] && src < dst) {
		is_conflict = 1;
		coloring[src] = MAXCOLOR;
	}
}

__device__ __forceinline__ bool conflictDetect2(int src, int dst, int *coloring, int *degree, bool &is_conflict) {
	if (coloring[src] == coloring[dst]) {
		bool is_victim;
		if (degree[src] == degree[dst])
			is_victim = (src < dst) ? true : false;
		else is_victim = (degree[src] < degree[dst]) ? true : false;
		if (is_victim) {
			is_conflict = 1;
			coloring[src] = MAXCOLOR;
		}
	}
}
/*
__device__ __forceinline__ unsigned LaneId() {
	unsigned ret;
	asm("mov.u32 %0, %laneid;" : "=r"(ret));
	return ret;
}

__device__ __forceinline__ void resolveByCta(int tid, int m, const int *row_offsets, const int *column_indices, const int * __restrict__ degree, Worklist2 &inwl, int *coloring, bool *is_conflict) {
	__shared__ int owner;
	__shared__ int sh_vertex;
	owner = -1;
	int owner_tid;
	int id = tid;
	int vertex, size = 0;
	if(inwl.pop_id(id, vertex)) {
		//size = degree[vertex];
		size = __ldg(degree + vertex);
	}
	while(true) {
		if(size >= BLKSIZE)
			owner = threadIdx.x;
		__syncthreads();
		owner_tid = owner;
		if(owner_tid == -1)
			break;
		__syncthreads();
		if(owner == threadIdx.x) {
			sh_vertex = vertex;
			// mark this vertex as processed already
			//inwl.dwl[id] = -1;
			cub::ThreadStore<cub::STORE_CG>(inwl.dwl + id, -1);
			owner = -1;
			size = 0;
		}
		__syncthreads();
		int row_begin = row_offsets[sh_vertex];
		//row_begin = __ldg(row_offsets + sh_vertex);
		//row_end = row_offsets[sh_vertex + 1];
		//row_end = __ldg(row_offsets + sh_vertex + 1);
		//int neighbor_size = row_end - row_begin;
		//int neighbor_size = degree[sh_vertex];
		int neighbor_size = __ldg(degree + sh_vertex);
		int num = ((neighbor_size - 1) / BLKSIZE + 1) * BLKSIZE;
		for(int i = threadIdx.x; i < num; i += BLKSIZE) {
			if(i < neighbor_size) {
				int edge = row_begin + i;
				int dst = column_indices[edge];
				//int dst = __ldg(column_indices + edge);
				if(conflictDetect(sh_vertex, dst, coloring, is_conflict[owner_tid]))
					break;
			}
		}
	}
}

#define WARP_SIZE 32
#define LOG_WARP_SIZE 5
#define NUM_WARPS (BLKSIZE / WARP_SIZE)
__device__ __forceinline__ void resolveByWarp(int id, int m, const int *row_offsets, const int *column_indices, const int* __restrict__ degree, Worklist2 &inwl, int *coloring, bool *is_conflict) {
	unsigned warp_id = threadIdx.x >> LOG_WARP_SIZE;
	unsigned lane_id = LaneId();
	__shared__ int owner[NUM_WARPS];
	__shared__ int sh_vertex[NUM_WARPS];
	owner[warp_id] = -1;
	int size = 0;
	int vertex;
	if(inwl.pop_id(id, vertex)) {
		if (vertex != -1) {
			//row_begin = row_offsets[vertex];
			//row_begin = __ldg(row_offsets + vertex);
			//row_end = row_offsets[vertex + 1];
			//row_end = __ldg(row_offsets + vertex + 1);
			//size = row_end - row_begin;
			//size = degree[vertex];
			size = __ldg(degree + vertex);
		}
	}
	while(__any(size) >= WARP_SIZE) {
		if(size >= WARP_SIZE)
			owner[warp_id] = lane_id;
		if(owner[warp_id] == lane_id) {
			sh_vertex[warp_id] = vertex;
			// mark this vertex as processed already
			//inwl.dwl[id] = -1;
			cub::ThreadStore<cub::STORE_CG>(inwl.dwl + id, -1);
			owner[warp_id] = -1;
			size = 0;
		}
		int winner = sh_vertex[warp_id];
		int winner_tid = warp_id * WARP_SIZE + owner[warp_id];
		int row_begin = row_offsets[winner];
		//row_begin = __ldg(row_offsets + sh_vertex[warp_id]);
		//row_end = row_offsets[sh_vertex[warp_id] + 1];
		//row_end = __ldg(row_offsets + sh_vertex[warp_id] + 1);
		//int neighbor_size = row_end - row_begin;
		//int neighbor_size = degree[winner];
		int neighbor_size = __ldg(degree + winner);
		int num = ((neighbor_size + WARP_SIZE - 1) / WARP_SIZE) * WARP_SIZE;
		for(int i = lane_id; i < num; i+= WARP_SIZE) {
			if(i < neighbor_size) {
				int edge = row_begin + i;
				int dst = column_indices[edge];
				//int dst = __ldg(column_indices + edge);
				if(conflictDetect(winner, dst, coloring, is_conflict[winner_tid])) break;
			}
		}
	}
}

__global__ void conflictResolve(int m, const int *row_offsets, const int *column_indices, const int * __restrict__ degree, Worklist2 inwl, Worklist2 outwl, int *coloring) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int tx = threadIdx.x;
	__shared__ bool is_conflict[BLKSIZE];

	typedef cub::BlockScan<int, BLKSIZE> BlockScan;
	__shared__ BlockScan::TempStorage temp_storage;
	__shared__ int gather_offsets[SCRATCHSIZE];
	__shared__ int src[BLKSIZE];
	__shared__ short srcIndex[BLKSIZE];
	int id = tid;
	//int total_inputs = (*inwl.dindex - 1) / (gridDim.x * blockDim.x) + 1;
	//for (int id = tid; total_inputs > 0; id += blockDim.x * gridDim.x, total_inputs--) {
		is_conflict[tx] = 0;
		__syncthreads();
		int old_vertex;
		inwl.pop_id(id, old_vertex);
		resolveByCta(id, m, row_offsets, column_indices, degree, inwl, coloring, is_conflict);
		//resolveByWarp(id, m, row_offsets, column_indices, degree, inwl, coloring, is_conflict);

		gather_offsets[tx] = 0;
		src[tx] = 0;
		srcIndex[tx] = 0;
		__syncthreads();
		int vertex;
		int neighbor_size = 0;
		int row_begin = 0;
		int row_end = 0;
		int scratch_offset = 0;
		int total_edges = 0;
		if (inwl.pop_id(id, vertex)) {
			if (vertex != -1) {
				row_begin = row_offsets[vertex];
				//row_begin = __ldg(row_offsets + vertex);
				row_end = row_offsets[vertex + 1];
				//row_end = __ldg(row_offsets + vertex + 1);
				neighbor_size = row_end - row_begin;
			}
		}
		BlockScan(temp_storage).ExclusiveSum(neighbor_size, scratch_offset, total_edges);
		int done = 0;
		int neighborsdone = 0;

		while (total_edges > 0) {
			__syncthreads();
			int i;
			for (i = 0; (neighborsdone + i) < neighbor_size && (scratch_offset + i - done) < SCRATCHSIZE; i++) {
				int ii = scratch_offset + i - done;
				gather_offsets[ii] = row_begin + neighborsdone + i;
				src[ii] = vertex;
				srcIndex[ii] = tx;
			}
			neighborsdone += i;
			scratch_offset += i;
			__syncthreads();

			if (tx < total_edges) {
				int edge = gather_offsets[tx];
				int dst = column_indices[edge];
				//int dst = __ldg(column_indices + edge);
				int index = srcIndex[tx];
				int srcsrc = src[tx];
				conflictDetect(srcsrc, dst, coloring, is_conflict[index]);
			}
			total_edges -= BLKSIZE;
			done += BLKSIZE;
		}
		__syncthreads();
		outwl.push_1item<BlockScan>((int)is_conflict[tx], old_vertex, BLKSIZE);
	//}
}
*/

__device__ void conflictResolve(int m, const int* __restrict__ csrRowPtr, const int* __restrict__ csrColInd, Worklist2 &inwl, Worklist2 &outwl, int * degree, int *coloring) {
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int total_inputs = (*inwl.dindex + gridDim.x * blockDim.x - 1)/(gridDim.x * blockDim.x);
	for (int id = tid; total_inputs > 0; id += blockDim.x * gridDim.x, total_inputs--) {
		bool is_conflict = 0;
		int vertex;
		if (inwl.pop_id(id, vertex)) {
			int row_begin = __ldg(csrRowPtr + vertex);
			int row_end= __ldg(csrRowPtr + vertex + 1);
			for (int offset = row_begin; offset < row_end; offset ++) {
				int neighbor = __ldg(csrColInd + offset);
				//conflictDetect1(vertex, neighbor, coloring, is_conflict);
				conflictDetect2(vertex, neighbor, coloring, degree, is_conflict);
				if(is_conflict) break;
			}
		}   
		outwl.push_1item<BlockScan>((int)is_conflict, vertex, BLKSIZE);
	}
}

#ifdef FUSION
__global__ void color_kernel(int m, int *csrRowPtr, int *csrColInd, int *degree, Worklist2 inwl, Worklist2 outwl, int *coloring, GlobalBarrier gb) {
	Worklist2 *in;
	Worklist2 *out;
	Worklist2 *tmp; 
	in = &inwl; out = &outwl;
	while (*in->dindex > 0) {
		firstFit(m, csrRowPtr, csrColInd, *in, coloring);
		gb.Sync();
		conflictResolve(m, csrRowPtr, csrColInd, *in, *out, degree, coloring);
		gb.Sync();
		tmp = in; 
		in = out;
		out = tmp;
		*out->dindex = 0;
	}
}
#endif

void color(int m, int nnz, int *csrRowPtr, int *csrColInd, int *coloring, int num_SMs) {
	double starttime, endtime;
  cudaEvent_t cudaStart, cudaEnd;
	double runtime[ITERATIONS];
	int colors[ITERATIONS];
	int iterations[ITERATIONS];
	int *d_csrRowPtr, *d_csrColInd, *d_coloring, *d_degree;
//	printf("Graph coloring data-driven Best version\n");

	int *degree = (int *)malloc(m * sizeof(int));
	for(int i = 0; i < m; i ++) {
		degree[i] = csrRowPtr[i + 1] - csrRowPtr[i];
	}
	for (int i = 0; i < ITERATIONS; i++) {
  cudaEventCreate(&cudaStart);
  cudaEventCreate(&cudaEnd);
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrRowPtr, (m + 1) * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_csrColInd, nnz * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_degree, m * sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void **)&d_coloring, m * sizeof(int)));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrRowPtr, csrRowPtr, (m + 1) * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_csrColInd, csrColInd, nnz * sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL(cudaMemcpy(d_degree, degree, m * sizeof(int), cudaMemcpyHostToDevice));
	cudaDeviceSynchronize();
	cudaDeviceProp deviceProp;
	cudaGetDeviceProperties(&deviceProp, 0);
	int nSM = deviceProp.multiProcessorCount;
	//int nSM = num_SMs;
#ifdef FUSION
	const size_t max_blocks = maximum_residency(color_kernel, BLKSIZE, 0);
	//printf("max_blocks=%d\n", max_blocks);
	GlobalBarrierLifetime gb;
	gb.Setup(nSM * max_blocks);
#else
	const size_t max_blocks_1 = maximum_residency(firstFit, BLKSIZE, 0); 
	const size_t max_blocks_2 = maximum_residency(conflictResolve, BLKSIZE, 0);
	//printf("max_blocks_1=%d, max_blocks_2=%d\n", max_blocks_1, max_blocks_2);
#endif
		Worklist2 inwl(m), outwl(m);
		Worklist2 *inwlptr = &inwl, *outwlptr = &outwl;
		CUDA_SAFE_CALL(cudaMemcpy(inwl.dindex, &m, sizeof(int), cudaMemcpyHostToDevice));
    cudaEventRecord(cudaStart, 0);
		initialize <<<((m - 1) / BLKSIZE + 1), BLKSIZE>>> (d_coloring, m);
		iterations[i] = 0;
		starttime = rtclock();
		int nitems = m;
		thrust::sequence(thrust::device, inwl.dwl, inwl.dwl + m);
#ifdef FUSION
		color_kernel<<<nSM * max_blocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, d_degree, inwl, outwl, d_coloring, gb);
#else
		while (nitems > 0) {
			iterations[i] ++;
			//printf("in_nitems[%d]=%d\n", iteration, nitems);
			int nblocks = (nitems - 1) / BLKSIZE + 1;
			int nblocks_1 = nSM * max_blocks_1;
			int nblocks_2 = nSM * max_blocks_2;
			if(nblocks < nblocks_1) nblocks_1 = nblocks;
			if(nblocks < nblocks_2) nblocks_2 = nblocks;
			//printf("Iteration %d (nitems=%d, nblocks=%d, blocksize=%d).\n", iteration, nitems, nblocks, BLKSIZE);
			firstFit<<<nblocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, d_degree, *inwlptr, d_coloring);
			conflictResolve<<<nblocks, BLKSIZE>>>(m, d_csrRowPtr, d_csrColInd, d_degree, *inwlptr, *outwlptr, d_coloring);
			nitems = outwlptr->nitems();
			Worklist2 * tmp = inwlptr;
			inwlptr = outwlptr;
			outwlptr = tmp;
			outwlptr->reset();
		}
#endif
		CUDA_SAFE_CALL(cudaDeviceSynchronize());
		endtime = rtclock();
    cudaEventRecord(cudaEnd, 0);
		colors[i] = thrust::reduce(thrust::device, d_coloring, d_coloring + m, 0, thrust::maximum<int>()) + 1;
	  CUDA_SAFE_CALL(cudaMemcpy(coloring, d_coloring, m * sizeof(int), cudaMemcpyDeviceToHost));
	  CUDA_SAFE_CALL(cudaFree(d_csrRowPtr));
	  CUDA_SAFE_CALL(cudaFree(d_csrColInd));
	  CUDA_SAFE_CALL(cudaFree(d_coloring));
    CUDA_SAFE_CALL(cudaFree(d_degree));
    cudaDeviceSynchronize();
    float time;
    cudaEventElapsedTime(&time, cudaStart, cudaEnd);
		runtime[i] = time/*1000.0f * (endtime - starttime);*/ ;
    cudaEventDestroy(cudaStart);
    cudaEventDestroy(cudaEnd);
	}
	double total_runtime = 0.0;
	int total_colors = 0;
	int total_iterations = 0;
	for (int i = 0; i < ITERATIONS; i++) {
		total_runtime += runtime[i];
		total_colors += colors[i];
		total_iterations += iterations[i];
//		printf("[%d %.2f %d] ", colors[i], runtime[i], iterations[i]);
	}
  double avg_runtime = (double)total_runtime / ITERATIONS;
	double avg_colors = (double)total_colors / ITERATIONS;
	double avg_iterations = (double)total_iterations / ITERATIONS;
  cout << avg_runtime << " " << avg_colors << endl;
//	printf("\navg_runtime %f ms, avg_colors %.2f avg_iterations %.2f\n", avg_runtime, avg_colors, avg_iterations);
}
