/******************************************************************************
 * 
 * Copyright 2010-2011 Duane Merrill
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. 
 * 
 * For more information, see our Google Code project site: 
 * http://code.google.com/p/back40computing/
 * 
 * Thanks!
 * 
 ******************************************************************************/


/******************************************************************************
 * Downsweep CTA processing abstraction
 ******************************************************************************/

#pragma once

#include <b40c/util/basic_utils.cuh>
#include <b40c/partition/downsweep/cta.cuh>

#include <b40c/graph/bfs/partition_compact/downsweep/tile.cuh>

namespace b40c {
namespace graph {
namespace bfs {
namespace partition_compact {
namespace downsweep {


/**
 * CTA
 *
 * Derives from partition::downsweep::Cta
 */
template <typename KernelPolicy>
struct Cta :
	partition::downsweep::Cta<
		KernelPolicy,
		Cta<KernelPolicy>,			// This class
		Tile>						// bfs::partition_compact::downsweep::Tile
{
	//---------------------------------------------------------------------
	// Typedefs and Constants
	//---------------------------------------------------------------------

	// Base class type
	typedef partition::downsweep::Cta<KernelPolicy, Cta, Tile> Base;

	typedef typename KernelPolicy::VertexId 				VertexId;
	typedef typename KernelPolicy::ValidFlag				ValidFlag;
	typedef typename KernelPolicy::SizeT 					SizeT;
	typedef typename KernelPolicy::SmemStorage				SmemStorage;
	typedef typename KernelPolicy::Grid::LanePartial		LanePartial;

	typedef typename KernelPolicy::KeyType 					KeyType;
	typedef typename KernelPolicy::ValueType				ValueType;


	//---------------------------------------------------------------------
	// Members
	//---------------------------------------------------------------------

	// Validity flags
	ValidFlag 			*&d_flags_in;

	// Number of GPUs to partition the frontier into
	int num_gpus;

	//---------------------------------------------------------------------
	// Methods
	//---------------------------------------------------------------------

	/**
	 * Constructor
	 */
	__device__ __forceinline__ Cta(
		SmemStorage 	&smem_storage,
		int 			num_gpus,
		VertexId 		*&d_in,
		VertexId 		*&d_out,
		VertexId 		*&d_parent_in,
		VertexId 		*&d_parent_out,
		ValidFlag		*&d_flags_in,
		SizeT 			*&d_spine,
		LanePartial		base_composite_counter,
		int				*raking_segment) :
			Base(
				smem_storage,
				d_in,							// d_in_keys
				d_out,							// d_out_keys
				(ValueType *&) d_parent_in,		// d_in_values
				(ValueType *&) d_parent_out,	// d_out_values
				d_spine,
				base_composite_counter,
				raking_segment),
			d_flags_in(d_flags_in),
			num_gpus(num_gpus)
	{}
};


} // namespace downsweep
} // namespace partition_compact
} // namespace bfs
} // namespace graph
} // namespace b40c

