/******************************************************************************
 * 
 * Copyright 2010 Duane Merrill
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
 * Common B40C Routines 
 ******************************************************************************/

#pragma once

namespace b40c {
namespace util {


/******************************************************************************
 * Macro utilities
 ******************************************************************************/

/**
 * Select maximum
 */
#define B40C_MAX(a, b) ((a > b) ? a : b)


/**
 * Select maximum
 */
#define B40C_MIN(a, b) ((a < b) ? a : b)

/**
 * Return the size in quad-words of a number of bytes
 */
#define B40C_QUADS(bytes) (((bytes + sizeof(uint4) - 1) / sizeof(uint4)))

/******************************************************************************
 * Simple templated utilities
 ******************************************************************************/

/**
 * Supress warnings for unused constants
 */
template <typename T>
__host__ __device__ __forceinline__ void SuppressUnusedConstantWarning(const T) {}


/**
 * Perform a swap
 */
template <typename T> 
void __host__ __device__ __forceinline__ Swap(T &a, T &b) {
	T temp = a;
	a = b;
	b = temp;
}


template <typename K, int magnitude, bool shift_left> struct MagnitudeShiftOp;

/**
 * MagnitudeShift().  Allows you to shift left for positive magnitude values, 
 * right for negative.   
 * 
 * N.B. This code is a little strange; we are using this meta-programming 
 * pattern of partial template specialization for structures in order to 
 * decide whether to shift left or right.  Normally we would just use a 
 * conditional to decide if something was negative or not and then shift 
 * accordingly, knowing that the compiler will elide the untaken branch, 
 * i.e., the out-of-bounds shift during dead code elimination. However, 
 * the pass for bounds-checking shifts seems to happen before the DCE 
 * phase, which results in a an unsightly number of compiler warnings, so 
 * we force the issue earlier using structural template specialization.
 */
template <typename K, int magnitude> 
__device__ __forceinline__ K MagnitudeShift(K key)
{
	return MagnitudeShiftOp<K, (magnitude > 0) ? magnitude : magnitude * -1, (magnitude > 0)>::Shift(key);
}

template <typename K, int magnitude>
struct MagnitudeShiftOp<K, magnitude, true>
{
	__device__ __forceinline__ static K Shift(K key)
	{
		return key << magnitude;
	}
};

template <typename K, int magnitude>
struct MagnitudeShiftOp<K, magnitude, false>
{
	__device__ __forceinline__ static K Shift(K key)
	{
		return key >> magnitude;
	}
};


/******************************************************************************
 * Metaprogramming Utilities
 ******************************************************************************/

/**
 * Null type
 */
struct NullType {};


/**
 * Int2Type
 */
template <int N>
struct Int2Type
{
	enum {VALUE = N};
};


/**
 * Statically determine log2(N), rounded up, e.g.,
 * 		Log2<8>::VALUE == 3
 * 		Log2<3>::VALUE == 2
 */
template <int N, int CURRENT_VAL = N, int COUNT = 0>
struct Log2
{
	// Inductive case
	static const int VALUE = Log2<N, (CURRENT_VAL >> 1), COUNT + 1>::VALUE;
};

template <int N, int COUNT>
struct Log2<N, 0, COUNT>
{
	// Base case
	static const int VALUE = (1 << (COUNT - 1) < N) ?
		COUNT :
		COUNT - 1;
};


/**
 * If/Then/Else
 */
template <bool IF, typename ThenType, typename ElseType>
struct If
{
	// true
	typedef ThenType Type;
};

template <typename ThenType, typename ElseType>
struct If<false, ThenType, ElseType>
{
	// false
	typedef ElseType Type;
};


/**
 * Equals 
 */
template <typename A, typename B>
struct Equals
{
	enum { VALUE = 0 };
};

template <typename A>
struct Equals <A, A>
{
	enum { VALUE = 1 };
};


/**
 * Is volatile
 */
template <typename Tp>
struct IsVolatile
{
	enum { VALUE = 0 };
};
template <typename Tp>
struct IsVolatile<Tp volatile>
{
	enum { VALUE = 1 };
};


/**
 * Removes pointers
 */
template <typename Tp, typename Up>
struct RemovePointersHelper
{
	typedef Tp Type;
};
template <typename Tp, typename Up>
struct RemovePointersHelper<Tp, Up*>
{
	typedef typename RemovePointersHelper<Up, Up>::Type Type;
};
template <typename Tp>
struct RemovePointers : RemovePointersHelper<Tp, Tp> {};



} // namespace util
} // namespace b40c

