/*
 * Copyright (c) <2008 - 2020>, University of Washington, Simon Fraser University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *   
 * Redistributions of source code must retain the above copyright notice, this list
 * of conditions and the following disclaimer.
 * - Redistributions in binary form must reproduce the above copyright notice, this
 *   list of conditions and the following disclaimer in the documentation and/or other
 *   materials provided with the distribution.
 * - Neither the name of the <ORGANIZATION> nor the names of its contributors may be
 *   used to endorse or promote products derived from this software without specific
 *   prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
 * Author: 
 *        Faraz Hach (fhach AT cs DOT sfu DOT ca)
 *        Iman Sarrafi (isarrafi AT cs DOT sfu DOT ca)
 *        Ermin Hodzic (ermin_hodzic AT sfu DOT ca)
*/

#include <stdio.h>
#include "Common.h"
#include "Reads.h"
/**********************************************/
static inline int lg( int x ) {
	// floor( log2( x ) )
	int r = 0;
	while ( x > 1 ) r++, x >>= 1;
	return r;
}
/**********************************************/
void heapSortGI( GeneralIndex * A, int N ) {
	// N - size of the array to be sorted
	if ( N <= 1 ) return;
	int i, j, maxi;
	// We make the array a max heap using "bubble-down" operation
	for ( i = N/2 - 1; i >= 0; i-- ) {
		j = i;
		// bubble-down
		while ( 2 * j + 2 < N ) {
			if ( A[ 2 * j + 1 ].checksum > A[ 2 * j + 2 ].checksum ) maxi = 2 * j + 1;
			else maxi = 2 * j + 2;
			if ( A[ maxi ].checksum > A[ j ].checksum ) {
				GeneralIndex temp = A[ j ];
				A[ j ] = A[ maxi ];
				A[ maxi ] = temp;
				j = maxi;
			}
			else break;
		}
		if ( 2 * j + 1 < N && A[ 2 * j + 1 ].checksum > A[ j ].checksum ) {
			GeneralIndex temp = A[ j ];
			A[ j ] = A[ 2 * j + 1 ];
			A[ 2 * j + 1 ] = temp;
		}
	}
	// Extracting maximum and moving it to the end of the array as long as size is greater than 1
	while ( --N ) {
		GeneralIndex temp = A[ 0 ];
		A[ 0 ] = A[ N ];
		A[ N ] = temp;
		j = 0;
		// bubble-down
		while ( 2 * j + 2 < N ) {
			if ( A[ 2 * j + 1 ].checksum > A[ 2 * j + 2 ].checksum ) maxi = 2 * j + 1;
			else maxi = 2 * j + 2;
			if ( A[ maxi ].checksum > A[ j ].checksum ) {
				GeneralIndex temp = A[ j ];
				A[ j ] = A[ maxi ];
				A[ maxi ] = temp;
				j = maxi;
			}
			else break;
		}
		if ( 2 * j + 1 < N && A[ 2 * j + 1 ].checksum > A[ j ].checksum ) {
			GeneralIndex temp = A[ j ];
			A[ j ] = A[ 2 * j + 1 ];
			A[ 2 * j + 1 ] = temp;
		}
	}
}
/**********************************************/
void insertionSortGI( GeneralIndex * A, const int left, const int right ) {
	// left, right - starting and ending index of the interval to be sorted
	int i, j;
	for ( i = left + 1; i <= right; i++ ) {
		j = i;
		GeneralIndex temp = A[ i ];
		while ( j > left && A[ j - 1 ].checksum > temp.checksum ) {
			A[ j ] = A[ j - 1 ];
			j--;
		}
		A[ j ] = temp;
	}
}
/**********************************************/
void quickSortGI( GeneralIndex * A, const int left, const int right, int depth ) {
	// left, right - starting and ending index of the interval to be sorted
	// depth - max recursion depth allowed
	if ( right - left <= 16 ) return;	// leave it for insertionSort
	
	if ( depth == 0 ) heapSortGI( A + left, right - left + 1 );	// too many recursive calls, switch to heapsort to ensure O(nlogn) time
	else {
		int mid = ( left + right ) / 2;
		int small = left - 1;
		int i;
		GeneralIndex temp;
		// median of A[ left ], A[ mid ], A[ right ], to be stored in A[ right ]
		if ( A[ mid ].checksum < A[ right ].checksum ) {
			temp = A[ right ];
			A[ right ] = A[ mid ];
			A[ mid ] = temp;
		}
		if ( A[ left ].checksum < A[ right ].checksum ) {
			temp = A[ right ];
			A[ right ] = A[ left ];
			A[ left ] = temp;
		}
		if ( A[ left ].checksum < A[ mid ].checksum ) {
			temp = A[ right ];
			A[ right ] = A[ left ];
			A[ left ] = temp;
		}
		else {
			temp = A[ right ];
			A[ right ] = A[ mid ];
			A[ mid ] = temp;
		}
		// partitioning the array with respect to A[ right ]
		for ( i = left; i < right; i++ ) {
			if ( A[ i ].checksum < A[ right ].checksum ) {
				temp = A[ ++small ];
				A[ small ] = A[ i ];
				A[ i ] = temp;
			}
		}
		temp = A[ small + 1 ];
		A[ small + 1 ] = A[ right ];
		A[ right ] = temp;
		// recursion
		quickSortGI( A, left, small, depth - 1 );
		quickSortGI( A, small + 2, right, depth - 1 );
	}
}
/**********************************************/
void introSortGI( GeneralIndex * A, const int left, const int right ) {
	// main sort, limiting recursion depth to 2 times the logarithm of its length
	quickSortGI( A, left, right, 2 * lg( right - left + 1 ) );
	// finish-up small unsorted pieces
	insertionSortGI( A, left, right );
}
/**********************************************/
void heapSortPair( Pair * A, int N ) {
	// N - size of the array to be sorted
	if ( N <= 1 ) return;
	int i, j, maxi;
	// We make the array a max heap using "bubble-down" operation
	for ( i = N/2 - 1; i >= 0; i-- ) {
		j = i;
		// bubble-down
		while ( 2 * j + 2 < N ) {
			if ( A[ 2 * j + 1 ].hv > A[ 2 * j + 2 ].hv || (A[ 2 * j + 1 ].hv == A[ 2 * j + 2 ].hv && A[ 2 * j + 1 ].checksum > A[ 2 * j + 2 ].checksum)  ) maxi = 2 * j + 1;
			else maxi = 2 * j + 2;
			if ( A[ maxi ].hv > A[ j ].hv ||  (A[ maxi ].hv == A[ j ].hv && A[ maxi ].checksum > A[ j ].checksum)  ) {
				Pair temp = A[ j ];
				A[ j ] = A[ maxi ];
				A[ maxi ] = temp;
				j = maxi;
			}
			else break;
		}
		if ( 2 * j + 1 < N && ( A[ 2 * j + 1 ].hv > A[ j ].hv || (A[ 2 * j + 1 ].hv == A[ j ].hv && A[ 2 * j + 1 ].checksum > A[ j ].checksum) ) ) {
			Pair temp = A[ j ];
			A[ j ] = A[ 2 * j + 1 ];
			A[ 2 * j + 1 ] = temp;
		}
	}
	// Extracting maximum and moving it to the end of the array as long as size is greater than 1
	while ( --N ) {
		Pair temp = A[ 0 ];
		A[ 0 ] = A[ N ];
		A[ N ] = temp;
		j = 0;
		// bubble-down
		while ( 2 * j + 2 < N ) {
			if ( A[ 2 * j + 1 ].hv > A[ 2 * j + 2 ].hv || (A[ 2 * j + 1 ].hv == A[ 2 * j + 2 ].hv && A[ 2 * j + 1 ].checksum > A[ 2 * j + 2 ].checksum) ) maxi = 2 * j + 1;
			else maxi = 2 * j + 2;
			if ( A[ maxi ].hv > A[ j ].hv || (A[ maxi ].hv == A[ j ].hv && A[ maxi ].checksum > A[ j ].checksum) ) {
				Pair temp = A[ j ];
				A[ j ] = A[ maxi ];
				A[ maxi ] = temp;
				j = maxi;
			}
			else break;
		}
		if ( 2 * j + 1 < N && ( A[ 2 * j + 1 ].hv > A[ j ].hv || (A[ 2 * j + 1 ].hv == A[ j ].hv && A[ 2 * j + 1 ].checksum > A[ j ].checksum) ) ) {
			Pair temp = A[ j ];
			A[ j ] = A[ 2 * j + 1 ];
			A[ 2 * j + 1 ] = temp;
		}
	}
}
/**********************************************/
void insertionSortPair( Pair * A, const int left, const int right ) {
	// left, right - starting and ending index of the interval to be sorted
	int i, j;
	for ( i = left + 1; i <= right; i++ ) {
		j = i;
		Pair temp = A[ i ];
		while ( j > left && ( A[ j - 1 ].hv > temp.hv || (A[ j - 1 ].hv == temp.hv && A[ j - 1 ].checksum > temp.checksum) ) ) {
			A[ j ] = A[ j - 1 ];
			j--;
		}
		A[ j ] = temp;
	}
}
/**********************************************/
void quickSortPair( Pair * A, const int left, const int right, int depth ) {
	// left, right - starting and ending index of the interval to be sorted
	// depth - max recursion depth allowed
	if ( right - left <= 16 ) return;	// leave it for insertionSort
	
	if ( depth == 0 ) heapSortPair( A + left, right - left + 1 );	// too many recursive calls, switch to heapsort to ensure O(nlogn) time
	else {
		int mid = ( left + right ) / 2;
		int small = left - 1;
		int i;
		Pair temp;
		// median of A[ left ], A[ mid ], A[ right ], to be stored in A[ right ]
		if ( A[ mid ].hv < A[ right ].hv || (A[ mid ].hv == A[ right ].hv && A[ mid ].checksum < A[ right ].checksum) ) {
			temp = A[ right ];
			A[ right ] = A[ mid ];
			A[ mid ] = temp;
		}
		if ( A[ left ].hv < A[ right ].hv || (A[ left ].hv == A[ right ].hv && A[ left ].checksum < A[ right ].checksum) ) {
			temp = A[ right ];
			A[ right ] = A[ left ];
			A[ left ] = temp;
		}
		if ( A[ left ].hv < A[ mid ].hv || (A[ left ].hv == A[ mid ].hv && A[ left ].checksum < A[ mid ].checksum) ) {
			temp = A[ right ];
			A[ right ] = A[ left ];
			A[ left ] = temp;
		}
		else {
			temp = A[ right ];
			A[ right ] = A[ mid ];
			A[ mid ] = temp;
		}
		// partitioning the array with respect to A[ right ]
		for ( i = left; i < right; i++ ) {
			if ( A[ i ].hv < A[ right ].hv || (A[ i ].hv == A[ right ].hv && A[ i ].checksum < A[ right ].checksum) ) {
				temp = A[ ++small ];
				A[ small ] = A[ i ];
				A[ i ] = temp;
			}
		}
		temp = A[ small + 1 ];
		A[ small + 1 ] = A[ right ];
		A[ right ] = temp;
		// recursion
		quickSortPair( A, left, small, depth - 1 );
		quickSortPair( A, small + 2, right, depth - 1 );
	}
}
/**********************************************/
void introSortPair( Pair * A, const int left, const int right ) {
	// main sort, limiting recursion depth to 2 times the logarithm of its length
	quickSortPair( A, left, right, 2 * lg( right - left + 1 ) );
	// finish-up small unsorted pieces
	insertionSortPair( A, left, right );
}

