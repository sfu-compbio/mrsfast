/*
 * Copyright (c) <2008 - 2009>, University of Washington, Simon Fraser University
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


#ifndef __MRS_FAST__
#define __MRS_FAST__

#include "Reads.h"

#define MAP_CHUNKS 15

// Pair is used to pre-processing and making the read index table
typedef struct
{
	int hv;
	int seqInfo;
} Pair;

typedef struct
{
	int hv;
	unsigned int *seqInfo;
} ReadIndexTable;


typedef struct mn
{
	int loc;
	char dir;
	char err;
	float score;
	char md[10];
	char chr[10];
} FullMappingInfo;

typedef struct lc
{
	int loc[MAP_CHUNKS];
	struct lc *next;
} MappingLocations;

typedef struct inf
{
	int size;
	MappingLocations *next;
} MappingInfo;

typedef struct 
{
	FILE * fp;
	char name[400];
} FILE_STRUCT;

extern long long			verificationCnt;
extern long long			mappingCnt;
extern long long			mappedSeqCnt;
extern long long			completedSeqCnt;

void initFAST(	Read *seqList,
				int seqListSize,
				int *samplingLocs,
				int samplingLocsSize, 
				char *fileName);

void finalizeFAST();

int mapSingleEndSeq();
int mapPaiedEndSeq();
void outputPairedEnd();
void outputPairedEndDiscPP();
#endif
