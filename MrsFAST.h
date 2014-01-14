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
 */

#ifndef __MRS_FAST__
#define __MRS_FAST__

#include "Reads.h"

#define MAP_CHUNKS 15

enum BestMappingStatus { unset = 0, first_mate, second_mate, trans_loc, improper, proper };

typedef struct mn
{
	int loc;
	char dir;
	char err;
	char mderr;
	float score;
	int hits;
	int secondBestHits;
	int secondBestErrors;
	char md[40];
	char chr[40];
} FullMappingInfo;


typedef struct mnp
{
	int loc1;
	char dir1;
	char err1;
	char mderr1;
	char md1[40];
	char chr1[40];
	int loc2;
	char dir2;
	char err2;
	char mderr2;
	char md2[40];
	char chr2[40];
	enum BestMappingStatus status;
} BestMappingInfoPE;



typedef struct lc
{
	int loc[MAP_CHUNKS];
	char err[MAP_CHUNKS];
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

void initializeFAST(int seqListSize);
void initFASTChunk(Read *seqList, int seqListSize);
void initFASTContig();
void finalizeFAST();

void mapSeq(unsigned char contigFlag);
#endif
