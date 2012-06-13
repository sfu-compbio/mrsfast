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

extern long long			verificationCnt;
extern long long			mappingCnt;
extern long long			mappedSeqCnt;
extern long long			completedSeqCnt;
extern long long			regTotal;
extern long long			metTotal;


void initFAST(int *samplingLocs, int samplingLocsSize);
void finalizeFAST();
int mapSingleEndSeq(char *seqName, char *seq, char* seqQual, unsigned char seqHits, int seqNo, short mappingDirection);
int mapPairedEndSeq(char *seq1Name, char *seq1, char* seq1Qual, unsigned int seq1Hits, int seq1No, short mappingDirection1,
					char *seq2Name, char *seq2, char* seq2Qual, unsigned int seq2Hits, int seq2No, short mappingDirection2);
int mapSingleEndSeqBS(  char *seqName, char *seq, char* seqQual, unsigned int seqHits, int seqNo, short mappingDirection, int type);
int mapPairedEndSeqBS(	char *seq1Name, char *seq1, char* seq1Qual, unsigned int seq1Hits, int seq1No, short mappingDirection1, int type1,
						char *seq2Name, char *seq2, char* seq2Qual, unsigned int seq2Hits, int seq2No, short mappingDirection2, int type2);
#endif
