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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"
#include "Reads.h"
#include "HashTable.h"
#include "Output.h"
#include "MrsFAST.h"
#include "RefGenome.h"

#ifdef MRSFAST_SSE4
#include <smmintrin.h>
#include <nmmintrin.h>
#endif


float calculateScore(int index, CompressedSeq *cmpSeq, char *qual, int *err);
unsigned char		mrFAST = 0;
char				*versionNumberF="0.0";

long long			verificationCnt = 0;
long long 			mappingCnt = 0;
long long			mappedSeqCnt = 0;
long long			completedSeqCnt = 0;
char				*mappingOutput;
/**********************************************/
int					_msf_refGenLength = 0;	
CompressedSeq		*_msf_crefGen = NULL;
int					_msf_crefGenLen = 0;
int					_msf_refGenOffset = 0;
char				*_msf_refGenName = NULL;
char				*_msf_tempSeq = NULL;

int					_msf_refGenBeg;
int					_msf_refGenEnd;

IHashTable			*_msf_hashTable = NULL;

int					*_msf_samplingLocs;
int					*_msf_samplingLocsEnds;
int					_msf_samplingLocsSize;
int					*_msf_samplingLocsSeg;
int					*_msf_samplingLocsOffset;
int 				*_msf_samplingLocsLen;
int					*_msf_samplingLocsLenFull;

Read				*_msf_seqList;
int					_msf_seqListSize;

ReadIndexTable		*_msf_rIndex = NULL;
int					_msf_rIndexSize;
int					_msf_rIndexMax;

SAM					_msf_output;

OPT_FIELDS			*_msf_optionalFields;

char				*_msf_op;

char				_msf_numbers[200][3];
char				_msf_cigar[5];

MappingInfo			*_msf_mappingInfo;

int					*_msf_seqHits;
int					_msf_openFiles = 0;
int					_msf_maxLSize=0;
int					_msf_maxRSize=0;
int 				typeSize = 63; //32
char 				*errCnt;
/**********************************************/
int compare (const void *a, const void *b)
{
	return ((Pair *)a)->hv - ((Pair *)b)->hv;
}
/**********************************************/
void initializeMasksErrCntShrink()
{
	errCnt = (char *)getMem(1 << 24);

	int i,x;

	for (i = 0; i<(1<<24); i++)
	{
		errCnt[i] = 0;
		for (x = 0; x < 8; x++)
		{
			if (i & (7 << 3*x))
			{
				errCnt[i]++;
			}
		}
	}
}
/**********************************************/
void preProcessReads()
{
	int i=0;
	int j=0;
	int pos = 0;
	char rseq[SEQ_LENGTH+1];

	_msf_rIndexMax = -1;

	int tmpSize = _msf_seqListSize*_msf_samplingLocsSize*2;
	Pair *tmp = getMem(sizeof(Pair)*tmpSize);

	for (i=0; i<_msf_seqListSize; i++)
	{
		for (j=0; j< _msf_samplingLocsSize; j++)
		{
			tmp[pos].hv = hashVal(_msf_seqList[i].seq+_msf_samplingLocs[j]);
			tmp[pos].seqInfo = pos;
			pos++;
		}
		for (j=0; j<_msf_samplingLocsSize; j++)
		{
			reverseComplete(_msf_seqList[i].seq, rseq, SEQ_LENGTH);
			tmp[pos].hv = hashVal(rseq+_msf_samplingLocs[j]);
			tmp[pos].seqInfo = pos;
			pos++;
		}
	}

	qsort(tmp, tmpSize, sizeof(Pair), compare);


	int uniq = 0;
	int prev = -2;
	int beg = -1;
	int end = -1;

	for (i=0; i<tmpSize; i++)
	{
		if (prev != tmp[i].hv)
		{
			uniq ++;
			prev = tmp[i].hv;
		}
	}

	_msf_rIndexSize = uniq;
	_msf_rIndex = getMem(sizeof(ReadIndexTable)*_msf_rIndexSize);
	prev = -2;

	j=0;
	beg =0;
	while (beg < tmpSize)
	{
		end = beg;
		while (end+1<tmpSize && tmp[end+1].hv==tmp[beg].hv)
			end++;

		_msf_rIndex[j].hv = tmp[beg].hv;
		_msf_rIndex[j].seqInfo = getMem(sizeof(int)*(end-beg+2));
		_msf_rIndex[j].seqInfo[0] = end-beg+1;
		if ((end-beg+1) > _msf_rIndexMax)
			_msf_rIndexMax = end-beg+1;

		for (i=1; i<=_msf_rIndex[j].seqInfo[0]; i++)
		{
			_msf_rIndex[j].seqInfo[i]=tmp[beg+i-1].seqInfo;
		}
		j++;
		beg = end+1;
	}
	freeMem(tmp, sizeof(Pair)*tmpSize);
}
/**********************************************/
void initFAST(Read *seqList, int seqListSize, int *samplingLocs, int samplingLocsSize, char *genFileName)
{
	if (_msf_tempSeq == NULL)
	{
		_msf_tempSeq = getMem(SEQ_LENGTH+1);
	}

	if (_msf_optionalFields == NULL)
	{	
		_msf_op = getMem(SEQ_LENGTH);
		if (pairedEndMode)
		{
			_msf_optionalFields = getMem(4*sizeof(OPT_FIELDS));
		}
		else
		{
			_msf_optionalFields = getMem(2*sizeof(OPT_FIELDS));
		}

		int i;
		for (i=0; i<200;i++)
		{
			sprintf(_msf_numbers[i],"%d%c",i, '\0');
		}
		sprintf(_msf_cigar, "%dM", SEQ_LENGTH);
	}

	if (_msf_samplingLocsEnds == NULL)
	{
		int i;
		_msf_samplingLocs = samplingLocs;
		_msf_samplingLocsSize = samplingLocsSize;

		_msf_samplingLocsEnds = malloc(sizeof(int)*_msf_samplingLocsSize);
		_msf_samplingLocsSeg = malloc(sizeof(int)*_msf_samplingLocsSize);
		_msf_samplingLocsOffset = malloc(sizeof(int)*_msf_samplingLocsSize);
		_msf_samplingLocsLen = malloc(sizeof(int)*_msf_samplingLocsSize);
		_msf_samplingLocsLenFull = malloc(sizeof(int)*_msf_samplingLocsSize);
		for (i=0; i<_msf_samplingLocsSize; i++)
		{
			_msf_samplingLocsEnds[i]	=_msf_samplingLocs[i]+WINDOW_SIZE-1;
			_msf_samplingLocsSeg[i]		= _msf_samplingLocs[i]/(sizeof(CompressedSeq)*8/3);
			_msf_samplingLocsOffset[i]	=_msf_samplingLocs[i]%(sizeof(CompressedSeq)*8/3);
			_msf_samplingLocsLen[i]		= _msf_samplingLocs[i+1]-_msf_samplingLocs[i];
			_msf_samplingLocsLenFull[i]	= SEQ_LENGTH - _msf_samplingLocs[i];
		}

		_msf_seqList = seqList;
		_msf_seqListSize = seqListSize;

		initializeMasksErrCntShrink();		// TODO: ifdef
		preProcessReads();
	}
	if (_msf_refGenName == NULL)
	{
		_msf_refGenName = getMem(SEQ_LENGTH);
	}
	_msf_refGenLength = getRefGenLength();
	_msf_crefGen = getCmpRefGenome();
	_msf_crefGenLen = getCmpRefGenLen();
	_msf_refGenOffset = getRefGenomeOffset();
	sprintf(_msf_refGenName,"%s%c", getRefGenomeName(), '\0');

	if (pairedEndMode && _msf_seqHits == NULL)
	{
		_msf_mappingInfo  = getMem(seqListSize * sizeof (MappingInfo));

		int i=0;
		for (i=0; i<seqListSize; i++)
		{
			//_msf_mappingInfo[i].next = getMem(sizeof(MappingLocations));	DEL
			_msf_mappingInfo[i].next = NULL;
			_msf_mappingInfo[i].size = 0;
		}

		_msf_seqHits = getMem((_msf_seqListSize/2) * sizeof(int));


		for (i=0; i<_msf_seqListSize/2; i++)
		{
			_msf_seqHits[i] = 0;
		}

//		initLoadingRefGenome(genFileName, NULL, &i);						DEL
	}

	if (_msf_refGenOffset == 0)
	{
		_msf_refGenBeg = 1;
	}
	else
	{
		_msf_refGenBeg = CONTIG_OVERLAP - SEQ_LENGTH + 2;
	}
	_msf_refGenEnd = _msf_refGenLength - SEQ_LENGTH + 1;

}
/**********************************************/
void finalizeFAST()
{
	freeMem(_msf_seqHits, (_msf_seqListSize/2) * sizeof(int));
	freeMem(_msf_refGenName, SEQ_LENGTH);
	int i;
	for (i=0; i<_msf_rIndexSize; i++)
	{
		freeMem(_msf_rIndex[i].seqInfo, _msf_rIndex[i].seqInfo[0]+1);
	}
	freeMem(_msf_rIndex, _msf_rIndexSize);
	freeMem(_msf_tempSeq, SEQ_LENGTH+1);
	freeMem(errCnt, 1<<24);
}
/**********************************************/
inline int countErrors(CompressedSeq *ref, int refOff, CompressedSeq *seq, int seqOff, int len, int *errSamp)
{
	CompressedSeq NMASK = 0x4924924924924924;

	int refALS = refOff * 3;
	int segALS = seqOff * 3;


	int refARS = typeSize - refALS;
	int segARS = typeSize - segALS;	
	
	CompressedSeq tmpref, tmpseq, diff;		// k for ref, and k for seg		DEL comment
	int err = 0;


	while(len >= 21)
	{
		tmpref = (*ref << refALS) | (*(++ref) >> refARS);
		tmpseq = (*seq << segALS) | (*(++seq) >> segARS);
		diff = (tmpref ^ tmpseq) & 0x7fffffffffffffff;

		*errSamp |= (tmpseq & NMASK);

#ifdef MRSFAST_SSE4
		err += _mm_popcnt_u64(((diff >> 1) | (diff >> 2) | diff ) &  0x9249249249249249);
#else
		err += errCnt[diff & 0xffffff] + errCnt[(diff>>24)&0xffffff] + errCnt[(diff>>48)&0xfffff];
#endif

		if (err > errThreshold)
			return errThreshold+1;
		len -= 21;
	}

	if (len)
	{
		tmpref = (*ref << refALS) | (*(++ref) >> refARS);
		tmpseq = (*seq << segALS) | (*(++seq) >> segARS);
		diff = (tmpref ^ tmpseq) & 0x7fffffffffffffff;
		
		diff >>= (typeSize - len*3);
		tmpseq  >>= (typeSize - len*3);
		
		*errSamp |= (tmpseq & NMASK);

#ifdef MRSFAST_SSE4
		err += _mm_popcnt_u64(((diff >> 1) | (diff >> 2) | diff ) &  0x9249249249249249);
#else
		err += errCnt[diff & 0xffffff] + errCnt[(diff>>24)&0xffffff] + errCnt[(diff>>48)&0xfffff];
#endif

		
		if (err > errThreshold)
			return errThreshold+1;
	}

	*errSamp |= err;
	return err;
}
/**********************IMAN************************/
inline int verifySingleEnd(int index, CompressedSeq *seq, int offset)
{
	int segLen, cmpSegLen, curOff, sampleErrors=0, err = 0, refLoc, segLoc;
	index--;

	CompressedSeq *refSeg = _msf_crefGen+index/21;
	int refOff = index % 21;

	CompressedSeq *refCurSeg = refSeg+_msf_samplingLocsSeg[offset];

	int refCurOff=refOff+_msf_samplingLocsOffset[offset];
	if (refCurOff>=21)
	{
		refCurSeg++;
		refCurOff-=21;
	}

	err = countErrors(refCurSeg, refCurOff, seq+_msf_samplingLocsSeg[offset], _msf_samplingLocsOffset[offset], _msf_samplingLocsLen[offset], &sampleErrors);
	if (sampleErrors || err)
		return -1;

	verificationCnt++;
	err = 0; 

	refCurSeg = refSeg;
	refCurOff = refOff;

	for (curOff = 0; curOff < offset; curOff++)
	{
		sampleErrors=0;

		refCurSeg = refSeg+_msf_samplingLocsSeg[curOff];
		refCurOff = refOff+_msf_samplingLocsOffset[curOff];
		if(refCurOff>=21)
		{
			refCurSeg++;
			refCurOff-=21;
		}

		err += countErrors(refCurSeg, refCurOff, seq+_msf_samplingLocsSeg[curOff], _msf_samplingLocsOffset[curOff], _msf_samplingLocsLen[curOff], &sampleErrors);
		if (err > errThreshold || sampleErrors==0)
			return -1;
	}

	if (offset != _msf_samplingLocsSize-1)
	{
		offset++;
		refCurSeg = refSeg+_msf_samplingLocsSeg[offset];
		refCurOff = refOff+_msf_samplingLocsOffset[offset];
		if(refCurOff>=21)
		{
			refCurSeg++;
			refCurOff-=21;
		}
		err += countErrors(refCurSeg, refCurOff, seq+_msf_samplingLocsSeg[offset], _msf_samplingLocsOffset[offset], _msf_samplingLocsLenFull[offset], &sampleErrors);
		if (err > errThreshold)
			return -1;
	}
	return err;
}
/**********************************************/
int calculateMD(int index, CompressedSeq *cmpSeq, int err, char **opSeq)
{
	index--;
	int i;
	short matchCnt = 0;
	char *op = *opSeq;
	int pp = 0;

	if (err>0 || err == -1 )
	{
		int mod = index % 21;
		int refALS = mod * 3;
		int refARS = typeSize - refALS;
		CompressedSeq tmpref, *refPos = _msf_crefGen + index/21;
		CompressedSeq *ref = refPos;

		CompressedSeq diffMask = 7;
		int shifts = (20 - mod) * 3;
		CompressedSeq diff;

		err = 0;

		for (i=0; i < SEQ_LENGTH; i++)
		{
			if (diffMask == 7)
			{
				diffMask = 0x7000000000000000;
				tmpref = (*ref << refALS) | (*(++ref) >> refARS);
				diff = (tmpref ^ *(cmpSeq++));
			}
			else
				diffMask >>= 3;

			if (diff & diffMask)		// ref[index + i - 1 ] != ver[i]
			{
				err++;
				if (matchCnt)
				{
					if (matchCnt < 10)
					{
						op[pp++]=_msf_numbers[matchCnt][0];
					}
					else if (matchCnt < 100)
					{
						op[pp++]=_msf_numbers[matchCnt][0];
						op[pp++]=_msf_numbers[matchCnt][1];
					}
					else
					{
						op[pp++]=_msf_numbers[matchCnt][0];
						op[pp++]=_msf_numbers[matchCnt][1];
						op[pp++]=_msf_numbers[matchCnt][2];
					}

					matchCnt = 0;
				}
				op[pp++] = alphabet[ (*refPos >> shifts) & 7 ];
			}
			else
			{
				matchCnt++;
			}

			if (shifts == 0)
			{
				refPos++;
				shifts = 60;
			}
			else
				shifts -= 3;

		}
	}
	if (err == 0)
	{
		matchCnt = SEQ_LENGTH;
	}

	if (matchCnt>0)
	{
		if (matchCnt < 10)
		{
			op[pp++]=_msf_numbers[matchCnt][0];
		}
		else if (matchCnt < 100)
		{
			op[pp++]=_msf_numbers[matchCnt][0];
			op[pp++]=_msf_numbers[matchCnt][1];
		}
		else
		{
			op[pp++]=_msf_numbers[matchCnt][0];
			op[pp++]=_msf_numbers[matchCnt][1];
			op[pp++]=_msf_numbers[matchCnt][2];
		}
	}
	op[pp]='\0';
	
	return err;
}

/**********************************************/
void mapSingleEndSeqListBal(unsigned int *l1, int s1, unsigned int *l2, int s2, int dir)
{

	if (s1 == 0 || s2 == 0)
	{
		return;
	}
	else if (s1 == s2 && s1 <= 50)
	{

		int j = 0;
		int z = 0;
		int *locs;
		int *seqInfo;
		CompressedSeq *_tmpCmpSeq;
		char *_tmpQual, *_tmpSeq;
		char rqual[QUAL_LENGTH+1];
		rqual[QUAL_LENGTH]='\0';
		char rseq[SEQ_LENGTH+21];		// 20 more bytes should be allocated because the decompress function needs extra space for decompressing			DEL +21 and comment
		rseq[SEQ_LENGTH]='\0';

		if (dir > 0)
		{
			locs		= (int *) l1;
			seqInfo		= (int *) l2;
		}
		else
		{
			locs		= (int *) l2;
			seqInfo		= (int *) l1;
		}


		for (j=0; j<s2; j++)
		{
			int re = _msf_samplingLocsSize * 2;
			int r = seqInfo[j]/re;
			if (maxHits!=0 && _msf_seqList[r].hits[0] == maxHits)
			{
				continue;
			}

			int x = seqInfo[j] % re;
			int o = x % _msf_samplingLocsSize;
			char d = (x/_msf_samplingLocsSize)?1:0;


			if (d)
			{
				reverse(_msf_seqList[r].qual, rqual, QUAL_LENGTH);
				_tmpQual = rqual;
				_tmpSeq = _msf_seqList[r].rseq;
			//	decompressSequence(_msf_seqList[r].crseq, rseq);		DEL
			//	reverseComplete(_msf_seqList[r].seq, rseq, SEQ_LENGTH);
			//	_tmpSeq = rseq;
				_tmpCmpSeq = _msf_seqList[r].crseq;
			}
			else
			{
				_tmpQual = _msf_seqList[r].qual;
				_tmpSeq = _msf_seqList[r].seq;
				_tmpCmpSeq = _msf_seqList[r].cseq;
			}


			for (z=0; z<s1; z++)
			{
				int genLoc = locs[z]-_msf_samplingLocs[o];


				if ( genLoc < _msf_refGenBeg || genLoc > _msf_refGenEnd )
					continue;

				int err = -1;



				err = verifySingleEnd(genLoc, _tmpCmpSeq, o);

				if (err != -1)
				{
					calculateMD(genLoc, _tmpCmpSeq, err, &_msf_op);
				
					mappingCnt++;
					_msf_seqList[r].hits[0]++;

					_msf_output.QNAME		= _msf_seqList[r].name;
					_msf_output.FLAG		= 16 * d;
					_msf_output.RNAME		= _msf_refGenName;
					_msf_output.POS			= genLoc + _msf_refGenOffset;
					_msf_output.MAPQ		= 255;
					_msf_output.CIGAR		= _msf_cigar;
					_msf_output.MRNAME		= "*";
					_msf_output.MPOS		= 0;
					_msf_output.ISIZE		= 0;
					_msf_output.SEQ			= _tmpSeq;
					_msf_output.QUAL		= _tmpQual;

					_msf_output.optSize		= 2;
					_msf_output.optFields	= _msf_optionalFields;

					_msf_optionalFields[0].tag = "NM";
					_msf_optionalFields[0].type = 'i';
					_msf_optionalFields[0].iVal = err;

					_msf_optionalFields[1].tag = "MD";
					_msf_optionalFields[1].type = 'Z';
					_msf_optionalFields[1].sVal = _msf_op;

					output(_msf_output);

					if (_msf_seqList[r].hits[0] == 1)
					{
						mappedSeqCnt++;
					}

					if ( maxHits == 0 )
					{
						_msf_seqList[r].hits[0] = 2;
					}


					if ( maxHits!=0 && _msf_seqList[r].hits[0] == maxHits)
					{
						completedSeqCnt++;
						break;
					}
				}

			}
		}
	}
	else
	{
		int tmp1=s1/2, tmp2= s2/2;
		if (tmp1 != 0)
			mapSingleEndSeqListBal(l1, tmp1, l2+tmp2, s2-tmp2, dir);
		mapSingleEndSeqListBal(l2+tmp2, s2-tmp2, l1+tmp1, s1-tmp1, -dir);
		if (tmp2 !=0)
			mapSingleEndSeqListBal(l1+tmp1, s1-tmp1, l2, tmp2, dir);
		if (tmp1 + tmp2 != 0)
			mapSingleEndSeqListBal(l2, tmp2, l1, tmp1, -dir);
	}
}


/**********************************************/
void mapSingleEndSeqListTOP(unsigned int *l1, int s1, unsigned int *l2, int s2)
{
	if (s1 < s2)
	{
		mapSingleEndSeqListBal(l1, s1, l2, s1,1);
		mapSingleEndSeqListTOP(l1, s1, l2+s1, s2-s1);		
	}
	else if (s1 > s2)
	{
		mapSingleEndSeqListBal(l1, s2, l2, s2,1);
		mapSingleEndSeqListTOP(l1+s2, s1-s2, l2, s2);
	}
	else
	{
		mapSingleEndSeqListBal(l1, s1, l2, s2,1);
	}
}


/**********************************************/
void mapSingleEndSeqList(unsigned int *l1, int s1, unsigned int *l2, int s2)
{
	if ( s2/s1 <= 2)
	{
		int j = 0;
		int z = 0;
		int *locs = (int *) l1;
		int *seqInfo = (int *) l2;
		CompressedSeq *_tmpCmpSeq;
		char *_tmpQual, *_tmpSeq;
		char rqual[QUAL_LENGTH+1];
		rqual[QUAL_LENGTH]='\0';

		for (j=0; j<s2; j++)
		{
			int re = _msf_samplingLocsSize * 2;
			int r = seqInfo[j]/re;
			if (maxHits!=0 && _msf_seqList[r].hits[0] == maxHits)
			{
				continue;
			}

			int x = seqInfo[j] % re;
			int o = x % _msf_samplingLocsSize;
			char d = (x/_msf_samplingLocsSize)?1:0;


			if (d)
			{
				reverse(_msf_seqList[r].qual, rqual, QUAL_LENGTH);
				_tmpQual = rqual;
				_tmpCmpSeq = _msf_seqList[r].crseq;
				_tmpSeq = _msf_seqList[r].seq;
			}
			else
			{
				_tmpQual = _msf_seqList[r].qual;
				_tmpSeq = _msf_seqList[r].seq;
				_tmpCmpSeq = _msf_seqList[r].cseq;
			}


			for (z=0; z<s1; z++)
			{
				int genLoc = locs[z]-_msf_samplingLocs[o];


				if ( genLoc < _msf_refGenBeg || genLoc > _msf_refGenEnd )
					continue;

				int err = -1;



				err = verifySingleEnd(genLoc, _tmpCmpSeq, o);



				if (err != -1)
				{
					calculateMD(genLoc, _tmpCmpSeq, err, &_msf_op);
					mappingCnt++;
					_msf_seqList[r].hits[0]++;

					_msf_output.QNAME		= _msf_seqList[r].name;
					_msf_output.FLAG		= 16 * d;
					_msf_output.RNAME		= _msf_refGenName;
					_msf_output.POS			= genLoc + _msf_refGenOffset;
					_msf_output.MAPQ		= 255;
					_msf_output.CIGAR		= _msf_cigar;
					_msf_output.MRNAME		= "*";
					_msf_output.MPOS		= 0;
					_msf_output.ISIZE		= 0;
					_msf_output.SEQ			= _tmpSeq;
					_msf_output.QUAL		= _tmpQual;

					_msf_output.optSize		= 2;
					_msf_output.optFields	= _msf_optionalFields;

					_msf_optionalFields[0].tag = "NM";
					_msf_optionalFields[0].type = 'i';
					_msf_optionalFields[0].iVal = err;

					_msf_optionalFields[1].tag = "MD";
					_msf_optionalFields[1].type = 'Z';
					_msf_optionalFields[1].sVal = _msf_op;

					output(_msf_output);


					if (_msf_seqList[r].hits[0] == 1)
					{
						mappedSeqCnt++;
					}

					if ( maxHits == 0 )
					{
						_msf_seqList[r].hits[0] = 2;
					}


					if ( maxHits!=0 && _msf_seqList[r].hits[0] == maxHits)
					{
						completedSeqCnt++;
						break;
					}
				}

			}
		}
	}
	else if (s1 == 1)
	{
		int tmp = s2/2;
		mapSingleEndSeqList(l1, s1, l2, tmp);
		mapSingleEndSeqList(l1, s1, l2+tmp, s2-tmp);
	}
	else if (s2 == 1)
	{
		int tmp = s1/2;
		mapSingleEndSeqList(l1, tmp, l2, s2);
		mapSingleEndSeqList(l1+tmp, s1-tmp, l2, s2);
	}
	else
	{
		int tmp1=s1/2, tmp2= s2/2;
		mapSingleEndSeqList(l1, tmp1, l2, tmp2);
		mapSingleEndSeqList(l1+tmp1, s1-tmp1, l2, tmp2);
		mapSingleEndSeqList(l1+tmp1, s1-tmp1, l2+tmp2, s2-tmp2);
		mapSingleEndSeqList(l1, tmp1, l2+tmp2, s2-tmp2);
	}
}
/**********************************************/
int	 mapSingleEndSeq()
{
	int i = 0;
	unsigned int *locs = NULL;
	unsigned int *seqInfo = NULL;
	while ( i < _msf_rIndexSize )
	{
		locs = getCandidates (_msf_rIndex[i].hv);
		if ( locs != NULL)
		{
			seqInfo  = _msf_rIndex[i].seqInfo;
			mapSingleEndSeqListTOP (locs+1, locs[0], seqInfo+1, seqInfo[0]);			
		}
		i++;
	}
	return 1;
}


/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
int compareOut (const void *a, const void *b)
{
	FullMappingInfo *aInfo = (FullMappingInfo *)a;
	FullMappingInfo *bInfo = (FullMappingInfo *)b;
	return aInfo->loc - bInfo->loc;
}

/**********************************************/
void mapPairedEndSeqList(unsigned int *l1, int s1, unsigned int *l2, int s2)
{
	if ( s2/s1 <= 2)
	{
		int j = 0;
		int z = 0;
		int *locs = (int *) l1;
		int *seqInfo = (int *) l2;
		CompressedSeq *_tmpCmpSeq;
		char *_tmpQual;// *_tmpSeq;
		char rqual[QUAL_LENGTH+1];
		rqual[QUAL_LENGTH]='\0';

		for (j=0; j<s2; j++)
		{
			int re = _msf_samplingLocsSize * 2;
			int r = seqInfo[j]/re;

			if (pairedEndDiscordantMode && (_msf_seqList[r].hits[0] == 1 || (_msf_seqHits[r/2] > DISCORDANT_CUT_OFF) ))
			{
				continue;
			}

			int x = seqInfo[j] % re;
			int o = x % _msf_samplingLocsSize;
			char d = (x/_msf_samplingLocsSize)?-1:1;


			if (d==-1)
			{
				_tmpCmpSeq = _msf_seqList[r].crseq;
			}
			else
			{
				_tmpCmpSeq = _msf_seqList[r].cseq;
			}


			for (z=0; z<s1; z++)
			{
				int genLoc = locs[z]-_msf_samplingLocs[o];

				if ( genLoc < _msf_refGenBeg || genLoc > _msf_refGenEnd )
					continue;

				int err = -1;



				err = verifySingleEnd(genLoc, _tmpCmpSeq, o);


				if (err != -1)
				{
					MappingLocations *parent = NULL;
					MappingLocations *child = _msf_mappingInfo[r].next;

					genLoc+= _msf_refGenOffset;
					int i = 0;
					for (i=0; i<(_msf_mappingInfo[r].size/MAP_CHUNKS); i++)
					{
						parent = child;
						child = child->next;
					}

					if (child==NULL)
					{
						MappingLocations *tmp = getMem(sizeof(MappingLocations));
						tmp->next = NULL;
						tmp->loc[0]=genLoc * d;
						if (parent == NULL)
							_msf_mappingInfo[r].next = tmp;
						else
							parent->next = tmp;
					}
					else
					{
						child->loc[_msf_mappingInfo[r].size % MAP_CHUNKS] = genLoc * d;
					}



					_msf_mappingInfo[r].size++;

				}

			}
		}
	}
	else if (s1 == 1)
	{
		int tmp = s2/2;
		mapPairedEndSeqList(l1, s1, l2, tmp);
		mapPairedEndSeqList(l1, s1, l2+tmp, s2-tmp);
	}
	else if (s2 == 1)
	{
		int tmp = s1/2;
		mapPairedEndSeqList(l1, tmp, l2, s2);
		mapPairedEndSeqList(l1+tmp, s1-tmp, l2, s2);
	}
	else
	{
		int tmp1=s1/2, tmp2= s2/2;
		mapPairedEndSeqList(l1, tmp1, l2, tmp2);
		mapPairedEndSeqList(l1+tmp1, s1-tmp1, l2, tmp2);
		mapPairedEndSeqList(l1+tmp1, s1-tmp1, l2+tmp2, s2-tmp2);
		mapPairedEndSeqList(l1, tmp1, l2+tmp2, s2-tmp2);
	}
}

/**********************************************/
int	 mapPairedEndSeq()
{
	int i = 0;
	unsigned int *locs = NULL;
	unsigned int *seqInfo = NULL;
	while ( i < _msf_rIndexSize )
	{
		locs = getCandidates (_msf_rIndex[i].hv);
		if ( locs != NULL)
		{
			seqInfo  = _msf_rIndex[i].seqInfo;
			mapPairedEndSeqList(locs+1, locs[0], seqInfo+1, seqInfo[0]);

		}
		i++;
	}


	char fname1[FILE_NAME_LENGTH];
	char fname2[FILE_NAME_LENGTH];
	MappingLocations *cur, *tmp;
	int tmpOut;
	int j;
	int lmax=0, rmax=0;

	sprintf(fname1, "%s__%s__%d__1",mappingOutputPath, mappingOutput, _msf_openFiles);
	sprintf(fname2, "%s__%s__%d__2",mappingOutputPath, mappingOutput, _msf_openFiles);

	FILE* out;
	FILE* out1 = fileOpen(fname1, "w");
	FILE* out2 = fileOpen(fname2, "w");

	_msf_openFiles++;

	for (i=0; i<_msf_seqListSize; i++)
	{

		if (i%2==0)
		{
			out = out1;

			if (lmax <  _msf_mappingInfo[i].size)
			{
				lmax = _msf_mappingInfo[i].size;
			}
		}
		else
		{
			out = out2;
			if (rmax < _msf_mappingInfo[i].size)
			{	
				rmax = _msf_mappingInfo[i].size;
			}
		}

		tmpOut = fwrite(&(_msf_mappingInfo[i].size), sizeof(int), 1, out);					
		if (_msf_mappingInfo[i].size > 0)
		{
			cur = _msf_mappingInfo[i].next;
			for (j=0; j < _msf_mappingInfo[i].size; j++)
			{
				if ( j>0  && j%MAP_CHUNKS==0)
				{
					cur = cur->next;
				}
				tmpOut = fwrite(&(cur->loc[j % MAP_CHUNKS]), sizeof(int), 1, out);
			}
			_msf_mappingInfo[i].size = 0;
			//	_msf_mappingInfo[i].next = NULL;			DEL
		}
	}

	_msf_maxLSize += lmax;
	_msf_maxRSize += rmax;

	fclose(out1);
	fclose(out2);

	//fprintf(stdout, "%d %d\n", _msf_maxLSize, _msf_maxRSize);		DEL

}

/**********************************************/
void outputPairedEnd()
{

	char *curGen;
	char *curGenName;
	int tmpOut;

	//	loadRefGenome(&_msf_refGen, &_msf_refGenName, &tmpOut);		DEL
	CompressedSeq * tmpCrefgen = _msf_crefGen;		//				DEL
	_msf_crefGen = getCmpRefGenOrigin();

	FILE* in1[_msf_openFiles];
	FILE* in2[_msf_openFiles];

	char fname1[_msf_openFiles][FILE_NAME_LENGTH];	
	char fname2[_msf_openFiles][FILE_NAME_LENGTH];	

	// discordant
	FILE *out, *out1, *out2;

	char fname3[FILE_NAME_LENGTH];
	char fname4[FILE_NAME_LENGTH];
	char fname5[FILE_NAME_LENGTH];

	if (pairedEndDiscordantMode)
	{
		sprintf(fname3, "%s__%s__disc", mappingOutputPath, mappingOutput);
		sprintf(fname4, "%s__%s__oea1", mappingOutputPath, mappingOutput);
		sprintf(fname5, "%s__%s__oea2", mappingOutputPath, mappingOutput);
		out = fileOpen(fname3, "a");
		out1 = fileOpen(fname4, "a");
		out2 = fileOpen(fname5, "a");
	}



	int i;

	FullMappingInfo *mi1 = getMem(sizeof(FullMappingInfo) * _msf_maxLSize);
	FullMappingInfo *mi2 = getMem(sizeof(FullMappingInfo) * _msf_maxRSize);


	for (i=0; i<_msf_openFiles; i++)
	{
		sprintf(fname1[i], "%s__%s__%d__1", mappingOutputPath, mappingOutput, i);
		sprintf(fname2[i], "%s__%s__%d__2", mappingOutputPath, mappingOutput, i);
		in1[i] = fileOpen(fname1[i], "r");
		in2[i] = fileOpen(fname2[i], "r");
	}


	int size;
	int j, k;
	int size1, size2;

	for (i=0; i<_msf_seqListSize/2; i++)
	{
		size1 = size2 = 0;
		for (j=0; j<_msf_openFiles; j++)
		{
			tmpOut = fread(&size, sizeof(int), 1, in1[j]);
			if ( size > 0 )
			{
				for (k=0; k<size; k++)
				{

					mi1[size1+k].dir = 1;
					tmpOut = fread (&(mi1[size1+k].loc), sizeof(int), 1, in1[j]);
					if (mi1[size1+k].loc<1)
					{	
						mi1[size1+k].loc *= -1;
						mi1[size1+k].dir = -1;
					}
				}
				qsort(mi1+size1, size, sizeof(FullMappingInfo), compareOut);
				size1+=size;
			}
		}

		for (j=0; j<_msf_openFiles; j++)
		{
			tmpOut = fread(&size, sizeof(int), 1, in2[j]);
			if ( size > 0 )
			{
				for (k=0; k<size; k++)
				{

					mi2[size2+k].dir = 1;
					tmpOut = fread (&(mi2[size2+k].loc), sizeof(int), 1, in2[j]);

					if (mi2[size2+k].loc<1)
					{	
						mi2[size2+k].loc *= -1;
						mi2[size2+k].dir = -1;
					}
				}
				qsort(mi2+size2, size, sizeof(FullMappingInfo), compareOut);
				size2+=size;
			}
		}

		int lm, ll, rl, rm;
		int pos = 0;

		if (pairedEndDiscordantMode)
		{

			for (j=0; j<size1; j++)
			{
				lm = mi1[j].loc - maxPairEndedDiscordantDistance + 1;
				ll = mi1[j].loc - minPairEndedDiscordantDistance + 1;
				rl = mi1[j].loc + minPairEndedDiscordantDistance - 1;
				rm = mi1[j].loc + maxPairEndedDiscordantDistance - 1;

				while (pos<size2 && mi2[pos].loc < lm)
				{
					pos++;
				}

				k = pos;
				while (k<size2 && mi2[k].loc<=rm)
				{
					if ( mi2[k].loc <= ll || mi2[k].loc >= rl)
					{
						if ( (mi1[j].loc < mi2[k].loc && mi1[j].dir==1 && mi2[k].dir == -1)  ||  
								(mi1[j].loc > mi2[k].loc && mi1[j].dir==-1 && mi2[k].dir == 1) )
						{
							_msf_seqList[i*2].hits[0]=1;
							_msf_seqList[i*2+1].hits[0]=1;
							size1=0;
							size2=0;
							break;
						}
					}
					k++;
				}
			}

			_msf_seqHits[i]+= size1*size2;
			if (_msf_seqHits[i]> DISCORDANT_CUT_OFF)
			{
				_msf_seqList[i*2].hits[0]=1;
				_msf_seqList[i*2+1].hits[0]=1;
				size1=0;
				size2=0;
			}
		}
		//if (i == 6615)										DEL
		//	fprintf(stdout, "%d %d\n", size1, size2);	

		CompressedSeq *cseq1, *cseq2, *crseq1, *crseq2;
		char *seq1, *seq2, *qual1, *qual2, *rseq1, *rseq2;
		char rqual1[QUAL_LENGTH+1], rqual2[QUAL_LENGTH+1];
		rqual1[QUAL_LENGTH] = rqual2[QUAL_LENGTH] = '\0';
		seq1 = _msf_seqList[i*2].seq;
		rseq1 = _msf_seqList[i*2].rseq;
		cseq1 = _msf_seqList[i*2].cseq;
		crseq1 = _msf_seqList[i*2].crseq;
		qual1 = _msf_seqList[i*2].qual;
		reverse(_msf_seqList[i*2].qual, rqual1, QUAL_LENGTH);

		seq2 = _msf_seqList[i*2+1].seq;
		rseq2 = _msf_seqList[i*2+1].rseq;	
		cseq2 = _msf_seqList[i*2+1].cseq;
		crseq2 = _msf_seqList[i*2+1].crseq;

		qual2 = _msf_seqList[i*2+1].qual;
		reverse(_msf_seqList[i*2+1].qual, rqual2, QUAL_LENGTH);


		if (pairedEndDiscordantMode)
		{
			for (k=0; k<size1; k++)
			{
				int tm = -1;
				mi1[k].score = calculateScore(mi1[k].loc, (mi1[k].dir==-1)?crseq1:cseq1, (mi1[k].dir==-1)?rqual1:qual1, &tm);
				mi1[k].err = tm;
			}

			for (k=0; k<size2; k++)
			{
				int tm = -1;
				mi2[k].score = calculateScore(mi2[k].loc, (mi2[k].dir==-1)?crseq2:cseq2, (mi2[k].dir==-1)?rqual2:qual2, &tm);
				mi2[k].err = tm;
			}

		}
		else
		{
			for (k=0; k<size1; k++)
			{
				mi1[k].err = calculateMD(mi1[k].loc, (mi1[k].dir==-1)?crseq1:cseq1, -1, &_msf_op);
				sprintf(mi1[k].md, "%s", _msf_op);
			}

			for (k=0; k<size2; k++)
			{
				mi2[k].err = calculateMD(mi2[k].loc, (mi2[k].dir==-1)?crseq2:cseq2, -1, &_msf_op);
				sprintf(mi2[k].md, "%s", _msf_op);
			}
		}
		pos = 0;

		for (j=0; j<size1; j++)
		{
			lm = mi1[j].loc - maxPairEndedDistance + 1;
			ll = mi1[j].loc - minPairEndedDistance + 1;
			rl = mi1[j].loc + minPairEndedDistance - 1;
			rm = mi1[j].loc + maxPairEndedDistance - 1;

			//fprintf(stdout, "%d %d %d %d %d\n",lm, ll,mi1[j].loc ,rl, rm); 		DEL

			while (pos<size2 && mi2[pos].loc < lm)
			{
				pos++;
			}

			//fprintf(stdout, "POS: %d %d \n", pos, mi2[pos].loc);					DEL

			k = pos;
			while (k<size2 && mi2[k].loc <= rm)
			{
				if (mi2[k].loc <= ll || mi2[k].loc >= rl)
				{
					if (pairedEndDiscordantMode)
					{
						int tmp;
						int rNo = i;
						int loc = mi1[j].loc*mi1[j].dir;
						int err = mi1[j].err;
						float sc = mi1[j].score;

						char l = strlen(_msf_refGenName);

						tmp = fwrite(&rNo, sizeof(int), 1, out);

						tmp = fwrite(&l, sizeof(char), 1, out);
						tmp = fwrite(_msf_refGenName, sizeof(char), l, out);

						tmp = fwrite(&loc, sizeof(int), 1, out);
						tmp = fwrite(&err, sizeof(char), 1, out);
						tmp = fwrite(&sc, sizeof(float), 1, out);


						loc = mi2[k].loc*mi2[k].dir;
						err = mi2[k].err;
						sc = mi2[k].score;

						tmp = fwrite(&loc, sizeof(int), 1, out);
						tmp = fwrite(&err, sizeof(char), 1, out);
						tmp = fwrite(&sc, sizeof(float), 1, out);
					} // end discordant
					else
					{ //start sampe
						CompressedSeq *cmpSeq, *cmpRseq;
						char *seq;
						char *qual;
						char d1;
						char d2;
						int isize;
						int proper=0;
						// ISIZE CALCULATION
						// The distance between outer edges								
						isize = abs(mi1[j].loc - mi2[k].loc)+SEQ_LENGTH-1;												
						if (mi1[j].loc - mi2[k].loc > 0)
						{
							isize *= -1;
						}

						d1 = (mi1[j].dir == -1)?1:0;
						d2 = (mi2[k].dir == -1)?1:0;

						if ( d1 )
						{
							seq = rseq1;
							qual = rqual1;
						}
						else
						{
							seq = seq1;
							qual = qual1;
						}

						if ( (mi1[j].loc < mi2[k].loc && !d1 && d2) ||
								(mi1[j].loc > mi2[k].loc && d1 && !d2) )
						{
							proper = 2;
						}
						else
						{
							proper = 0;
						}


						_msf_output.POS			= mi1[j].loc;
						_msf_output.MPOS		= mi2[k].loc;
						_msf_output.FLAG		= 1+proper+16*d1+32*d2+64;
						_msf_output.ISIZE		= isize;
						_msf_output.SEQ			= seq,
							_msf_output.QUAL		= qual;
						_msf_output.QNAME		= _msf_seqList[i*2].name;
						_msf_output.RNAME		= _msf_refGenName;
						_msf_output.MAPQ		= 255;
						_msf_output.CIGAR		= _msf_cigar;
						_msf_output.MRNAME		= "=";

						_msf_output.optSize	= 2;
						_msf_output.optFields	= _msf_optionalFields;

						_msf_optionalFields[0].tag = "NM";
						_msf_optionalFields[0].type = 'i';
						_msf_optionalFields[0].iVal = mi1[j].err;

						_msf_optionalFields[1].tag = "MD";
						_msf_optionalFields[1].type = 'Z';
						_msf_optionalFields[1].sVal = mi1[j].md;


						output(_msf_output);

						if ( d2 )
						{
							seq = rseq2;
							qual = rqual2;
						}
						else
						{
							seq = seq2;
							qual = qual2;
						}

						_msf_output.POS			= mi2[k].loc;
						_msf_output.MPOS		= mi1[j].loc;
						_msf_output.FLAG		= 1+proper+16*d2+32*d1+128;
						_msf_output.ISIZE		= -isize;
						_msf_output.SEQ			= seq,
							_msf_output.QUAL		= qual;
						_msf_output.QNAME		= _msf_seqList[i*2].name;
						_msf_output.RNAME		= _msf_refGenName;
						_msf_output.MAPQ		= 255;
						_msf_output.CIGAR		= _msf_cigar;
						_msf_output.MRNAME		= "=";

						_msf_output.optSize	= 2;
						_msf_output.optFields	= _msf_optionalFields;

						_msf_optionalFields[0].tag = "NM";
						_msf_optionalFields[0].type = 'i';
						_msf_optionalFields[0].iVal = mi2[k].err;;

						_msf_optionalFields[1].tag = "MD";
						_msf_optionalFields[1].type = 'Z';
						_msf_optionalFields[1].sVal = mi2[k].md;

						output(_msf_output);
					} //end sampe
				}
				k++;
			}
		}
	}

	if (pairedEndDiscordantMode)
	{
		fclose(out);
	}


	freeMem(mi1, sizeof(FullMappingInfo)*_msf_maxLSize);
	freeMem(mi2, sizeof(FullMappingInfo)*_msf_maxRSize);

	for (i=0; i<_msf_openFiles; i++)
	{
		fclose(in1[i]);
		fclose(in2[i]);
		//fprintf(stdout, "%s %s \n", fname1[i], fname2[i]);
		unlink(fname1[i]);
		unlink(fname2[i]);
	}
	_msf_openFiles = 0;

	_msf_crefGen = tmpCrefgen;		///////////////////
}

/**********************************************/
/**********************************************/
/**********************************************/
/**********************************************/
float calculateScore(int index, CompressedSeq *cmpSeq, char *qual, int *err)
{
	index--;
	float score = 1;
	int i;
	int mod = index % 21;
	int refALS = mod * 3;
	int refARS = typeSize - refALS;
	CompressedSeq tmpref, *refPos = _msf_crefGen + index/21;
	CompressedSeq *ref = refPos;								// DEL refPos

	CompressedSeq diffMask = 7;
	int shifts = (20 - mod) * 3;								// DEL shifts
	CompressedSeq diff;

	*err = 0;

	for (i=0; i < SEQ_LENGTH; i++)
	{
		if (diffMask == 7)
		{
			diffMask = 0x7000000000000000;
			tmpref = (*ref << refALS) | (*(++ref) >> refARS);
			diff = (tmpref ^ *(cmpSeq++));
		}
		else
			diffMask >>= 3;

		if (diff & diffMask)		// ref[index + i - 1 ] != ver[i]
		{
			(*err)++;
			score *= 0.001 + 1/pow( 10, ((qual[i]-33)/10.0) );
		}
	}

	return score;
}

/**********************************************/
void outputPairedEndDiscPP()
{
	char genName[SEQ_LENGTH];
	char fname1[FILE_NAME_LENGTH];
	char fname2[FILE_NAME_LENGTH];
	char fname3[FILE_NAME_LENGTH];
	char fname4[FILE_NAME_LENGTH];
	char fname5[FILE_NAME_LENGTH];
	char fname6[FILE_NAME_LENGTH];
	char l;
	int loc1, loc2;
	char err1, err2;
	char dir1, dir2;
	float sc1, sc2, lsc=0;
	int flag = 0;
	int rNo,lrNo = -1;
	int tmp;
	FILE *in, *in1, *in2, *out, *out1, *out2;

	sprintf(fname1, "%s__%s__disc", mappingOutputPath, mappingOutput);
	sprintf(fname2, "%s__%s__oea1", mappingOutputPath, mappingOutput);
	sprintf(fname3, "%s__%s__oea2", mappingOutputPath, mappingOutput);
	sprintf(fname4, "%s%s_DIVET.vh", mappingOutputPath, mappingOutput);
	sprintf(fname5, "%s%s_OEA1.vh", mappingOutputPath, mappingOutput);
	sprintf(fname6, "%s%s_OEA2.vh", mappingOutputPath, mappingOutput);

	in   = fileOpen(fname1, "r");
	in1  = fileOpen(fname2, "r");
	in2  = fileOpen(fname3, "r");
	out  = fileOpen(fname4, "w");
	out1 = fileOpen(fname5, "w");
	out2 = fileOpen(fname6, "w");
	if (in != NULL)
	{
		flag = fread(&rNo, sizeof(int), 1, in);
	}
	else
	{
		flag  = 0;
	}


	while (flag)
	{

		tmp = fread(&l, sizeof(char), 1, in);
		tmp = fread(genName, sizeof(char), l, in);
		genName[l]='\0';
		tmp = fread(&loc1, sizeof(int), 1, in);
		tmp = fread(&err1, sizeof(char), 1, in);
		tmp = fread(&sc1, sizeof(float), 1, in);
		tmp = fread(&loc2, sizeof(int), 1, in);
		tmp = fread(&err2, sizeof(char), 1, in);
		tmp = fread(&sc2, sizeof(float), 1, in);

		//if (rNo ==6615)								DEL
		//	fprintf(stdout, "%s %d: %d %0.20f %d %d %0.20f\n", genName, loc1, err1, sc1, loc2, err2, sc2);

		if (_msf_seqList[rNo*2].hits[0] % 2 == 0 && _msf_seqHits[rNo] < DISCORDANT_CUT_OFF)
		{
			dir1 = dir2 = 'F';

			if (loc1 < 0)
			{
				dir1 = 'R';
				loc1 = -loc1;
			}

			if (loc2 < 0)
			{
				dir2 = 'R';
				loc2 = -loc2;
			}

			if (rNo != lrNo)
			{
				int j;
				for (j=0; j<SEQ_LENGTH; j++)
				{
					lsc += _msf_seqList[rNo*2].qual[j]+_msf_seqList[rNo*2+1].qual[j];
				}
				lsc /= 2*SEQ_LENGTH;
				lsc -= 33;
				lrNo = rNo;
			}

			int inv = 0;
			int eve = 0;
			int dist = 0;
			char event;

			//fprintf(stdout, "%c %c ", dir1, dir2);

			if ( dir1 == dir2 )
			{
				event = 'V';
				//fprintf(stdout, "Inverstion \n");
			}
			else
			{
				if (loc1 < loc2)
				{

					//fprintf(stdout, "< %d ", loc2-loc1-SEQ_LENGTH);

					if (dir1 == 'R' && dir2 == 'F')
					{
						event = 'E';

						//fprintf(stdout, "Everted \n");
					}
					else if ( loc2 - loc1 >= maxPairEndedDiscordantDistance )
					{
						event = 'D';
						//fprintf(stdout, "Deletion \n");
					}
					else
					{
						event = 'I';
						//fprintf(stdout, "Insertion \n");
					}
				}
				else if (loc2 < loc1)
				{
					//fprintf(stdout, "> %d ", loc1-loc2-SEQ_LENGTH);
					if (dir2 == 'R' && dir1 == 'F')
					{
						event = 'E';
						//fprintf(stdout, "Everted \n");
					}
					else if ( loc1 - loc2 >= maxPairEndedDiscordantDistance )
					{
						event = 'D';
						//fprintf(stdout, "Deletion \n");
					}
					else
					{
						event = 'I';
						//fprintf(stdout, "Insertion \n");
					}
				}
			}
			_msf_seqList[rNo*2].hits[0] = 2;
			fprintf(out, "%s\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%c\t%d\t%0.0f\t%0.20f\n",
					_msf_seqList[rNo*2].name, genName, loc1, (loc1+SEQ_LENGTH-1), dir1, genName, loc2, (loc2+SEQ_LENGTH-1), dir2, event, (err1+err2), lsc, sc1*sc2);
		}
		flag = fread(&rNo, sizeof(int), 1, in);

	}

	/*
	   MappingInfoNode *lr[_msf_seqListSize/2];
	   MappingInfoNode *rr[_msf_seqListSize/2];
	   MappingInfoNode *cur, *tmpDel, *cur2;


	   int ls[_msf_seqListSize/2];
	   int rs[_msf_seqListSize/2];


	   int i=0;

	   for (i = 0; i<_msf_seqListSize/2; i++)
	   {
	   lr[i] = rr[i] = NULL;
	   ls[i] = rs[i] = 0;
	   }



	   if (in1 != NULL)
	   {
	   flag = fread(&rNo, sizeof(int), 1, in1);
	   }
	   else
	   {
	   flag  = 0;
	   }


	   while (flag)
	   {
	   tmp = fread(&loc1, sizeof(int), 1, in1);
	   tmp = fread(&err1, sizeof(char), 1, in1);
	   tmp = fread(&sc1, sizeof(float), 1, in1);
	   tmp = fread(&l, sizeof(char), 1, in1);
	   tmp = fread(genName, sizeof(char), l, in1);
	   genName[l]='\0';

	   if (_msf_seqList[rNo*2].hits[0] == 0)
	   {

	   if ( ls[rNo] < DISCORDANT_CUT_OFF )
	   {
	   ls[rNo]++;

	   cur = lr[rNo];

	   if (cur !=NULL)
	   {
	   if (err1 == cur->err)
	   {
	   MappingInfoNode *nr = getMem(sizeof(MappingInfoNode));

	   nr->loc = loc1;
	   nr->err = err1;
	   nr->score = sc1;
	   nr->next = lr[rNo];
	   sprintf(nr->chr,"%s", genName);
	   lr[rNo] = nr;
	   }
	   else if (err1 < cur->err)
	   {
	   MappingInfoNode *nr = getMem(sizeof(MappingInfoNode));

	   nr->loc = loc1;
	   nr->err = err1;
	   nr->score = sc1;
	   sprintf(nr->chr,"%s", genName);
	   nr->next = NULL;
	   lr[rNo] = nr;
	while (cur!=NULL)
	{
		tmpDel = cur;
		cur = cur->next;
		freeMem(tmpDel, sizeof(MappingInfoNode));
	}
}
}
else
{

	MappingInfoNode *nr = getMem(sizeof(MappingInfoNode));

	nr->loc = loc1;
	nr->err = err1;
	nr->score = sc1;
	sprintf(nr->chr,"%s", genName);
	nr->next = NULL;
	lr[rNo] = nr;
}

if (ls[rNo] > DISCORDANT_CUT_OFF)
{
	cur = lr[rNo];
	while (cur!=NULL)
	{
		tmpDel = cur;
		cur = cur->next;
		freeMem(tmpDel, sizeof(MappingInfoNode));
	}
}
}

}
flag = fread(&rNo, sizeof(int), 1, in1);

}


if (in2 != NULL)
{
	flag = fread(&rNo, sizeof(int), 1, in2);
}
else
{
	flag  = 0;
}


while (flag)
{
	tmp = fread(&loc1, sizeof(int), 1, in2);
	tmp = fread(&err1, sizeof(char), 1, in2);
	tmp = fread(&sc1, sizeof(float), 1, in2);
	tmp = fread(&l, sizeof(char), 1, in2);
	tmp = fread(genName, sizeof(char), l, in2);
	genName[l]='\0';

	if (_msf_seqList[rNo*2].hits[0] == 0)
	{

		if ( rs[rNo] < DISCORDANT_CUT_OFF )
		{
			rs[rNo]++;

			cur = rr[rNo];

			if (cur !=NULL)
			{
				if (err1 == cur->err)
				{
					MappingInfoNode *nr = getMem(sizeof(MappingInfoNode));

					nr->loc = loc1;
					nr->err = err1;
					nr->score = sc1;
					nr->next = rr[rNo];
					sprintf(nr->chr,"%s", genName);
					rr[rNo] = nr;
				}
				else if (err1 < cur->err)
				{
					MappingInfoNode *nr = getMem(sizeof(MappingInfoNode));

					nr->loc = loc1;
					nr->err = err1;
					nr->score = sc1;
					sprintf(nr->chr,"%s", genName);
					nr->next = NULL;
					rr[rNo] = nr;
					while (cur!=NULL)
					{
						tmpDel = cur;
						cur = cur->next;
						freeMem(tmpDel, sizeof(MappingInfoNode));
					}
				}
			}
			else
			{

				MappingInfoNode *nr = getMem(sizeof(MappingInfoNode));

				nr->loc = loc1;
				nr->err = err1;
				nr->score = sc1;
				sprintf(nr->chr,"%s", genName);
				nr->next = NULL;
				rr[rNo] = nr;
			}

			if (rs[rNo] > DISCORDANT_CUT_OFF)
			{
				cur = rr[rNo];
				while (cur!=NULL)
				{
					tmpDel = cur;
					cur = cur->next;
					freeMem(tmpDel, sizeof(MappingInfoNode));
				}
			}
		}
	}
	flag = fread(&rNo, sizeof(int), 1, in2);

}


for (i=0; i<_msf_seqListSize/2; i++)
{
	int j;
	for (j=0; j<SEQ_LENGTH; j++)
	{
		lsc += _msf_seqList[i*2].qual[j]+_msf_seqList[i*2+1].qual[j];
	}
	lsc /= 2*SEQ_LENGTH;
	lsc -= 33;
	if (ls[i] * rs[i] < DISCORDANT_CUT_OFF && ls[i] & rs[i] > 0)
	{
		cur = lr[i];
		while (cur != NULL)
		{
			cur2 = rr[i];
			if (cur->loc < 0)
			{
				dir1 = 'R';
				loc1 = -cur->loc;
			}
			else
			{
				dir1 = 'F';
				loc1 = cur->loc;
			}
			while (cur2 != NULL)
			{

				if (cur2->loc < 0)
				{
					dir2 = 'R';
					loc2 = -cur2->loc;
				}
				else
				{
					dir2 = 'F';
					loc2 = cur2->loc;
				}

				fprintf(out, "%s\t%s\t%d\t%d\t%c\t%s\t%d\t%d\t%c\t%c\t%d\t%0.0f\t%0.20f\n",
						_msf_seqList[i*2].name, cur->chr, loc1, (loc1+SEQ_LENGTH-1), dir1, cur2->chr, loc2, (loc2+SEQ_LENGTH-1), dir2, 'T', (cur->err+cur2->err), lsc, cur->score*cur2->score);
				cur2 = cur2->next;
			}
			cur = cur->next;
		}
	}

}*/


fclose(in);
fclose(in1);
fclose(in2);
fclose(out);
fclose(out1);
fclose(out2);

unlink(fname1);
unlink(fname2);
unlink(fname3);
unlink(fname5);
unlink(fname6);
}
