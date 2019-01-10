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
#include <ctype.h>
#include <zlib.h>
#include <pthread.h>
#include "Common.h"
#include "Reads.h"
#include "MrsFAST.h"
#include "Sort.h"

FILE *_r_fp1;
FILE *_r_fp2;
FILE *_r_umfp;
gzFile _r_gzumfp;
gzFile _r_gzfp1;
gzFile _r_gzfp2;
Read *_r_seq;
int _r_uncroppedSeqLength;
int _r_seqCnt;
int _r_samplingLocsSize;
int *_r_samplingLocs;
int *_r_samplingLocsSeg;
int *_r_samplingLocsOffset;
int *_r_samplingLocsLen;
int *_r_samplingLocsLenFull;
int *_r_indexSize;
pthread_t *_r_threads;
unsigned char _r_fastq;
ReadIndexTable **_r_readIndex;
int *_r_readIndexSize;
int _r_maxSeqCnt;
int	_r_firstIteration = 1;
long long _r_readMemUsage = 0;
char *_r_alphIndex = NULL;
char checkSumLength = 0;
char *_r_buf1;
char *_r_buf2;
int *_r_buf1_size;
int *_r_buf2_size;
int *_r_buf1_pos;
int *_r_buf2_pos;
/**********************************************/
void (*readBuffer1)();
void (*readBuffer2)();
/**********************************************/
void readBufferGZ1()
{
	(*_r_buf1_size) = gzread(_r_gzfp1, _r_buf1, 10000000);
	(*_r_buf1_pos) = 0;
}
/**********************************************/
void readBufferGZ2()
{
	if (_r_buf1 == _r_buf2)
		readBuffer1();
	else
	{
		(*_r_buf2_size) = gzread(_r_gzfp2, _r_buf2, 10000000);
		(*_r_buf2_pos) = 0;		
	}
}
/**********************************************/
void readBufferTxT1()
{
	(*_r_buf1_size) = fread(_r_buf1, 1, 10000000, _r_fp1);
	(*_r_buf1_pos) = 0;
}
/**********************************************/
void readBufferTxT2()
{
	if (_r_buf1 == _r_buf2)
		readBuffer1();
	else
	{
		(*_r_buf2_size) = fread(_r_buf2, 1, 10000000, _r_fp2);
		(*_r_buf2_pos) = 0;		
	}
}
/**********************************************/
int readFirstSeq(char *seq , int line)
{
	int i=0, l=0, j = 0;
	char cur;

	while (1)
	{
		if (*_r_buf1_pos ==*_r_buf1_size)
		{
			readBuffer1();
			if (*_r_buf1_size == 0)
				return 0;
		}

		cur = _r_buf1[*_r_buf1_pos];		
		(*_r_buf1_pos)++;
		j++;

		if ( cur == '\n')
		{
			if (l>0) i=l;
			seq[i]='\0';
			return i;
		}
		
		if (l==0 && cur == ' ')
		{
			l = i;
		}

		if (cropSize>0 && line%2==0 && i==cropSize)
			continue;
		if (tailCropSize>0 && line%2==0 && j <= _r_uncroppedSeqLength - tailCropSize)
			continue;

		seq[i++]=cur;

	}
}
/**********************************************/
int readSecondSeq( char *seq, int line )
{
	int i=0, l=0, j = 0;;
	char cur;
	while (1)
	{
		if (*_r_buf2_pos ==*_r_buf2_size)
		{
			readBuffer2();
			if (*_r_buf2_size == 0)
				return 0;
		}

		cur = _r_buf2[*_r_buf2_pos];		
		(*_r_buf2_pos)++;
		j++;
		if ( cur == '\n')
		{
			if (l>0) i=l;
			seq[i]='\0';
			return i;
		}
		
		if (cur == ' '&& l==0)
		{
			l = i;
		}

		if (cropSize>0 && line%2==0 && i==cropSize)
			continue;
		if (tailCropSize>0 && line%2==0 && j <= _r_uncroppedSeqLength - tailCropSize)
			continue;

		seq[i++]=cur;

	}
}
/**********************************************/
int compare (const void *a, const void *b)
{
	Pair *x=(Pair *)a;
	Pair *y=(Pair *)b;

	if (x->hv == y->hv)
		return x->checksum - y->checksum;
	else
		return x->hv - y->hv;
}
/**********************************************/
void getReadIndex(ReadIndexTable ***rIndex, int **rIndexSize)
{
	*rIndex = _r_readIndex;
	*rIndexSize = _r_readIndexSize;
}

/**********************************************/
void *preProcessReads(int *idp)
{
	int id = *idp;
	int i=0, j=0, pos=0, tmpSize=0;

	int32_t	hvtmp, cstmp;

	int div = _r_seqCnt / THREAD_COUNT;
	div += (_r_seqCnt % THREAD_COUNT)?1:0;
	Pair *tmp = getMem(sizeof(Pair)*(div * _r_samplingLocsSize*2));
	char alphCnt[5];
	char *a, *b;

	for (i=id*div; i<div*(id+1) && i<_r_seqCnt; i++)
	{

		alphCnt[0]=alphCnt[1]=alphCnt[2]=alphCnt[3]=alphCnt[4]=0;
		a = _r_seq[i].seq;
		b = _r_seq[i].rseq+SEQ_LENGTH;
		*(b--)='\0';
		for (j = 0; j<SEQ_LENGTH; j++)
		{
			*a = toupper(*a);
			*b = reverseCompleteChar(*a);
			alphCnt[_r_alphIndex[*a]]++;
			a++;
			b--;			
		}

		if (alphCnt[4]>errThreshold)
			_r_seq[i].hits[0]=1;

		_r_seq[i].alphCnt[0] = alphCnt[0];
		_r_seq[i].alphCnt[1] = alphCnt[1];
		_r_seq[i].alphCnt[2] = alphCnt[2];
		_r_seq[i].alphCnt[3] = alphCnt[3];

		compressSequence(_r_seq[i].seq, SEQ_LENGTH, _r_seq[i].cseq);
		compressSequence(_r_seq[i].rseq, SEQ_LENGTH, _r_seq[i].crseq);

		if (_r_seq[i].hits[0] == 1)			// marked reads are not indexed
		{
			_r_seq[i].hits[0] = 0;
			for (j=0; j< 2*_r_samplingLocsSize; j++)
			{
				tmp[pos].hv = -1;
				tmp[pos].checksum = 0;
				tmp[pos].seqInfo = pos +(div*id*2*_r_samplingLocsSize);
				pos++;
			}
		}
		else
		{
			for (j=0; j< _r_samplingLocsSize; j++)
			{
				hvtmp = hashVal(_r_seq[i].seq+_r_samplingLocs[j]);
				cstmp = checkSumVal(_r_seq[i].seq+_r_samplingLocs[j]+WINDOW_SIZE);
				if (hvtmp == -1 || cstmp == -1)
				{
					tmp[pos].hv = -1;
					tmp[pos].checksum = 0;
				}
				else
				{
					tmp[pos].hv = hvtmp;
					tmp[pos].checksum = cstmp;
				}
				tmp[pos].seqInfo = pos +(div*id*2*_r_samplingLocsSize);
				pos++;
			}

			for (j=0; j<_r_samplingLocsSize; j++)
			{
				hvtmp = hashVal(_r_seq[i].rseq+_r_samplingLocs[j]);
				cstmp = checkSumVal(_r_seq[i].rseq+_r_samplingLocs[j]+WINDOW_SIZE);
				
				if (hvtmp == -1  || cstmp == -1)
				{
					tmp[pos].hv = -1;
					tmp[pos].checksum = 0;
				}
				else
				{
					tmp[pos].hv = hvtmp;
					tmp[pos].checksum = cstmp;
				}
				tmp[pos].seqInfo = pos+(div*id*2*_r_samplingLocsSize);
				pos++;
			}

		}
		tmpSize+=2*_r_samplingLocsSize;
	}
	
	introSortPair( tmp, 0, tmpSize-1);

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

	_r_readIndexSize[id] = uniq;
	_r_readIndex[id] = getMem(sizeof(ReadIndexTable)*_r_readIndexSize[id]);


	prev = -2;
	j=0;
	beg =0;
	while (beg < tmpSize)
	{
		end = beg;
		while (end+1<tmpSize && tmp[end+1].hv==tmp[beg].hv)
			end++;

		_r_readIndex[id][j].hv = tmp[beg].hv;

		_r_readIndex[id][j].list = getMem(sizeof(GeneralIndex)*(end-beg+2));
		_r_readIndex[id][j].list[0].info = end-beg+1;

		for (i=1; i <= _r_readIndex[id][j].list[0].info; i++)
		{
			_r_readIndex[id][j].list[i].info=tmp[beg+i-1].seqInfo;
			_r_readIndex[id][j].list[i].checksum=tmp[beg+i-1].checksum;
		}

		j++;
		beg = end+1;
	}
	freeMem(tmp, sizeof(Pair)*(div*_r_samplingLocsSize*2));
	return NULL;
}
/**********************************************/
void preProcessReadsMT()
{
	_r_readIndexSize = getMem(sizeof(int)*THREAD_COUNT);
	_r_readIndex = getMem(sizeof(ReadIndexTable*)*THREAD_COUNT);

	_r_threads = getMem(sizeof(pthread_t)*THREAD_COUNT);
	int i=0; 
	for (i=0; i<THREAD_COUNT; i++)
		pthread_create(_r_threads+i, NULL, (void*)preProcessReads, THREAD_ID+i);
	
	for (i=0; i<THREAD_COUNT; i++)
		pthread_join(_r_threads[i], NULL);
	freeMem(_r_threads, sizeof(pthread_t)*THREAD_COUNT);	
}

/**********************************************/
void calculateSamplingLocations()
{
	int i;
	_r_samplingLocsSize = errThreshold + 1;
	_r_samplingLocs = getMem(sizeof(int)*(_r_samplingLocsSize+1));
	for (i=0; i<_r_samplingLocsSize; i++)
	{
		_r_samplingLocs[i] = (SEQ_LENGTH / _r_samplingLocsSize) *i;
		if ( _r_samplingLocs[i] + WINDOW_SIZE > SEQ_LENGTH)
			_r_samplingLocs[i] = SEQ_LENGTH - WINDOW_SIZE;
	}
	_r_samplingLocs[_r_samplingLocsSize]=SEQ_LENGTH;
	
	int size = sizeof(int)*_r_samplingLocsSize;
	_r_samplingLocsSeg = getMem(size);
	_r_samplingLocsOffset = getMem(size);
	_r_samplingLocsLen = getMem(size);
	_r_samplingLocsLenFull = getMem(size);
	for (i=0; i<_r_samplingLocsSize; i++)
	{
		_r_samplingLocsSeg[i]		= _r_samplingLocs[i] / (sizeof(CompressedSeq)*8/3);
		_r_samplingLocsOffset[i]	= _r_samplingLocs[i] % (sizeof(CompressedSeq)*8/3);
		_r_samplingLocsLen[i]		= _r_samplingLocs[i+1] - _r_samplingLocs[i];
		_r_samplingLocsLenFull[i]	= SEQ_LENGTH - _r_samplingLocs[i];
	}
	

	// Outputing the sampling locations
	/*int j;
 	for (i=0; i<SEQ_LENGTH; i++)
	{
		fprintf(stderr, "-");
	}
	fprintf(stderr, "\n");

	for ( i=0; i<_r_samplingLocsSize; i++ )
	{
		for ( j=0; j<_r_samplingLocs[i]; j++ )
			fprintf(stderr," ");
		for (j=0; j<WINDOW_SIZE; j++)
			fprintf(stderr,"+");
		fprintf(stderr, "\n");
		fflush(stderr);
	}
	for ( i=0; i<SEQ_LENGTH; i++ )
	{
		fprintf(stderr, "-");
	}
	fprintf(stderr, "\n"); */
}

/**********************************************/
int initRead(char *fileName1, char *fileName2)
{
	char dummy[SEQ_MAX_LENGTH];
	char ch;
	int i, maxCnt=0;

	_r_buf1 = getMem(10000000);
	_r_buf1_pos = getMem(sizeof(int));
	_r_buf1_size = getMem(sizeof(int));
	*_r_buf1_size = *_r_buf1_pos = 0; 
	if ( pairedEndMode && fileName2 != NULL )
	{
		_r_buf2 = getMem(10000000);
		_r_buf2_pos = getMem(sizeof(int));
		_r_buf2_size = getMem(sizeof(int));
	}
	else
	{
		_r_buf2 = _r_buf1;
		_r_buf2_pos = _r_buf1_pos;
		_r_buf2_size = _r_buf1_size;
	}


	if (!seqCompressed)
	{
		_r_fp1 = fileOpen( fileName1, "r");

		if (_r_fp1 == NULL)
			return 0;

		ch = fgetc(_r_fp1);

		if ( pairedEndMode) 
		{
			if ( fileName2 == NULL )
			{
				_r_fp2 = _r_fp1;
			}
			else
			{
				_r_fp2 = fileOpen ( fileName2, "r" );
				if (_r_fp2 == NULL)
					return 0;
			}
		}

		readBuffer1 = &readBufferTxT1;
		readBuffer2 = &readBufferTxT2;
	}
	else
	{

		_r_gzfp1 = fileOpenGZ (fileName1, "r");

		if (_r_gzfp1 == NULL)
		{
			return 0;
		}

		ch = gzgetc(_r_gzfp1);

		if ( pairedEndMode && fileName2 != NULL )
		{
			_r_gzfp2 = fileOpenGZ ( fileName2, "r" );
			if (_r_gzfp2 == NULL)
			{
				return 0;
			}
		}
		else
		{
			_r_gzfp2 = _r_gzfp1;
		}

		readBuffer1 = &readBufferGZ1;
		readBuffer2 = &readBufferGZ2;
	}

	if (!seqCompressed)
		rewind(_r_fp1);
	else
		gzrewind(_r_gzfp1);

	if (ch == '>')
		_r_fastq = 0;
	else
		_r_fastq = 1;
	
	readFirstSeq(dummy,1);
	int nameLen = strlen(dummy);
	readFirstSeq(dummy,2);
	*_r_buf1_pos = 0;
	int seqLen = strlen(dummy);
	SEQ_LENGTH = 0;
	i = 0;
	while (i<seqLen && !isspace(dummy[i]))
	{
		i++;
		SEQ_LENGTH++;
	}
	
	_r_uncroppedSeqLength = SEQ_LENGTH;

	if (cropSize > 0)
		SEQ_LENGTH = cropSize;
	if (tailCropSize > 0)
		SEQ_LENGTH = tailCropSize;

	if ( SEQ_LENGTH >= SEQ_MAX_LENGTH )
	{
		fprintf(stderr, "ERR: Read Length is greater than the MAX length we can process (Current Max: %d).\n", SEQ_MAX_LENGTH);
		exit(EXIT_FAILURE);
	}

	if (_r_fastq)
	{
		QUAL_LENGTH = SEQ_LENGTH;
	}
	else
	{
		QUAL_LENGTH = 1;
	}

	CMP_SEQ_LENGTH = calculateCompressedLen(SEQ_LENGTH);

	//TODO MEMORY CALCULATION FIX
	double readMem = sizeof(Read) + (2 + (SEQ_LENGTH * 2) + QUAL_LENGTH + 3 + (CMP_SEQ_LENGTH * 2 * 8) + (nameLen+10) + 4);
	readMem += ((bestMappingMode) ?(sizeof(FullMappingInfo)) :0);
	if (pairedEndMode)
		readMem += sizeof(MappingInfo) + sizeof(MappingLocations);

	_r_maxSeqCnt = (int)(((MAX_MEMORY-1.2) * (1 << 30))/readMem);
	if ( pairedEndMode && _r_maxSeqCnt % 2 )
		_r_maxSeqCnt ++;
	_r_maxSeqCnt -= _r_maxSeqCnt % THREAD_COUNT;

//_r_maxSeqCnt = 500000;

	_r_seq = getMem(sizeof(Read)*_r_maxSeqCnt);

	int maxErrThreshold = (SEQ_LENGTH/WINDOW_SIZE) - 1;
	if (errThreshold == -1)
	{
		errThreshold = SEQ_LENGTH*6/100;
		fprintf(stderr, "# Errors: %d\n", errThreshold);
	}
	if (errThreshold > maxErrThreshold && SEQ_LENGTH>0)
	{
		errThreshold = maxErrThreshold;
		fprintf(stderr, "# Error: %d (full sensitivity)\n", errThreshold);
	}


	checkSumLength = (SEQ_LENGTH / (errThreshold+1)) - WINDOW_SIZE;
	if (checkSumLength > sizeof(CheckSumType)*4)
		checkSumLength = sizeof(CheckSumType)*4;

	calculateSamplingLocations();


	if (!nohitDisabled)
	{
	  if (!outCompressed)
	    _r_umfp = fileOpen(unmappedOutput, "w");
	  else 
	    _r_gzumfp = fileOpenGZ(unmappedOutput, "w");
	}

	_r_alphIndex = getMem(128);		// used in readChunk()
	_r_alphIndex['A'] = 0;
	_r_alphIndex['C'] = 1;
	_r_alphIndex['G'] = 2;
	_r_alphIndex['T'] = 3;
	_r_alphIndex['N'] = 4;

	return 1;
}

/**********************************************/
int readChunk(Read **seqList, unsigned int *seqListSize)
{
	double startTime=getTime();

	char seq1[SEQ_MAX_LENGTH];
	char name1[SEQ_MAX_LENGTH];
	char qual1[SEQ_MAX_LENGTH];

	char seq2[SEQ_MAX_LENGTH];
	char name2[SEQ_MAX_LENGTH];
	char qual2[SEQ_MAX_LENGTH];


	char dummy[SEQ_MAX_LENGTH];
	int  size;

	int maxCnt = 0;
	_r_seqCnt = 0;
	_r_readMemUsage = 0;
	
	int i;//, len;

	int namelen;
	int eliminate = 0;
	
	while( (namelen = readFirstSeq(name1,1)) )
	{

		if (pairedEndMode)
		{
			if (name1[namelen-2]=='/' && name1[namelen-1]=='1')
			{
				namelen -= 2;
				name1[namelen]='\0';
			}
		}
		size = sizeof(uint16_t) + (SEQ_LENGTH * 2) + QUAL_LENGTH + 3 + (CMP_SEQ_LENGTH << 4) + namelen +/* 1 +*/ 4;	
		_r_seq[_r_seqCnt].hits	= getMem(size);		
		_r_readMemUsage += size;		
		_r_seq[_r_seqCnt].seq	= (char *)(_r_seq[_r_seqCnt].hits + 1);
		_r_seq[_r_seqCnt].rseq	= (char *)(_r_seq[_r_seqCnt].seq + SEQ_LENGTH + 1);
		_r_seq[_r_seqCnt].qual	= (char *)(_r_seq[_r_seqCnt].rseq + SEQ_LENGTH + 1);
		_r_seq[_r_seqCnt].cseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].qual + QUAL_LENGTH + 1);
		_r_seq[_r_seqCnt].crseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].cseq + CMP_SEQ_LENGTH);
		_r_seq[_r_seqCnt].name	= (char *)(_r_seq[_r_seqCnt].crseq + CMP_SEQ_LENGTH);
		_r_seq[_r_seqCnt].alphCnt = (unsigned char *)(_r_seq[_r_seqCnt].name + namelen);// + 1);
		_r_seq[_r_seqCnt].hits[0] = 0;

		for (i=1; i<namelen+1; i++)
			_r_seq[_r_seqCnt].name[i-1] = name1[i];

		eliminate = 0;
		if ( readFirstSeq(_r_seq[_r_seqCnt].seq,2) != SEQ_LENGTH)
		{
		        if (cropSize > 0)
			  eliminate = 1;
			else
			  {
			    fprintf(stderr, "ERR: Inconsistent read length for %s\n", name1);
			    exit(EXIT_FAILURE);
			  }
		} 

		if ( _r_fastq )
		{
			readFirstSeq(dummy,3);
			readFirstSeq(_r_seq[_r_seqCnt].qual,4);
		}
		else
		{
			_r_seq[_r_seqCnt].qual = "*";
		}
		
		if (eliminate)
		        freeMem(_r_seq[_r_seqCnt].hits, size);
		else
		        _r_seqCnt++;

		if (pairedEndMode)
		{
			_r_seq[_r_seqCnt].hits	= getMem(size);		
			_r_readMemUsage += size;		
			_r_seq[_r_seqCnt].seq	= (char *) (_r_seq[_r_seqCnt].hits + 1);
			_r_seq[_r_seqCnt].rseq	= (char *)(_r_seq[_r_seqCnt].seq + SEQ_LENGTH + 1);
			_r_seq[_r_seqCnt].qual	= (char *)(_r_seq[_r_seqCnt].rseq + SEQ_LENGTH + 1);
			_r_seq[_r_seqCnt].cseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].qual + QUAL_LENGTH + 1);
			_r_seq[_r_seqCnt].crseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].cseq + CMP_SEQ_LENGTH);
			_r_seq[_r_seqCnt].name	= (char *)(_r_seq[_r_seqCnt].crseq + CMP_SEQ_LENGTH);
			_r_seq[_r_seqCnt].alphCnt = (unsigned char *)(_r_seq[_r_seqCnt].name + namelen);// + 1);
			_r_seq[_r_seqCnt].hits[0] = 0;
			
			readSecondSeq(name2, 1);

			for (i=1; i<namelen+1; i++)
				_r_seq[_r_seqCnt].name[i-1] = name1[i];

			if ( readSecondSeq(_r_seq[_r_seqCnt].seq,2) != SEQ_LENGTH)
			{
				fprintf(stderr, "ERR: Inconsistent read length for %s\n", name1);
				exit(EXIT_FAILURE);			
			} 

			if ( _r_fastq )
			{
				readSecondSeq(dummy,3);
				readSecondSeq(_r_seq[_r_seqCnt].qual,4);
			}
			else
			{
				_r_seq[_r_seqCnt].qual = "*";
			}
			_r_seqCnt++;
		}

		if (_r_seqCnt >= _r_maxSeqCnt)
			break;
	}
	*seqList = _r_seq;
	*seqListSize = _r_seqCnt;

	if (_r_seqCnt > 0)
	{
		preProcessReadsMT();
		fprintf(stderr, "| *Reading Input* | %15.2f | XXXXXXXXXXXXXXX | %15.2f | XXXXXXXXXXXXXXX %15d |\n", (getTime()-startTime), getMemUsage(), _r_seqCnt );
		_r_firstIteration = 0;
	}
	else if (_r_firstIteration)
	{
		fprintf(stderr, "ERR: No reads for mapping\n");
		exit(EXIT_FAILURE);
	}

	if (_r_seqCnt < _r_maxSeqCnt)		// reached end of file
		return 0;
	else
		return 1;
}
/**********************************************/
void outputUnmapped()
{
	if (nohitDisabled)
		return;

	if (pairedEndMode)
		_r_seqCnt /=2;

	int i = 0;
	for (i = 0; i < _r_seqCnt; i++)
	{
		if (pairedEndMode)
		{
			if (_r_seq[2*i].hits[0] == 0 && _r_fastq)
			{
			        if (!outCompressed)
				  fprintf(_r_umfp,"@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].qual, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[i*2+1].qual);
				else
				  gzprintf(_r_gzumfp,"@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].qual, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[i*2+1].qual);
			}
			else if (_r_seq[2*i].hits[0] == 0)
			{
			        if (!outCompressed)
				  fprintf(_r_umfp, ">%s/1\n%s\n>%s/2\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].name, _r_seq[i*2+1].seq);
				else
				  gzprintf(_r_gzumfp, ">%s/1\n%s\n>%s/2\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].name, _r_seq[i*2+1].seq);
			}
		}
		else
		{
			if (_r_seq[i].hits[0] == 0 && _r_fastq)
			{
			        if (!outCompressed)
				  fprintf(_r_umfp,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);
				else
				  gzprintf(_r_gzumfp,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);				
			}
			else if (_r_seq[i].hits[0] == 0)
			{
			        if (!outCompressed)
				  fprintf(_r_umfp,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
				else 
				  gzprintf(_r_gzumfp,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
			}
		}
	}

	if (pairedEndMode)
		_r_seqCnt *= 2;
}
/**********************************************/
void releaseChunk()
{
	outputUnmapped();

	int i, j;
	for (i = 0; i < _r_seqCnt; i++)
		freeMem(_r_seq[i].hits, 0);
	memUsage -= _r_readMemUsage;
	_r_readMemUsage = 0;

	for (i = 0; i < THREAD_COUNT; i++)
	{
		for (j = 0; j < _r_readIndexSize[i]; j++)
			freeMem(_r_readIndex[i][j].list, (_r_readIndex[i][j].list[0].info+1)*sizeof(GeneralIndex));
		freeMem(_r_readIndex[i], sizeof(ReadIndexTable)*_r_readIndexSize[i]);
	}
	freeMem(_r_readIndex, sizeof(ReadIndexTable*)*THREAD_COUNT);
	freeMem(_r_readIndexSize, sizeof(int)*THREAD_COUNT);
}
/**********************************************/
void getSamplingLocsInfo(int **samplingLocs, int **samplingLocsSeg, int **samplingLocsOffset, int **samplingLocsLen, int **samplingLocsLenFull, int *samplingLocsSize)
{
	*samplingLocs = _r_samplingLocs;
	*samplingLocsSeg = _r_samplingLocsSeg;
	*samplingLocsOffset = _r_samplingLocsOffset;
	*samplingLocsLen = _r_samplingLocsLen;
	*samplingLocsLenFull = _r_samplingLocsLenFull;
	*samplingLocsSize = _r_samplingLocsSize;
}

/**********************************************/
void finalizeReads()
{
	if (!seqCompressed)
	{
		fclose(_r_fp1);
		if ( pairedEndMode && _r_fp2 != _r_fp1 )
		{
			fclose(_r_fp2);
		}
	}
	else
	{
		gzclose(_r_gzfp1);
		if ( pairedEndMode && _r_gzfp2 != _r_gzfp1)
		{
			gzclose(_r_gzfp2);
		}
	}
	freeMem(_r_seq, sizeof(Read)*_r_maxSeqCnt);
	freeMem(_r_samplingLocs, sizeof(int)*(_r_samplingLocsSize+1));
	int size = sizeof(int)*_r_samplingLocsSize;
	freeMem(_r_samplingLocsSeg, size);
	freeMem(_r_samplingLocsOffset, size);
	freeMem(_r_samplingLocsLen, size);
	freeMem(_r_samplingLocsLenFull, size);
	freeMem(_r_alphIndex, 128);

	if (pairedEndMode && _r_buf1 != _r_buf2)
	{
		freeMem(_r_buf2, 10000000);
		freeMem(_r_buf2_pos, sizeof(int));
		freeMem(_r_buf2_size, sizeof(int));
	}
	freeMem(_r_buf1, 10000000);
	freeMem(_r_buf1_pos, sizeof(int));
	freeMem(_r_buf1_size, sizeof(int));

	if (!nohitDisabled)
	{
	        if (!outCompressed)
		  fclose(_r_umfp);
		else
		  gzclose(_r_gzumfp);	    
	}
}

/**********************************************/
/*int checkAllReads()
{
	char seq[SEQ_MAX_LENGTH];
	char name[SEQ_MAX_LENGTH];
	FILE *fp = _r_fp1;
	int flag, firstIteration = 1;
	int i = 0, seqCnt[2];
	seqCnt[0] = seqCnt[1] = 0;
	
	do {
		while (readSeq(name, fp))		// name
		{
			readSeq(seq, fp);			// seq
			seqCnt[i]++;
			
			if ( strlen(seq)-1 != SEQ_LENGTH )
			{
				rewind(fp);
				fprintf(stderr, "ERR: Inconsistent read length %s", name);
				return 0;
			}
			if (_r_fastq)
			{
				readSeq(seq, fp);		// 3rd line
				readSeq(seq, fp);		// qual
			}
		}

		if (firstIteration && pairedEndMode && _r_fp2 != _r_fp1)
		{
			flag = 1;
			rewind(_r_fp1);
			fp = _r_fp2;
			firstIteration = 0;
			i++;
		}
		else
		{
			flag = 0;
		}

	} while (flag);
	
	rewind(fp);

	if (pairedEndMode)
	{
		if (_r_fp1 == _r_fp2)
		{
			if (seqCnt[0] & 1)
			{
				fprintf(stderr, "ERR: In paired-end mode, number of reads must be divisible by 2\n");
				return 0;
			}
		}
		else
		{
			if (seqCnt[0] != seqCnt[1])
			{
				fprintf(stderr, "ERR: Number of reads must be equal in the input files\n");
				return 0;
			}
		}
	}
	
	fprintf(stdout, "Input check: OK\n");
	return 1;
}
*/
