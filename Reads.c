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

FILE *_r_fp1;
FILE *_r_fp2;
FILE *_r_umfp;
gzFile _r_gzfp1;
gzFile _r_gzfp2;
Read *_r_seq;
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

/**********************************************/
char *(*readFirstSeq)(char *);
char *(*readSecondSeq)(char *);
/**********************************************/
char *readFirstSeqTXT( char *seq )
{
	return fgets(seq, SEQ_MAX_LENGTH, _r_fp1);
}

/**********************************************/
char *readSecondSeqTXT( char *seq )
{
	return fgets(seq, SEQ_MAX_LENGTH, _r_fp2);
}
/**********************************************/
	
char *readFirstSeqGZ( char *seq )
{
	return gzgets(_r_gzfp1, seq, SEQ_MAX_LENGTH);
}

/**********************************************/
char *readSecondSeqGZ( char *seq )
{
	return gzgets(_r_gzfp2, seq, SEQ_MAX_LENGTH);
}
/**********************************************/
int compare (const void *a, const void *b)
{
	return ((Pair *)a)->hv - ((Pair *)b)->hv;
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
	int i=0; 
	int j=0;
	int pos = 0;
	char rseq[SEQ_LENGTH+1];
	int tmpSize;

	int div = _r_seqCnt / THREAD_COUNT ; // Should be fixed for PE
	div += (_r_seqCnt % THREAD_COUNT)?1:0;
	Pair *tmp = getMem(sizeof(Pair)*(div * _r_samplingLocsSize*2));

	tmpSize = 0;
	pos = 0;
	for (i=id*div; i<div*(id+1) && i<_r_seqCnt; i++)
	{
		for (j=0; j< _r_samplingLocsSize; j++)
		{
			tmp[pos].hv = hashVal(_r_seq[i].seq+_r_samplingLocs[j]);
			tmp[pos].seqInfo = pos +(div*id*2*_r_samplingLocsSize);
			pos++;
		}
		for (j=0; j<_r_samplingLocsSize; j++)
		{
			reverseComplete(_r_seq[i].seq, rseq, SEQ_LENGTH);
			tmp[pos].hv = hashVal(rseq+_r_samplingLocs[j]);
			tmp[pos].seqInfo = pos+(div*id*2*_r_samplingLocsSize);
			pos++;
		}
		tmpSize+=2*_r_samplingLocsSize;
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
		_r_readIndex[id][j].seqInfo = getMem(sizeof(int)*(end-beg+2));
		_r_readIndex[id][j].seqInfo[0] = end-beg+1;

		for (i=1; i <= _r_readIndex[id][j].seqInfo[0]; i++)
		{
			_r_readIndex[id][j].seqInfo[i]=tmp[beg+i-1].seqInfo;
		}
		j++;
		beg = end+1;
	}
	freeMem(tmp, sizeof(Pair)*(div*_r_samplingLocsSize*2));
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
/*	int j;
 	for (i=0; i<SEQ_LENGTH; i++)
	{
		fprintf(stdout, "-");
	}
	fprintf(stdout, "\n");

	for ( i=0; i<_r_samplingLocsSize; i++ )
	{
		for ( j=0; j<samLocs[i]; j++ )
			fprintf(stdout," ");
		for (j=0; j<WINDOW_SIZE; j++)
			fprintf(stdout,"+");
		fprintf(stdout, "\n");
		fflush(stdout);
	}
	for ( i=0; i<SEQ_LENGTH; i++ )
	{
		fprintf(stdout, "-");
	}
	fprintf(stdout, "\n");*/
}
/**********************************************/

int initRead(char *fileName1, char *fileName2)
{
	char dummy[SEQ_MAX_LENGTH];
	char ch;
	int maxCnt=0;

	if (!seqCompressed)
	{
		_r_fp1 = fileOpen( fileName1, "r");

		if (_r_fp1 == NULL)
			return 0;

		ch = fgetc(_r_fp1);
		if ( pairedEndMode && fileName2 != NULL )
		{
			_r_fp2 = fileOpen ( fileName2, "r" );
			if (_r_fp2 == NULL)
			{
				return 0;
			}
		}
		else
		{
			_r_fp2 = _r_fp1;
		}

		readFirstSeq = &readFirstSeqTXT;
		readSecondSeq = &readSecondSeqTXT;
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
			_r_fp2 = fileOpenGZ ( fileName2, "r" );
			if (_r_fp2 == NULL)
			{
				return 0;
			}
		}
		else
		{
			_r_fp2 = _r_fp1;
		}

		readFirstSeq = &readFirstSeqGZ;
		readSecondSeq = &readSecondSeqGZ;
	}

	if (ch == '>')
		_r_fastq = 0;
	else
		_r_fastq = 1;
	
	readFirstSeq(dummy);
	int nameLen = strlen(dummy);
	readFirstSeq(dummy);
	int seqLen = strlen(dummy);
	SEQ_LENGTH = seqLen - 1;
	int cmpLen = calculateCompressedLen(seqLen);
	double readMem = (2 + (seqLen * 3) + 3 + (cmpLen << 4) + nameLen + 1 + 4);
	if (bestMappingMode)
		readMem += ((bestMappingMode) ?(sizeof(FullMappingInfo)) :0);
	if (pairedEndMode)
		readMem += sizeof(MappingInfo) + 2*sizeof(MappingLocations);

	_r_maxSeqCnt = (int)(((MAX_MEMORY-1) * (1 << 30))/readMem);
	if ( pairedEndMode && _r_maxSeqCnt % 2 )
		_r_maxSeqCnt ++;
	_r_maxSeqCnt -= _r_maxSeqCnt % THREAD_COUNT;

//_r_maxSeqCnt =500000;

	if (!seqCompressed)
	{
		rewind(_r_fp1);
	}
	else
	{
		gzrewind(_r_gzfp1);
	}
	
	_r_seq = getMem(sizeof(Read)*_r_maxSeqCnt);

	int maxErrThreshold = (SEQ_LENGTH/WINDOW_SIZE) - 1;
	if (errThreshold == -1)
	{
		errThreshold = SEQ_LENGTH*6/100;
		fprintf(stdout, "# Errors: %d\n", errThreshold);
	}
	if (errThreshold > maxErrThreshold && SEQ_LENGTH>0)
	{
		errThreshold = maxErrThreshold;
		fprintf(stdout, "# Error: %d (full sensitivity)\n", errThreshold);
	}




	calculateSamplingLocations();



	fprintf(stdout, "%d max num of reads %d\n", SEQ_LENGTH, _r_maxSeqCnt);
	

	char *t = unmappedOutput;
	unmappedOutput = getMem(100);
	strcpy(unmappedOutput, t);
	if (strcmp(unmappedOutput, "") == 0)
		sprintf(unmappedOutput, "%s%s.nohit", mappingOutputPath, mappingOutput );
	_r_umfp = fopen(unmappedOutput, "w");

	return 1;
}


/**********************************************/
int readChunk(Read **seqList, unsigned int *seqListSize)
{
	double startTime=getTime();

	char seq1[SEQ_MAX_LENGTH];
	CompressedSeq cseq1[CMP_SEQ_MAX_LENGTH];		// cmp of seq
	char name1[SEQ_MAX_LENGTH];
	char qual1[SEQ_MAX_LENGTH];

	char seq2[SEQ_MAX_LENGTH];
	CompressedSeq cseq2[CMP_SEQ_MAX_LENGTH];		// cmp of seq2
	char name2[SEQ_MAX_LENGTH];
	char qual2[SEQ_MAX_LENGTH];

	char rseq[SEQ_MAX_LENGTH];
	CompressedSeq crseq1[CMP_SEQ_MAX_LENGTH];		// cmp of rseq1

	char dummy[SEQ_MAX_LENGTH];
	unsigned char alphCnt[5], alphCnt2[5];
	char ch;
	int err1, err2, size;
	//int nCnt;
	unsigned char *nCnt;
	int discarded = 0;
	int maxCnt = 0;
	_r_seqCnt = 0;
	_r_readMemUsage = 0;
	
	int i, len;
	unsigned int *setZero;
	char alphIndex[128];
	alphIndex['A'] = 0;
	alphIndex['C'] = 1;
	alphIndex['G'] = 2;
	alphIndex['T'] = 3;
	alphIndex['N'] = 4;


	while( readFirstSeq(name1) )
	{
		err1 = 0;
		err2 = 0;
		readFirstSeq(seq1);
		name1[strlen(name1)-1] = '\0';
		for (i=0; i<strlen(name1);i++)
		{
			if (name1[i] == ' ')
			{
				name1[i] = '\0';
				break;
			}

		}

		if ( _r_fastq )
		{
			readFirstSeq(dummy);
			readFirstSeq(qual1);
			qual1[strlen(qual1)-1] = '\0';

		}
		else
		{
			sprintf(qual1, "*");
		}

		// Cropping
		if (cropSize > 0)
		{
			seq1[cropSize] = '\0';
			if ( _r_fastq )
				qual1[cropSize] = '\0';
		}


		nCnt = alphCnt + 4;
		*nCnt = 0;
		setZero = (unsigned int *)alphCnt;
		*setZero = 0;
		len = strlen(seq1);

		for (i = 0; i < len; i++)
		{
			seq1[i] = toupper (seq1[i]);
			
			if (isspace(seq1[i]))
			{
				seq1[i] = '\0';
				break;
			}
			else
			{
				alphCnt[alphIndex[seq1[i]]]++;
			}
		}

		if (*nCnt > errThreshold)
		{
			err1 = 1;
		}

		// Reading the second seq of pair-ends
		if (pairedEndMode)
		{
			readSecondSeq(name2);
			readSecondSeq(seq2);
			name2[strlen(name2)-1] = '\0';
			for (i=0; i<strlen(name2);i++)
			{
				if (name2[i] == ' ')
				{
					name2[i] = '\0';
					break;
				}

			}

			if ( _r_fastq )
			{
				readSecondSeq(dummy);
				readSecondSeq(qual2);

				qual2[strlen(qual2)-1] = '\0';
			}
			else
			{
				sprintf(qual2, "*");
			}


			// Cropping
			if (cropSize > 0)
			{
				seq2[cropSize] = '\0';
				if ( _r_fastq )
					qual2[cropSize] = '\0';
			}

			nCnt = alphCnt2 + 4;
			*nCnt = 0;
			setZero = (unsigned int *)alphCnt2;
			*setZero = 0;
			len = strlen(seq2);

			for (i = 0; i < len; i++)
			{
				seq2[i] = toupper (seq2[i]);
				
				if (isspace(seq2[i]))
				{
					seq2[i] = '\0';
					break;
				}
				else
				{
					alphCnt2[alphIndex[seq2[i]]]++;
				}
			}

			if (*nCnt > errThreshold)
			{
				err2 = 1;
			}
		}

		if (!pairedEndMode && !err1)
		{
			int _mtmp = strlen(seq1);
			int cmpLen = calculateCompressedLen(_mtmp);
			int namelen = strlen(name1);
			size = 2 + (_mtmp * 3) + 3 + (cmpLen << 4) + namelen + 1 + 4;
			_r_seq[_r_seqCnt].hits	= getMem(size);
			_r_readMemUsage += size;
			_r_seq[_r_seqCnt].seq	= (char *) (_r_seq[_r_seqCnt].hits + 1);
			_r_seq[_r_seqCnt].rseq	= (char *)(_r_seq[_r_seqCnt].seq + _mtmp + 1);
			_r_seq[_r_seqCnt].qual	= (char *)(_r_seq[_r_seqCnt].rseq + _mtmp + 1);
			_r_seq[_r_seqCnt].cseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].qual + _mtmp + 1);
			_r_seq[_r_seqCnt].crseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].cseq + cmpLen);
			_r_seq[_r_seqCnt].name	= (char *)(_r_seq[_r_seqCnt].crseq + cmpLen);
			_r_seq[_r_seqCnt].alphCnt = (char *)(_r_seq[_r_seqCnt].name + namelen + 1);

			for (i = 0; i < 4; i++)
				_r_seq[_r_seqCnt].alphCnt[i] = alphCnt[i];

			reverseComplete(seq1, rseq, _mtmp);
			rseq[_mtmp] = '\0';
			compressSequence(seq1, _mtmp, _r_seq[_r_seqCnt].cseq);
			compressSequence(rseq, _mtmp, _r_seq[_r_seqCnt].crseq);

			int i;

			_r_seq[_r_seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				_r_seq[_r_seqCnt].seq[i] = seq1[i];
				_r_seq[_r_seqCnt].rseq[i] = rseq[i];
				_r_seq[_r_seqCnt].qual[i] = qual1[i];
			}
			_r_seq[_r_seqCnt].qual[_mtmp] = '\0';
			sprintf(_r_seq[_r_seqCnt].name,"%s%c", ((char*)name1)+1,'\0');

			_r_seqCnt++;
		}
		else if (pairedEndMode && !err1 && !err2)
		{
			// Naming Conventions X/1, X/2 OR X
			int tmplen = strlen(name1);
			if (strcmp(name1, name2) != 0)
			{
				tmplen = strlen(name1)-2;
			}
		
			//first seq
			int _mtmp = strlen(seq1);	
			int cmpLen = calculateCompressedLen(_mtmp);
			size = 2 + (_mtmp * 3) + 3 + (cmpLen << 4) + tmplen + 1 + 4;
			_r_seq[_r_seqCnt].hits	= getMem(size);
			_r_readMemUsage += size;
			_r_seq[_r_seqCnt].seq	= (char *) (_r_seq[_r_seqCnt].hits + 1);
			_r_seq[_r_seqCnt].rseq	= (char *)_r_seq[_r_seqCnt].seq + _mtmp + 1;
			_r_seq[_r_seqCnt].qual	= (char *)_r_seq[_r_seqCnt].rseq + _mtmp + 1;
			_r_seq[_r_seqCnt].cseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].qual + _mtmp + 1);
			_r_seq[_r_seqCnt].crseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].cseq + cmpLen);
			_r_seq[_r_seqCnt].name	= (char *)(_r_seq[_r_seqCnt].crseq + cmpLen);
			_r_seq[_r_seqCnt].alphCnt = (char *)(_r_seq[_r_seqCnt].name + tmplen + 1);

			for (i = 0; i < 4; i++)
				_r_seq[_r_seqCnt].alphCnt[i] = alphCnt[i];

			reverseComplete(seq1, rseq, _mtmp);
			rseq[_mtmp] = '\0';
			compressSequence(seq1, _mtmp, _r_seq[_r_seqCnt].cseq);
			compressSequence(rseq, _mtmp, _r_seq[_r_seqCnt].crseq);

			int i;

			_r_seq[_r_seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				_r_seq[_r_seqCnt].seq[i] = seq1[i];
				_r_seq[_r_seqCnt].rseq[i] = rseq[i];
				_r_seq[_r_seqCnt].qual[i] = qual1[i];
			}
			
			name1[tmplen]='\0';
			sprintf(_r_seq[_r_seqCnt].name,"%s%c", ((char*)name1)+1,'\0');


			_r_seqCnt++;

			//second seq
			size = 2 + (_mtmp * 3) + 3 + (cmpLen << 4) + tmplen + 1 + 4;
			_r_seq[_r_seqCnt].hits	= getMem(size);
			_r_readMemUsage += size;
			_r_seq[_r_seqCnt].seq	= (char *) (_r_seq[_r_seqCnt].hits + 1);
			_r_seq[_r_seqCnt].rseq	= (char *)_r_seq[_r_seqCnt].seq + _mtmp + 1;
			_r_seq[_r_seqCnt].qual	= (char *)_r_seq[_r_seqCnt].rseq + _mtmp + 1;
			_r_seq[_r_seqCnt].cseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].qual + _mtmp + 1);
			_r_seq[_r_seqCnt].crseq	= (CompressedSeq *)(_r_seq[_r_seqCnt].cseq + cmpLen);
			_r_seq[_r_seqCnt].name	= (char *)(_r_seq[_r_seqCnt].crseq + cmpLen);
			_r_seq[_r_seqCnt].alphCnt = (char *)(_r_seq[_r_seqCnt].name + tmplen + 1);

			for (i = 0; i < 4; i++)
				_r_seq[_r_seqCnt].alphCnt[i] = alphCnt2[i];
			
			reverseComplete(seq2, rseq, _mtmp);
			rseq[_mtmp] = '\0';
			compressSequence(seq2, _mtmp, _r_seq[_r_seqCnt].cseq);
			compressSequence(rseq, _mtmp, _r_seq[_r_seqCnt].crseq);

			_r_seq[_r_seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				_r_seq[_r_seqCnt].seq[i] = seq2[i];
				_r_seq[_r_seqCnt].rseq[i] = rseq[i];
				_r_seq[_r_seqCnt].qual[i] = qual2[i];
			}

			name2[tmplen]='\0';
			sprintf(_r_seq[_r_seqCnt].name,"%s%c", ((char*)name2)+1,'\0');

			_r_seqCnt++;

		}
		else
		{
			discarded++;
		}

		if (_r_seqCnt == _r_maxSeqCnt)
			break;
	}

	if (_r_seqCnt > 0)
	{
		QUAL_LENGTH = SEQ_LENGTH = strlen(_r_seq[0].seq);
		CMP_SEQ_LENGTH = calculateCompressedLen(SEQ_LENGTH);
		if (! _r_fastq)
		{
			QUAL_LENGTH = 1;
		}
	}

/*	if (pairedEndMode)
	{
		_r_seqCnt /= 2;
	}*/


	// Closing Files
	
	*seqList = _r_seq;
	*seqListSize = _r_seqCnt;


	if (_r_seqCnt > 0)
	{
		preProcessReadsMT();
		fprintf(stdout, "%d sequences are read in %0.2f. (%d discarded) [Mem:%0.2f M]\n", _r_seqCnt, (getTime()-startTime), discarded, getMemUsage());
		_r_firstIteration = 0;
	}
	else if (_r_firstIteration)
	{
		fprintf(stdout, "ERR: No reads can be found for mapping\n");
	}

	if (_r_seqCnt < _r_maxSeqCnt)		// reached end of file
		return 0;
	else
		return 1;
}
/**********************************************/
void releaseChunk()
{
	if (pairedEndMode)
		_r_seqCnt /=2;

	int i = 0, j = 0;
	for (i = 0; i < _r_seqCnt; i++)
	{
		if (pairedEndMode && _r_seq[2*i].hits[0] == 0 && _r_fastq)
		{
			fprintf(_r_umfp,"@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].qual, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[i*2+1].qual);
		}
		else if (pairedEndMode && _r_seq[2*i].hits[0] == 0)
		{
			fprintf(_r_umfp, ">%s/1\n%s\n>%s/2\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].name, _r_seq[i*2+1].seq);
		}
		else if (_r_seq[i].hits[0] == 0 && _r_fastq)
		{
			fprintf(_r_umfp,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);
		}
		else if (_r_seq[i].hits[0] == 0)
		{
			fprintf(_r_umfp,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
		}
	}

	if (pairedEndMode)
		_r_seqCnt *= 2;

	for (i = 0; i < _r_seqCnt; i++)
		freeMem(_r_seq[i].hits,0);
	memUsage -= _r_readMemUsage;
	_r_readMemUsage = 0;


	for (i = 0; i < THREAD_COUNT; i++)
	{
		for (j = 0; j < _r_readIndexSize[i]; j++)
			freeMem(_r_readIndex[i][j].seqInfo, (_r_readIndex[i][j].seqInfo[0]+1)*sizeof(int));
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
		if ( pairedEndMode && _r_fp2 != NULL )
		{
			fclose(_r_fp2);
		}
	}
	else
	{
		gzclose(_r_gzfp1);
		if ( pairedEndMode && _r_fp2 != NULL)
		{
			gzclose(_r_fp2);
		}
	}

	freeMem(_r_seq, sizeof(Read)*_r_maxSeqCnt);
	freeMem(_r_samplingLocs, sizeof(int)*(_r_samplingLocsSize+1));
	int size = sizeof(int)*_r_samplingLocsSize;
	freeMem(_r_samplingLocsSeg, size);
	freeMem(_r_samplingLocsOffset, size);
	freeMem(_r_samplingLocsLen, size);
	freeMem(_r_samplingLocsLenFull, size);
	freeMem(unmappedOutput, 100);
	fclose(_r_umfp);
}

