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
#include <pthread.h>
#include "Common.h"
#include "RefGenome.h"
#include "HashTable.h"
#include "Output.h"
#include "Sort.h"

/**********************************************/
FILE			*_ih_fp					= NULL;
IHashTable		*_ih_hashTable			= NULL;
int				_ih_maxHashTableSize	= 0;
unsigned int	_ih_hashTableMemSize	= 0;
GeneralIndex	*_ih_hashTableMem		= NULL;
int				_ih_refGenLen			= 0;
CompressedSeq	*_ih_crefGen			= NULL;
int				_ih_crefGenLen			= 0;
char			*_ih_refGenName			= NULL;
long long		_ih_memUsage			= 0;
int				_ih_refGenOff			= 0;
unsigned char	*_ih_IOBuffer			= NULL;
unsigned int	_ih_IOBufferSize		= (1 << 24);
int				MAX_GENOME_INFO_SIZE	= 10000000;
int				_ih_maxChrLength		= 0;
CompressedSeq	*_ih_crefGenOrigin		= NULL;		// only used in pairedEndMode
unsigned char	*_ih_alphCnt			= NULL;
long int		_ih_contigStartPos		= 0;
int				_ih_chrCnt				= 0;
char			**_ih_chrNames			= 0;
pthread_t		*_ih_threads = NULL;
pthread_mutex_t	_ih_writeLock;



/**********************************************/
static inline int encodeVariableByte(unsigned char *buffer, unsigned int value)		// returns number of bytes written to buffer
{
	int t = 0;
	do {
		buffer[t++] = (unsigned char) (value & 127);
		value /= 128;
	} while (value != 0);
	buffer[t-1] |= 128;
	return t;
}
/**********************************************/
static inline unsigned int decodeVariableByte(unsigned char *buffer, unsigned int *result)		// returns number of bytes read from the buffer
{
	int i = 0;
	char t;
	*result = 0;
	do {
		t = buffer[i];
		*result |= ((t&127) <<(7*i));
		i++;
	} while ((t & 128) == 0);
	return i;
}
/**********************************************/
int compareCheckSumHT (const void *a, const void *b)
{
	return ((GeneralIndex *)a)->checksum - ((GeneralIndex *)b)->checksum;
}
/**********************************************/
void initSavingIHashTable(char *fileName, char *genomeMetaInfo, int genomeMetaInfoLength)
{
	// file header:
	// 1 byte (magicNumber): Magic number of HashTable (0: <v3, 1: bisulfite <v1.26.4, 2: >v3)
	// 1 byte (WINDOW_SIZE): Windows Size of indexing
	// 4 bytes (_ih_hsahTableMemSize): HashTbleMemSize: maximum number of elements that can be saved.
	// 4 bytes (_ih_IOBufferSize): memory required for reading hash table. In case the value is changed for loading.
	// 4 bytes (CONTIG_MAX_SIZE): maximum number of characters that can be in a contig. In case the value is changed for loading
	// n bytes (genomeMetaInfo): number of chromosomes, their names and lengths

	_ih_fp = fileOpen(fileName, "w");
	unsigned char magicNumber = 2;

	int tmp;
	tmp = fwrite(&magicNumber, sizeof(magicNumber), 1, _ih_fp);
	tmp = fwrite(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ih_fp);
	tmp = fwrite(&_ih_hashTableMemSize, sizeof(_ih_hashTableMemSize), 1, _ih_fp);
	tmp = fwrite(&_ih_IOBufferSize, sizeof(_ih_IOBufferSize), 1, _ih_fp);
	tmp = fwrite(&CONTIG_MAX_SIZE, sizeof(CONTIG_MAX_SIZE), 1, _ih_fp);
	tmp = fwrite(genomeMetaInfo, sizeof(char), genomeMetaInfoLength, _ih_fp);

	_ih_IOBuffer = getMem(_ih_IOBufferSize);
}
/**********************************************/
void finalizeSavingIHashTable()
{
	// seeking back to hashTableMemSize to update the value
	fseek(_ih_fp, 2, SEEK_SET);
	fwrite(&_ih_hashTableMemSize, sizeof(_ih_hashTableMemSize), 1, _ih_fp);

	freeMem(_ih_IOBuffer,_ih_IOBufferSize);
	fclose(_ih_fp);
}
/**********************************************/
void saveHashTable(unsigned int *hashTable,  unsigned int size, unsigned int maxSize, char *refGen, char *refGenName, int refGenOffset, unsigned char lastContig)
{
	// 1 byte (extraInfo): Reserved; in case the contig has extra information
	// 2 bytes (len): Length of the reference genome name
	// n bytes (refGenName): Reference genome name
	// 4 bytes (refGenOfsset): Offset of the contig from the beginning of the chromosome
	// 4 bytes (refGenLength): Length of reference genome
	// n bytes (crefGen): compressed reference genome
	// 4 bytes (size): number of hashValues in hashTable with more than 0 locations
	// n bytes (bufferSize and buffer): array of bufferSize/buffer which includes encoded values of hashValue, count of locations 

	int tmp, i;

	unsigned char extraInfo = lastContig;
	tmp = fwrite (&extraInfo, sizeof(extraInfo), 1, _ih_fp);

	short len = strlen(refGenName);
	tmp = fwrite(&len, sizeof(len), 1, _ih_fp);
	tmp = fwrite(refGenName, sizeof(char), len, _ih_fp);

	tmp = fwrite(&refGenOffset, sizeof(refGenOffset), 1, _ih_fp);

	unsigned int refGenLength = strlen(refGen);
	tmp = fwrite(&refGenLength, sizeof(refGenLength), 1, _ih_fp);
	
	unsigned int crefGenLength = calculateCompressedLen(refGenLength);
	CompressedSeq *crefGen = getMem(crefGenLength * sizeof(CompressedSeq));
	compressSequence(refGen, refGenLength, crefGen);
	tmp = fwrite(crefGen, sizeof(CompressedSeq), crefGenLength, _ih_fp);
	freeMem(crefGen, crefGenLength * sizeof(CompressedSeq));

	unsigned int memSize = 0;
	for (i=0; i<maxSize; i++)
	{
		if (hashTable[i])
			memSize += hashTable[i] + 1;
	}

	if (_ih_hashTableMemSize < memSize)
		_ih_hashTableMemSize = memSize;

	tmp = fwrite(&size, sizeof(size), 1, _ih_fp);

	int j=0, k = 0, prevHV = 0;
	for (i=0; i<maxSize; i++)
	{
		if (hashTable[i])
		{
			int hvDiff = i - prevHV;	// save hvDiff
			prevHV = i;
			k += encodeVariableByte(_ih_IOBuffer + k, hvDiff);
			k += encodeVariableByte(_ih_IOBuffer + k, hashTable[i]);
			if (k > _ih_IOBufferSize - 10 )
			{
				fwrite(&k, sizeof(int), 1, _ih_fp);
				fwrite(_ih_IOBuffer, sizeof(unsigned char), k, _ih_fp);
				k = 0;
			}
		}
	}
	if (k)
	{
		fwrite(&k, sizeof(int), 1, _ih_fp);
		fwrite(_ih_IOBuffer, sizeof(unsigned char), k, _ih_fp);
	}
}
/**********************************************/
int generateHashTable(char *fileName, char *indexName)
{
	double          startTime           = getTime();
	unsigned int	hashTableSize		= 0;
	int				refGenOff			= 0;
	unsigned int 	hashTableMaxSize	= (1 << 2*WINDOW_SIZE);		// 4^WINDOW_SIZE
	unsigned int	*hashTable			= getMem(hashTableMaxSize * sizeof(unsigned int));
	char			*genomeMetaInfo		= getMem(MAX_GENOME_INFO_SIZE);
	int				genomeMetaInfoLength;
	char 			*refGenName			= NULL;
	char			*refGen				= NULL;
	char			*c					= NULL;
	char			*prev = getMem (CONTIG_NAME_SIZE);
	int				i, hv, l, flag, stack , val, loc;
	char			lookup[128];
	unsigned int	windowMask			= 0xffffffff >> (sizeof(unsigned int)*8 - WINDOW_SIZE*2);

	memset(lookup, 4, 128);
	lookup['A'] = 0;
	lookup['C'] = 1;
	lookup['G'] = 2;
	lookup['T'] = 3;
	lookup['N'] = 4;


	//Loading Fasta File
	prev[0]='\0';

	if (!initLoadingRefGenome(fileName, genomeMetaInfo, &genomeMetaInfoLength))
		return 0;		
	initSavingIHashTable(indexName, genomeMetaInfo, genomeMetaInfoLength);
	fprintf(stdout, "Generating Index from %s", fileName);
	fflush(stdout);

	do
	{
		flag = 	 loadRefGenome (&refGen, &refGenName, &refGenOff, &_ih_refGenLen);

		memset(hashTable, 0, hashTableMaxSize * sizeof(unsigned int));
		hashTableSize = 0;

		if ( strcmp(prev, refGenName) != 0)
		{
			fprintf(stdout, "\n - %s ", refGenName);
			fflush(stdout);
			sprintf(prev, "%s", refGenName);
		}
		else
		{
			fprintf(stdout, ".");
			fflush(stdout);
		}
		
		c = refGen;
		i = hv = val = 0;
		stack = 1;
		loc = -WINDOW_SIZE+1;

		while (i++ < _ih_refGenLen)
		{
			loc++;
			val = lookup[*(c++)];

			if (val != 4 && stack == WINDOW_SIZE)
			{
				hv = ((hv << 2)|val)&windowMask;
				if (hashTable[hv]++ == 0)
					hashTableSize++;
			}
			else
			{
				if (val == 4)
				{
					stack = 1;
					hv = 0;
				}
				else
				{
					stack ++;
					hv = (hv <<2)|val;
				}

			}
		}

		saveHashTable(hashTable, hashTableSize, hashTableMaxSize, refGen, refGenName, refGenOff, flag);
	} while (flag);

	freeMem(prev, CONTIG_NAME_SIZE);
	freeMem(hashTable, sizeof(unsigned int)*hashTableMaxSize);
	freeMem(genomeMetaInfo, MAX_GENOME_INFO_SIZE);

	finalizeLoadingRefGenome();
	finalizeSavingIHashTable();

	fprintf(stdout, "\nDONE in %0.2fs!\n", (getTime()-startTime));
	return 1;
}
/**********************************************/
void rewindHashTable()
{
	fseek(_ih_fp, _ih_contigStartPos, SEEK_SET);
}
/**********************************************/
int checkHashTable(char *fileName)
{
	_ih_fp = fileOpen(fileName, "r");

	unsigned char magicNumber;
	int tmp;

	tmp = fread(&magicNumber, sizeof(magicNumber), 1, _ih_fp);
	if (magicNumber == 1)
	{
		fprintf(stdout, "Error: Please use version 1.2.6.4 in bisulfite mode.\n");
		return 0;
	}
	else if (magicNumber == 0)
	{
		fprintf(stdout, "Error: Please use version 2.x.x.x or upgrade your index.\n");
		return 0;
	}

	tmp = fread(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ih_fp);
	tmp = fread(&_ih_hashTableMemSize, sizeof(_ih_hashTableMemSize), 1, _ih_fp);
	tmp = fread(&_ih_IOBufferSize, sizeof(_ih_IOBufferSize), 1, _ih_fp);
	tmp = fread(&CONTIG_MAX_SIZE, sizeof(CONTIG_MAX_SIZE), 1, _ih_fp);
	fclose(_ih_fp);
	_ih_fp = NULL;
	return 1;
}

/**********************************************/
int initLoadingHashTable(char *fileName)
{
	// file header:
	// 1 byte (magicNumber): Magic number of HashTable (0: <v3, 1: bisulfite <v1.26.4, 2: >v3)
	// 1 byte (WINDOW_SIZE): Windows Size of indexing
	// 4 bytes (_ih_hsahTableMemSize): HashTbleMemSize: maximum number of elements that can be saved.
	// 4 bytes (_ih_IOBufferSize): memory required for reading hash table. In case the value is changed for loading.
	// 4 bytes (CONTIG_MAX_SIZE): maximum number of characters that can be in a contig. In case the value is changed for loading
	// n bytes (genomeMetaInfo): number of chromosomes, their names and lengths

	if (_ih_fp == NULL)		// first time
		_ih_fp = fileOpen(fileName, "r");
	else
		rewind(_ih_fp);			// rewind for mapping the next chunk of reads

	int i, numOfChrs, nameLen;
	unsigned char magicNumber;
	int tmp;

	_ih_threads = getMem(sizeof(pthread_t) * THREAD_COUNT);

	tmp = fread(&magicNumber, sizeof(magicNumber), 1, _ih_fp);
	tmp = fread(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ih_fp);
	tmp = fread(&_ih_hashTableMemSize, sizeof(_ih_hashTableMemSize), 1, _ih_fp);

	_ih_hashTableMem = getMem(_ih_hashTableMemSize*sizeof(GeneralIndex));

	tmp = fread(&_ih_IOBufferSize, sizeof(_ih_IOBufferSize), 1, _ih_fp);
	_ih_IOBuffer = getMem(_ih_IOBufferSize);

	tmp = fread(&CONTIG_MAX_SIZE, sizeof(CONTIG_MAX_SIZE), 1, _ih_fp);

	// Reading Meta
	char *strtmp = getMem(2*CONTIG_NAME_SIZE);

	tmp = fread(&_ih_chrCnt, sizeof(int), 1, _ih_fp);

	_ih_chrNames = getMem(_ih_chrCnt * sizeof(char *));
	for (i = 0; i < _ih_chrCnt; i++)
	{
		_ih_chrNames[i] = getMem(CONTIG_NAME_SIZE);
		tmp = fread(&nameLen, sizeof(int), 1, _ih_fp);
		tmp = fread(_ih_chrNames[i], sizeof(char), nameLen, _ih_fp);
		_ih_chrNames[i][nameLen] = '\0';
		tmp = fread(&_ih_refGenLen, sizeof(int), 1, _ih_fp);
		
		sprintf(strtmp,"@SQ\tSN:%s\tLN:%d%c", _ih_chrNames[i], _ih_refGenLen, '\0');
		outputMeta(strtmp);
		
		if (_ih_refGenLen > _ih_maxChrLength)
			_ih_maxChrLength = _ih_refGenLen;
	}
	freeMem(strtmp, 2*CONTIG_NAME_SIZE);
	// Reading Meta End

	if (pairedEndMode)
	{
		_ih_crefGenOrigin = getMem((calculateCompressedLen(_ih_maxChrLength)+1) * sizeof(CompressedSeq));
		_ih_crefGen = _ih_crefGenOrigin; 
	}
	else
	{
		_ih_crefGen = getMem((calculateCompressedLen(CONTIG_MAX_SIZE)+1) * sizeof(CompressedSeq));
	}

	_ih_maxHashTableSize = pow(4, WINDOW_SIZE);

	_ih_hashTable = getMem (sizeof(IHashTable) * _ih_maxHashTableSize);
	memset(_ih_hashTable, 0, _ih_maxHashTableSize * sizeof(IHashTable));
	_ih_refGenName = getMem(CONTIG_NAME_SIZE);
	_ih_refGenName[0] = '\0';
	if (!SNPMode)
		_ih_alphCnt = getMem(CONTIG_MAX_SIZE * 4);

	_ih_contigStartPos = ftell(_ih_fp);

	return 1;
}
/**********************************************/
void finalizeLoadingHashTable()
{
	int i;
	freeMem(_ih_threads, sizeof(pthread_t) * THREAD_COUNT);

	freeMem(_ih_hashTableMem, _ih_hashTableMemSize * sizeof(GeneralIndex));
	freeMem(_ih_IOBuffer, _ih_IOBufferSize);
	if (pairedEndMode)
		freeMem(_ih_crefGenOrigin, (calculateCompressedLen(_ih_maxChrLength)+1) * sizeof(CompressedSeq));
	else
		freeMem(_ih_crefGen, (calculateCompressedLen(CONTIG_MAX_SIZE)+1) * sizeof(CompressedSeq));
	freeMem(_ih_hashTable, sizeof(IHashTable)* _ih_maxHashTableSize);
	freeMem(_ih_refGenName, CONTIG_NAME_SIZE);	
	if (!SNPMode)
		freeMem(_ih_alphCnt, CONTIG_MAX_SIZE * 4);
	for (i = 0; i < _ih_chrCnt; i++)
		freeMem(_ih_chrNames[i], CONTIG_NAME_SIZE);
	freeMem(_ih_chrNames, _ih_chrCnt * sizeof(char *));
	fclose(_ih_fp);
}
/**********************************************/
void *calculateHashTableOnFly(int *idp)
{
	int id = *idp;

	int windowMaskSize = WINDOW_SIZE + checkSumLength;
	unsigned long long windowMask =   0xffffffffffffffff >> (sizeof(unsigned long long)*8 - windowMaskSize*2);
	unsigned long long checkSumMask = 0xffffffffffffffff >> (sizeof(unsigned long long)*8 - (checkSumLength)*2);
	if (checkSumLength == 0)
		checkSumMask = 0;

	CompressedSeq *cnext = (_ih_crefGen);
	CompressedSeq cdata = *(cnext++);

	int i = 0;
	unsigned long long hv = 0;
	unsigned long long hvtemp;
	int pos, val, t = 0, stack = 1;
	int loc = -WINDOW_SIZE - checkSumLength + 1 ;
	int x;
	// calculate refGen hashValues
	while (i++ < _ih_refGenLen ) // BORDER LINE CHECK
	{
		loc++;
		val = (cdata >> 60) & 7;
		if (++t == 21)
		{
			t = 0;
			cdata = *(cnext++);
		}
		else
		{
			cdata <<= 3;
		}

		if (val != 4 && stack == windowMaskSize)
		{
			hv = ((hv << 2)|val)&windowMask;
			hvtemp = hv >> (checkSumLength<<1);

			if (hvtemp % THREAD_COUNT == id)
			{
				++_ih_hashTable[hvtemp].list;
				_ih_hashTable[hvtemp].list->info= loc;
				_ih_hashTable[hvtemp].list->checksum= hv & checkSumMask;
			}
		}
		else
		{
			if (val == 4)		// N
			{
				stack = 1;
				hv = 0;
			}
			else
			{
				stack ++;
				hv = (hv <<2)|val;
			}

		}
	}
	return NULL;
}
/**********************************************/
void *sortHashTable(int *id)
{
	int cnt;
	int i;
	for (i=*id; i<_ih_maxHashTableSize;i+=THREAD_COUNT)
	{
		if (_ih_hashTable[i].list == NULL) continue;
		cnt = 0;
		while (_ih_hashTable[i].list->info != _ih_refGenLen+1)
		{
			_ih_hashTable[i].list--;
			cnt++;
		}
		_ih_hashTable[i].list[0].info=cnt;
		if (cnt)
			introSortGI(_ih_hashTable[i].list, 1 , _ih_hashTable[i].list[0].info);
	}
	return NULL;
}
/**********************************************/
void *countQGrams(int *idp)
{
	int id = *idp;

	CompressedSeq *cnext, cdata;
	int i, t, val;

	int rgBlockSize = _ih_crefGenLen / THREAD_COUNT;
	int rgBlockStart = (rgBlockSize * id * 21);
	int rgBlockLen = rgBlockSize * 21;
	int rgBlockIt = rgBlockLen + SEQ_LENGTH - 1;
	if (id == THREAD_COUNT - 1)
	{
		rgBlockLen = _ih_refGenLen - id*rgBlockSize*21;
		rgBlockIt = rgBlockLen;
	}

	cnext = _ih_crefGen+(id*rgBlockSize);
	cdata = *(cnext++);
	t = 0;
	char outgoingChar[SEQ_LENGTH];
	unsigned int *copy = (unsigned int *)(_ih_alphCnt+4*rgBlockStart);
	unsigned char *cur = (unsigned char *)copy;		// current loc
	*copy = 0;

	for (i = 0; i < SEQ_LENGTH; i++)
	{
		val = (cdata >> 60) & 7;
		outgoingChar[i] = val;

		if (++t == 21)
		{
			t = 0;
			cdata = *(cnext++);
		}
		else
		{
			cdata <<= 3;
		}
		if (val != 4)
			(*(cur+val)) ++;
	}

	int o = 0;

	while (i++ < rgBlockIt) // BORDER LINE CHECK
	{
		cur = (unsigned char *)++copy;
		val = (cdata >> 60) & 7;
		if (++t == 21)
		{
			t = 0;
			cdata = *(cnext++);
		}
		else
		{
			cdata <<= 3;
		}

		*copy = *(copy-1);	// copies all 4 bytes at once
		if (val != 4)
			(*(cur + val)) ++;
		if (outgoingChar[o]!= 4)
			(*(cur + outgoingChar[o])) --;
		outgoingChar[o] = val;
		o = (++o == SEQ_LENGTH) ?0 :o;
	}
	return NULL;
}
/**********************************************/
int  loadHashTable(double *loadTime)
{
	// 1 byte (extraInfo): Reserved; in case the contig has extra information
	// 2 bytes (len): Length of the reference genome name
	// n bytes (refGenName): Reference genome name
	// 4 bytes (refGenOfsset): Offset of the contig from the beginning of the chromosome
	// 4 bytes (refGenLength): Length of reference genome
	// n bytes (crefGen): compressed reference genome
	// 4 bytes (size): number of hashValues in hashTable with more than 0 locations
	// n bytes (bufferSize and buffer): array of bufferSize/buffer which includes encoded values of hashValue, count of locations 

	int tmp;
	double startTime = getTime();
	unsigned char extraInfo = 0;
	short len;
	unsigned int hashTableSize;
	unsigned int tmpSize;
	int i = 0, j;

	if ( fread(&extraInfo, sizeof(extraInfo), 1, _ih_fp) != sizeof(extraInfo) )
	{
		return 0;
	}

	memset(_ih_hashTable, 0, _ih_maxHashTableSize * sizeof(IHashTable));

	// Reading Chr Name
	tmp = fread(&len, sizeof(len), 1, _ih_fp);
	tmp = fread(_ih_refGenName, sizeof(char), len, _ih_fp);
	_ih_refGenName [len] ='\0';

	tmp = fread(&_ih_refGenOff, sizeof (_ih_refGenOff), 1, _ih_fp);

	// Reading Size and Content of Ref Genome
	tmp = fread(&_ih_refGenLen, sizeof(_ih_refGenLen), 1, _ih_fp);

	_ih_crefGenLen = calculateCompressedLen(_ih_refGenLen);
	if (pairedEndMode)
	{
		_ih_crefGen = _ih_crefGenOrigin + _ih_refGenOff/21;
	}
	tmp = fread(_ih_crefGen, sizeof(CompressedSeq), _ih_crefGenLen, _ih_fp);


	//Reading Hashtable Size and Content
	GeneralIndex *mem =_ih_hashTableMem;

	tmp = fread(&hashTableSize, sizeof(hashTableSize), 1, _ih_fp);

	int index = 0, bytesToRead;
	unsigned int diff;
	unsigned long long hv=0;
	i = 0;
	while (i < hashTableSize)
	{
		fread(&bytesToRead, sizeof(int), 1, _ih_fp);
		fread(_ih_IOBuffer, sizeof(unsigned char), bytesToRead, _ih_fp);
		index = 0;
		while (index < bytesToRead)
		{
			index += decodeVariableByte(_ih_IOBuffer + index, &diff);
			index += decodeVariableByte(_ih_IOBuffer + index, &tmpSize);
			hv += diff;
			_ih_hashTable[hv].list = mem;
			mem->info = _ih_refGenLen+1;
			mem += (tmpSize + 1);
			i++;
		}
	}

	// // creating hash table
	for (i = 0; i < THREAD_COUNT; i++)
		pthread_create(_ih_threads + i, NULL, (void*)calculateHashTableOnFly, THREAD_ID + i);
	for (i = 0; i < THREAD_COUNT; i++)
		pthread_join(_ih_threads[i], NULL);

	// sorting based on checksum
	for (i = 0; i < THREAD_COUNT; i++)
		pthread_create(_ih_threads + i, NULL, (void*)sortHashTable, THREAD_ID + i);
	for (i = 0; i < THREAD_COUNT; i++)
		pthread_join(_ih_threads[i], NULL);

	// calculate alphabet count for each location in genome
	if (!SNPMode)
	{
		for (i = 0; i < THREAD_COUNT; i++)
			pthread_create(_ih_threads + i, NULL, (void*)countQGrams, THREAD_ID + i);
		for (i = 0; i < THREAD_COUNT; i++)
			pthread_join(_ih_threads[i], NULL);
	}

	*loadTime = getTime()-startTime;
	return extraInfo;
}
/**********************************************/
GeneralIndex *getCandidates(int hv)
{
	if ( hv != -1 && _ih_hashTable[hv].list != NULL && _ih_hashTable[hv].list[0].info != 0)
		return _ih_hashTable[hv].list;
	else 
		return NULL;
}
/**********************************************/
char *getRefGenomeName()
{
	return _ih_refGenName;
}
/**********************************************/
int getRefGenomeOffset()
{
	return _ih_refGenOff;
}
/**********************************************/
HashTable *getHashTable()
{
	return NULL;
}
/**********************************************/
CompressedSeq *getCmpRefGenome()
{
	return _ih_crefGen;
}
/**********************************************/
int getRefGenLength()
{
	return _ih_refGenLen;
}
/**********************************************/
int getCmpRefGenLength()
{
	return _ih_crefGenLen;
}
/**********************************************/
unsigned char *getAlphabetCount()
{
	return _ih_alphCnt;
}
/**********************************************/
CompressedSeq *getCmpRefGenOrigin()
{
	return _ih_crefGenOrigin;
}
/**********************************************/
int getChrCnt()
{
	return _ih_chrCnt;
}
/**********************************************/
char **getChrNames()
{
	return _ih_chrNames;
}
/**********************************************/
int getMaxChrLength()
{
	return _ih_maxChrLength;
}
