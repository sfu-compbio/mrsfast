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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"
#include "RefGenome.h"
#include "HashTableBS.h"
/**********************************************/
FILE			*_ihbs_fp					= NULL;
IHashTableBS	*_ihbs_hashTable			= NULL;
int				_ihbs_maxHashTableSize		= 0;
char			*_ihbs_refGen				= NULL;
char			*_ihbs_refGenName			= NULL;
long long		_ihbs_memUsage				= 0;
int				_ihbs_refGenOff				= 0;
/**********************************************/

int hashValBS(char *seq, short type)
{
	/*
	 * 0 ==> A=0  C=1 G=2 T=1 ( CT considered as one character )
	 * 1 ==> A=0  C=1 G=0 T=2 ( AG considered as one character )
	 */
	int i=0;
	int val=0, numericVal=0;

	if (!type)
	{

		while(i<WINDOW_SIZE)
		{
			switch (seq[i])
			{
			    case 'A':
				    numericVal = 0;
				    break;
			    case 'C':
				case 'T':
				    numericVal = 1;
				    break;
			    case 'G' :
				    numericVal = 2;
			    	break;
				default:
					return -1;
					break;
			}
			val = (val << 2)|numericVal;
			i++;
		}
	}

	else
	{
		while(i<WINDOW_SIZE)
		{
			switch (seq[i])
			{
			    case 'A':
				case 'G':
				    numericVal = 0;
				    break;
			    case 'C':
					numericVal = 1;
				    break;
			    case 'T' :
				    numericVal = 2;
					break;
				default:
					return -1;
					break;
			}
			val = (val << 2)|numericVal;
			i++;
		}
	}
    return val;
}

/**********************************************/
void freeIHashTableContentBS(IHashTableBS *hashTable, unsigned int maxSize)
{
	int i=0, j=0;
	for (i=0; i<maxSize; i++)
	{
		for (j=0; j<2; j++)
		{
			if (hashTable[i].locs[j] != NULL)
			{

				freeMem(hashTable[i].locs[j], (hashTable[i].locs[j][0]+1)*(sizeof(unsigned int)));
				hashTable[i].locs[j] = NULL;
			}
		}
	}
}
/**********************************************/
void initSavingIHashTableBS(char *fileName)
{
	int tmp;

	_ihbs_fp = fileOpen(fileName, "w");
	unsigned char bIndex = 1;				// Bisulfite Index
	
	// First Two bytes are indicating the type of the index & window size
	tmp = fwrite(&bIndex, sizeof(bIndex), 1, _ihbs_fp);
	tmp = fwrite(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ihbs_fp);
	
}
/**********************************************/
void finalizeSavingIHashTableBS()
{
	fclose(_ihbs_fp);
}
/**********************************************/
void saveIHashTableBS(IHashTableBS *hashTable,  unsigned int * size, unsigned int maxSize, char *refGen, char *refGenName, int refGenOffset)
{

	int tmp;

	// Every Chunk starts with a byte indicating whether it has extra info;
	unsigned char extraInfo = 0;
	tmp = fwrite (&extraInfo, sizeof(extraInfo), 1, _ihbs_fp);

	short len = strlen(refGenName);
	tmp = fwrite(&len, sizeof(len), 1, _ihbs_fp);
	tmp = fwrite(refGenName, sizeof(char), len, _ihbs_fp);

	tmp = fwrite(&refGenOffset, sizeof(refGenOffset), 1, _ihbs_fp);
	
	unsigned int refGenLength = strlen(refGen);
	tmp = fwrite(&refGenLength, sizeof(refGenLength), 1, _ihbs_fp);
	tmp = fwrite(refGen, sizeof(char), refGenLength, _ihbs_fp);



	int i=0,j=0;
	unsigned char cnt=0;
	int k=0;
	for (k=0; k<2; k++)
	{
		tmp = fwrite(&(size[k]), sizeof(size[k]), 1, _ihbs_fp);
		for (i=0; i<maxSize; i++)
		{
			if (hashTable[i].locs[k] != NULL)
			{
				tmp = fwrite(&i, sizeof(i), 1, _ihbs_fp);

				if (hashTable[i].locs[k][0] < 250)
				{
					cnt = hashTable[i].locs[k][0];
					tmp = fwrite(&cnt, sizeof(cnt), 1, _ihbs_fp);
				}
				else
				{
					cnt =0;
					tmp = fwrite (&cnt, sizeof(cnt), 1, _ihbs_fp);
					tmp = fwrite (&(hashTable[i].locs[k][0]), sizeof(hashTable[i].locs[k][0]), 1, _ihbs_fp);
				}

				for (j=1; j<=hashTable[i].locs[k][0]; j++)
					tmp = fwrite(&(hashTable[i].locs[k][j]), sizeof(hashTable[i].locs[k][j]), 1, _ihbs_fp);
			}
		}
	}
}
/**********************************************/
unsigned int addIHashTableLocationBS(IHashTableBS *hashTable, int hv, int location, short type)
{
	unsigned int	sizeInc				= 0;

	if (hashTable[hv].locs[type] == NULL)
	{
		sizeInc = 1;
		hashTable[hv].locs[type] = getMem (sizeof(unsigned int)*2);
		hashTable[hv].locs[type][0]=1;
		hashTable[hv].locs[type][1]=location;
	}
	else
	{
		int size = hashTable[hv].locs[type][0];
		int i;
		unsigned int *tmp = getMem( (size + 2) * sizeof(unsigned int) );

		for (i = 0; i <= size; i++)
		{
			tmp[i] = hashTable[hv].locs[type][i];
		}
		size++;
		tmp[0] = size;
		tmp[size] = location;
		
		freeMem(hashTable[hv].locs[type], (hashTable[hv].locs[type][0]*(sizeof(unsigned int))));
		hashTable[hv].locs[type] = tmp;
	}
	return sizeInc;
}
/**********************************************/
void generateIHashTableBS(char *fileName, char *indexName)
{
	double          startTime           = getTime();
	unsigned int	hashTableSize[2]	= {0,0};
	unsigned int 	hashTableMaxSize	= pow(4, WINDOW_SIZE);
	IHashTableBS	*hashTable			= getMem(sizeof(IHashTableBS)*hashTableMaxSize);
	char 			*refGenName;
	char			*refGen;
	int				refGenOff			= 0;
	int i, hv, l, flag,k;

	
	for (k=0; k<2; k++)
	{
		for ( i = 0; i < hashTableMaxSize; i++)
		{
			hashTable[i].locs[k] = NULL;
		}
	}

	//Loading Fasta File
	if (!initLoadingRefGenome(fileName))
		return;		
	initSavingIHashTableBS(indexName);
	
	fprintf(stdout, "Generating Index from %s", fileName);
	fflush(stdout);

	char *prev = getMem (CONTIG_NAME_SIZE);
	prev[0]='\0';

	do
	{
		flag = 	 loadRefGenome (&refGen, &refGenName, &refGenOff);	

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
		
		l = strlen(refGen) - WINDOW_SIZE;

		for (i=0; i < l; i++)
		{
			for (k=0; k<2; k++)
			{
				hv = hashValBS(refGen + i, k);
				if (hv != -1)
				{
					hashTableSize[k] += addIHashTableLocationBS(hashTable, hv, i+1, k);
				}
			}
		}

		saveIHashTableBS(hashTable, hashTableSize, hashTableMaxSize, refGen, refGenName, refGenOff);
		freeIHashTableContentBS(hashTable, hashTableMaxSize);
		hashTableSize[0] = hashTableSize[1] = 0 ;
	} while (flag);

	freeMem(prev, CONTIG_NAME_SIZE);
	freeMem(hashTable, sizeof(IHashTableBS)*hashTableMaxSize);

	finalizeLoadingRefGenome();
	finalizeSavingIHashTableBS();

	fprintf(stdout, "\nDONE in %0.2fs!\n", (getTime()-startTime));
}
/**********************************************/
void finalizeLoadingIHashTableBS()
{
	freeIHashTableContentBS(_ihbs_hashTable, _ihbs_maxHashTableSize);
	freeMem(_ihbs_hashTable, sizeof(IHashTableBS)* _ihbs_maxHashTableSize);
	freeMem(_ihbs_refGen, strlen(_ihbs_refGen)+1) ;
	freeMem(_ihbs_refGenName, strlen(_ihbs_refGenName)+1);
	fclose(_ihbs_fp);
}
/**********************************************/
int  loadIHashTableBS(double *loadTime)
{
	double startTime = getTime();
	unsigned char extraInfo = 0;
	short len;
	unsigned int refGenLength;
	unsigned int hashTableSize[2];
	unsigned int tmpSize;
	int i=0,j=0, k=0;
	int tmp;

	if ( fread(&extraInfo, sizeof(extraInfo), 1, _ihbs_fp) != sizeof(extraInfo) )
		return 0;

	freeIHashTableContentBS(_ihbs_hashTable, _ihbs_maxHashTableSize);
	freeMem(_ihbs_refGen, strlen(_ihbs_refGen)+1) ;
	freeMem(_ihbs_refGenName, strlen(_ihbs_refGenName)+1);

	// Reading Chr Name
	tmp = fread(&len, sizeof(len), 1, _ihbs_fp);
	_ihbs_refGenName = getMem(sizeof(char)* (len+1));
	tmp = fread(_ihbs_refGenName, sizeof(char), len, _ihbs_fp);
	_ihbs_refGenName [len] ='\0';

	tmp = fread(&_ihbs_refGenOff, sizeof (_ihbs_refGenOff), 1, _ihbs_fp);

	// Reading Size and Content of Ref Genome
	tmp = fread(&refGenLength, sizeof(refGenLength), 1, _ihbs_fp);
	_ihbs_refGen = getMem(sizeof(char)*(refGenLength+1));
	tmp = fread(_ihbs_refGen, sizeof(char), refGenLength, _ihbs_fp);
	_ihbs_refGen[refGenLength]='\0';

	//Reading Hashtable Size and Content
	for (k=0; k<2; k++)
	{
		tmp = fread(&(hashTableSize[k]), sizeof(hashTableSize[k]), 1, _ihbs_fp);

		unsigned int hv;
		unsigned char cnt=0;
		for (i=0; i<hashTableSize[k]; i++)
		{
			tmp = fread(&hv, sizeof(hv), 1, _ihbs_fp);
			tmp = fread(&cnt, sizeof(cnt), 1, _ihbs_fp);

			if (cnt>0)
			{
				tmpSize = cnt;
			}
			else
			{
				tmp = fread(&tmpSize, sizeof(tmpSize), 1, _ihbs_fp);
			}

			_ihbs_hashTable[hv].locs[k] = getMem( sizeof(unsigned int)* (tmpSize+1) );
			_ihbs_hashTable[hv].locs[k][0] = tmpSize;

			for (j=1; j<=tmpSize; j++)
				tmp = fread(&(_ihbs_hashTable[hv].locs[k][j]), sizeof(_ihbs_hashTable[hv].locs[k][j]), 1, _ihbs_fp);
		}
	}

	*loadTime = getTime()-startTime;

	return 1;
}
/**********************************************/
unsigned int *getIHashTableCandidatesBS(int hv, short type)
{
	if ( hv != -1 )
		return _ihbs_hashTable[hv].locs[type];
	else 
		return NULL;
}
/**********************************************/
void configHashTableBS()
{
	if (WINDOW_SIZE <= 14)
	{
		generateHashTableBS = &generateIHashTableBS;
		loadHashTableBS = &loadIHashTableBS;
		finalizeLoadingHashTableBS = &finalizeLoadingIHashTableBS;
		getCandidatesBS = &getIHashTableCandidatesBS;
	}
}
/**********************************************/
int initLoadingHashTableBS(char *fileName)
{
	int i;	
	int tmp;
	unsigned char bsIndex;
	_ihbs_fp = fileOpen(fileName, "r");	

	if (_ihbs_fp == NULL)
		return 0;

	tmp = fread(&bsIndex, sizeof(bsIndex), 1, _ihbs_fp);
	if (!bsIndex)
	{
		fprintf(stdout, "Error: Wrong Type of Index indicated");
		return 0;
	}
	
	tmp = fread(&WINDOW_SIZE, sizeof(WINDOW_SIZE), 1, _ihbs_fp);

	configHashTableBS();

	if (_ihbs_maxHashTableSize != pow(4, WINDOW_SIZE))
	{

		if (_ihbs_hashTable != NULL)
		{
			freeIHashTableContentBS(_ihbs_hashTable, _ihbs_maxHashTableSize);
			freeMem(_ihbs_hashTable, sizeof(IHashTableBS)* _ihbs_maxHashTableSize);
			freeMem(_ihbs_refGen, strlen(_ihbs_refGen)+1) ;
			freeMem(_ihbs_refGenName, strlen(_ihbs_refGenName)+1);
		}

		_ihbs_maxHashTableSize = pow(4, WINDOW_SIZE);

		_ihbs_hashTable = getMem (sizeof(IHashTableBS) * _ihbs_maxHashTableSize);

		int k=0;

		for (k=0; k<2; k++)
		{
			for (i=0; i<_ihbs_maxHashTableSize; i++)
				_ihbs_hashTable[i].locs[k] = NULL;
		}
		_ihbs_refGen = getMem(1);
		_ihbs_refGen[0]='\0';
		_ihbs_refGenName = getMem(1);
		_ihbs_refGenName[0] = '\0';
	}

	return 1;
}
/**********************************************/
char *getRefGenomeBS()
{
	return _ihbs_refGen;

}
/**********************************************/
char *getRefGenomeNameBS()
{
	return _ihbs_refGenName;

}
/**********************************************/
int getRefGenomeOffsetBS()
{
	return _ihbs_refGenOff;

}
/**********************************************/
HashTableBS *getHashTableBS()
{
	return NULL;
}
