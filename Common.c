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

/*
 * Author         : Faraz Hach
 * Email          : fhach AT cs DOT sfu
 * Last Update    : 2009-02-01
 */


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <zlib.h>
#include <string.h>
#include <math.h>
#include "Common.h"


unsigned short 			SEQ_LENGTH = 0;
unsigned short 			QUAL_LENGTH = 0;
unsigned short			CMP_SEQ_LENGTH = 0;
long long				memUsage = 0;
char					*alphabet = "ACGTN";
/**********************************************/
FILE *fileOpen(char *fileName, char *mode)
{
	FILE *fp;
	fp = fopen (fileName, mode);
	if (fp == NULL)
	{
		fprintf(stdout, "Error: Cannot Open the file %s\n", fileName);
		fflush(stdout);
		exit(0);
	}
	return fp;
}
/**********************************************/
gzFile fileOpenGZ(char *fileName, char *mode)
{
	gzFile gzfp;
	gzfp = gzopen (fileName, mode);
	if (gzfp == Z_NULL)
	{
		fprintf(stdout, "Error: Cannot Open the file %s\n", fileName);
		fflush(stdout);
		exit(0);
	}
	return gzfp;
}
/**********************************************/
double getTime(void)
{
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec+t.tv_usec/1000000.0;
}

/**********************************************/
inline char reverseCompleteChar(char c)
{
	/*char rc[26];
	rc['A']='T';
	rc['T']='A';
	rc['C']='G';
	rc['G']='C';
	return rc[c];*/
	char ret;
	switch (c)
	{
		case 'A': 
					ret = 'T';
					break;
		case 'T':
					ret = 'A';
					break;
		case 'C':	
					ret = 'G';
					break;
		case 'G':
					ret = 'C';
					break;
		default:
					ret = 'N';
					break;
	}
	return ret;
}
/**********************************************/
inline void reverseComplete (char *seq, char *rcSeq , int length)
{
	char rc[100];
	memset(rc, 'N', 100);
	rc['A']='T';
	rc['T']='A';
	rc['C']='G';
	rc['G']='C';
	rc['N']='N';
	int i;
	seq+=length-1;
	for (i=0; i<length; i++)
	{
		rcSeq[i]=rc[*(seq--)];
		//rcSeq[i]=reverseCompleteChar (seq[length-1-i]) ;
	}
}
/**********************************************/
void * getMem(size_t size)
{
	memUsage+=size;
	return malloc(size);
}
/**********************************************/
void freeMem(void *ptr, size_t size)
{
	memUsage-=size;
	free(ptr);
}
/**********************************************/
double getMemUsage()
{
	return memUsage/1048576.0;
}
/**********************************************/
inline void reverse (char *seq, char *rcSeq , int length)
{
	int i;
	seq += length-1;
	for (i=0; i<length; i++)
	{
		//rcSeq[i]=seq[length-1-i];
		rcSeq[i]=*(seq--);
	}
}
/**********************************************/
void stripPath(char *full, char **path, char **fileName)
{
	int i;
	int pos = -1;

	for (i=strlen(full)-1; i>=0; i--)
	{
		if (full[i]=='/')
		{
			pos = i;
			break;
		}

	}

	if (pos != -1)
	{
		sprintf(*fileName, "%s%c", (full+pos+1), '\0');
		full[pos+1]='\0';
		sprintf(*path,"%s%c", full, '\0');
	}
	else
	{
		sprintf(*fileName, "%s%c", full, '\0');
		sprintf(*path,"%c", '\0');
	}
}
/**********************************************/
inline int calculateCompressedLen(int normalLen)
{
	return (normalLen / 21) + ((normalLen%21)?1:0);
	//return (normalLen >> 4) + ((normalLen & 15)?1:0);
}
/**********************************************/
void compressSequence(char *seq, int seqLen, CompressedSeq *cseq)
{
	CompressedSeq val = 0;
	int i = 0, pos = 0;
	
	*cseq = 0;
	while (pos < seqLen)
	{
		if (i == 0)
			*cseq = 0;
		*cseq <<= 3;
		switch (seq[pos++])
		{
			case 'A':
				break;
			case 'C':
				*cseq |= 1;				
				break;
			case 'G':
				*cseq |= 2;
				break;
			case 'T':
				*cseq |= 3;
				break;
			case 'N':
				*cseq |= 4;
				break;
			default:
				*cseq |= 4;
				break;
		}

		if (++i == 21)
		{
			i = 0;
			cseq++;
		}
	}
	if (i > 0)
	{
		*cseq <<= (3*(21-i));
	}
}

/**********************************************/
char 		**decompressLookup		= NULL;
void initDecompressLookup()
{
	int i, j, t;

	decompressLookup = getMem((1 << 15) * sizeof(char *));
	for (i = 0; i < (1 << 15); i++)
	{
		decompressLookup[i] = getMem(5*sizeof(char));
		t = i;
		for (j = 4; j >= 0; j--)
		{
			decompressLookup[i][j] = alphabet[t&7];
			t >>= 3;
		}
	}
}
/**********************************************/
void decompressSequence(CompressedSeq *cseq, char *seq)
{
	int i;
	char *seqq = seq;
	
	for (i = 0; i < CMP_SEQ_LENGTH; i++)
	{
		memcpy(seq, decompressLookup[(*cseq >> 48) & 0x7fff], 5);
		memcpy(seq + 5, decompressLookup[(*cseq >> 33) & 0x7fff], 5);
		memcpy(seq + 10, decompressLookup[(*cseq >> 18) & 0x7fff], 5);
		memcpy(seq + 15, decompressLookup[(*cseq >> 3) & 0x7fff], 5);
		seq += 20;
		*(seq++) = alphabet[*cseq & 7];
		cseq++;
	}
	seqq[SEQ_LENGTH] = '\0';
}
/**********************************************/
finalizeDecompressLookup()
{
	int i;
	for (i = 0; i < (1 << 15); i++)
		freeMem(decompressLookup[i], 5);
	freeMem(decompressLookup, (1 << 15)*sizeof(char *));
	decompressLookup = NULL;
}
