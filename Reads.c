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
 * Last Update    : 2009-12-08
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include "Common.h"
#include "Reads.h"

FILE *_r_fp1;
FILE *_r_fp2;
gzFile _r_gzfp1;
gzFile _r_gzfp2;
Read *_r_seq;
int _r_seqCnt;
int *_r_samplingLocs;

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
int readAllReads(char *fileName1,
						char *fileName2,
						int compressed,
						unsigned char *fastq,
						unsigned char pairedEnd,
						Read **seqList,
						unsigned int *seqListSize)
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
	char ch;
	int err1, err2;
	int nCnt;
	int discarded = 0;
	int seqCnt = 0;
	int maxCnt = 0;
	int i;
	Read *list = NULL;


	if (!compressed)
	{
		_r_fp1 = fileOpen( fileName1, "r");

		if (_r_fp1 == NULL)
		{
			return 0;
		}

		ch = fgetc(_r_fp1);

		if ( pairedEnd && fileName2 != NULL )
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

		if ( pairedEnd && fileName2 != NULL )
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
		*fastq = 0;
	else
		*fastq = 1;
	
	// Counting the number of lines in the file
	while (readFirstSeq(dummy)) maxCnt++;

	if (!compressed)
	{
		rewind(_r_fp1);
	}
	else
	{
		gzrewind(_r_gzfp1);
	}

	// Calculating the Maximum # of sequences
	if (*fastq)
	{
		maxCnt /= 4;
	}
	else
	{
		maxCnt /= 2;
	}



	if (pairedEnd && fileName2 != NULL )
		maxCnt *= 2;

	list = getMem(sizeof(Read)*maxCnt);


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

		if ( *fastq )
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
			if ( *fastq )
				qual1[cropSize] = '\0';
		}


		nCnt = 0;
		for (i=0; i<strlen(seq1); i++)
		{
			seq1[i] = toupper (seq1[i]);
			if (seq1[i] == 'N')
			{
				nCnt++;
			}
			else if (isspace(seq1[i]))
			{

				seq1[i] = '\0';
				break;
			}
		}

		if (nCnt > errThreshold)
		{
			err1 = 1;
		}

		// Reading the second seq of pair-ends
		if (pairedEnd)
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

			if ( *fastq )
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
				if ( *fastq )
					qual2[cropSize] = '\0';
			}


			nCnt = 0;
			for (i=0; i<strlen(seq2); i++)
			{
				seq2[i] = toupper (seq2[i]);
				if (seq2[i] == 'N')
				{
					nCnt++;

				}
				else if (isspace(seq2[i]))
				{
					seq2[i] = '\0';
				}
			}
			if (nCnt > errThreshold)
			{
				err2 = 1;
			}
		}

		if (!pairedEnd && !err1)
		{
			int _mtmp = strlen(seq1);
			int cmpLen = calculateCompressedLen(_mtmp);
			list[seqCnt].hits	= getMem (2 + (_mtmp * 3) + 3 + (cmpLen << 4) + strlen(name1) + 1);
			list[seqCnt].seq	= (char *) (list[seqCnt].hits + 1);
			list[seqCnt].rseq	= (char *)list[seqCnt].seq + _mtmp + 1;
			list[seqCnt].qual	= (char *)list[seqCnt].rseq + _mtmp + 1;
			list[seqCnt].cseq	= (CompressedSeq *)(list[seqCnt].qual + _mtmp + 1);
			list[seqCnt].crseq	= (CompressedSeq *)(list[seqCnt].cseq + cmpLen);
			list[seqCnt].name	= (char *)(list[seqCnt].crseq + cmpLen);
			//
			/*fprintf(stderr, "%d %d %d %d %d %d\n", 
				(list[seqCnt].hits), 
				(list[seqCnt].seq), 
				list[seqCnt].qual, 
				list[seqCnt].cseq, 
				list[seqCnt].crseq, 
				list[seqCnt].name);
			exit (0);*/

			reverseComplete(seq1, rseq, _mtmp);
			rseq[_mtmp] = '\0';
			compressSequence(seq1, _mtmp, list[seqCnt].cseq);
			compressSequence(rseq, _mtmp, list[seqCnt].crseq);

			int i;

			list[seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				list[seqCnt].seq[i] = seq1[i];
				list[seqCnt].rseq[i] = rseq[i];
				list[seqCnt].qual[i] = qual1[i];
			}
			
			sprintf(list[seqCnt].name,"%s%c", ((char*)name1)+1,'\0');

			seqCnt++;
		}
		else if (pairedEnd && !err1 && !err2)
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
			list[seqCnt].hits = getMem (2 + (_mtmp * 3) + 3 + (cmpLen << 4) + tmplen + 1);
			list[seqCnt].seq	= (char *) (list[seqCnt].hits + 1);
			list[seqCnt].rseq	= (char *)list[seqCnt].seq + _mtmp + 1;
			list[seqCnt].qual	= (char *)list[seqCnt].rseq + _mtmp + 1;
			list[seqCnt].cseq	= (CompressedSeq *)(list[seqCnt].qual + _mtmp + 1);
			list[seqCnt].crseq	= (CompressedSeq *)(list[seqCnt].cseq + cmpLen);
			list[seqCnt].name	= (char *)(list[seqCnt].crseq + cmpLen);

			reverseComplete(seq1, rseq, _mtmp);
			rseq[_mtmp] = '\0';
			compressSequence(seq1, _mtmp, list[seqCnt].cseq);
			compressSequence(rseq, _mtmp, list[seqCnt].crseq);

			int i;

			list[seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				list[seqCnt].seq[i] = seq1[i];
				list[seqCnt].rseq[i] = rseq[i];
				list[seqCnt].qual[i] = qual1[i];
			}
			
			name1[tmplen]='\0';
			sprintf(list[seqCnt].name,"%s%c", ((char*)name1)+1,'\0');


			seqCnt++;

			//second seq
			list[seqCnt].hits = getMem (2 + (_mtmp * 3) + 3 + (cmpLen << 4) + tmplen + 1);
			list[seqCnt].seq	= (char *) (list[seqCnt].hits + 1);
			list[seqCnt].rseq	= (char *)list[seqCnt].seq + _mtmp + 1;
			list[seqCnt].qual	= (char *)list[seqCnt].rseq + _mtmp + 1;
			list[seqCnt].cseq	= (CompressedSeq *)(list[seqCnt].qual + _mtmp + 1);
			list[seqCnt].crseq	= (CompressedSeq *)(list[seqCnt].cseq + cmpLen);
			list[seqCnt].name	= (char *)(list[seqCnt].crseq + cmpLen);

			reverseComplete(seq2, rseq, _mtmp);
			rseq[_mtmp] = '\0';
			compressSequence(seq2, _mtmp, list[seqCnt].cseq);
			compressSequence(rseq, _mtmp, list[seqCnt].crseq);

			list[seqCnt].hits[0] = 0;

			for (i=0; i<=_mtmp; i++)
			{
				list[seqCnt].seq[i] = seq2[i];
				list[seqCnt].rseq[i] = rseq[i];
				list[seqCnt].qual[i] = qual2[i];
			}

			name2[tmplen]='\0';
			sprintf(list[seqCnt].name,"%s%c", ((char*)name2)+1,'\0');

			seqCnt++;

		}
		else
		{
			discarded++;
		}

	}

	if (seqCnt > 0)
	{
		QUAL_LENGTH = SEQ_LENGTH = strlen(list[0].seq);
		CMP_SEQ_LENGTH = calculateCompressedLen(SEQ_LENGTH);
		if (! *fastq)
		{
			QUAL_LENGTH = 1;
		}
	}
	else
	{
		fprintf(stdout, "ERR: No reads can be found for mapping\n");
		return 0;
	}


	if (pairedEnd)
	{
//		seqCnt /= 2;
	}


	// Closing Files
	if (!compressed)
	{
		fclose(_r_fp1);
		if ( pairedEnd && fileName2 != NULL )
		{
			fclose(_r_fp2);
		}
	}
	else
	{
		gzclose(_r_gzfp1);
		if ( pairedEnd && fileName2 != NULL)
		{
			gzclose(_r_fp2);
		}
	}

	*seqList = list;
	*seqListSize = seqCnt;

	_r_seq = list;
	_r_seqCnt = seqCnt;

	fprintf(stdout, "%d sequences are read in %0.2f. (%d discarded) [Mem:%0.2f M]\n", seqCnt, (getTime()-startTime), discarded, getMemUsage());
	return 1;
}
/**********************************************/
void loadSamplingLocations(int **samplingLocs, int *samplingLocsSize)
{
	int i;
	int samLocsSize = errThreshold + 1;
	int *samLocs = getMem(sizeof(int)*(samLocsSize+1));
	for (i=0; i<samLocsSize; i++)
	{
		samLocs[i] = (SEQ_LENGTH / samLocsSize) *i;
		if ( samLocs[i] + WINDOW_SIZE > SEQ_LENGTH)
			samLocs[i] = SEQ_LENGTH - WINDOW_SIZE;
	}
	samLocs[samLocsSize]=SEQ_LENGTH;
	
	// Outputing the sampling locations
/*	int j;
 	for (i=0; i<SEQ_LENGTH; i++)
	{
		fprintf(stdout, "-");
	}
	fprintf(stdout, "\n");

	for ( i=0; i<samLocsSize; i++ )
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
	*samplingLocs = samLocs;
	*samplingLocsSize = samLocsSize;
	_r_samplingLocs = samLocs;
}

void finalizeReads(char *fileName)
{
	FILE *fp1=NULL;

	if (fileName != NULL)
	{
		fp1 = fileOpen(fileName, "w");
	}
	if (pairedEndMode)
		_r_seqCnt /=2;

	int i=0;
	for (i = 0; i < _r_seqCnt; i++)
	{
		if (pairedEndMode && _r_seq[2*i].hits[0] == 0 &&  strcmp(_r_seq[2*i].qual,"*")!=0)
		{
			fprintf(fp1,"@%s/1\n%s\n+\n%s\n@%s/2\n%s\n+\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].qual, _r_seq[i*2].name, _r_seq[i*2+1].seq, _r_seq[i*2+1].qual);
		}
		else if (pairedEndMode && _r_seq[2*i].hits[0] == 0)
		{
			fprintf(fp1, ">%s/1\n%s\n>%s/2\n%s\n", _r_seq[i*2].name, _r_seq[i*2].seq, _r_seq[i*2].name, _r_seq[i*2+1].seq);
		}
		else if (_r_seq[i].hits[0] == 0 && strcmp(_r_seq[i].qual, "*")!=0)
		{
			fprintf(fp1,"@%s\n%s\n+\n%s\n", _r_seq[i].name, _r_seq[i].seq, _r_seq[i].qual);
		}
		else if (_r_seq[i].hits[0] == 0)
		{
			fprintf(fp1,">%s\n%s\n", _r_seq[i].name, _r_seq[i].seq);
		}
	}

	fclose(fp1);
	if (pairedEndMode)
		_r_seqCnt *= 2;

	for (i = 0; i < _r_seqCnt; i++)
	{
		freeMem(_r_seq[i].hits,0);
	}


	freeMem(_r_seq,0);
	freeMem(_r_samplingLocs,0);
}

