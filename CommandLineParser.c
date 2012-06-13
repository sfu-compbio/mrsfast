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
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include "Common.h"
#include "CommandLineParser.h"

int						indexingMode;
int						searchingMode;
int						bisulfiteMode;
int						pairedEndMode;
int						seqCompressed;
int						outCompressed;
int						cropSize = 0;
int						progressRep = 0;
char					*seqFile1;
char					*seqFile2;
char					*mappingOutput = "output";
char					*unmappedOutput = "unmapped";
short					minPairEndedDistance=-1;
short					maxPairEndedDistance=-1;
char					fileName[1000][2][FILE_NAME_LENGTH];
int						fileCnt;
void printHelp();

int parseCommandLine (int argc, char *argv[])
{
	
	int o;
	int index;
	char *fastaFile = NULL;
	char *fastaOutputFile = NULL;
	char *indexFile = NULL;
	char *batchFile = NULL ;

	static struct option longOptions[] = 
	{
		{"index",			no_argument,		&indexingMode,		1},
		{"search",			no_argument, 		&searchingMode,		1},
		{"bs",				no_argument,		&bisulfiteMode,		1},
		{"pe",				no_argument,		&pairedEndMode,		1},
		{"seqcomp",			no_argument,		&seqCompressed,		1},
		{"outcomp",			no_argument,		&outCompressed,		1},
		{"progress",		no_argument,		&progressRep,		1},
		{"help",			no_argument,		0,					'h'},
		{"version",			no_argument,		0,					'V'},
		{"seq",				required_argument,	0,					'x'},
		{"seq1",			required_argument,	0,					'x'},
		{"seq2",			required_argument,	0,					'y'},
		{"ws",				required_argument,  0,					'w'},
		{"min",				required_argument,  0,					'l'},
		{"max",				required_argument,  0,					'm'},
		{"crop",			required_argument,  0,					'c'}
	};

		

	while ( (o = getopt_long ( argc, argv, "f:b:i:u:o:s:e:n:hV", longOptions, &index))!= -1 )
	{
		switch (o)
		{
			case 'c': 
					cropSize = atoi(optarg);
					break;
			case 'w':
					WINDOW_SIZE = atoi(optarg);
					break;
			case 'x':
					seqFile1 = optarg;
					break;
			case 'y':
					seqFile2 = optarg;
					break;
			case 'f':
					fastaFile = optarg;
					break;
			case 'i':
					indexFile = optarg;
					break;
			case 'b':
					batchFile = optarg;
					break;
			case 'u':
					unmappedOutput = optarg;
					break;
			case 's':
					fastaOutputFile = optarg;
					break;
			case 'o':
					mappingOutput = optarg;
					break;
			case 'n':
					maxHits = atoi(optarg);
					break;
			case 'e':
					errThreshold = atoi(optarg);
					break;
			case 'l':
					minPairEndedDistance = atoi(optarg);
					break;
			case 'm':
					maxPairEndedDistance = atoi(optarg);
					break;					
			case 'h':
					printHelp();
					return 0;
					break;
			case 'V':
					fprintf(stdout, "%s.%s\n", versionNumber, versionNumberF);
					return 0;
					break;
		}

	}

	if (indexingMode + searchingMode != 1)
	{
		fprintf(stdout, "ERROR: Indexing / Searching mode should be selected\n");
		return 0;
	}

	if (WINDOW_SIZE > 14 || WINDOW_SIZE < 11)
	{
		fprintf(stdout, "ERROR: Window size should be in [12..14]\n");
		return 0;
	}


	if ( indexingMode )
	{
		if (batchFile == NULL && fastaFile == NULL)
		{
			fprintf(stdout, "ERROR: Reference(s) should be indicated for indexing\n");
			return 0;
		}

		if (batchFile != NULL && fastaFile != NULL)
		{
			fprintf(stdout, "ERROR: -b cannot be used with -f. \n");
			return 0;
		}

		if (batchFile != NULL && fastaOutputFile != NULL)
		{
			fprintf(stdout, "ERROR: -s cannot be used with -b. \n");
			return 0;
		}
	}


	if ( searchingMode )
	{
		if (batchFile == NULL && indexFile == NULL)
		{
			fprintf(stdout, "ERROR: Index File(s) should be indiciated for searching\n");
			return 0;
		}

		if (batchFile != NULL && indexFile != NULL)
		{
			fprintf(stdout, "ERROR: -b cannot be used with -i.\n");
			return 0;
		}

		if (seqFile1 == NULL && seqFile2 == NULL)
		{
			fprintf(stdout, "ERROR: Please indicate a sequence file for searching.\n");
			return 0;
		}


		if (!pairedEndMode && !bisulfiteMode && seqFile2 != NULL)
		{
			fprintf(stdout, "ERROR: Second File can be indicated in pairedend/bisulfite mode\n");
			return 0;
		}

		if (pairedEndMode && (minPairEndedDistance <0 || maxPairEndedDistance < 0 || minPairEndedDistance > maxPairEndedDistance))
		{
			fprintf(stdout, "ERROR: Please enter a valid range for pairedend sequences.\n");
			return 0;
		}

		if (pairedEndMode && seqFile1 == NULL)
		{
			fprintf(stdout, "ERROR: Please indicate the first file for pairedend search.\n");
			return 0;
		}
	}

	int i = 0;


	if (batchFile != NULL)
	{
		FILE *fp = fileOpen(batchFile, "r");

		if (fp == NULL)
			return 0;

		fileCnt  = 0;

		while ( fgets(fileName[fileCnt][0], FILE_NAME_LENGTH, fp))
		{
			for (i = strlen(fileName[fileCnt][0])-1; i>=0; i--)
				if ( !isspace(fileName[fileCnt][0][i]))
					break;
			fileName[fileCnt][0][i+1] = '\0';

			if (strcmp(fileName[fileCnt][0], "") != 0)
			{
				fileCnt++;
			}
		}
	}
	else
	{
		if (indexingMode)
		{
			sprintf(fileName[fileCnt][0], "%s", fastaFile);
		}
		else
		{
			sprintf(fileName[fileCnt][0], "%s", indexFile);
		}
		fileCnt++;
	}


	if (indexingMode)
	{
		for (i = 0;  i <fileCnt; i++)
		{
			if (bisulfiteMode)
				sprintf(fileName[i][1], "%s.bsindex", fileName[i][0]); 
			else
				sprintf(fileName[i][1], "%s.index", fileName[i][0]); 
		}
	}

	if (fastaOutputFile)
	{
		sprintf(fileName[0][1],"%s",fastaOutputFile);
	}

	return 1;
}


void printHelp()
{
	char *errorType;
	if (mrFAST)
	{
		fprintf(stdout,"mrFAST : Micro-Read Fast Alignment Search Tool.\n\n");
		fprintf(stdout,"Usage: mrFAST [options]\n\n");
		errorType="edit distance";
	}
	else
	{
		fprintf(stdout,"mrsFAST : Micro-Read Substitutions (only) Fast Alignment Search Tool.\n\n");
		fprintf(stdout,"Usage: mrsFAST [options]\n\n");
		errorType="hamming distance";
	}
	
	fprintf(stdout,"Indexing Options:\n");
	fprintf(stdout," --index\t\tMake index from a fasta file. \n");
	fprintf(stdout," -f [fastafile]\t\tInput fasta file for indexing. The output will be saved into '[fastafile].index' unless it is specified with -f. \n");
	fprintf(stdout," -b [file]\t\tIndex a set of fasta files listed in [file].\n");
	fprintf(stdout," -ws [int]\t\tSet window size for indexing (default:12-max:14).\n");
	fprintf(stdout," -s [file]\t\tSave index in the specified file.\n");
	fprintf(stdout," -bs \t\t\tBisulfite mode.");
	fprintf(stdout,"\n\n");
	fprintf(stdout,"Searching Options:\n");
	fprintf(stdout," --search\t\tSearch an index file\n");
	fprintf(stdout," --pe \t\t\tSearch will be done in paired-end mode.\n");
	fprintf(stdout," --bs \t\t\tSearch will be done in bisulfite mode.\n");
	fprintf(stdout," -i [indexfile]\t\tIndex file for search.\n");
	fprintf(stdout," -b [file]\t\tSearch a set of index files listed in [file].\n");
	fprintf(stdout," --seq [file]\t\tInput sequences in fasta/fastq format [file]. If paired-end sequences are interleaved, use this option.\n");
	fprintf(stdout," --seq1 [file]\t\tInput sequences in fasta/fastq format [file] (First file). Use this option to indicate the first file of pair-end sequences. You can use this option alone in bisulfite mode. \n");
	fprintf(stdout," --seq2 [file]\t\tInput sequences in fasta/fastq format [file] (Second file). Use this option to indicate the second file of pair-end sequences. You can use this option alone in bisulfite mode. \n");
	fprintf(stdout," -o [file]\t\tOutput of the mapped sequences. The default is output\n");
	fprintf(stdout," --seqcomp \t\tIndicates that the input sequences are compressed(gz).\n");
	fprintf(stdout," --outcomp \t\tIndicates that output file should be compressed(gz).\n");
	fprintf(stdout," -u [file]\t\tSave unmapped sequences to the [file] in fasta format.\n");
	fprintf(stdout," -n [int]\t\tMaximum number of locations a sequence can map to (default 1). To report all use 0.\n");
	fprintf(stdout," -e [int]\t\t%s (default 2).\n", errorType);
	fprintf(stdout," --min [int]\t\tMin distance allowed between two pairend sequences.(Should be used with --pe)\n");
	fprintf(stdout," --max [int]\t\tMax distance allowed between two pairend sequences.(Should be used with --pe)\n");
	fprintf(stdout," --crop [int]\t\tCrops the input sequences after the specified number of base pairs.\n");
}
