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
#include <getopt.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include "Common.h"
#include "CommandLineParser.h"

int						uniqueMode=1;
int						indexingMode;
int						searchingMode;
int						pairedEndMode;
int						pairedEndDiscordantMode;
int						pairedEndProfilingMode = 0;
int						bestMappingMode = 0;
int						SNPMode = 0;
int						seqCompressed;
int						outCompressed;
int						isCloud;
int						cropSize = 0;
int						tailCropSize = 0;
int						progressRep = 0;
int						nohitDisabled = 0;
int						noSamHeader = 0;
int						minPairEndedDistance=-1;
int						maxPairEndedDistance=-1;
int						minPairEndedDiscordantDistance=-1;
int						maxPairEndedDiscordantDistance=-1;
int						errThreshold = -1;
char					*seqFile1;
char					*seqFile2;
char					fileName[3][FILE_NAME_LENGTH];
char					*unmappedOutput;
char					*mappingOutputPath;
char					*concordantStatOutput;
short					maxHits = 0;
unsigned char			WINDOW_SIZE = 12;
unsigned int			CONTIG_SIZE;
unsigned int			CONTIG_MAX_SIZE;
unsigned int			THREAD_COUNT = 1;
unsigned short			DISCORDANT_CUT_OFF = 300;
double					MAX_MEMORY = 4;// GB
int						THREAD_ID[255];
int						SNP_QUAL_THRESHOLD = 53;



/*
#if (defined(__MACH__) && defined(__APPLE__))
#include <mach-o/getsect.h>
#else
extern char _binary_HELP_start;
extern char _binary_HELP_end;
#endif
*/

/*
void printHelp()
{
#if (defined(__MACH__) && defined(__APPLE__))
	size_t i, sz = getsectbyname("binary", "HELP")->size;
	const uint8_t *c =  (const uint8_t*) getsectbyname("binary", "HELP")->addr;
	for (i = 0; i < sz; i++) 
		putchar(c[i]); 
#else
	char *c;
	for (c = &_binary_HELP_start; c != &_binary_HELP_end; c++)
		putchar(*c);
#endif
	exit(EXIT_SUCCESS);
	}*/

int parseCommandLine (int argc, char *argv[])
{
	int index, len, o;
	char *fastaFile = NULL;
	char *fastaOutputFile = NULL;
	char *indexFile = NULL;
	char *SNPFile = NULL;

	mappingOutput = getMem(FILE_NAME_LENGTH);
	mappingOutputPath = getMem(FILE_NAME_LENGTH);
	unmappedOutput = getMem(FILE_NAME_LENGTH);
	concordantStatOutput = getMem(FILE_NAME_LENGTH);
	strcpy(mappingOutput, "output");
	strcpy(unmappedOutput, "output.nohit");
	strcpy(concordantStatOutput, "concordant.statistic");
	mappingOutputPath[0] = '\0';
	isCloud = 0;

	static struct option longOptions[] = 
	{
		{"pe",						no_argument,		&pairedEndMode,		1},
		{"discordant-vh",			no_argument,		&pairedEndDiscordantMode,	1},
		{"seqcomp",					no_argument,		&seqCompressed,		1},
		{"outcomp",					no_argument,		&outCompressed,		1},
		{"progress",				no_argument,		&progressRep,		1},
		{"best",					no_argument,		&bestMappingMode,	1},
		{"cloud",					no_argument,		&isCloud,	1},
		{"disable-nohits",			no_argument,		&nohitDisabled,		1},
		{"disable-sam-header",		no_argument,		&noSamHeader,		1},
		{"index",					required_argument,	0, 					'i'},
		{"search",					required_argument,	0,					's'},
		{"help",					no_argument,		0,					'h'},
		{"version",					no_argument,		0,					'v'},
		{"seq",						required_argument,	0,					'x'},
		{"seq1",					required_argument,	0,					'x'},
		{"seq2",					required_argument,	0,					'y'},
		{"ws",						required_argument,  0,					'w'},
		{"min",						required_argument,  0,					'l'},
		{"max",						required_argument,  0,					'm'},
		{"crop",					required_argument,  0,					'c'},
		{"tail-crop",				required_argument,  0,					'f'},
		{"threads",					required_argument,  0,					't'},
		{"mem",						required_argument,  0,					'z'},
		{"snp",						required_argument,  0,					'p'},
		{"max-discordant-cutoff",	required_argument,  0,                  'd'},
		{"snp-qual",				required_argument,  0,                  'q'},
		{0,0,0,0}
	};



	while ( (o = getopt_long ( argc, argv, "f:i:u:o:s:e:n:t:bhv", longOptions, &index))!= -1 )
	{
		switch (o)
		{
			case 'i':
				indexingMode = 1;
				fastaFile = optarg;
				break;
			case 's':
				searchingMode = 1;
				fastaFile = optarg;
				break;
			case 'c': 
				cropSize = atoi(optarg);
				break;
			case 'f': 
				tailCropSize = atoi(optarg);
				break;
			case 'w':
				if (searchingMode == 1)
				{
					fprintf(stderr, "Error: Window size can only be set in indexing mode.\n");
					return 0;
				}
				WINDOW_SIZE = atoi(optarg);
				break;
			case 'x':
				seqFile1 = optarg;
				break;
			case 'y':
				seqFile2 = optarg;
				break;
			case 'u':
				strcpy(unmappedOutput, optarg);
				break;
			case 'o':
				stripPath (optarg, &mappingOutputPath, &mappingOutput);
				sprintf(unmappedOutput, "%s%s.nohit", mappingOutputPath, mappingOutput );
				break;
			case 'n':
				maxHits = atoi(optarg);
				break;
			case 'e':
				errThreshold = atoi(optarg);
				break;
			case 'q':
				SNP_QUAL_THRESHOLD = atoi(optarg);
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
			case 'v':
				fprintf(stderr, "Version: %s\nBuild Date: %s\n", MRSFAST_VERSION, BUILD_DATE);
				return 0;
				break;
			case 't':
				THREAD_COUNT = atoi(optarg);
				if (THREAD_COUNT == 0 || THREAD_COUNT > sysconf( _SC_NPROCESSORS_ONLN ))
					THREAD_COUNT = sysconf( _SC_NPROCESSORS_ONLN );
				break;
			case 'z':
				MAX_MEMORY = atoi(optarg);
				break;
			case 'd':
				DISCORDANT_CUT_OFF = atoi(optarg);
				break;
			case 'p':
				SNPMode = 1;
				SNPFile = optarg;
				break;
		}

	}

#ifndef MRSFAST_SSE4
	if (searchingMode)
		fprintf(stderr, "==> This version is compiled without any SSE4 optimization <==\n");
#endif
	if (bestMappingMode)
	{
		nohitDisabled = 1;
	}

	if (indexingMode + searchingMode != 1)
	{
		fprintf(stderr, "ERROR: Indexing / Searching mode should be selected\n");
		return 0;
	}

	if (WINDOW_SIZE > 14 || WINDOW_SIZE < 8)
	{
		fprintf(stderr, "ERROR: Window size should be in [8..14]\n");
		return 0;
	}

	if (MAX_MEMORY < 2)
		fprintf(stderr, "ERROR: At least 2 GB of memory is required for running mrsFAST\n");


	if ( indexingMode )
	{
		CONTIG_SIZE		= 80000000;
		CONTIG_MAX_SIZE	= 120000000;

		if (fastaFile == NULL)
		{
			fprintf(stderr, "ERROR: Reference(s) should be indicated for indexing\n");
			return 0;
		}
	}

	if (maxHits)
	{
		if (maxHits < 0)
		{
			fprintf(stderr, "ERROR: Number of maximum hits must be greater than 0\n");
			return 0;
		}

		if (bestMappingMode)
		{
			fprintf(stderr, "ERROR: Maximum number of mappings could not be set in best mapping mode. Maximum mappings input ignored\n");
			maxHits = 0;
		}
	}

	if ( searchingMode )
	{
		CONTIG_SIZE		= 300000000;
		CONTIG_MAX_SIZE	= 300000000;

		if ( cropSize && tailCropSize)
		{
			fprintf(stderr, "ERROR: Sequences can be cropped from only one side\n");
			return 0;
		}

		if (pairedEndDiscordantMode)
		{
			pairedEndDiscordantMode = pairedEndMode = 1;
		}
		
		if (fastaFile == NULL)
		{
			fprintf(stderr, "ERROR: Index File(s) should be indiciated for searching\n");
			return 0;
		}

		if (seqFile1 == NULL && seqFile2 == NULL)
		{
			fprintf(stderr, "ERROR: Please indicate a sequence file for searching.\n");
			return 0;
		}
		
		if (!pairedEndMode && seqFile2 != NULL)
		{
			fprintf(stderr, "ERROR: Second File can be indicated in pairedend mode\n");
			return 0;
		}

		if (pairedEndMode)
		{
			if (minPairEndedDistance < 0 && maxPairEndedDistance < 0)
			{
				pairedEndProfilingMode = 1;
			}
			else if ( minPairEndedDistance <0 || maxPairEndedDistance < 0 || minPairEndedDistance > maxPairEndedDistance )
			{
				fprintf(stderr, "ERROR: Please enter a valid range for pairedend sequences.\n");
				return 0;
			}

			if (seqFile1 == NULL)
			{
				fprintf(stderr, "ERROR: Please indicate the first file for pairedend search.\n");
				return 0;
			}
		}
	}

	int i = 0;

	sprintf(fileName[0], "%s", fastaFile);
	sprintf(fileName[1], "%s.index", fileName[0]);
	if (SNPMode)
		sprintf(fileName[2], "%s", SNPFile);

	if (!indexingMode)
	{
		fprintf(stderr, "# Threads: %d\n", THREAD_COUNT);
		for (i = 0; i < 255; i++)
			THREAD_ID[i] = i;
	}

	char fname1[FILE_NAME_LENGTH];
	char fname2[FILE_NAME_LENGTH];
	char fname3[FILE_NAME_LENGTH];
	char fname4[FILE_NAME_LENGTH];
	char fname5[FILE_NAME_LENGTH];

	/* change defaut output filenames */
	if (!strcmp(mappingOutput, "output") && searchingMode)
	  {
	    sprintf(mappingOutput, "%s-output", seqFile1);
	    fprintf(stderr, "seqFile1 is %s\noutput file is %s\n", seqFile1, mappingOutput);
	    if (!outCompressed)
	      sprintf(unmappedOutput, "%s-output.nohit.fastq", seqFile1 );
	    else
	      sprintf(unmappedOutput, "%s-output.nohit.fastq.gz", seqFile1 );
	    char tmp_fname[FILE_NAME_LENGTH];
	    strcpy(tmp_fname, seqFile1);
	    stripPath (tmp_fname, &mappingOutputPath, &mappingOutput);
	  }
	
	// Why is this one here?
	if (pairedEndMode)
	{
		sprintf(fname1, "%s__%s__1", mappingOutputPath, mappingOutput);
		sprintf(fname2, "%s__%s__2", mappingOutputPath, mappingOutput);
		sprintf(fname3, "%s__%s__disc", mappingOutputPath, mappingOutput);
		sprintf(fname4, "%s__%s__oea1", mappingOutputPath, mappingOutput);
		sprintf(fname5, "%s__%s__oea2", mappingOutputPath, mappingOutput);
		unlink(fname1);
		unlink(fname2);
		unlink(fname3);
		unlink(fname4);
		unlink(fname5);
	}
	initCommon();
	return 1;
}
/**********************************************/
void finalizeCommandParser()
{
	freeMem(mappingOutput, FILE_NAME_LENGTH);
	freeMem(unmappedOutput, FILE_NAME_LENGTH);
	freeMem(mappingOutputPath, FILE_NAME_LENGTH);
	freeMem(concordantStatOutput, FILE_NAME_LENGTH);
}

void printHelp()
{
  fprintf(stderr, "mrsFAST-Ultra(1)             mrsfast-Ultra Manual             mrsFAST-Ultra(1)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "NNAAMMEE\n");
  fprintf(stderr, "       mrsfast-ultra\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "DDEESSCCRRIIPPTTIIOONN\n");
  fprintf(stderr, "       mrsFAST is a cache oblivious read mapping tool. mrsFAST capable of map-\n");
  fprintf(stderr, "       ping single and paired end reads to  the  reference  genome.  Bisulfite\n");
  fprintf(stderr, "       treated  sequences are not supported in this version.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "IINNSSTTAALLLLAATTIIOONN\n");
  fprintf(stderr, "       To  install  mrsFAST-ultra, please download the source zip package from\n");
  fprintf(stderr, "       http://sourceforge.net/projects/mrsfast/.  After  unzipping  the  down-\n");
  fprintf(stderr, "       loaded  file \"mrsfast-ultra-3.X.X.zip\", change the current directory to\n");
  fprintf(stderr, "       the source directory \"mrsfast-ultra-3.X.X\", and run \"make\" in the  ter-\n");
  fprintf(stderr, "       minal.  The  binary  file  \"mrsfast\" will be created, which is ready to\n");
  fprintf(stderr, "       use.\n");
  fprintf(stderr, "       $ unzip mrsfast-ultra-3.X.X.zip\n");
  fprintf(stderr, "       $ cd mrsfast-ultra-3.X.X\n");
  fprintf(stderr, "       $ make\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "SSYYNNOOPPSSIISS\n");
  fprintf(stderr, "       mrsfast --index [file] [OPTIONS]\n");
  fprintf(stderr, "       mrsfast --search [index] --seq [file] [OPTIONS]\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "OOPPTTIIOONNSS\n");
  fprintf(stderr, "   GGEENNEERRAALL OOPPTTIIOONNSS\n");
  fprintf(stderr, "       --hh     Prints this help file.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       --vv,, ----vveerrssiioonn\n");
  fprintf(stderr, "              Prints the version of software.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   IINNDDEEXXIINNGG OOPPTTIIOONNSS\n");
  fprintf(stderr, "       ----wwss _w_i_n_d_o_w___s_i_z_e\n");
  fprintf(stderr, "              Index the reference genome with sliding a window  of  size  _w_i_n_-\n");
  fprintf(stderr, "              _d_o_w___s_i_z_e (default: 12).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "   MMAAPPPPIINNGG OOPPTTIIOONNSS\n");
  fprintf(stderr, "       ----mmeemm _m\n");
  fprintf(stderr, "              Use maximum _m GB of memory (default: 4).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----tthhrreeaaddss _t\n");
  fprintf(stderr, "              Use  _t  number  of cores for mapping the sequences (default: 1).\n");
  fprintf(stderr, "              Use _0 to use all the available cores in the system.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----sseeqq _f_i_l_e\n");
  fprintf(stderr, "              Set the input sequence to _f_i_l_e_.  In paired-end mode, _f_i_l_e should\n");
  fprintf(stderr, "              be used if the read sequences are interleaved.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----sseeqq11 _f_i_l_e\n");
  fprintf(stderr, "              Set  the  input sequence (left mate) to _f_i_l_e_.  Paired-end option\n");
  fprintf(stderr, "              only.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----sseeqq22 _f_i_l_e\n");
  fprintf(stderr, "              Set the input sequence (right mate) to _f_i_l_e_.  Paired-end  option\n");
  fprintf(stderr, "              only.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----sseeqqccoommpp\n");
  fprintf(stderr, "              Input file is compressed through gzip.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       --oo _f_i_l_e\n");
  fprintf(stderr, "              Output the mapping record into _f_i_l_e (default: output.sam)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----ddiissaabbllee--ssaamm--hheeaaddeerr\n");
  fprintf(stderr, "              Do not generate SAM header.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       --uu _f_i_l_e\n");
  fprintf(stderr, "              Output unmapped reads in _f_i_l_e (default: output.nohit). This file\n");
  fprintf(stderr, "              will be generated in all mapping mode.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----ddiissaabbllee--nnoohhiittss\n");
  fprintf(stderr, "              Do not output unmapped reads.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----oouuttccoommpp\n");
  fprintf(stderr, "              Compress the output _f_i_l_e by gzip.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       --nn _c_u_t_-_o_f_f\n");
  fprintf(stderr, "              Output the mapping for the read sequences that  have  less  than\n");
  fprintf(stderr, "              _c_u_t_-_o_f_f number of mappings. Cannot be used with ----bbeesstt or ----ddiiss--\n");
  fprintf(stderr, "              ccoorrddaanntt--vvhh options.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----ccrroopp _n\n");
  fprintf(stderr, "              Trim the reads to _n base pairs from begining of the read.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----ttaaiill--ccrroopp _n\n");
  fprintf(stderr, "              Trim the reads to _n base pairs from end of the read.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       --ee _e_r_r_o_r_-_t_h_r_e_s_h_o_l_d\n");
  fprintf(stderr, "              Allow up to _e_r_r_o_r_-_t_h_r_e_s_h_o_l_d mismatches in the mappings.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----bbeesstt Find the best mapping location of given reads.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----ppee   Map the reads in Paired-End  mode.   mmiinn  and  mmaaxx  of  template\n");
  fprintf(stderr, "              length  will  be  calculated  if  not  provided by corresponding\n");
  fprintf(stderr, "              options.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----mmiinn _m_i_n_-_d_i_s_c_o_r_d_a_n_t_-_l_e_n_g_t_h\n");
  fprintf(stderr, "              Use _m_i_n_-_d_i_s_c_o_r_d_a_n_t_-_l_e_n_g_t_h for minimum length of concordant  map-\n");
  fprintf(stderr, "              ping. Paired-end option only.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----mmaaxx _m_a_x_-_d_i_s_c_o_r_d_a_n_t_-_l_e_n_g_t_h\n");
  fprintf(stderr, "              Use  _m_a_x_-_d_i_s_c_o_r_d_a_n_t_-_l_e_n_g_t_h for maximum length of concordant map-\n");
  fprintf(stderr, "              ping. Paired-end option only.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----ddiissccoorrddaanntt--vvhh\n");
  fprintf(stderr, "              Map the reads in discordant fashion that  can  be  processed  by\n");
  fprintf(stderr, "              Variation  Hunter / Common Law. Output will be generate in DIVET\n");
  fprintf(stderr, "              format.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----mmaaxx--ddiissccoorrddaanntt--ccuuttooffff _m\n");
  fprintf(stderr, "              Allow _m discordant mappings per read. Should be only  used  with\n");
  fprintf(stderr, "              ----ddiissccoorrddaanntt--vvhh option.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----ssnnpp _s_n_p_-_f_i_l_e\n");
  fprintf(stderr, "              Map the reads in SNP aware mode. In this mode mrsFAST-Ultra tol-\n");
  fprintf(stderr, "              erates the mismatches in known SNP  locations  reported  by  the\n");
  fprintf(stderr, "              provided SNP database. The SNP index _s_n_p_-_f_i_l_e\n");
  fprintf(stderr, "               should  be  created  from  the  dbSNP  (.vcf)  file  using  the\n");
  fprintf(stderr, "              snp_indexer binary.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       ----ssnnpp--qquuaall _q_u_a_l_i_t_y_-_t_h_r_e_s_h_o_l_d\n");
  fprintf(stderr, "              In SNP-aware mode, a mismatch at a reported SNP location will be\n");
  fprintf(stderr, "              ignored  only  if  the corresponding read location has a quality\n");
  fprintf(stderr, "              higher than or equal to the _q_u_a_l_i_t_y_-_t_h_r_e_s_h_o_l_d  _q_u_a_l_i_t_y_-_t_h_r_e_s_h_o_l_d\n");
  fprintf(stderr, "              is  a  Phred-Value  base  33. The default is 53 (base call error\n");
  fprintf(stderr, "              0.01).\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "EEXXAAMMPPLLEESS\n");
  fprintf(stderr, "       Indexing reference genome:\n");
  fprintf(stderr, "       $ ./mrsfast --index refgen.fasta\n");
  fprintf(stderr, "       $ ./mrsfast --index refgen.fasta --ws 14\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Single-end mapping:\n");
  fprintf(stderr, "       $ ./mrsfast --search refgen.fa --seq reads.fastq\n");
  fprintf(stderr, "       $ ./mrsfast --search refgen.fa --seq reads.fastq -e 3 -n 10 --threads 4\n");
  fprintf(stderr, "       $ ./mrsfast --search refgen.fa --seq reads.fastq -e 3 --best -o output\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Paired-end mapping:\n");
  fprintf(stderr, "       $ ./mrsfast --search refgen.fasta --pe --seq pe-reads.fastq  --min  100\n");
  fprintf(stderr, "       --max 400\n");
  fprintf(stderr, "       $  ./mrsfast --search refgen.fasta --pe --seq1 first-mates.fastq --seq2\n");
  fprintf(stderr, "       second-mates.fastq -e 3 --threads 4\n");
  fprintf(stderr, "       $ ./mrsfast --search refgen.fasta --pe --seq1 first-mates.fastq  --seq2\n");
  fprintf(stderr, "       second-mates.fastq --min 100 --max 400 --best -o output\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Discordant mapping:\n");
  fprintf(stderr, "       $   ./mrsfast   --search   refgen.fasta   --pe   --discordant-vh  --seq\n");
  fprintf(stderr, "       reads.fastq --min 100 --max 400\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "BBUUGGSS\n");
  fprintf(stderr, "       Please report the  bugs  through  mrsfast's  web  page  at  http://mrs-\n");
  fprintf(stderr, "       fast.sourceforge.net\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "AAuutthhoorrss\n");
  fprintf(stderr, "       Faraz Hach (fhach@sfu.ca)\n");
  fprintf(stderr, "       Iman Sarrafi (isarrafi@sfu.ca)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "RReeffeerreennccee\n");
  fprintf(stderr, "       Please cite the following paper for publications if using mrsFAST:\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Faraz Hach, Fereydoun Hormozdiari, Can Alkan, Farhad Hormozdiari, Inanc\n");
  fprintf(stderr, "       Birol, Evan E Eichler and S Cenk Sahinalp, \"mrsFAST: a  cache-oblivious\n");
  fprintf(stderr, "       algorithm for short-read mapping\", Nature Methods 7, 576-577 (2010)\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Please  cite  the  following  paper  for publications if using mrsFAST-\n");
  fprintf(stderr, "       Ultra:\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Faraz Hach, Iman Sarrafi, Farhad Hormozdiari, Can Alkan, Evan E.  Eich-\n");
  fprintf(stderr, "       ler,  S. Cenk Sahinalp, \"mrsFAST-Ultra: a compact, SNP-aware mapper for\n");
  fprintf(stderr, "       high performance sequencing applications\", Nucl.  Acids  Res.  (1  July\n");
  fprintf(stderr, "       2014) 42 (W1): W494-W500.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "CCOOPPYYRRIIGGHHTT\n");
  fprintf(stderr, "       Copyright (c) <2012-2020>, Simon Fraser University All rights reserved.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       Redistribution and use in source and binary forms, with or without mod-\n");
  fprintf(stderr, "       ification, are permitted provided that  the  following  conditions  are\n");
  fprintf(stderr, "       met:\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       1      Redistributions  of  source code must retain the above copyright\n");
  fprintf(stderr, "              notice, this list of conditions and the following disclaimer.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       2      Redistributions in binary form must reproduce  the  above  copy-\n");
  fprintf(stderr, "              right  notice,  thislist  of  conditions  and the following dis-\n");
  fprintf(stderr, "              claimer in the documentation  and/or  other  materials  provided\n");
  fprintf(stderr, "              with the distribution.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       3      Neither the name of the Simon Fraser University nor the names of\n");
  fprintf(stderr, "              its contributors may be used  to  endorse  or  promote  products\n");
  fprintf(stderr, "              derived  from  this software without specific prior written per-\n");
  fprintf(stderr, "              mission.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "       THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS\n");
  fprintf(stderr, "       IS\"  AND  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED\n");
  fprintf(stderr, "       TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTIC-\n");
  fprintf(stderr, "       ULAR  PURPOSE  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR\n");
  fprintf(stderr, "       CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,  INCIDENTAL,  SPECIAL,\n");
  fprintf(stderr, "       EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING, BUT NOT LIMITED TO,\n");
  fprintf(stderr, "       PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;  LOSS  OF  USE,  DATA,  OR\n");
  fprintf(stderr, "       PROFITS;  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF\n");
  fprintf(stderr, "       LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,  OR  TORT  (INCLUDING\n");
  fprintf(stderr, "       NEGLIGENCE  OR  OTHERWISE)  ARISING  IN  ANY WAY OUT OF THE USE OF THIS\n");
  fprintf(stderr, "       SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "mrsFAST-Ultra             Last Updated: December 3, 2018          mrsFAST-Ultra(1)\n");
}
