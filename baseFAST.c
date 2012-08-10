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
#include "CommandLineParser.h"
#include "Reads.h"
#include "Output.h"
#include "HashTable.h"
#include "MrsFAST.h"
char 				*versionNumber = "3.0";			// Current Version
unsigned char		seqFastq;
void				printStat();

int main(int argc, char *argv[])
{
	if (!parseCommandLine(argc, argv))
		return 1;

	/****************************************************
	 * INDEXING
	 ***************************************************/
	if (indexingMode)
	{
		int i;
		for (i = 0; i < fileCnt; i++)
		{
			if (!generateHashTable(fileName[i][0], fileName[i][1]))
				return 1;
		}
	}
	/****************************************************
	 * SEARCHING
	 ***************************************************/
	else
	{
		Read *seqList;
		unsigned int seqListSize;
		int fc;
		int samplingLocsSize;
		int *samplingLocs;
		double totalLoadingTime = 0;
		double totalMappingTime = 0;
		double startTime;
		double loadingTime;
		double mappingTime;
		double lstartTime;
		double tmpTime;;
		int	flag;
		double maxMem=0;
		char fname1[FILE_NAME_LENGTH];
		char fname2[FILE_NAME_LENGTH];
		char fname3[FILE_NAME_LENGTH];
		char fname4[FILE_NAME_LENGTH];
		char fname5[FILE_NAME_LENGTH];

		// Why is this one here?
		if (pairedEndMode)
		{
			//Switching to Inferred Size 
			minPairEndedDistance = minPairEndedDistance - SEQ_LENGTH + 2;
			maxPairEndedDistance = maxPairEndedDistance - SEQ_LENGTH + 2;
			if (pairedEndDiscordantMode)
			{
				maxPairEndedDiscordantDistance = maxPairEndedDiscordantDistance - SEQ_LENGTH + 2;
				minPairEndedDiscordantDistance = minPairEndedDiscordantDistance - SEQ_LENGTH + 2;
			}

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

		// Loading Sequences & Sampling Locations
		startTime = getTime();

		if (!checkHashTable(fileName[fc][1]))
		{
			return 1;
		}


		if (!readAllReads(seqFile1, seqFile2, seqCompressed, &seqFastq, pairedEndMode, &seqList, &seqListSize))
		{
			return 1;
		}
		
		//loadSamplingLocations(&samplingLocs, &samplingLocsSize);
		totalLoadingTime += getTime()-startTime;
		
		// Preparing output
		initOutput(mappingOutput, outCompressed);

		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stdout, "| %15s | %15s | %15s | %15s | %15s %15s |\n","Genome Name","Loading Time", "Mapping Time", "Memory Usage(M)","Total Mappings","Mapped reads");
		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");

		for (fc = 0; fc <fileCnt; fc++)
		{
			if (!initLoadingHashTable(fileName[fc][1]))
			{
				return 1;
			}

			loadSamplingLocations(&samplingLocs, &samplingLocsSize);

			mappingTime = 0;
			loadingTime = 0;
			flag = 1;
			do
			{
				flag = loadHashTable ( &tmpTime );  			// Reading a fragment

				loadingTime += tmpTime;

				lstartTime = getTime();
				initFAST(seqList, seqListSize);
				mapSeq(flag);
				mappingTime += getTime() - lstartTime;

				if (maxMem < getMemUsage())
					maxMem = getMemUsage();

				if (flag == 0 || flag == 2)
				{
					totalMappingTime += mappingTime;
					totalLoadingTime += loadingTime;


					fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
							getRefGenomeName(),loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
					fflush(stdout);

					loadingTime = 0;
					mappingTime = 0;
					maxMem = 0;
				}
				else if (progressRep)
				{
					fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n",
							getRefGenomeName(),loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
					fflush(stdout);
				}
			} while (flag);
		} // end for;


		finalizeFAST();
		finalizeLoadingHashTable();

		finalizeOutput();

		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stdout, "%19s%16.2f%18.2f\n\n", "Total:",totalLoadingTime, totalMappingTime);
		fprintf(stdout, "%-30s%10.2f\n","Total Time:", totalMappingTime+totalLoadingTime);
		fprintf(stdout, "%-30s%10d\n","Total No. of Reads:", seqListSize);
		fprintf(stdout, "%-30s%10lld\n","Total No. of Mappings:", mappingCnt);
		fprintf(stdout, "%-30s%10.0f\n\n","Avg No. of locations verified:", ceil((float)verificationCnt/seqListSize));

		if (strcmp(unmappedOutput, "") == 0)
		{
			char fn[strlen(mappingOutputPath) + strlen(mappingOutput) + 6 ];
			sprintf(fn, "%s%s.nohit", mappingOutputPath, mappingOutput );
			finalizeReads(fn);
		}
		else
		{
			finalizeReads(unmappedOutput);
		}
	}

	return 0;
}

