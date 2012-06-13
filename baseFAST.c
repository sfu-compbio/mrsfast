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
#include "CommandLineParser.h"
#include "Reads.h"
#include "Output.h"
#include "HashTable.h"
#include "HashTableBS.h"
#include "MrsFAST.h"


char 				*versionNumber = "1.1";			// Current Version
unsigned char		seqFastq;

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
		/********************************
		 * Regulard Mode
		 ********************************/
		if (!bisulfiteMode)
		{
			configHashTable();
			for (i = 0; i < fileCnt; i++)
			{
				generateHashTable(fileName[i][0], fileName[i][1]);
			}
		}
		else
		/********************************
		 * Bisulfite Mode
		 ********************************/
		{
			for (i = 0; i < fileCnt; i++)
			{
				configHashTableBS();
				generateHashTableBS(fileName[i][0], fileName[i][1]);
			}
		}
	}
	/****************************************************
	 * SEARCHING
	 ***************************************************/
	else
	{
		Read *seqList;
		unsigned int seqListSize;
		int fc, sc;
		int samplingLocsSize;
		int *samplingLocs;
		double totalLoadingTime = 0;
		double totalMappingTime = 0;
		double startTime;
		double loadingTime;
		double mappingTime;
		double lstartTime;
		double tmpTime;;
		char *prevGen = getMem(CONTIG_NAME_SIZE);
		prevGen[0]='\0';
		char *curGen;
		int	flag;
		double maxMem=0;

		// Loading Sequences & Sampling Locations
		startTime = getTime();
		if (bisulfiteMode && !pairedEndMode && seqFile1 == NULL)
		{
			if (!readAllReads(seqFile2, seqFile1, seqCompressed, &seqFastq, pairedEndMode, &seqList, &seqListSize))
			{
				return 1;
			}
		}
		else
		{
			if (!readAllReads(seqFile1, seqFile2, seqCompressed, &seqFastq, pairedEndMode, &seqList, &seqListSize))
			{
				return 1;			
			}
		}

			

		loadSamplingLocations(&samplingLocs, &samplingLocsSize);
		totalLoadingTime += getTime()-startTime;


		if (pairedEndMode)
		{
				minPairEndedDistance = minPairEndedDistance + SEQ_LENGTH + 1;
				maxPairEndedDistance = maxPairEndedDistance + SEQ_LENGTH + 1;
		}

		// Preparing output
		initOutput(mappingOutput, outCompressed);

		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stdout, "| %15s | %15s | %15s | %15s | %15s %15s |\n","Genome Name","Loading Time", "Mapping Time", "Memory Usage(M)","Total Mappings","Mapped reads");
		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");

		/********************************
		 * Regular Mode
		 ********************************/
		if (!bisulfiteMode)
		{
			Read *cseq, *cseqP;
			char rq1[SEQ_LENGTH+1];
			char rq2[SEQ_LENGTH+1];
			rq1[SEQ_LENGTH] = rq2[SEQ_LENGTH] = '\0';
			if (!pairedEndMode)
			{
				for (fc = 0; fc <fileCnt; fc++)
				{

					if (!initLoadingHashTable(fileName[fc][0]))
					{
						return 1;
					}
					mappingTime = 0;
					loadingTime = 0;
					prevGen[0] = '\0';
					flag = 1;
					int seqNo = 0;

					do 
					{
						flag = loadHashTable ( &tmpTime );  			// Reading a fragment
						curGen = getRefGenomeName();
						
						// First Time
						if (flag && prevGen[0]== '\0')
						{
							sprintf(prevGen, "%s", curGen);
						}

						if ( !flag || strcmp(prevGen, curGen)!=0)						
						{

							fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n", 
									prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
							fflush(stdout);

							totalMappingTime += mappingTime;
							totalLoadingTime += loadingTime;

							loadingTime = 0;
							mappingTime = 0;
							maxMem = 0;

							if (!flag)
							{
								break;
							}
						}
						else if (progressRep && mappingTime != 0)
						{
							fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n", 
									prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
							fflush(stdout);
						}

						sprintf(prevGen, "%s", curGen);

						loadingTime += tmpTime;
						lstartTime = getTime();

						initFAST(samplingLocs, samplingLocsSize);
						seqNo = 1;
						for (sc= 0; sc < seqListSize; sc++)
						{
							cseq = &(seqList[sc]);
							reverse(cseq->qual, rq1, SEQ_LENGTH);
							cseq->hits += mapSingleEndSeq ( cseq->name, cseq->seq, cseq->qual, cseq->hits, seqNo++, FORWARD);
							cseq->hits += mapSingleEndSeq ( cseq->name, cseq->rseq, rq1, cseq->hits, seqNo++, REVERSE);
						}
						mappingTime += getTime() - lstartTime;
						if (maxMem < getMemUsage())
						{
							maxMem = getMemUsage();
						}
					} while (flag); 

				} // end for;
				finalizeFAST();
				finalizeLoadingHashTable();
			}			
			// Pairedend Mapping Mode
			else
			{
				for (fc = 0; fc <fileCnt; fc++)
				{
					if (!initLoadingHashTable(fileName[fc][0]))
					{
						return 1;
					}
					mappingTime = 0;
					loadingTime = 0;
					prevGen[0] = '\0';
					flag = 1;
					int seqNo = 1;
					do 
					{
						flag = loadHashTable ( &tmpTime );  			// Reading a fragment
						
						curGen = getRefGenomeName();
						
						// First Time
						if (flag && prevGen[0]== '\0')
						{
							sprintf(prevGen, "%s", curGen);
						}

						if ( !flag || strcmp(prevGen, curGen)!=0)						
						{
							fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n", 
									prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
							fflush(stdout);

							totalMappingTime += mappingTime;
							totalLoadingTime += loadingTime;

							loadingTime = 0;
							mappingTime = 0;
							maxMem = 0;

							if (!flag)
							{
								break;
							}
						}
						else if (progressRep && mappingTime != 0)
						{
							fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n", 
									prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
							fflush(stdout);
						}

						sprintf(prevGen, "%s", curGen);

						loadingTime += tmpTime;
						lstartTime = getTime();

						initFAST(samplingLocs, samplingLocsSize);
						seqNo = 1;
						for (sc= 0; sc < seqListSize; sc++)
						{
							cseq = &(seqList[sc*2]);
							cseqP = &(seqList[sc*2+1]);
							
							reverse(cseq->qual, rq1, SEQ_LENGTH);
							reverse(cseqP->qual, rq2, SEQ_LENGTH);
	
                            cseq->hits += mapPairedEndSeq ( cseq->name, cseq->seq, cseq->qual, cseq->hits, seqNo, FORWARD,
															  cseqP->name, cseqP->seq, cseqP->qual, cseqP->hits, seqNo+2, FORWARD);

							cseq->hits += mapPairedEndSeq ( cseq->name, cseq->seq, cseq->qual, cseq->hits, seqNo, FORWARD, 
															  cseqP->name, cseqP->rseq, rq2, cseqP->hits, seqNo+3, REVERSE);

							cseq->hits += mapPairedEndSeq ( cseq->name, cseq->rseq, rq1, cseq->hits, seqNo+1, REVERSE,
															  cseqP->name, cseqP->rseq, rq2, cseqP->hits, seqNo+3, REVERSE);

							cseq->hits += mapPairedEndSeq ( cseq->name, cseq->rseq, rq1, cseq->hits, seqNo+1, REVERSE,
															  cseqP->name, cseqP->seq, cseqP->qual, cseqP->hits, seqNo+2, FORWARD);

							seqNo += 4; 
						}
						mappingTime += getTime() - lstartTime;
						if (maxMem < getMemUsage())
						{
							maxMem = getMemUsage();
						}
					} while (flag); 
				} // end for;
				finalizeFAST();
				finalizeLoadingHashTable();
			}
		}
		/********************************
		 * Bisulfite Mode
		 ********************************/
		else {
			int BSTableType;
			if (seqFile1!=NULL)
			{
				BSTableType = 0;
			}
			else
			{
				BSTableType = 1;
			}
			Read *cseq;
			Read *cseqP;
			char rq1[SEQ_LENGTH+1];
			char rq2[SEQ_LENGTH+1];
			rq1[SEQ_LENGTH] = rq2[SEQ_LENGTH] = '\0';
			/********************************
			 * Single End
			 ********************************/
			if (!pairedEndMode)
			{
				for (fc = 0; fc <fileCnt; fc++)
				{
					if (!initLoadingHashTableBS(fileName[fc][0]))
					{
						return 1;
					}
					mappingTime = 0;
					loadingTime = 0;
					prevGen[0] = '\0';
					flag = 1;
					int seqNo = 0;
					do 
					{
						flag = loadHashTableBS ( &tmpTime );  			// Reading a fragment
						
						curGen = getRefGenomeNameBS();
						
						// First Time
						if (flag && prevGen[0]== '\0')
						{
							sprintf(prevGen, "%s", curGen);
						}

						if ( !flag || strcmp(prevGen, curGen)!=0)						
						{
							fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n", 
									prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
							fflush(stdout);

							totalMappingTime += mappingTime;
							totalLoadingTime += loadingTime;

							loadingTime = 0;
							mappingTime = 0;
							maxMem = 0;

							if (!flag)
							{
								break;
							}
						}
						else if (progressRep && mappingTime != 0)
						{
							fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n", 
									prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
							fflush(stdout);
						}

						sprintf(prevGen, "%s", curGen);

						loadingTime += tmpTime;
						lstartTime = getTime();

						initFAST(samplingLocs, samplingLocsSize);
						seqNo = 1;
						for (sc= 0; sc < seqListSize; sc++)
						{
							cseq = &(seqList[sc]);
							reverse(cseq->qual, rq1, SEQ_LENGTH);
							cseq->hits += mapSingleEndSeqBS( cseq->name, cseq->seq, cseq->qual, cseq->hits, seqNo++, FORWARD, BSTableType);
							cseq->hits += mapSingleEndSeqBS( cseq->name, cseq->rseq, rq1, cseq->hits, seqNo++, REVERSE, !BSTableType);
						}
						mappingTime += getTime() - lstartTime;
						if (maxMem < getMemUsage())
						{
							maxMem = getMemUsage();
						}
					} while (flag); 

				} // end for;
				finalizeFAST();
				finalizeLoadingHashTableBS();
			}
			/********************************
			 * Paired End
			 ********************************/
			else
			{



				for (fc = 0; fc <fileCnt; fc++)
				{
					if (!initLoadingHashTableBS(fileName[fc][0]))
					{
						return 1;
					}
					mappingTime = 0;
					loadingTime = 0;
					prevGen[0] = '\0';
					flag = 1;
					int seqNo = 1;
					do 
					{
						flag = loadHashTableBS ( &tmpTime );  			// Reading a fragment
						
						curGen = getRefGenomeNameBS();
						
						// First Time
						if (flag && prevGen[0]== '\0')
						{
							sprintf(prevGen, "%s", curGen);
						}

						if ( !flag || strcmp(prevGen, curGen)!=0)						
						{
							fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n", 
									prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
							fflush(stdout);

							totalMappingTime += mappingTime;
							totalLoadingTime += loadingTime;

							loadingTime = 0;
							mappingTime = 0;
							maxMem = 0;

							if (!flag)
							{
								break;
							}
						}
						else if (progressRep && mappingTime != 0)
						{
							fprintf(stdout, "| %15s | %15.2f | %15.2f | %15.2f | %15lld %15lld |\n", 
									prevGen,loadingTime, mappingTime, maxMem, mappingCnt , mappedSeqCnt);
							fflush(stdout);
						}

						sprintf(prevGen, "%s", curGen);

						loadingTime += tmpTime;
						lstartTime = getTime();

						initFAST(samplingLocs, samplingLocsSize);
						seqNo = 1;
						for (sc= 0; sc < seqListSize; sc++)
						{
							cseq = &(seqList[sc*2]);
							cseqP = &(seqList[sc*2+1]);
	
							reverse(cseq->qual, rq1, SEQ_LENGTH);
							reverse(cseqP->qual, rq2, SEQ_LENGTH);

                            cseq->hits += mapPairedEndSeqBS ( cseq->name, cseq->seq, cseq->qual, cseq->hits, seqNo, FORWARD, BSTableType,
															  cseqP->name, cseqP->seq, cseqP->qual, cseqP->hits, seqNo+2, FORWARD, !BSTableType);
							cseqP->hits = cseq->hits;

							cseq->hits += mapPairedEndSeqBS ( cseq->name, cseq->seq, cseq->qual, cseq->hits, seqNo, FORWARD, BSTableType,
															  cseqP->name, cseqP->rseq, rq2, cseqP->hits, seqNo+3, REVERSE, BSTableType);
							cseqP->hits = cseq->hits;

							cseq->hits += mapPairedEndSeqBS ( cseq->name, cseq->rseq, rq1, cseq->hits, seqNo+1, REVERSE, !BSTableType,
															  cseqP->name, cseqP->rseq, rq2, cseqP->hits, seqNo+3, REVERSE, BSTableType);
							cseqP->hits = cseq->hits;

							cseq->hits += mapPairedEndSeqBS ( cseq->name, cseq->rseq, rq1, cseq->hits, seqNo+1, REVERSE, !BSTableType,
															  cseqP->name, cseqP->seq, cseqP->qual, cseqP->hits, seqNo+2, FORWARD, !BSTableType);
							cseqP->hits = cseq->hits;

							seqNo += 4; 
						}
						mappingTime += getTime() - lstartTime;
						if (maxMem < getMemUsage())
						{
							maxMem = getMemUsage();
						}
					} while (flag); 

				} // end for;
				finalizeFAST();
				finalizeLoadingHashTableBS();


			}
		}
		////////////////////////
		finalizeOutput();

		fprintf(stdout, "-----------------------------------------------------------------------------------------------------------\n");
		fprintf(stdout, "%19s%16.2f%18.2f\n\n", "Total:",totalLoadingTime, totalMappingTime);
		fprintf(stdout, "%-30s%10.2f\n","Total Time:", totalMappingTime+totalLoadingTime);
		fprintf(stdout, "%-30s%10d\n","Total No. of Reads:", seqListSize);
		fprintf(stdout, "%-30s%10lld\n","Total No. of Mappings:", mappingCnt);
		fprintf(stdout, "%-30s%10.0f\n\n","Avg No. of locations verified:", ceil((float)verificationCnt/seqListSize));
		
		if (bisulfiteMode)
		{
			fprintf(stdout, "%-30s%10.2f\n\n","Bisulfite Conversion Rate:", (float)regTotal/(regTotal+metTotal)*100);			
		}

		int cof = (pairedEndMode)?2:1;

		if (progressRep && maxHits != 0)
		{
			int frequency[maxHits+1];
			int i;
			for ( i=0 ; i <= maxHits; i++)
			{
				frequency[i] = 0;
			}


			for (fc = 0; fc < seqListSize; fc++)
			{
				frequency[seqList[fc*cof].hits]++;
			}
			frequency[maxHits] = completedSeqCnt;
			for ( i=0 ; i <= maxHits; i++)
			{
				fprintf(stdout, "%-30s%10d%10d%10.2f%%\n","Reads Mapped to ", i, frequency[i], 100*(float)frequency[i]/(float)seqListSize);
			}
		}

		
		finalizeReads(unmappedOutput);
		freeMem(prevGen, CONTIG_NAME_SIZE);
	}



	return 1;
}
