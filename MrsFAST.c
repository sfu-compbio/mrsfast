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
#include "Common.h"
#include "HashTable.h"
#include "HashTableBS.h"
#include "Output.h"

unsigned char		mrFAST = 0;
char				*versionNumberF="6.2";


long long			verificationCnt = 0;
long long 			mappingCnt = 0;
long long			completedSeqCnt = 0;
long long			mappedSeqCnt = 0;	
long long			metTotal = 0;
long long			regTotal = 0;

char				*_msf_refGen = NULL;
int					_msf_refGenLength = 0;
int					_msf_refGenOffset = 0;
char				*_msf_refGenName = NULL;

HashTable			*_msf_hashTable = NULL;
HashTableBS			*_msf_hashTableBS = NULL;

int					*_msf_verifiedLocs = NULL;
int					*_msf_verifiedLocsP = NULL;

int					*_msf_samplingLocs;
int					_msf_samplingLocsSize;

SAM					_msf_output;
SAM					_msf_outputP;

OPT_FIELDS			*_msf_optionalFields;

char				*_msf_op;
char				*_msf_opP;

/**********************************************/
void initFAST(int *samplingLocs, int samplingLocsSize)
{
	if (_msf_optionalFields == NULL)	
	{
		_msf_op = getMem(SEQ_LENGTH);
		if (pairedEndMode)
		{
			_msf_opP = getMem(SEQ_LENGTH);
			_msf_optionalFields = getMem(4*sizeof(OPT_FIELDS));
		}	
		else
		{
			_msf_optionalFields = getMem(2*sizeof(OPT_FIELDS));
		}
	}

	if (_msf_verifiedLocs != NULL)
		freeMem(_msf_verifiedLocs, sizeof(int) * (_msf_refGenLength+1));
	
	if (_msf_verifiedLocsP != NULL )
		freeMem(_msf_verifiedLocsP, sizeof(int) * (_msf_refGenLength+1));

	_msf_samplingLocs = samplingLocs;
	_msf_samplingLocsSize = samplingLocsSize;

	if (!bisulfiteMode)
	{
		_msf_refGen =  getRefGenome();
		_msf_refGenLength = strlen(_msf_refGen);
		_msf_refGenOffset = getRefGenomeOffset();
		_msf_refGenName = getRefGenomeName();

		_msf_hashTable = getHashTable();
	}
	else
	{
		_msf_refGen =  getRefGenomeBS();
		_msf_refGenLength = strlen(_msf_refGen);
		_msf_refGenOffset = getRefGenomeOffsetBS();
		_msf_refGenName = getRefGenomeNameBS();

		_msf_hashTableBS = getHashTableBS();
	}

	
	_msf_verifiedLocs = getMem(sizeof(int)*(_msf_refGenLength+1));
	int i;
	for (i=0; i<=_msf_refGenLength; i++)
		_msf_verifiedLocs[i] = 0;

	if (pairedEndMode)
	{
		_msf_verifiedLocsP = getMem(sizeof(int)*(_msf_refGenLength+1));
		int i;
		for (i=0; i<=_msf_refGenLength; i++)
			_msf_verifiedLocsP[i] = 0;
	}
}
/**********************************************/
void finalizeFAST()
{
	freeMem(_msf_op, SEQ_LENGTH);

	if (pairedEndMode)
	{
		freeMem(_msf_opP, SEQ_LENGTH);
		freeMem(_msf_optionalFields, 4*sizeof(OPT_FIELDS));
	}
	else
	{
		freeMem(_msf_optionalFields, 1*sizeof(OPT_FIELDS));
	}

	if (_msf_verifiedLocs != NULL)
		freeMem(_msf_verifiedLocs, sizeof(int) * (_msf_refGenLength+1));
}
/**********************************************/
int verify(int index, char* seq, char **opSeq)
{
	int i;
	int err, errCnt;
	char *ref;
	char *ver;
	short matchCnt = 0;
	char *op = *opSeq;
	int pp = 0;

	verificationCnt ++;

	errCnt = 0;
	
	ref = _msf_refGen + index-1;
	ver = seq;
	
	for (i=0; i < SEQ_LENGTH; i++)
	{	
		err = * ref != *ver;
		errCnt += err;

		if (errCnt > errThreshold)
		{	
			errCnt = -1;
			break;
		}
		else
		{
			if (err && matchCnt)
			{
				if (matchCnt < 10)
				{
					op[pp++]=48+matchCnt;
					op[pp++]=*ref;
				}
				else if (matchCnt < 100)
				{
					op[pp++]=48+(matchCnt/10);
					op[pp++]=48+(matchCnt%10);
					op[pp++] = *ref;
				}
				else
				{
					op[pp++] = 48+matchCnt/100;
					matchCnt %= 100;
					op[pp++] = 48+matchCnt/10;
					op[pp++] = 48+matchCnt%10;
					op[pp++] = *ref;
				}

				matchCnt = 0;
			}
			else if (err && matchCnt == 0)
			{
				op[pp++]=*ref;
			}
			else
			{
				matchCnt++;
			}
		}
		ref++;
		ver++;
	}
	if (matchCnt)
	{
		if (matchCnt < 10)
		{
			op[pp++]=48+matchCnt;
		}
		else if (matchCnt < 100)
		{
			op[pp++]=48+matchCnt/10;
			op[pp++]=48+matchCnt%10;
		}
		else
		{
			op[pp++]=48+matchCnt/100;
			matchCnt %= 100;
			op[pp++]=48+matchCnt/10;
			op[pp++]=48+matchCnt%10;
		}
	}
	op[pp]='\0';
	return errCnt;
}
/**********************************************/
int mapSingleEndSeq(char *seqName, char *seq, char* seqQual, unsigned char seqHits, int seqNo, short mappingDirection)
{

	if (maxHits !=0 && maxHits == seqHits)
	{
		return 0;
	}


	int i = 0;
	int j = 0;
	int flag = 1;
	int offset;
	int genLoc;
	int err;
	int size;
	int hv;
	int hits = 0;
	unsigned int *locs;
	char cigar[5];

	sprintf(cigar, "%dM", SEQ_LENGTH);

	while ( flag && i < _msf_samplingLocsSize )
	{
		offset = _msf_samplingLocs[i];
		hv = hashVal(seq + offset);
		j = 1;
		locs = getCandidates(hv);
		if (locs != NULL)
		{
			size = locs[0];
		}
		else
		{
			size = 0;
		}

		while (j <= size)
		{
			genLoc = locs[j] - offset;

			if ( genLoc <  1 ||
					genLoc > (_msf_refGenLength- SEQ_LENGTH + 1) ||
					_msf_verifiedLocs[genLoc] ==  seqNo ||
					_msf_verifiedLocs[genLoc] == -seqNo )
			{
				j++;
			}
			else
			{
				err = verify(genLoc, seq, &_msf_op);
				if ( err != -1 )
				{
					_msf_verifiedLocs[genLoc] = seqNo;
					mappingCnt++;

					_msf_output.QNAME		= seqName;
					_msf_output.FLAG		= 16 * mappingDirection; 
					_msf_output.RNAME		= _msf_refGenName;
					_msf_output.POS			= genLoc + _msf_refGenOffset;
					_msf_output.MAPQ		= 255;
					_msf_output.CIGAR		= cigar;
					_msf_output.MRNAME		= "*";
					_msf_output.MPOS		= 0;
					_msf_output.ISIZE		= 0;
					_msf_output.SEQ			= seq;					
					_msf_output.QUAL		= seqQual;

					_msf_output.optSize		= 2;
					_msf_output.optFields	= _msf_optionalFields;

					_msf_optionalFields[0].tag = "NM";
					_msf_optionalFields[0].type = 'i';
					_msf_optionalFields[0].iVal = err;


					_msf_optionalFields[1].tag = "MD";
					_msf_optionalFields[1].type = 'Z';
					_msf_optionalFields[1].sVal = _msf_op;

					output(_msf_output);

					if (hits+seqHits == 0)
					{
						mappedSeqCnt ++;
					}

					if ( maxHits == 0 && hits+seqHits == 0)
					{
						hits=1;
					}
					else if (maxHits != 0 )
					{
						hits++;
						if (hits+seqHits == maxHits)
						{
							completedSeqCnt++;
							flag =0;
							break;
						}
					}
				}
				else
				{
					_msf_verifiedLocs[genLoc] = -seqNo;
				}
				j++;

			}
		}
		i++;
	}
	return hits;
}

/**********************************************/
int mapPairedEndSeq(	char *seqName, char *seq1, char* seq1Qual, unsigned int seq1Hits, int seq1No, short mappingDirection1,
						char *seq2Name, char *seq2, char* seq2Qual, unsigned int seq2Hits, int seq2No, short mappingDirection2)
{
	
	if (maxHits !=0 && maxHits == seq1Hits)
	{
		return 0;
	}


	int i = 0;
	int j = 0;

	int flag = 1;
	int offset1;
	int offset2;
	int genLoc1;
	int genLoc2;
	int size1;
	int size2;
	int hv1;
	int hv2;
	unsigned int *locs1;
	unsigned int *locs2;

	int p1, p2, p3;

	int err1, err2;
	int hits = 0;

	int nv = 0;
	int k;

	int lm, ll, rl, rm; // leftmost, leftleast, rightleast, rightmost;

	char cigar[5];

	sprintf(cigar, "%dM", SEQ_LENGTH);

	while (flag && i < _msf_samplingLocsSize)
	{
		offset1 = _msf_samplingLocs[i];
		hv1 = hashVal(seq1 + offset1);
		locs1 = getCandidates(hv1);
		if (locs1 != NULL)
		{
			size1 = locs1[0];
		}
		else
		{
			size1 = 0;
		}

        j = 0;
        while (flag && j < _msf_samplingLocsSize)
        {
            offset2 = _msf_samplingLocs[j];
            hv2 = hashVal(seq2 + offset2);
            locs2 = getCandidates(hv2);
			if ( locs2 != NULL)
			{
				size2 = locs2[0];
			}
			else
			{
				size2 = 0;
			}
            if (size2>0)
            {
                p1=1;
                p2=1;
                p3=1;
				// Main Loop
                while (flag && p1<=size1 && p2 <=size2)
                {
					genLoc1 = locs1[p1] - offset1;
                    nv = 0;
                    if ( genLoc1 <  1 ||
                         genLoc1 > (_msf_refGenLength - SEQ_LENGTH + 1) ||
                         _msf_verifiedLocs[genLoc1] == -seq1No )
                    {
                        p1++;
                        continue;
                    }
                    else
                    {
                        if (_msf_verifiedLocs[genLoc1] == seq1No)
                        {
                            err1 = verify(genLoc1, seq1, &_msf_op);
                        }
                        else
                        {
                            err1 = verify(genLoc1, seq1, &_msf_op);

                            if ( err1 != -1 )
                            {
                                _msf_verifiedLocs[genLoc1] = seq1No;
                                nv = 1;
                            }
                            else
                            {
                                _msf_verifiedLocs[genLoc1] = -seq1No;
                                p1++;
                                continue;
                            }
                        }

                    }

					lm = genLoc1 - maxPairEndedDistance + 1;
					ll = genLoc1 - minPairEndedDistance + 1;
					rl = genLoc1 + minPairEndedDistance - 1;
					rm = genLoc1 + maxPairEndedDistance - 1;

//					fprintf(stdout, "%d %d [%d] %d %d\n", lm, ll,genLoc1, rl, rm);

                    while (p2<=size2 && (locs2[p2]-offset2) < lm)
                    {
                        p2++;
                    }
                    while (p3<=size2 && (locs2[p3]-offset2) <= rm)
                    {
                        p3++;
                    }

                    k=p2;
                    while (flag && k<=p3)
                    {
						genLoc2 = locs2[k] - offset2;
                        if ( genLoc2 >  ll && genLoc2 < rl )
                        {
                            k++;
                            continue;
                        }

                        if (genLoc2 < 1 ||
                            genLoc2 > (_msf_refGenLength - SEQ_LENGTH + 1)||
                            _msf_verifiedLocsP[genLoc2] == -seq2No)
                        {
                            k++;
                        }
                        else
                        {

                            err2 = -1;
                            if (_msf_verifiedLocsP[genLoc2] == seq2No)
                            {
                                err2 = verify(genLoc2, seq2, &_msf_opP);
                            }
                            else if (_msf_verifiedLocsP[genLoc2] != -seq2No)
                            {
                                err2 = verify(genLoc2, seq2, &_msf_opP);
                            }

                            if ( (err2 != -1 && _msf_verifiedLocsP[genLoc2] != seq2No) ||
								 (_msf_verifiedLocsP[genLoc2] == seq2No && nv) )
                            {
								mappingCnt++;
                                _msf_verifiedLocsP[genLoc2] = seq2No;
								
								_msf_output.POS			= genLoc1 + _msf_refGenOffset;
								_msf_output.MPOS		= genLoc2 + _msf_refGenOffset;
								_msf_output.FLAG		= 1+2+16*mappingDirection1+32*mappingDirection2+64;
								_msf_output.ISIZE		= abs(genLoc2 - genLoc1) + SEQ_LENGTH;
								_msf_output.SEQ			= seq1,
								_msf_output.QUAL		= seq1Qual;
								_msf_output.QNAME		= seqName;
								_msf_output.RNAME		= _msf_refGenName;
								_msf_output.MAPQ		= 255;
								_msf_output.CIGAR		= cigar;
								_msf_output.MRNAME		= "=";
								
								_msf_output.optSize	= 2;
								_msf_output.optFields	= _msf_optionalFields;

								_msf_optionalFields[0].tag = "NM";
								_msf_optionalFields[0].type = 'i';
								_msf_optionalFields[0].iVal = err1;

								_msf_optionalFields[1].tag = "MD";
								_msf_optionalFields[1].type = 'Z';
								_msf_optionalFields[1].sVal = _msf_op;


								output(_msf_output);
								
								_msf_output.POS			= genLoc2 + _msf_refGenOffset;
								_msf_output.MPOS		= genLoc1 + _msf_refGenOffset;
								_msf_output.FLAG		= 1+2+16*mappingDirection2+32*mappingDirection1+128;
								_msf_output.ISIZE		= abs(genLoc2 - genLoc1) + SEQ_LENGTH;
								_msf_output.SEQ			= seq2,
								_msf_output.QUAL		= seq2Qual;
								_msf_output.QNAME		= seqName;
								_msf_output.RNAME		= _msf_refGenName;
								_msf_output.MAPQ		= 255;
								_msf_output.CIGAR		= cigar;
								_msf_output.MRNAME		= "=";
								
								_msf_output.optSize	= 2;
								_msf_output.optFields	= _msf_optionalFields;

								_msf_optionalFields[0].tag = "NM";
								_msf_optionalFields[0].type = 'i';
								_msf_optionalFields[0].iVal = err2;

								_msf_optionalFields[1].tag = "MD";
								_msf_optionalFields[1].type = 'Z';
								_msf_optionalFields[1].sVal = _msf_opP;

								output(_msf_output);
                                
								if (hits+seq1Hits == 0)
								{
									mappedSeqCnt ++;
								}

								if ( maxHits == 0 && hits+seq1Hits == 0)
								{
									hits=1;
								}
								else if (maxHits != 0 )
								{
									hits++;
									if (hits+seq1Hits == maxHits)
									{
										completedSeqCnt++;
										flag =0;
										break;
									}
								}
                            }
                            else if (err2 == -1)
                            {
                                _msf_verifiedLocsP[genLoc2] = -seq2No;
                            }
                            k++;
                        }
                    }
                    p1++;
                }// END of Mail Loop;
            }
            j++;
		}
		i++;
	}

	return hits;
}


/**********************************************/
int verifyBS(int index, char *seq, char **opSeq, int type, int *met, int *reg)
{
	int i;
	int errCnt, err;
	char *ref;
	char *ver;
	short matchCnt=0;

	verificationCnt ++;
	errCnt = 0; 
	ref = _msf_refGen +index - 1;
	ver = seq;

	*met = 0;
	*reg = 0;

	int pp = 0;
	char *op=*opSeq;

	for (i=0; i < SEQ_LENGTH; i++)
	{
		err = 0;
	    if (!type && *ver == 'T' && *ref == 'C')
        {
            (*met)++;
        }
        else if (!type && *ver == 'C' && *ref == 'C')
        {
            (*reg)++;
        }
        else if (type && *ver == 'A' && *ref == 'G')
        {
            (*met)++;
        }
        else if (type && *ver == 'A' && *ref == 'A')
        {
            (*reg)++;
        }
        else
		{
            err =(*ref != *ver);
		}

		errCnt += err;

		if (errCnt > errThreshold)
		{
			errCnt = -1;
			break;
		}
		else
		{
			if (err && matchCnt)
			{
				if (matchCnt < 10)
				{
					op[pp++]=48+matchCnt;
					op[pp++]=*ref;
				}
				else if (matchCnt < 100)
				{
					op[pp++]=48+(matchCnt/10);
					op[pp++]=48+(matchCnt%10);
					op[pp++] = *ref;
				}
				else
				{
					op[pp++] = 48+matchCnt/100;
					matchCnt %= 100;
					op[pp++] = 48+matchCnt/10;
					op[pp++] = 48+matchCnt%10;
					op[pp++] = *ref;
				}

				matchCnt = 0;
			}
			else if (err && matchCnt == 0)
			{
				op[pp++]=*ref;
			}
			else
			{
				matchCnt++;
			}
		}
		ref++;
		ver++;
	}
	if (matchCnt)
	{
		if (matchCnt < 10)
		{
			op[pp++]=48+matchCnt;
		}
		else if (matchCnt < 100)
		{
			op[pp++]=48+matchCnt/10;
			op[pp++]=48+matchCnt%10;
		}
		else
		{
			op[pp++]=48+matchCnt/100;
			matchCnt %= 100;
			op[pp++]=48+matchCnt/10;
			op[pp++]=48+matchCnt%10;
		}
	}
	op[pp]='\0';

	return errCnt;
}
/**********************************************/
int mapSingleEndSeqBS( char *seqName, char *seq, char* seqQual, unsigned int seqHits, int seqNo, short mappingDirection, int type)
{

	if (maxHits !=0 && maxHits == seqHits)
	{
		return 0;
	}


	int i = 0;
	int j = 0;
	int flag = 1;
	int offset;
	int genLoc;
	int err;
	int size;
	int hv;
	int hits = 0;
	int reg = 0;
	int met = 0;
	unsigned int *locs;
	char cigar[5];

	sprintf(cigar, "%dM", SEQ_LENGTH);
	
	while (flag && i <  _msf_samplingLocsSize )
	{
		offset = _msf_samplingLocs[i];
		hv = hashValBS(seq + offset, type);
		j = 1;
		locs = getCandidatesBS(hv, type);
		if (locs != NULL)
		{
			size = locs[0];
		}
		else
		{
			size = 0;
		}

		while (j <= size)
		{
			genLoc = locs[j] - offset;
			if ( genLoc <  1 ||
				 genLoc > ( _msf_refGenLength - SEQ_LENGTH + 1) ||
				 _msf_verifiedLocs[genLoc] ==  seqNo ||
				 _msf_verifiedLocs[genLoc] == -seqNo )
			{
				j++;
			}
			else
			{
				err = verifyBS(genLoc, seq, &_msf_op, type, &met, &reg);
				if ( err != -1 )
				{

					metTotal += met;
					regTotal += reg;					
					mappingCnt++;

					_msf_verifiedLocs[genLoc] = seqNo;

					_msf_output.QNAME		= seqName;
					_msf_output.FLAG		= 16 * mappingDirection; 
					_msf_output.RNAME		= _msf_refGenName;
					_msf_output.POS			= genLoc + _msf_refGenOffset;
					_msf_output.MAPQ		= 255;
					_msf_output.CIGAR		= cigar;
					_msf_output.MRNAME		= "*";
					_msf_output.MPOS		= 0;
					_msf_output.ISIZE		= 0;
					_msf_output.SEQ			= seq;					
					_msf_output.QUAL		= seqQual;

					_msf_output.optSize		= 2;
					_msf_output.optFields	= _msf_optionalFields;

					_msf_optionalFields[0].tag = "NM";
					_msf_optionalFields[0].type = 'i';
					_msf_optionalFields[0].iVal = err;


					_msf_optionalFields[1].tag = "MD";
					_msf_optionalFields[1].type = 'Z';
					_msf_optionalFields[1].sVal = _msf_op;

					output(_msf_output);



					if (hits+seqHits == 0)
					{
						mappedSeqCnt ++;
					}

					if ( maxHits == 0 && hits+seqHits == 0)
					{
						hits=1;
					}
					else if (maxHits != 0 )
					{
						hits++;
						if (hits+seqHits == maxHits)
						{
							completedSeqCnt++;
							flag =0;
							break;
						}
					}
				}
				else
				{
					_msf_verifiedLocs[genLoc] = -seqNo;
				}
				j++;

			}
		}
		i++;
	}
	return hits;
}
/**********************************************/
int mapPairedEndSeqBS(	char *seqName, char *seq1, char* seq1Qual, unsigned int seq1Hits, int seq1No, short mappingDirection1, int type1,
						char *seq2Name, char *seq2, char* seq2Qual, unsigned int seq2Hits, int seq2No, short mappingDirection2, int type2)
{
	char tmp1[2*SEQ_LENGTH+2];
	char tmp2[2*SEQ_LENGTH+2];
	tmp1[2*SEQ_LENGTH+2]=tmp2[2*SEQ_LENGTH+2] = '\0';

	if (maxHits !=0 && maxHits == seq1Hits)
	{
		return 0;
	}


	int i = 0;
	int j = 0;

	int flag = 1;
	int offset1;
	int offset2;
	int genLoc1;
	int genLoc2;
	int size1;
	int size2;
	int hv1;
	int hv2;
	unsigned int *locs1;
	unsigned int *locs2;

	int p1, p2, p3;

	int err1, err2;
	int hits = 0;

	int nv = 0;
	int k;

	int lm, ll, rl, rm; // leftmost, leftleast, rightleast, rightmost;

	int reg1=0, met1=0, reg2=0, met2=0;

	char cigar[5];

	sprintf(cigar, "%dM", SEQ_LENGTH);
	
	while (flag && i < _msf_samplingLocsSize)
	{
		offset1 = _msf_samplingLocs[i];
		hv1 = hashValBS(seq1 + offset1, type1);
		locs1 = getCandidatesBS(hv1, type1);
		if (locs1 != NULL)
		{
			size1 = locs1[0];
		}
		else
		{
			size1 = 0;
		}
        j = 0;
        while (flag && j < _msf_samplingLocsSize)
        {
            offset2 = _msf_samplingLocs[j];
            hv2 = hashValBS(seq2 + offset2, type2);
			locs2 = getCandidatesBS(hv2, type2);
			if (locs2 != NULL)
			{
				size2 = locs2[0];
			}
			else
			{
				size2 = 0;
			}
			if (size2>0)
            {
                p1=1;
                p2=1;
                p3=1;
				// Main Loop
                while (flag && p1<=size1 && p2 <=size2)
                {
					genLoc1 = locs1[p1] - offset1;
                    nv = 0;
                    if ( genLoc1 <  1 ||
                         genLoc1 > (_msf_refGenLength - SEQ_LENGTH + 1) ||
                         _msf_verifiedLocs[genLoc1] == -seq1No )
                    {
                        p1++;
                        continue;
                    }
                    else
                    {
                        if (_msf_verifiedLocs[genLoc1] == seq1No)
                        {
                            err1 = verifyBS(genLoc1, seq1, &_msf_op, type1, &met1, &reg1);
                        }
                        else
                        {
                            err1 = verifyBS(genLoc1, seq1, &_msf_op, type1, &met1, &reg1);

                            if ( err1 != -1 )
                            {
                                _msf_verifiedLocs[genLoc1] = seq1No;
                                nv = 1;
                            }
                            else
                            {
                                _msf_verifiedLocs[genLoc1] = -seq1No;
                                p1++;
                                continue;
                            }
                        }

                    }

					lm = genLoc1-maxPairEndedDistance+1;
					ll = genLoc1-minPairEndedDistance+1;
					rl = genLoc1+minPairEndedDistance-1;
					rm = genLoc1+maxPairEndedDistance-1;

					//fprintf(stderr, "%d %d [%d] %d %d\n", lm, ll,genLoc1, rl, rm);

                    while (p2<=size2 && (locs2[p2]-offset2) < lm)
                    {
                        p2++;
                    }
                    while (p3<=size2 && (locs2[p3]-offset2) <= rm)
                    {
                        p3++;
                    }

                    k=p2;
                    while (flag && k<p3)
                    {
						genLoc2 = locs2[k] - offset2;
                        if ( genLoc2 >  ll && genLoc2 < rl )
                        {
                            k++;
                            continue;
                        }

                        if (genLoc2 < 1 ||
                            genLoc2 > (_msf_refGenLength - SEQ_LENGTH + 1)||
                            _msf_verifiedLocsP[genLoc2] == -seq2No)
                        {
                            k++;
                        }
                        else
                        {

                            err2 = -1;
                            if (_msf_verifiedLocsP[genLoc2] == seq2No)
                            {
                                err2=verifyBS(genLoc2, seq2, &_msf_opP, type2, &met2, &reg2);
                            }
                            else if (_msf_verifiedLocsP[genLoc2] != -seq2No)
                            {
                                err2 = verifyBS(genLoc2, seq2, &_msf_opP, type2, &met2, &reg2);
                            }

                            if ( (err2 != -1 && _msf_verifiedLocsP[genLoc2] != seq2No) ||
								 (_msf_verifiedLocsP[genLoc2] == seq2No && nv) )
                            {

								metTotal += met1+met2;
								regTotal += reg1+reg2;
                                mappingCnt++;

                                _msf_verifiedLocsP[genLoc2] = seq2No;

								_msf_output.POS			= genLoc1 + _msf_refGenOffset;
								_msf_output.MPOS		= genLoc2 + _msf_refGenOffset;
								_msf_output.FLAG		= 1+2+16*mappingDirection1+32*mappingDirection2+64;
								_msf_output.ISIZE		= abs(genLoc2 - genLoc1) + SEQ_LENGTH;
								_msf_output.SEQ			= seq1,
								_msf_output.QUAL		= seq1Qual;
								_msf_output.QNAME		= seqName;
								_msf_output.RNAME		= _msf_refGenName;
								_msf_output.MAPQ		= 255;
								_msf_output.CIGAR		= cigar;
								_msf_output.MRNAME		= "=";
								
								_msf_output.optSize	= 2;
								_msf_output.optFields	= _msf_optionalFields;

								_msf_optionalFields[0].tag = "NM";
								_msf_optionalFields[0].type = 'i';
								_msf_optionalFields[0].iVal = err1;

								_msf_optionalFields[1].tag = "MD";
								_msf_optionalFields[1].type = 'Z';
								_msf_optionalFields[1].sVal = _msf_op;


								output(_msf_output);
								
								_msf_output.POS			= genLoc2 + _msf_refGenOffset;
								_msf_output.MPOS		= genLoc1 + _msf_refGenOffset;
								_msf_output.FLAG		= 1+2+16*mappingDirection2+32*mappingDirection1+128;
								_msf_output.ISIZE		= abs(genLoc2 - genLoc1) + SEQ_LENGTH;
								_msf_output.SEQ			= seq2,
								_msf_output.QUAL		= seq2Qual;
								_msf_output.QNAME		= seqName;
								_msf_output.RNAME		= _msf_refGenName;
								_msf_output.MAPQ		= 255;
								_msf_output.CIGAR		= cigar;
								_msf_output.MRNAME		= "=";
								
								_msf_output.optSize	= 2;
								_msf_output.optFields	= _msf_optionalFields;

								_msf_optionalFields[0].tag = "NM";
								_msf_optionalFields[0].type = 'i';
								_msf_optionalFields[0].iVal = err2;

								_msf_optionalFields[1].tag = "MD";
								_msf_optionalFields[1].type = 'Z';
								_msf_optionalFields[1].sVal = _msf_opP;

								output(_msf_output);
                                
								if (hits+seq1Hits == 0)
								{
									mappedSeqCnt ++;
								}

								if ( maxHits == 0 && hits+seq1Hits == 0)
								{
									hits=1;
								}
								else if (maxHits != 0 )
								{
									hits++;
									if (hits+seq1Hits == maxHits)
									{
										completedSeqCnt++;
										flag =0;
										break;
									}
								}
                            }
                            else if (err2 == -1)
                            {
                                _msf_verifiedLocsP[genLoc2] = -seq2No;
                            }
                            k++;
                        }
                    }
                    p1++;
                }// END of Mail Loop;
            }
            j++;
		}
		i++;
	}

	return hits;
}
