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
#include "HashTable.h"
#include "SNPReader.h"

CompressedSeq	*_snp_SNPMap		= NULL;
int				_snp_SNPMapLength	= 0;
ChrSNPs			*_snp_chrSNPs		= NULL;
int				_snp_chrCnt			= 0;
char 			**_snp_chrNames		= NULL;
int				_snp_currentChr		= 0;
int				_snp_currentLoc		= 0;

/**********************************************/
int findChrIndex(char *chrName)
{
	int i;
	char cname[CONTIG_NAME_SIZE];	// chr name in FASTA file

	for (i = 0; i < _snp_chrCnt; i++)
	{
		strcpy(cname, _snp_chrNames[i]);
		if (strlen(cname) > 3 && cname[0] == 'c' && cname[1] == 'h' && cname[2] == 'r')
			strcpy(cname, _snp_chrNames[i] + 3);		// get rid of the potential "chr" at the beginning
		if (! strcmp(cname, "MT"))						// change "MT" to "M" for consistency with dbSNP naming
			cname[1] = '\0';

		if (! strcmp(chrName, cname))
			return i;
	}
	return -1;
}
/**********************************************/
void initLoadingSNPs(char *fileName)
{
	int i, loc, t, locCnt, chrIndex, nameLen;
	char cname[CONTIG_NAME_SIZE];	// chromosome name from dbSNP
	int ccnt;						// number of chromosomes in dbSNP
	//int chrNameOffset = 0;			// used to trim the "chr" at the beginning of chromosome names
	SNPLoc *dummy = getMem(MAX_SNP_PER_CHR * sizeof(SNPLoc));

	_snp_chrCnt = getChrCnt();
	_snp_chrNames = getChrNames();
	_snp_chrSNPs = getMem(_snp_chrCnt * sizeof(ChrSNPs));
	
	for (i = 0; i < _snp_chrCnt; i++)		// FASTA chromosomes
	{
		//chrNameOffset = (strlen(_snp_chrNames[i]) > 3 && chrNames[i][0] == 'c' && chrNames[i][1] == 'h' && chrNames[i][2] == 'r') ?3 :0;
		_snp_chrSNPs[i].chrName = _snp_chrNames[i];// + chrNameOffset;

		_snp_chrSNPs[i].locCnt = 0;
		_snp_chrSNPs[i].snpLocs = NULL;	//getMem(MAX_SNP_PER_CHR * sizeof(SNPLoc));
	}

	_snp_SNPMapLength = (calculateCompressedLen(CONTIG_MAX_SIZE)+1) * sizeof(CompressedSeq);
	_snp_SNPMap = getMem(_snp_SNPMapLength);

	FILE *fp = fopen(fileName, "rt");
	t = fread(&ccnt, sizeof(int), 1, fp);		// ccnt = number of chromosomes in dbSNP

	// look for each dbSNP chromosome in the reference
	for (i = 0; i < ccnt; i++)
	{
		t = fread(&nameLen, sizeof(int), 1, fp);
		t = fread(cname, sizeof(char), nameLen, fp);
		t = fread(&locCnt, sizeof(int), 1, fp);

		cname[nameLen] = '\0';
		chrIndex = findChrIndex(cname);

		if (chrIndex != -1)		// found in FASTA chromosomes
		{
			_snp_chrSNPs[chrIndex].locCnt = locCnt;
			_snp_chrSNPs[chrIndex].snpLocs = getMem(locCnt * sizeof(SNPLoc));
			t = fread(_snp_chrSNPs[chrIndex].snpLocs, sizeof(SNPLoc), locCnt, fp);
		}
		else					// not found
		{
			t = fread(dummy, sizeof(SNPLoc), locCnt, fp);	// read dummy
			fprintf(stdout, "Warning: chromosome %s is present in the SNP database but not found in the reference genome\n", cname);
		}
	}

	fclose(fp);
	freeMem(dummy, MAX_SNP_PER_CHR * sizeof(SNPLoc));
}
/**********************************************/
void finalizeSNPs()
{
	int i;
	for (i = 0; i < _snp_chrCnt; i++)
		freeMem(_snp_chrSNPs[i].snpLocs, _snp_chrSNPs[i].locCnt * sizeof(SNPLoc));
	freeMem(_snp_chrSNPs, _snp_chrCnt * sizeof(ChrSNPs));
	freeMem(_snp_SNPMap, _snp_SNPMapLength);
}
/**********************************************/
CompressedSeq *loadSNPMap(char *chrName, int contigStartIndex, int contigLength, char *alt)
{
	//memset(_snp_SNPMap, -1, calculateCompressedLen(contigLength) * sizeof(CompressedSeq));
	memset(_snp_SNPMap, -1, _snp_SNPMapLength);
	int contigEnd = contigStartIndex + contigLength;
	int loc, offset;
	CompressedSeq *snp, mask;

	if ( strcmp(chrName, _snp_chrSNPs[_snp_currentChr].chrName) )		// new chr
	{
		_snp_currentChr++;
		_snp_currentLoc = 0;
	}

	if (_snp_chrSNPs[_snp_currentChr].locCnt)
	{
		int i = _snp_currentChr;		// just to make the code more readable
		int pos = _snp_currentLoc;

		while ( pos < _snp_chrSNPs[i].locCnt && _snp_chrSNPs[i].snpLocs[pos].loc < contigStartIndex )	// this should never happen!
			pos ++;

		while ( pos < _snp_chrSNPs[i].locCnt && _snp_chrSNPs[i].snpLocs[pos].loc < contigEnd )
		{
			loc = _snp_chrSNPs[i].snpLocs[pos].loc - contigStartIndex - 1;
			alt[loc] = _snp_chrSNPs[i].snpLocs[pos].alt;
			offset = loc % 21;
			mask = 0x7000000000000000;
			mask = ~(mask >> offset*3);
			
			snp = _snp_SNPMap + (loc/21);
			*snp &= mask;

			pos ++;
		}

		_snp_currentLoc = pos;
	}
	return _snp_SNPMap;
}
/**********************************************/
void rewindSNPIndex()
{
	_snp_currentChr = 0;
	_snp_currentLoc = 0;
}
