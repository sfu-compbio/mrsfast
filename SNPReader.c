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
int				_snp_currentChr		= 0;
int				_snp_currentLoc	= 0;

/**********************************************/
void initLoadingSNPs(char *fileName)
{
	int i, j, loc, t, found, count;
	char cname[CONTIG_NAME_SIZE];	// chromosome name
	int ccnt;						// chr count in the file
	char **chrNames;
	SNPLoc *dummy = getMem(MAX_SNP_PER_CHR * sizeof(SNPLoc));

	_snp_chrCnt = getChrCnt();
	chrNames = getChrNames();
	_snp_chrSNPs = getMem(_snp_chrCnt * sizeof(ChrSNPs));
	
	for (i = 0; i < _snp_chrCnt; i++)
	{
		_snp_chrSNPs[i].chrName = chrNames[i];
		_snp_chrSNPs[i].locCnt = 0;
		_snp_chrSNPs[i].snpLocs = getMem(MAX_SNP_PER_CHR * sizeof(SNPLoc));
	}

	_snp_SNPMapLength = (calculateCompressedLen(CONTIG_MAX_SIZE)+1) * sizeof(CompressedSeq);
	_snp_SNPMap = getMem(_snp_SNPMapLength);

	strcpy(cname, "chr");
	FILE *fp = fopen(fileName, "rt");
	t = fread(&ccnt, sizeof(int), 1, fp);	// must be 25 for homosapien

	for (i = 1; i <= ccnt; i++)
	{
		if (i < 23)
			sprintf(cname+3, "%d", i);
		else
		{
			switch (i)
			{
				case 23: cname[3] = 'X'; break;
				case 24: cname[3] = 'Y'; break;
				case 25: cname[3] = 'M'; break;
				default: cname[3] = 'T'; break;		// something invalid
			}
			cname[4] = '\0';
		}

		found = 0;
		for (j = 0; j < _snp_chrCnt; j++)
		{
			if (!strcmp(cname, _snp_chrSNPs[j].chrName))
			{
				found = 1;
				t = fread(&_snp_chrSNPs[j].locCnt, sizeof(int), 1, fp);
				t = fread(_snp_chrSNPs[j].snpLocs, sizeof(SNPLoc), _snp_chrSNPs[j].locCnt, fp);
				break;
			}
		}
		if (!found)
		{
			t = fread(&count, sizeof(int), 1, fp);
			t = fread(dummy, sizeof(SNPLoc), count, fp);
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
		freeMem(_snp_chrSNPs[i].snpLocs, MAX_SNP_PER_CHR * sizeof(SNPLoc));
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
