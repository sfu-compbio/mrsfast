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
#include "SnipReader.h"

int				MAX_SNIP_CNT		= 24000;	// 1/1000 * length of chr01
CompressedSeq	*_snp_snipMap		= NULL;
int				_snp_snipMapLength	= 0;
ChrSnips		*_snp_chrSnips		= NULL;
int				_snp_chrCnt			= 0;
int				_snp_currentChr		= 0;
int				_snp_currentLoc	= 0;

/**********************************************/
int cmp(const void *a, const void *b)
{
	int *x = (int *)a;
	int *y = (int *)b;
	return (*x - *y);
}
/**********************************************/
void initLoadingSnips(char *fileName)
{
	int maxLineLength = CONTIG_NAME_SIZE + 20;
	int i, loc;
	char cname[CONTIG_NAME_SIZE];	// chromosome name, read from file
	char line[maxLineLength];
	char **chrNames;

	_snp_chrCnt = getChrCnt();
	chrNames = getChrNames();
	_snp_chrSnips = getMem(_snp_chrCnt * sizeof(ChrSnips));
	
	for (i = 0; i < _snp_chrCnt; i++)
	{
		_snp_chrSnips[i].chrName = chrNames[i];
		_snp_chrSnips[i].locCnt = 0;
		_snp_chrSnips[i].locs = getMem(MAX_SNIP_CNT * sizeof(int));
	}

	_snp_snipMapLength = (calculateCompressedLen(CONTIG_MAX_SIZE)+1) * sizeof(CompressedSeq);
	_snp_snipMap = getMem(_snp_snipMapLength);

	FILE *fp = fopen(fileName, "rt");
	while ( fgets(line, maxLineLength, fp) )
	{
		sscanf(line, "%s %d", cname, &loc);
		// TODO: speed up lookup by sorting chrNames --> binary search
		for (i = 0; i < _snp_chrCnt; i++)
		{
			if (!strcmp(cname, _snp_chrSnips[i].chrName))
			{
				_snp_chrSnips[i].locs[_snp_chrSnips[i].locCnt++] = loc;
				break;
			}
		}
	}

	for (i = 0; i < _snp_chrCnt; i++)
		qsort(_snp_chrSnips[i].locs, _snp_chrSnips[i].locCnt, sizeof(int), cmp);

	fclose(fp);
}
/**********************************************/
void finalizeSnips()
{
	int i;
	for (i = 0; i < _snp_chrCnt; i++)
		freeMem(_snp_chrSnips[i].locs, MAX_SNIP_CNT * sizeof(int));
	freeMem(_snp_chrSnips, _snp_chrCnt * sizeof(ChrSnips));
	freeMem(_snp_snipMap, _snp_snipMapLength);
}
/**********************************************/
CompressedSeq *loadSnipMap(char *chrName, int contigStartIndex, int contigLength)
{
	memset(_snp_snipMap, -1, calculateCompressedLen(contigLength) * sizeof(CompressedSeq));
	//memset(_snp_snipMap, -1, _snp_snipMapLength);
	int contigEnd = contigStartIndex + contigLength;
	int loc, offset;
	CompressedSeq *snp, mask;

	if ( strcmp(chrName, _snp_chrSnips[_snp_currentChr].chrName) )		// new chr
	{
		_snp_currentChr++;
		_snp_currentLoc = 0;
	}

	if (_snp_chrSnips[_snp_currentChr].locCnt)
	{
		int i = _snp_currentChr;		// just to make the code more readable
		int pos = _snp_currentLoc;

		while ( pos < _snp_chrSnips[i].locCnt && _snp_chrSnips[i].locs[pos] < contigStartIndex )	// this should never happen!
			pos ++;

		while ( pos < _snp_chrSnips[i].locCnt && _snp_chrSnips[i].locs[pos] < contigEnd )
		{
			loc = _snp_chrSnips[i].locs[pos] - contigStartIndex - 1;
			offset = loc % 21;
			mask = 0x7000000000000000;
			mask = ~(mask >> offset*3);
			
			snp = _snp_snipMap + (loc/21);
			*snp &= mask;

			pos ++;
		}

		_snp_currentChr = i;
		_snp_currentLoc = pos;
	}
	return _snp_snipMap;
}
