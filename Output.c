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
#include <zlib.h>
#include <string.h>
#include "Common.h"
#include "Output.h"

FILE			*_out_fp;
gzFile			_out_gzfp;



void finalizeGZOutput()
{
	gzclose(_out_gzfp);
}

void finalizeTXOutput()
{
	fclose(_out_fp);
}


void gzOutputQ(SAM map)
{
	gzprintf(_out_gzfp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", 
			map.QNAME, 
			map.FLAG,
			map.RNAME, 
			map.POS,
			map.MAPQ,
			map.CIGAR,
			map.MRNAME,
			map.MPOS,
			map.ISIZE,
			map.SEQ,
			map.QUAL);
	
	int i;

	for ( i = 0; i < map.optSize; i++)
	{
		switch (map.optFields[i].type)
		{
			case 'A':
						gzprintf(_out_gzfp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
						break;
			case 'i':
						gzprintf(_out_gzfp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
						break;
			case 'f':
						gzprintf(_out_gzfp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
						break;
			case 'Z':
			case 'H':
						gzprintf(_out_gzfp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
						break;
		}
	}
	gzprintf(_out_gzfp, "\n");
}

void outputQ(SAM map)
{

	fprintf(_out_fp, "%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s", 
			map.QNAME, 
			map.FLAG,
			map.RNAME, 
			map.POS,
			map.MAPQ,
			map.CIGAR,
			map.MRNAME,
			map.MPOS,
			map.ISIZE,
			map.SEQ,
			map.QUAL);

	
	int i;
	for ( i = 0; i < map.optSize; i++)
	{
		switch (map.optFields[i].type)
		{
			case 'A':
						fprintf(_out_fp, "\t%s:%c:%c", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].cVal);
						break;
			case 'i':
						fprintf(_out_fp, "\t%s:%c:%d", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].iVal);
						break;
			case 'f':
						fprintf(_out_fp, "\t%s:%c:%f", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].fVal);
						break;
			case 'Z':
			case 'H':
						fprintf(_out_fp, "\t%s:%c:%s", map.optFields[i].tag, map.optFields[i].type, map.optFields[i].sVal);
						break;
		}
	}
	
	fprintf(_out_fp, "\n");
}

void outputBufferTxT(char *str, int size)
{
	fwrite(str, 1, size, _out_fp);
}

void outputBufferGZ(char *str, int size)
{
	gzwrite(_out_gzfp, str, size);
}

void outputMetaQ(char* str)
{
	fprintf(_out_fp, "%s\n", str);
}

void gzOutputMetaQ(char* str)
{
	gzprintf(_out_gzfp, "%s\n", str);
}

void noMetaOutput(char *str) {}

int initOutput ( char *fileName, int compressed)
{
	if (compressed)
	{
		char newFileName[strlen(mappingOutputPath)+strlen(fileName)+4];
		sprintf(newFileName, "%s%s.sam.gz", mappingOutputPath, fileName);
		_out_gzfp = fileOpenGZ(newFileName, "w1f");
		if (_out_gzfp == Z_NULL)
		{
			return 0;
		}
	
		finalizeOutput = &finalizeGZOutput;

		output = &gzOutputQ;
		outputMeta =&gzOutputMetaQ;
		outputBuffer = &outputBufferGZ;
	}
	else
	{
	
		char newFileName[strlen(mappingOutputPath)+strlen(fileName)+strlen(".sam")+1];
		if ( !strcmp(mappingOutputPath, "/dev/") && !strcmp(fileName, "null") )
		{
			sprintf(newFileName, "%s%s", mappingOutputPath, fileName);
			nohitDisabled = 1;
		}
		else
		{
			//sprintf(newFileName, "%s%s.sam", mappingOutputPath, fileName);
			sprintf(newFileName, "%s%s", mappingOutputPath, fileName);
		}

		_out_fp = fileOpen(newFileName, "w");
		if (_out_fp == NULL)
		{
			return 0;
		}

		finalizeOutput = &finalizeTXOutput;
		output = &outputQ;
		outputMeta = &outputMetaQ;
		outputBuffer = &outputBufferTxT;
	}
	
	if (noSamHeader)
		outputMeta = &noMetaOutput;

	outputMeta("@HD\tVN:1.4\tSO:unsorted");
	
	return 1;
}


