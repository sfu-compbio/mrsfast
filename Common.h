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

/*
 * Author         : Faraz Hach
 * Email          : fhach AT cs DOT sfu
 * Last Update    : 2009-12-08
 */


#ifndef __COMMON__
#define __COMMON__

#include <zlib.h>

#define SEQ_MAX_LENGTH		200			// Seq Max Length
#define CONTIG_OVERLAP		200 		// No. of characters overlapped between contings
#define CONTIG_NAME_SIZE	200			// Contig name max size
#define FILE_NAME_LENGTH	400			// Filename Max Length
#define DISCORDANT_CUT_OFF  800

extern unsigned int		CONTIG_SIZE;
extern unsigned int		CONTIG_MAX_SIZE;


extern unsigned char	WINDOW_SIZE				;		// WINDOW SIZE for indexing/searching
extern unsigned short	SEQ_LENGTH;						// Sequence(read) length
extern unsigned short	QUAL_LENGTH;

extern char				*versionNumber;
extern char				*versionNumberF;
extern unsigned char	mrFAST;


extern int				uniqueMode;
extern int				indexingMode;
extern int				searchingMode;
extern int				bisulfiteMode;
extern int				pairedEndMode;
extern int				pairedEndDiscordantMode;
extern int				pairedEndProfilingMode;
extern int				seqCompressed;
extern int				outCompressed;
extern int				cropSize;
extern int				progressRep;
extern char 			*seqFile1;
extern char				*seqFile2;
extern char				*seqUnmapped;
extern char				*mappingOutput;
extern char 			*mappingOutputPath;
extern char				*unmappedOutput;
extern unsigned char	seqFastq;
extern unsigned char	errThreshold;
extern unsigned char	maxHits;	
extern int				minPairEndedDiscordantDistance;
extern int				maxPairEndedDiscordantDistance;
extern int				minPairEndedDistance;
extern int				maxPairEndedDistance;
extern char				fileName[1000][2][FILE_NAME_LENGTH];
extern int				fileCnt;
extern long long		memUsage;

FILE	* fileOpen(char *fileName, char *mode);
gzFile	fileOpenGZ(char *fileName, char *mode);
double	getTime(void);
void	reverseComplete (char *seq, char *rcSeq , int length);
void	* getMem(size_t size);
void	freeMem(void * ptr, size_t size);
double	getMemUsage();
void 	reverse (char *seq, char *rcSeq , int length);
void 	stripPath(char *full, char **path, char **fileName);
#endif
