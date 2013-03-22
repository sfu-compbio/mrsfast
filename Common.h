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

#ifndef __COMMON__
#define __COMMON__

#if __GNUC__ >= 4 && __GNUC_MINOR__ >=4 && SSE4==1 
#define MRSFAST_SSE4
#endif

#if SSE4==1 && (__GNUC__ < 4 || __GNUC_MINOR < 4)
#define GCC_UPDATE
#endif

#include <zlib.h>
#include <stdint.h>

#define SEQ_MAX_LENGTH		500			// Seq Max Length
#define CMP_SEQ_MAX_LENGTH	10			// Compressed Seq Max Length
#define CONTIG_OVERLAP		1050 		// No. of characters overlapped between contings  --  equals 50 blocks of length 21
#define CONTIG_NAME_SIZE	200			// Contig name max size
#define FILE_NAME_LENGTH	400			// Filename Max Length
#define DISCORDANT_CUT_OFF  800


typedef uint64_t CompressedSeq;

extern unsigned int		CONTIG_SIZE;
extern unsigned int		CONTIG_MAX_SIZE;
extern unsigned int		THREAD_COUNT;
extern double			MAX_MEMORY;
extern int				THREAD_ID[255];

extern unsigned char	WINDOW_SIZE;					// WINDOW SIZE for indexing/searching
extern unsigned short	SEQ_LENGTH;						// Sequence(read) length
extern unsigned short	QUAL_LENGTH;
extern unsigned short	CMP_SEQ_LENGTH;

extern int				indexingMode;
extern int				searchingMode;
extern int				pairedEndMode;
extern int				pairedEndDiscordantMode;
extern int				pairedEndProfilingMode;
extern int				bestMappingMode;
extern int				SNPMode;
extern int				seqCompressed;
extern int				outCompressed;
extern int				cropSize;
extern int				progressRep;
extern int				nohitDisabled;
extern char 			*seqFile1;
extern char				*seqFile2;
extern char				*seqUnmapped;
extern char				*mappingOutput;
extern char 			*mappingOutputPath;
extern char				*unmappedOutput;
extern unsigned char	seqFastq;
extern int				errThreshold;
extern short			maxHits;	
extern int				minPairEndedDiscordantDistance;
extern int				maxPairEndedDiscordantDistance;
extern int				minPairEndedDistance;
extern int				maxPairEndedDistance;
extern char				fileName[3][FILE_NAME_LENGTH];
extern int				fileCnt;
extern long long		memUsage;
extern char				*alphabet;
extern char				checkSumLength;

#pragma pack(push, 1)
typedef struct
{
	uint8_t  checksum; 
	uint32_t info;				// ReadIndex => seqInfo | GenomeIndex ==> Loc
} GeneralIndex;
#pragma pack(pop)
typedef struct
{
	int hv;
	GeneralIndex *list;
} ReadIndexTable;

FILE	* fileOpen(char *fileName, char *mode);
gzFile	fileOpenGZ(char *fileName, char *mode);
double	getTime(void);
inline void	reverseComplete (char *seq, char *rcSeq , int length);
void	* getMem(size_t size);
void	freeMem(void * ptr, size_t size);
double	getMemUsage();
inline void 	reverse (char *seq, char *rcSeq , int length);
void 	stripPath(char *full, char **path, char **fileName);
void compressSequence(char *seq, int seqLen, CompressedSeq *cseq);
inline int 	calculateCompressedLen(int normalLen);
int	hashVal(char *seq);
int	checkSumVal(char *seq);
#endif
