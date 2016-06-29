#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "Common.h"

#define MAX_LINE_LENGTH		4000
#define MAX_NUM_OF_CHRS		100000
#define PROGRESS_METER_UNIT	4000000

/**********************************************/
FILE *fileOpen(char *fileName, char *mode)
{
	FILE *fp;
	fp = fopen (fileName, mode);
	if (fp == NULL)
	{
		fprintf(stdout, "Error: Cannot Open file \"%s\"\n", fileName);
		fflush(stdout);
		exit(EXIT_FAILURE);
	}
	return fp;
}

/**********************************************/
int cmp(const void *a, const void *b)
{
	SNPLoc *x = (SNPLoc *) a;
	SNPLoc *y = (SNPLoc *) b;
	return (x->loc - y->loc);
}

/**********************************************/
int findChrIndex(char *chr, ChrSNPs *chrInfo, int chrCount)
{
	int i;
	for (i = 0; i < chrCount; i++)
	{
		if (! strcmp(chr, chrInfo[i].chrName))
			return i;
	}
	return -1;
}
/**********************************************/
void freeMems(ChrSNPs *chrInfo, int chrCount)
{
	int i;
	for (i = 0; i < chrCount; i++)
	{
		free(chrInfo[i].chrName);
		free(chrInfo[i].snpLocs);
	}
	free(chrInfo);
}
/**********************************************/
void fixChromosomeName(char *cname)
{
	// any other unifying standard for chromosome names can be added here
	// this might be required to make sure different names for the same chromosome (like M and MT, or chr1 and 1) are unified
	if (! strcmp(cname, "MT"))
		cname[1] = '\0';
}
/**********************************************/

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Too few input arguments\nInputs must be:\n\t1. Input vcf (v4.0) file name\n\t2. Output SNP index file name\n");
		return 0;
	}
	FILE *inFile = 0, *outFile = 0;
	char *inFileName = argv[1];
	char *outFileName = argv[2];
	char line[MAX_LINE_LENGTH], chr[MAX_LINE_LENGTH], ref[MAX_LINE_LENGTH], alt[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
	int i, j, loc, chrIndex, nameLen;
	unsigned int snpCount = 0;
	int chrCount = 0;
	ChrSNPs *chrInfo = malloc(MAX_NUM_OF_CHRS * sizeof(ChrSNPs));
	
	// read file, count number of chromosomes and their locations
	inFile = fileOpen(inFileName, "r");
	fprintf(stdout, "Pre-processing VCF file ...\n");

	while ( fgets(line, MAX_LINE_LENGTH, inFile) )
	{
		if (line[0] == '#')		// comment line
			continue;

		sscanf(line, "%s%d%s%s%s", chr, &loc, dummy, ref, alt);		// read fields
		if (strlen(ref) != 1 || strlen(alt) != 1)					// only one bp variants
			continue;
		fixChromosomeName(chr);

		snpCount ++;
		chrIndex = findChrIndex(chr, chrInfo, chrCount);
		if (chrIndex == -1)		// new chr name
		{
			chrIndex = chrCount ++;
			// copy chr name
			nameLen = strlen(chr);
			chrInfo[chrIndex].chrName = malloc(nameLen + 1);
			strcpy(chrInfo[chrIndex].chrName, chr);
			// set number of locations
			chrInfo[chrIndex].locCnt = 1;
			chrInfo[chrIndex].snpLocs = NULL;
		}
		else
		{
			chrInfo[chrIndex].locCnt ++;
		}
	}

	fprintf(stdout, "Chromosomes: %d\n", chrCount);
	fprintf(stdout, "Valid SNP locations: %d\n", snpCount);

	// allocate SNPLocs for each chromosome
	for (i = 0; i < chrCount; i++)
	{
		chrInfo[i].snpLocs = malloc(chrInfo[i].locCnt * sizeof(SNPLoc));
		chrInfo[i].locCnt = 0;
		//fprintf(stdout, "%s\n", chrInfo[i].chrName);
	}

	// read file again, fill locations
	rewind(inFile);
	i = 0;
	fprintf(stdout, "Reading SNP locations ");
	fflush(stdout);

	while ( fgets(line, MAX_LINE_LENGTH, inFile) )
	{
		if (++i == PROGRESS_METER_UNIT)
		{
			fprintf(stdout, ".");
			fflush(stdout);
			i = 0;
		}
		if (line[0] == '#')		// comment line
			continue;

		sscanf(line, "%s%d%s%s%s", chr, &loc, dummy, ref, alt);		// read fields
		if (strlen(ref) != 1 || strlen(alt) != 1)					// only one bp variants
			continue;
		fixChromosomeName(chr);
		
		chrIndex = findChrIndex(chr, chrInfo, chrCount);
		assert(chrIndex != -1);
		j = chrInfo[chrIndex].locCnt;
		chrInfo[chrIndex].snpLocs[j].loc = loc;
		chrInfo[chrIndex].snpLocs[j].alt = alt[0];
		chrInfo[chrIndex].locCnt ++;
	}

	fclose(inFile);

	// sort locations for each chromosome
	fprintf(stdout, ".\nReformatting data ...\n");

	for (i = 0; i < chrCount; i++)
	{
		if (chrInfo[i].locCnt > 0)
			qsort(chrInfo[i].snpLocs, chrInfo[i].locCnt, sizeof(SNPLoc), cmp);
	}

	// write to output file
	fprintf(stdout, "Creating output in %s\n", outFileName);
	
	outFile = fileOpen(outFileName, "w");	
	fwrite(&chrCount, sizeof(int), 1, outFile);

	for (i = 0; i < chrCount; i++)
	{
		nameLen = strlen(chrInfo[i].chrName);
		fwrite(&nameLen, sizeof(int), 1, outFile);						// chr name length
		fwrite(chrInfo[i].chrName, sizeof(char), nameLen, outFile);		// chr name
		fwrite(&chrInfo[i].locCnt, sizeof(int), 1, outFile);			// num of locs
		fwrite(chrInfo[i].snpLocs, sizeof(SNPLoc), chrInfo[i].locCnt, outFile);		// all snp locations
	}

	fclose(outFile);
	fprintf(stdout, "%u SNP locations registered successfully\n", snpCount);

	freeMems(chrInfo, chrCount);
	return 0;
}
