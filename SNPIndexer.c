#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "Common.h"

#define MAX_LINE_LENGTH 2000
#define NUM_OF_CHRS 25

/**********************************************/
FILE *fileOpen(char *fileName, char *mode)
{
	FILE *fp;
	fp = fopen (fileName, mode);
	if (fp == NULL)
	{
		fprintf(stdout, "Error: Cannot Open the file %s\n", fileName);
		fflush(stdout);
		exit(0);
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
/**************************************/

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "Too few input arguments\nInputs must be:\n\t1. Input vcf (version 4) file name\n\t2. Output file name\n");
		return 0;
	}
	FILE *inFile = 0, *outFile = 0;
	char *inFileName = argv[1];
	char *outFileName = argv[2];
	char line[MAX_LINE_LENGTH], chr[MAX_LINE_LENGTH], ref[MAX_LINE_LENGTH], alt[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
	int i, loc, index;
	unsigned int snpCnt;
	int len[NUM_OF_CHRS];
	SNPLoc **chrSNPs = malloc(NUM_OF_CHRS * sizeof(SNPLoc *));
	
	for (i = 0; i < NUM_OF_CHRS; i++)
	{
		chrSNPs[i] = malloc(MAX_SNP_PER_CHR * sizeof(SNPLoc));
		len[i] = 0;
	}

	inFile = fileOpen(argv[1], "r");
	if (!inFile)
	{
		fprintf(stderr, "ERROR: Could not open VCF file %s\n", argv[1]);
		return 0;
	}

	fprintf(stdout, "Reading VCF file ");
	fflush(stdout);

	snpCnt = 0;
	i = 0;

	while ( fgets(line, MAX_LINE_LENGTH, inFile) )
	{
		if (++i == 2000000)
		{
			fprintf(stdout, ".");
			fflush(stdout);
			i = 0;
		}

		if (line[0] == '#')		// comment line
			continue;
		sscanf(line, "%s%d%s%s%s", chr, &loc, dummy, ref, alt);
		if (strlen(ref) != 1 || strlen(alt) != 1)
			continue;

		index = atoi(chr);
		if (index)		// a number
		{
			if (index > 22)
			{
				fprintf(stderr, "ERROR (line %d): %s is not a valid chromosome for homosapiens. Line Skipped\n", i, chr);
				continue;
			}
			index --;
		}
		else			// X, Y or M
		{
			chr[0] = toupper(chr[0]);
			chr[1] = toupper(chr[1]);
			if (!strcmp(chr, "MT"))
				chr[1] = '\0';
			if (strlen(chr) > 1)
			{
				fprintf(stderr, "ERROR (line %d): %s is not a valid chromosome for homosapiens. Line Skipped\n", i, chr);
				continue;
			}

			switch(chr[0])
			{
				case 'X': index = 22; break;
				case 'Y': index = 23; break;
				case 'M': index = 24; break;
				default : 
					fprintf(stderr, "ERROR (line %d): %s is not a valid chromosome for homosapiens. Line Skipped\n", i, chr);
					continue;
			}
		}

		//chrSNPs[index][++chrSNPs[index][0]] = loc;
		chrSNPs[index][len[index]].loc = loc;
		chrSNPs[index][len[index]].alt = alt[0];
		len[index] ++;
		snpCnt++;
	}

	fprintf(stdout, ".\nReformatting data\n");
	fclose(inFile);

	for (i = 0; i < NUM_OF_CHRS; i++)
	{
		if (len[i])
			qsort(chrSNPs[i], len[i], sizeof(SNPLoc), cmp);
	}

	// -------- write to output file ------- //
	outFile = fileOpen(outFileName, "w");
	
	fprintf(stdout, "Creating output in %s\n", outFileName);
	int n = NUM_OF_CHRS;
	fwrite(&n, sizeof(int), 1, outFile);

	for (i = 0; i < NUM_OF_CHRS; i++)
	{
		fwrite(len + i, sizeof(int), 1, outFile);		// length of list for this chr
		fwrite(chrSNPs[i], sizeof(SNPLoc), len[i], outFile);
	}

	fclose(outFile);
	fprintf(stdout, "%u SNP locations registered successfully\n", snpCnt);

	for (i = 0; i < NUM_OF_CHRS; i++)
		free(chrSNPs[i]);
	free(chrSNPs);

	return 0;
}
