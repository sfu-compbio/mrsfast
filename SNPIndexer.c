#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define FILE_NAME_LENGTH 400
#define MAX_LINE_LENGTH 2000
#define MAX_SNP_PER_CHR 4030000
#define NUM_OF_CHRS 25

int cmp(const void *a, const void *b)
{
	int *x = (int *)a;
	int *y = (int *)b;
	return (*x - *y);
}
/**************************************/

void stripPath(char *full, char **path, char **fileName)
{
	int i;
	int loc = -1;

	for (i=strlen(full)-1; i>=0; i--)
	{
		if (full[i]=='/')
		{
			loc = i;
			break;
		}
	}

	if (loc != -1)
	{
		sprintf(*fileName, "%s%c", (full+loc+1), '\0');
		full[loc+1]='\0';
		sprintf(*path,"%s%c", full, '\0');
	}
	else
	{
		sprintf(*fileName, "%s%c", full, '\0');
		sprintf(*path,"%c", '\0');
	}
}
/**************************************/

int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		fprintf(stderr, "ERROR: Too few input arguments\nInputs must be:\n\t1. Input vcf file name\n\t2. Output file name\n");
		return 0;
	}
	FILE *inFile = 0, *outFile = 0;
	char *filePath = malloc(FILE_NAME_LENGTH);
	char *inFileName = malloc(FILE_NAME_LENGTH);
	char *outFileName = argv[2];//malloc(FILE_NAME_LENGTH);
	char line[MAX_LINE_LENGTH], chr[MAX_LINE_LENGTH], ref[MAX_LINE_LENGTH], alt[MAX_LINE_LENGTH], dummy[MAX_LINE_LENGTH];
	int i, loc, index;
	unsigned int snpCnt;
	unsigned int **chrSNPs = malloc(NUM_OF_CHRS * sizeof(char *));
	
	for (i = 0; i < NUM_OF_CHRS; i++)
	{
		chrSNPs[i] = malloc(MAX_SNP_PER_CHR * sizeof(int));
		chrSNPs[i][0] = 0;		// num of locs
	}
	
	stripPath(argv[1], &filePath, &inFileName);
	//sprintf(outFileName, "%s%s", filePath, "out.snp");
	
	inFile = fopen(argv[1], "r");
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

		chrSNPs[index][++chrSNPs[index][0]] = loc;
		snpCnt++;
	}
	fprintf(stdout, ".\nReformatting data\n");
	fclose(inFile);

	for (i = 0; i < NUM_OF_CHRS; i++)
	{
		if (chrSNPs[i][0])
			qsort(chrSNPs[i]+1, chrSNPs[i][0], sizeof(int), cmp);
	}
	// -------- write to output file ------- //
	outFile = fopen(outFileName, "w");
	
	fprintf(stdout, "Creating output in %s\n", outFileName);
	int n = NUM_OF_CHRS;
	fwrite(&n, sizeof(int), 1, outFile);

	for (i = 0; i < NUM_OF_CHRS; i++)
		fwrite(chrSNPs[i], sizeof(int), chrSNPs[i][0]+1, outFile);

	fclose(outFile);
	
	for (i = 0; i < NUM_OF_CHRS; i++)
		free(chrSNPs[i]);
	free(chrSNPs);
	free(filePath);
	free(inFileName);
	free(outFileName);
	fprintf(stdout, "%lld SNP locations registered successfully\n", snpCnt);
}
