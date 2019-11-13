## 0. Installation

### From source
To install mrsFAST-Ultra, first you should fetch it from our git repository, or download one of the corresponding compressed zip/tar.gz packages. After downloading, change the current directory to the source directory `mrsfast`, and run `make` in the terminal. The `mrsfast` and `snp_indexer` binaries will be created, which are ready to use.

```bash
git clone https://github.com/sfu-compbio/mrsfast
cd mrsfast
make
```

If you are interested in a particular version, after downloading the git repo, checkout the that version and do `make`.

```bash
git clone https://github.com/sfu-compbio/mrsfast
cd mrsfast
git checkout v3.3.0
make
```
### From release
Alternatively, you can go to [releases page](https://github.com/sfu-compbio/mrsfast/releases) and click on the desired version and then click on download the zip or tar file, switch to directory and run `make`. 

To grab sample data and test `mrsfast`, please download it from our git repository or the compressed package.

```bash
git clone https://github.com/sfu-compbio/mrsfast mrsfast/sample-data -b sample-data
```
Note that the current version of mrsFAST-Ultra requires 64 bit architecture and does not support 32 bit machines.

### From `bioconda`
Install [conda](https://conda.io) and run:

```bash
conda install -c bioconda mrsfast
```

----

## 1. Indexing Reference Genome

In order to map read sequences to a reference genome, mrsFAST-Ultra first needs to creata an index from the genome fasta file. This command will create the file `genome.fa.index`.

```bash
$ ./mrsfast --index genome.fa
```


By default, the indexing window size (length of hash values stored in the index hash table) is 12. This value could be also determined with the `--ws` option. A maximum value of 14 would be logical for window size. Larger values could lead to excessive memory usage in the mapping stage.

```bash
$ ./mrsfast --index genome.fa --ws 14
```

----

## 2. Mapping Options

To perform read mapping, mrsFAST-Ultra should be executed with the `--search` option. By default, mrsFAST-Ultra is a all-mapper tool. This means that it finds and reports all the mappings for each input read. If no option is provided, mrsFAST-Ultra performs single-end mapping. This is an example of running mrsFAST-Ultra for mapping read sequences in the sample `reads.fq` file.

```bash
$ ./mrsfast --search genome.fa --seq reads.fq
```

By default, mrsFAST-Ultra sets the maximum error threshold to 6% so that %94 similarity between reads and genome is guaranteed. The error threshold could also be set using the `-e` option. Setting this value to zero means that only exact matches are desired.

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4
```

To trim input reads to a shorter length, `--crop` option can be used. For example if the input reads are of length 100 and `--crop 40` is used, only the first 40 base pairs of each read would be read and used for mapping.

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4 --crop 40
```

Alternatively, you may crop from tail of the reads by using `--tail-crop`.

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4 --tail-crop 40
```

mrsFAST-Ultra is able to perform mapping on multiple CPU cores in a parallel manner. The number of threads is determined using the `--threads` option. The default value is 1 which runs mrsFAST-Ultra on a single thread. If `--threads` is set to 0, mrsFAST-Ultra will use all the available cores in the system.

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4 --threads 4
```

When huge volumes of input reads are required to be mapped on a machine with limited memory resources, mrsFAST-Ultra is capable of adjusting itself with the specified memory limits. The total memory (in GB) available for running mrsFAST-Ultra should be specified with the `--mem` option. In this mode, mrsFAST-Ultra might perform mapping in several iterations and each time it only loads as many reads as allowed by the memory limit.

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4 --mem 8
```

In the limited mapping mode, mrsFAST-Ultra reports only up to a specified number of mappings for each read. The `-n` option sets the maximum number of mappings per read. Reads with mapping more than this value will not be printed in the output. This option is valid in both single-end and paired-end modes, but could not be used together with best mapping (`--best`) and paired-end discordant mode (`--discordant-vh`).

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4 -n 100
```

In best mapping mode, for each read mrsFAST-Ultra reports a single best location which has the smallest hamming distance among all of its possible mappings. In case of a tie, one of the mapping locations with the smallest hamming distance is reported at random. This option cannot be used in paired-end discordant mode (`--discordant-vh`).

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4 --best
```

Except the cases that are pointed out, obviously any combination of the above options can be used together in any of the single-end and paired-end modes.

```
$ ./mrsfast --search genome.fa --seq reads.fq --crop 60 -e 2 --mem 6 --threads 4 -n 100 -o mappings.sam --disable-nohits
```

----

## 3. Paired-end Mapping

As mentioned above, by default mrsFAST-Ultra will run in single-end mode. The paired-end mapping options is invoked using the `--pe` option. If the reads are in two different files, `--seq1` and `--seq2` should be used to indicate the input files. If the reads are interleaved in a single file, `--seq` is used to indicated the file.

```
$ ./mrsfast --search genome.fa --seq interleaved-reads.fq --pe -e 4
$ ./mrsfast --search genome.fa --seq1 mates1.fq --seq2 mates2.fq --pe -e 4
```

If the distance range between condordant pairs is not specified as above, mrsFAST-Ultra automatically decides about this range according to the mean and standard deviation of distances observed among the mates. The distance allowed between the paired-end reads should be specified with `--min` and `--max`. These values specify the minimum and maximum of the template length (the distance between outer edges of the mapping mates).

```
$ ./mrsfast --search genome.fa --seq1 mates1.fq --seq2 mates2.fq --pe -e 4 --min 100 --max 500
```

Again, any combination of the introduced mapping options could be used in the paired-end mode.

```
$ ./mrsfast --search genome.fa --seq1 mates1.fq --seq2 mates2.fq --pe -e 4 --min 100 --max 500 --threads 4 --best -o mappings.sam
```
mrsFAST-Ultra can report discordant paired-end mappings for structural variation detection using [Variation Hunter](http://variationhunter.sf.net). In this mode the `--min` and `--max` options will define the minimum and maximum inferred size for concordant mappings.

```
$ ./mrsfast --search genome.fa --seq1 mates1.fq --seq2 mates2.fq --discordant-vh -e 4 --min 100 --max 500
```

In paired-end discordant mode, an upper bound can be defined for maximum number of discordant mappings per read. This values is determined by `--max-discordant-cutoff` option. This option is only applicable in paired-end discordant mode.

```
$ ./mrsfast --search genome.fa --seq1 mates1.fq --seq2 mates2.fq --discordant-vh -e 4 \ 
  --min 100 --max 500 --max-discordant-cutoff 300
```

----

## 4. SNP-aware mode
mrsFAST-Ultra is able to do sequence mapping in SNP-aware mode. In this mode mrsFAST-Ultra tolerates the mismatches in known SNP locations provided by dbSNP database (see sample file `dbSNP.vcf`). To run in this mode, first, the snp_indexer binary should be used to create an index from the input dbSNP (vcf) file. The following command reads the sample `dbSNP.vcf` file and creates `snp.index` which is only readable by mrsFAST-Ultra. The current vcf format that is accepted by mrsFAST-Ultra is vcf version 4.

```
$ ./snp_indexer dbSNP.vcf snp.index
```

In the next step, using `--snp` option mrsFAST-Ultra accepts the snp.index file as an input, and therefore ignores all mismatches occurring in the known SNP locations. The following command line executes mrsFAST-Ultra in SNP-aware mode using the index file created in last step.

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4 --snp snp.index
```

To be able to distinguish mismatches occurring as a result of sequencing errors, and those caused by SNPs, mrsFAST-Ultra considers a quality threshold such that a mismatch at a reported SNP location will be ignored only if the corresponding read location has a quality higher than or equal to that quality threshold; otherwise the mismatch would affect the mapping as it is more probably caused by sequencing errors. The value of this threshold is set by the `--snp-qual` option. The default value is 53.

```
$ ./mrsfast --search genome.fa --seq reads.fq -e 4 --snp snp.index --snp-qual 60
```

The SNP-aware mode could be run together with any other combination of options both in single-end and paired-end modes (except `--discordant-vh`).

----

## 5. Input and Output Options

By default, mrsfast-Ultra outputs the mapping results in `output` which is written in standard SAM format. Also in single-end mode, the set of unmapped reads are printed in `output.nohit` file. The name of sam and nohit files can be set by the `-o` option.

```
$ ./mrsfast --search genome.fa --seq reads.fq -o mappings.sam
```

The name of the nohit file can be determined by the `-u` option.

```
$ ./mrsfast --search genome.fa --seq reads.fq -o mappings.sam -u unmapped
```

If the nohit file is not desired as output, it could be omitted by adding `--disable-nohits`.

```
$ ./mrsfast --search genome.fa --seq reads.fq -o mappings.sam --disable-nohits
```

The `--outcomp` option can be used to compress the mrsFAST-Ultra output file in gzip format.

```
$ ./mrsfast --search genome.fa --seq reads.fq -o mappings.sam --outcomp
```

mrsFAST-Ultra can process the input (FASTA/FASTQ) read files in compressed .gz format, using the `--seqcomp` option.

```
$ ./mrsfast --search genome.fa --seqcomp --seq reads.fq.gz -o mappings.sam
```

By default, mrsFAST-Ultra includes a SAM header in the output file. To make sure the SAM header does not appear in the output, `--disable-sam-header` can be used.

```
$ ./mrsfast --search genome.fa --seq reads.fq -o mappings.sam --disable-sam-header
```

----
## 6. mrsFAST-Ultra man page

To view the full list of mrsFAST-Ultra options and their descriptions, please run the following.

```
$ ./mrsfast -h
```

## 7. Citation
If you use mrsfast, please cite its [paper](https://www.nature.com/articles/nmeth0810-576):

> Hach, F. et al. mrsFAST: a cache-oblivious algorithm for short-read mapping. Nature Methods 7, 576–577 (2010)
