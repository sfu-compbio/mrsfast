---
title: Introduction
---

## Introduction

**mrsFAST** is designed to map short reads to reference genome assemblies; in a fast and memory-efficient manner. **mrsFAST** is a cache-oblivous short read mapper that optimizes cache usage to get higher performance.

Currently Supported Features:

- Mismatches, no indels
- Paired-end Mapping
- Discordant Paired-end Mapping (to be used in conjunction with [Variation Hunter](http://variationhunter.sourceforge.net))
- Best Mapping
- Limited Mapping
- SNP-aware Mapping
- Multi-threading

*Legacy Bisulfite Mapping Mode (see version 1.2.6.4): This is a legacy version and is not being supported at the moment. We are aware of a small bug in paired-end mode of mrsFAST-1.2.6.4 that causes mrsFAST to lose some of its mapping locations(<2%). This legacy version is not cache oblivious.*

### Publications

- [Nucleic Acid Research 2014](http://nar.oxfordjournals.org/content/42/W1/W494)
- [Nature Methods 2010](http://www.nature.com/nmeth/journal/v7/n8/full/nmeth0810-576.html)

### Usage

To download and use mrsFAST please use the side links.

### Support

Feel free to send your inquiries to [@isarrafi](http://github.com/isarrafi/) or [@fhach](http://github.com/fhach).

### Bug Reports

To report any bugs or issues please refer to the [issues](https://github.com/sfu-compbio/mrsfast/issues) page.

### Developers

mrsFAST-Ultra is brought to you by:

- [Faraz Hach](http://www.cs.sfu.ca/~fhach/personal/) [@fhach](http://github.com/fhach)
- Iman Sarrafi [@isarrafi](http://github.com/isarrafi/)
- Farhad Hormozdiari
- [Fereydoun Hormozdiari](http://www.gs.washington.edu/~fhormozd/)
- [Can Alkan](http://www.cs.bilkent.edu.tr/~calkan/)

From the [Lab for Computational Biology](http://compbio.cs.sfu.ca) at [Simon Fraser University](http://www.sfu.ca), [Eicher Lab](http://eichlerlab.gs.washington.edu/) at [University of Washington](http://www.washington.edu), and [Alkan Lab](http://www.cs.bilkent.edu.tr/~calkan/compgen/) at [Bilkent University](http://www.bilkent.edu.tr/).

---
