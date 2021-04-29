## NuclearPhaser: Phasing of dikaryotic genome assemblies with Hi-C data

  * [Background](#background)
  * [Overview of the method](#overview-of-the-method)
  * [Prerequisites for running NuclearPhaser](#prerequisites-for-running-nuclearphaser)
  * [High-level overview for running NuclearPhaser](#high-level-overview-for-running-nuclearphaser)
  * [Detailed instructions for running NuclearPhaser](#detailed-instructions-for-running-nuclearphaser)
  
### Background
Most animals and plants have more than one set of chromosomes and package these haplotypes into a single nucleus, usually most often diploid, within each cell. In contrast, many fungal species carry multiple haploid nuclei per cell. Rust fungi are such species with two nuclei (karyons) that contain a full set of haploid chromosomes each. This dikaryotic state has advantages for haplotype phase separation using Hi-C chromatin contact information as the two haplotypes are physically separated. It also means that, unlike in diploids, Hi-C chromatin contacts between haplotypes are false positive signals. 

### Overview of the method
NuclearPhaser is a method for phasing of dikaryotic genomes into the two haplotypes using Hi-C contact graphs. This is an overview of the phasing pipeline for dikaryons.

<img src="https://github.com/JanaSperschneider/NuclearPhaser/blob/main/PhasingGeneBins_v2.png" width="50%" height="50%">

### Prerequisites for running NuclearPhaser

* A clean genome assembly that has few collapsed regions and is cleaned from contaminant contigs

The following software needs to be installed:
* Python3 
* Networkx (https://networkx.org/) and community dectection (https://github.com/taynaud/python-louvain or https://python-louvain.readthedocs.io/en/latest/)
* BUSCO for generating a full table of BUSCO hits (https://busco.ezlab.org/)
* BioKanga for generating a table of gene hits (https://github.com/csiro-crop-informatics/biokanga)
* HiC-Pro for generating a Hi-C contact map (https://github.com/nservant/HiC-Pro)
* Minimap2 for aligning haplotigs (https://github.com/lh3/minimap2)

### High-level overview for running NuclearPhaser

#### Step1: Construct a highly confident subset of the two haplotypes that are expected to reside in separate nuclei

A gene binning step is used to find sets of homologous contigs which represent the two haplotypes. For this, genes that map exactly twice to the unphased assembly are used as phasing markers to assign homologous contigs into diploid scaffold bins *Bin_1,...,Bin_n*. Scaffold bins are constructed with a graph network approach where nodes are contigs and edges are the number of shared genes per Mb. Each strongly connected community in the graph is a diploid scaffold bin *Bin_x* and contains two subsets *Bin_xa* and *Bin_xb*. Thus, a scaffold bin is part of a chromosome where the two subsets represent the haplotypes. 

#### Step2: Separate the binned contigs into two haplotype sets representing their nuclear origin  

We use a graph based on Hi-C links between the scaffold bins, ignoring Hi-C links within scaffold bins for preliminary phasing. A graph network approach **should** return the two expected communities that represent a high proportion of the phased haplotypes, but might still include phase switches. Note that if you get more than two haplotypes in this step the pipeline will stop. It is likely you have either an assembly with too many phase switches or you have contaminant contigs in your assembly (e.g. plant or bacterial).

#### Step3: Fix phase switches in the two haplotype sets

This step requires some manual work at the moment. For each contig in the scaffold bins, we visualized the proportion of Hi-C contacts to haploypes A and B for each scaffold bin. As an example, see below contig tig00000828 from a HiCanu assembly and its associated haplotig alignments. Contig tig00000828 appears to switch phase at ~1.5-3.7 Mb, which overlaps with the corresponding haplotig alignment start and end points.

<img src="https://github.com/JanaSperschneider/NuclearPhaser/blob/main/tig00000828.png" width="50%" height="50%">

#### Step4: After fixing phase switches, run the pipeline again with the updated genome
After correction of phase switches, the input files (Gene mapping, BUSCO table & Hi-C contact map) need to be re-generated with the genome that has the phase switches corrected.

#### Step5: Obtain phased haplotypes
At the very end, NuclearPhaser will return the two haplotype sets in FASTA format.

### Detailed instructions for running NuclearPhaser

#### Step0: Generate all required input files
NuclearPhaser needs three input files and the clean genome assembly FASTA file. 

First, generate the gene hit table with biokanga blitz (https://github.com/csiro-crop-informatics/biokanga). You need a set of transcripts/genes for your species. This can come from previously published genome annotations, or it can be a gene set from a closely related species. You want to find genes that are highly conserved and occur exactly twice, think of it as housekeeping genes or 'haplotype phasing markers'. For example:

```
genome="Pt_Clean_Genome.fasta"
genes="Puctr1_GeneCatalog_transcripts_20131203.nt.fasta"

biokanga index --threads=4 -i ${genome} -o biokanga_index -r gene_mapping
biokanga blitz --sensitivity=2 --mismatchscore=1 --threads=4 -o Pt_Clean_Genome_GeneMapping.txt --in=${genes} --sfx=biokanga_index
```

The output should look like this:
```
head Pt_Clean_Genome_GeneMapping.txt

psLayout version 3
Generated by biokanga blitz, Version 4.4.2
match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts     tStarts
        match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
---------------------------------------------------------------------------------------------------------------------------------------------------------------
357     0       0       0       0       0       3       630     -       jgi|Puctr1|2|PTTG_25073T0       357     0       357     tig00000821     1756792 38612   39599   4       43,130,127,57,     0,43,173,300,   38612,38733,38967,39542,
356     1       0       0       0       0       3       630     -       jgi|Puctr1|2|PTTG_25073T0       357     0       357     tig00000656     4891113 3219391 3220378 4       43,130,127,57,     0,43,173,300,   3219391,3219512,3219746,3220321,
420     0       0       0       0       0       0       0       +       jgi|Puctr1|3|PTTG_04086T0       420     0       420     tig00000821     1756792 45261   45681   1       420,       0,      45261,
419     1       0       0       0       0       0       0       +       jgi|Puctr1|3|PTTG_04086T0       420     0       420     tig00000656     4891113 3226034 3226454 1       420,       0,      3226034,
317     46      0       0       1       5       1       4       +       jgi|Puctr1|3|PTTG_04086T0       420     0       368     tig00000533     5698764 1765183 1765550 2       271,92,    0,276,  1765183,1765458,
```

You can see here that gene PTTG_25073T0 will act like a phasing marker as it has two hits exactly to contigs tig00000821 and tig00000656.

Second, you need a full table of BUSCO hits (https://busco.ezlab.org/). Note that we have tested BUSCO v3 only and the latest BUSCO output format might be incompatible.
