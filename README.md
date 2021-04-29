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

Second, you need a full table of BUSCO hits (https://busco.ezlab.org/). Note that we have tested BUSCO v3 only and the latest BUSCO output format might be incompatible. For example:

```
run_BUSCO.py -i ${genome} -o buscov3_clean_assembly -l basidiomycota_odb9 -m geno -sp coprinus -c4
```

The output should look like this:
```
head full_table_buscov3_clean_assembly.tsv

# BUSCO version is: 3.1.0
# The lineage dataset is: basidiomycota_odb9 (Creation date: 2016-02-13, number of species: 25, number of BUSCOs: 1335)
# To reproduce this run: python /apps/busco/3.1.0/scripts/run_BUSCO.py -i Pt_Clean_Genome.fasta -o buscov3_clean_assembly -l basidiomycota_odb9/ -m genome -c 4 -sp coprinus
#
# Busco id      Status  Contig  Start   End     Score   Length
EOG092R000C     Duplicated      tig00000348     6772572 6790264 2620.1  2401
EOG092R000C     Duplicated      tig00001348     188359  205921  2621.3  2402
EOG092R000I     Duplicated      tig00000129     7090796 7107677 2028.3  1919
EOG092R000I     Duplicated      tig00001352     223668  240549  2028.4  1919
```
You can see here that duplicated BUSCOs will also act like a phasing marker as they have two hits exactly to two contigs.

Third, you need a Hi-C contact map in ginteractions format. We recommend following the HiC-Pro pipeline (https://github.com/nservant/HiC-Pro) and then converting the h5 format to ginteractions with hicexplorer (https://hicexplorer.readthedocs.io/en/latest/content/tools/hicConvertFormat.html).

We also recommend to set MAPQ=30 in the HiC-Pro config file and to obtain a matrix at resolution 20,000 bps.

Your Hi-C contact map output file in ginteractions format should look like this:

```
head HiC_MAPQ30.clean_assembly.20000.matrix.tsv
tig00000001     0       20000   tig00000001     0       20000   27.008251
tig00000001     0       20000   tig00000001     20000   40000   2.989856
tig00000001     0       20000   tig00000001     40000   60000   1.906241
tig00000001     0       20000   tig00000001     180000  200000  0.653958
tig00000001     0       20000   tig00000001     9280000 9300000 0.435236
tig00000001     0       20000   tig00000171     3740000 3760000 0.550744
tig00000001     0       20000   tig00000286     0       20000   0.415017
```

#### Step1: Run NuclearPhaser

Now you are ready to run the phasing pipeline like so:
```
python NuclearPhaser.py Pt_Clean_Genome_GeneMapping.txt full_table_buscov3_clean_assembly.tsv HiC_MAPQ30.clean_assembly.20000.matrix.tsv ${genome} /path_to_output_dir/Genome_Phasing_Output/
```

You should get an output like this:
```
Run minimap2
[M::mm_idx_gen::4.708*1.47] collected minimizers
[M::mm_idx_gen::5.214*1.69] sorted minimizers
[M::main::5.214*1.69] loaded/built the index for 600 target sequence(s)
[M::mm_mapopt_update::5.363*1.67] mid_occ = 343
[M::mm_idx_stat] kmer size: 19; skip: 19; is_hpc: 0; #seq: 600
[M::mm_idx_stat::5.467*1.66] distinct minimizers: 8631579 (16.13% are singletons); average occurrences: 2.975; average spacing: 9.988
[M::worker_pipeline::58.309*3.48] mapped 600 sequences
[M::main] Version: 2.16-r922
[M::main] CMD: minimap2 -k19 -w19 -m200 -DP -r1000 -t4 -o /path_to_output_dir/Genome_Phasing_Output/temp.paf genome.clean.fasta genome.clean.fasta
[M::main] Real time: 58.356 sec; CPU: 202.880 sec; Peak RSS: 3.928 GB
Minimap2 alignments are finished, now scan the paf file.
Done scanning the PAF alignment file.
Now find haplotigs
Mean alignment length is: 16528.201507812657
Find contig pairs that share genes.
Read in contact map
Hi-C contact matrix has resolution: 20000
Done
Construct graph from contig pairs to find well-connected groups
Construct community
Use 26 gene bins.

Contigs in bins: 193 with a total length of 240.48791 MB
Still unphased: 407 contigs with a total length of 15.975574 MB

Construct graph to phase the gene bins with Hi-C data

Construct community to phase gene bins into the haplotypes
----------
----------
Haplotype_0 ['Bin_10a', 'Bin_11b', 'Bin_12a', 'Bin_13a', 'Bin_14b', 'Bin_15a', 'Bin_16a', 'Bin_17a', 'Bin_18b', 'Bin_19b', 'Bin_1a', 'Bin_20a', 'Bin_21b', 'Bin_22a', 'Bin_23a', 'Bin_24a', 'Bin_25b', 'Bin_26a', 'Bin_2b', 'Bin_3b', 'Bin_4a', 'Bin_5a', 'Bin_6a', 'Bin_7b', 'Bin_8b', 'Bin_9a']
120.424197 MB
----------
----------
Haplotype_1 ['Bin_10b', 'Bin_11a', 'Bin_12b', 'Bin_13b', 'Bin_14a', 'Bin_15b', 'Bin_16b', 'Bin_17b', 'Bin_18a', 'Bin_19a', 'Bin_1b', 'Bin_20b', 'Bin_21a', 'Bin_22b', 'Bin_23b', 'Bin_24b', 'Bin_25a', 'Bin_26b', 'Bin_2a', 'Bin_3a', 'Bin_4b', 'Bin_5b', 'Bin_6b', 'Bin_7a', 'Bin_8a', 'Bin_9b']
120.063713 MB
----------
----------------------------------------
Recommended to do a DGenies dot-plot alignment of these two files at this stage to confirm that the gene binning & phasing went well:
/path_to_output_dir/Genome_Phasing_Output/Haplotype_0_genephasing.fasta
/path_to_output_dir/Genome_Phasing_Output/Haplotype_1_genephasing.fasta
----------------------------------------
```

In this example, NuclearPhaser found 26 gene/scaffold bins that contain 193 contigs with a total length of 240.5 MB. Also, in this example the Hi-C phasing signal is clear and the scaffold bins were phased into two haplotypes. It is a good idea to visualize the two scaffold bin haplotypes to see if the dot plot has nice synteny. Here we use Dgenies (http://dgenies.toulouse.inra.fr/) for visualization.

<img src="https://github.com/JanaSperschneider/NuclearPhaser/blob/main/map_Haplotype_1_genephasing_to_Haplotype_0_genephasing.png" width="25%" height="25%">

**Note that if you get more than two haplotypes at this stages, your data does not support that this is a dikaryon. This could be because of extensive phase switches, collapsed regions or the presence of contaminants. We are planning to add support for other ploidys in the future.** 

