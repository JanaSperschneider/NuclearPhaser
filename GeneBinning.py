import sys
import os
import networkx as nx
import community
import subprocess
import errno
import uuid
import shutil
import tempfile
import operator
from collections import defaultdict
#--------------------------------------
#--------------------------------------
THRESHOLD_GENE_DENSITY = 30.0
THRESHOLD_SHARED_GENES = 2.0
BUSCO_WEIGHT = 1.0
#--------------------------------------
#--------------------------------------
def genes_shared(BUSCO_PAIRS, CONTIG_PAIRS, LENGTHS):

  # Now merge the gene mapping information with the BUSCO information
  # to get the total numbers of genes shared per contig pair
  SHARED_GENES = {}

  # This is run first to capture contig pairs that only occur in the BUSCO set
  for pair in BUSCO_PAIRS:
    gene_density = 0.0
    if pair in CONTIG_PAIRS:
      total_genes_shared = BUSCO_PAIRS[pair] + CONTIG_PAIRS[pair]
    else:
      total_genes_shared = BUSCO_PAIRS[pair]

    if total_genes_shared > THRESHOLD_SHARED_GENES:
      # Shared genes per MB normalized by the shorter contig length
      shorter_contig_length = min(LENGTHS[pair[0]], LENGTHS[pair[1]])
      gene_density = total_genes_shared/float(shorter_contig_length) * 1000000.0

      if gene_density > THRESHOLD_GENE_DENSITY:
        SHARED_GENES[(pair[0], pair[1])] = gene_density
        SHARED_GENES[(pair[1], pair[0])] = gene_density

  for pair in CONTIG_PAIRS:
    gene_density = 0.0
    if pair in BUSCO_PAIRS:
      total_genes_shared = BUSCO_PAIRS[pair] + CONTIG_PAIRS[pair]
    else:
      total_genes_shared = CONTIG_PAIRS[pair]

    if total_genes_shared > THRESHOLD_SHARED_GENES:
      # Shared genes per MB normalized by the shorter contig length
      shorter_contig_length = min(LENGTHS[pair[0]], LENGTHS[pair[1]])
      gene_density = total_genes_shared/float(shorter_contig_length) * 1000000.0

      if gene_density > THRESHOLD_GENE_DENSITY:
        SHARED_GENES[(pair[0], pair[1])] = gene_density
        SHARED_GENES[(pair[1], pair[0])] = gene_density

  return SHARED_GENES
#--------------------------------------
def find_haplotigs(PATH_TO_GENOME_FASTA, OUTPUT_DIRECTORY_PATH, LENGTHS):
  #--------------------------------------
  print("Run minimap2")
  ParamList = ['minimap2', '-k19', '-w19', '-m200', '-DP', '-r1000', '-t4', PATH_TO_GENOME_FASTA, PATH_TO_GENOME_FASTA, '-o', OUTPUT_DIRECTORY_PATH + '/' + 'temp.paf']

  try:
      Process = subprocess.Popen(ParamList, shell=False)
      sts = Process.wait()
      cstdout, cstderr = Process.communicate()

      if Process.returncode:
          raise Exception("Calling minimap2 returned %s"%Process.returncode)
      if cstdout:
          pass
      elif cstderr:
          sys.exit()
  except:
      e = sys.exc_info()[1]
      print("Error calling minimap2: %s" % e)
      sys.exit(1)
  #--------------------------------------
  # Now go through the minimap2 output files
  #--------------------------------------
  print('Minimap2 alignments are finished, now scan the paf file.')

  HAPLOTIGS = defaultdict(list)
  ALL_HITS = defaultdict(list)

  size_distribution = []

  # Record all hits for each contig
  with open(OUTPUT_DIRECTORY_PATH + '/' + 'temp.paf') as infile:
    for line in infile:
      query_contig = line.split('\t')[0]
      hit_contig = line.split('\t')[5]
      query_start = int(line.split('\t')[2])
      query_end = int(line.split('\t')[3])
      hit_start = int(line.split('\t')[7])
      hit_end = int(line.split('\t')[8])
      alignment_block_length = int(line.split('\t')[10])

      if alignment_block_length > 5000:

        ALL_HITS[query_contig].append((hit_contig, query_start, query_end, hit_start, hit_end, alignment_block_length))
        size_distribution.append(alignment_block_length)

  print('Done scanning the PAF alignment file.')
  print("Now find haplotigs")
  mean = (sum(size_distribution)/len(size_distribution))
  print("Mean alignment length is:", mean)

  for query_contig, length in LENGTHS.items():
    HITS = defaultdict(list)

    if query_contig in ALL_HITS:
      # Go through all the hits for this contig and record the alignments
      for (hit_contig, query_start, query_end, hit_start, hit_end, alignment_block_length) in ALL_HITS[query_contig]:
        if alignment_block_length > mean:
          HITS[hit_contig].append((query_start, query_end))

      # Now merge the alignment coordinates allowing a gap
      for hit_contig, hits in HITS.items():

        merged_hits = merge_intervals(hits, mean)
        bases_aligned = sum([end-start for (start, end) in merged_hits])

        HAPLOTIGS[query_contig].append((bases_aligned, hit_contig, merged_hits))

  return HAPLOTIGS
#--------------------------------------
def write_shared_genes_to_file(OUTPUT_DIRECTORY_PATH, SHARED_GENES, BUSCO_PAIRS, CONTIG_PAIRS, LENGTHS):

  # Write the shared gene content to a file
  f = open(OUTPUT_DIRECTORY_PATH + '/shared_genes.txt', 'w')
  header = 'Contig_1\tContig_2\tLength contig_1 (bps)\tLength contig_2 (bps)\tNumber of shared genes\tShared genes per MB\n'
  f.writelines(header)

  for (contig1, contig2), gene_density in SHARED_GENES.items():
    total_genes_shared = 0
    if (contig1, contig2) in BUSCO_PAIRS:
      total_genes_shared += BUSCO_PAIRS[contig1, contig2]
    if (contig1, contig2) in CONTIG_PAIRS:
      total_genes_shared += CONTIG_PAIRS[contig1, contig2]

    f.writelines(contig1 + '\t' + contig2 + '\t' + str(LENGTHS[contig1]) + '\t' + str(LENGTHS[contig2]) + '\t' + str(total_genes_shared) + '\t' + str(gene_density) + '\n')

  f.close()

  return
#--------------------------------------
def merge_intervals(intervals):

  sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])

  merged = []

  for higher in sorted_by_lower_bound:
      if not merged:
          merged.append(higher)
      else:
          lower = merged[-1]
          # test for intersection between lower and higher:
          # we know via sorting that lower[0] <= higher[0]
          if higher[0] <= lower[1] + 1:
              upper_bound = max(lower[1], higher[1])
              merged[-1] = (lower[0], upper_bound)  # replace by merged interval

          else:
              merged.append(higher)
  return merged
#--------------------------------------
def merge_intervals(intervals, GAP):

  sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])

  merged = []

  for higher in sorted_by_lower_bound:
      if not merged:
          merged.append(higher)
      else:
          lower = merged[-1]
          # test for intersection between lower and higher:
          # we know via sorting that lower[0] <= higher[0]
          if higher[0] <= lower[1] + 1 + GAP:
              upper_bound = max(lower[1], higher[1])
              merged[-1] = (lower[0], upper_bound)  # replace by merged interval

          else:
              merged.append(higher)
  return merged
#--------------------------------------
def make_fasta_files(contigs, output_name, identifiers, sequences):
  # Now make combined FASTA files for contigs
  f = open(output_name, 'w')

  for contig in contigs:
    index = identifiers.index(contig)
    f.writelines('>' + identifiers[index] + '\n')
    f.writelines(sequences[index] + '\n')

  f.close()

  return
#--------------------------------------
def total_length(contigs, sequence_lengths):
  # Takes a list of contigs as input, returns total length

  total_length = sum([sequence_lengths[contig] for contig in contigs])

  return total_length
#--------------------------------------
def SimpleFastaParser(handle):
    for line in handle:
        if line[0] == ">":
            title = line[1:].rstrip()
            break

    lines = []
    for line in handle:
        if line[0] == ">":
            yield title, "".join(lines).replace(" ", "").replace("\r", "")
            lines = []
            title = line[1:].rstrip()
            continue
        lines.append(line.rstrip())

    yield title, "".join(lines).replace(" ", "").replace("\r", "")
#--------------------------------------
def read_gene_mapping(GENE_MAPPING):
  """ Function: read_gene_mapping()
      Purpose:  This function reads in a BioKanga blitz gene mapping results table and records contigs that share genes == likely haplotigs.

                psLayout version 3
                Generated by biokanga blitz, Version 4.4.2
                match   mis-    rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts      tStarts
                        match   match           count   bases   count   bases           name            size    start   end     name            size    start   end     count
                ---------------------------------------------------------------------------------------------------------------------------------------------------------------
                600     3       0       0       0       0       0       0       +       evm.TU.hcontig_014_192.4        1154    543     1146    tig00000148_arrow_arrow_pilon_pilon_1_2600000   2600000     371285  371888  1       603,    543,    371285,
                600     3       0       0       0       0       0       0       -       evm.TU.hcontig_014_192.4        1154    543     1146    tig00000338_arrow_arrow_pilon_pilon     204015  156237      156840  1       603,    8,      156237,
                551     0       0       0       0       0       0       0       +       evm.TU.hcontig_014_192.4        1154    0       551     tig00000148_arrow_arrow_pilon_pilon_1_2600000   2600000     370743  371294  1       551,    0,      370743,
                1401    2       0       0       1       13      2       17      +       evm.TU.hcontig_014_192.6        1419    0       1416    tig00000148_arrow_arrow_pilon_pilon_1_2600000   2600000     376787  378207  3       1119,79,205,    0,1119,1211,    376787,377907,378002,
                551     0       0       0       0       0       0       0       -       evm.TU.hcontig_014_192.4        1154    0       551     tig00000338_arrow_arrow_pilon_pilon     204015  156831      157382  1       551,    603,    156831,

      Input:    Biokanga blitz results table.

      Return:   A dictionary that holds contigs pairs and the number of genes they share.
  """
  SHARED_GENES = {}

  f = open(GENE_MAPPING, 'r')
  content = f.readlines()
  f.close()

  for line in content[5:]:
    gene = line.split('\t')[9]
    contig = line.split('\t')[13]
    match_score = int(line.split('\t')[0])

    if gene in SHARED_GENES:
      SHARED_GENES[gene] = list(set(SHARED_GENES[gene] + [(match_score, contig)]))
    else:
      SHARED_GENES[gene] = [(match_score, contig)]

  CONTIG_PAIRS = {}

  for gene, contig_list in SHARED_GENES.items():
    # The top matches will not appear first in the list by default, sorted by highest match score first
    contig_list = sorted(contig_list, reverse=True)

    # Collect the number of genes that two contigs share, use only genes that occur exactly twice
    if len(contig_list) == 2.0: # Genes occur exactly twice
      x = contig_list[0][1]
      y = contig_list[1][1]
      match_score_x = contig_list[0][0]
      match_score_y = contig_list[1][0]

      # Record pairs in both possible combinations
      pair = (x,y)
      if pair in CONTIG_PAIRS:
        CONTIG_PAIRS[pair] = CONTIG_PAIRS[pair] + 1
      else:
        CONTIG_PAIRS[pair] = 1

      pair = (y,x)
      if pair in CONTIG_PAIRS:
        CONTIG_PAIRS[pair] = CONTIG_PAIRS[pair] + 1
      else:
        CONTIG_PAIRS[pair] = 1

  return CONTIG_PAIRS
#--------------------------------------
def genes_shared_per_MB(contig1, contig2, GENE_PAIR_DIC):

  gene_density = 0.0

  shorter_contig_length = min(LENGTHS[contig1], LENGTHS[contig2])

  if (contig1, contig2) in GENE_PAIR_DIC:
    gene_density = GENE_PAIR_DIC[(contig1, contig2)]/float(shorter_contig_length) * 1000000.0

  return gene_density
#--------------------------------------
def gene_bins_to_file(OUTPUT_DIRECTORY_PATH, GENE_BINS, LENGTHS):

  f = open(OUTPUT_DIRECTORY_PATH + '/gene_binning.txt', 'w')
  header = 'Gene bin identifier\tTotal length (MB)\tContigs in haplotype 1\tTotal length haplotype 1 (MB)\tContigs in haplotype 2\tTotal length haplotype 2 (MB)\n'
  f.writelines(header)
  for key, values in GENE_BINS.items():

    haplotype1 = values[0]
    haplotype2 = values[1]

    length_of_partition = total_length(haplotype1 + haplotype2, LENGTHS)

    line = key + '\t' + str(length_of_partition/1000000.0) + '\t'
    line += ",".join(haplotype1) + '\t' + str(total_length(haplotype1, LENGTHS)/1000000.0) + '\t'
    line += ",".join(haplotype2) + '\t' + str(total_length(haplotype2, LENGTHS)/1000000.0) + '\n'
    f.writelines(line)
  f.close()

  return
#--------------------------------------
def partition_to_genebins(partition, LENGTHS, SHARED_GENES):

  GENE_BINS = {}

  count = 0
  for com in set(partition.values()) :
      list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]
      length_of_partition = total_length(list_nodes, LENGTHS)
      list_1, list_2 = [], []

      # Use partitions longer than 1 MB for good Hi-C signal later
      if length_of_partition/1000000.0 > 1.0 and len(list_nodes) >= 2.0:

        count += 1

        # If a group contains only two contigs, phasing is trivial
        if len(list_nodes) == 2.0:

          if LENGTHS[list_nodes[0]] > LENGTHS[list_nodes[1]]:
            list_1.append(list_nodes[0])
            list_2.append(list_nodes[1])
          else:
            list_1.append(list_nodes[1])
            list_2.append(list_nodes[0])
          GENE_BINS['Bin_' + str(count)] = [list_1, list_2]

        # If a group contains > two contigs, phasing is done with a graph community approach
        else:
          GENE_BINS = bin_to_haplotypes(list_nodes, SHARED_GENES, count, LENGTHS, GENE_BINS)

  return GENE_BINS
#--------------------------------------
def bin_to_haplotypes(list_nodes, SHARED_GENES, count, LENGTHS, GENE_BINS):

  G = nx.Graph()
  # Now split this contig set into the two haplotypes
  for node1 in sorted(list_nodes):
    for node2 in sorted(list_nodes):
      if node1 != node2:
        if (node1, node2) in SHARED_GENES:
          gene_density = SHARED_GENES[(node1, node2)]
        else:
          gene_density = 0.0

        # If contigs share genes, they are likely haplotigs.
        # Haplotigs should not be connected in phasing, so give them weight of zero
        if gene_density > 0.0:
          G.add_edge(node1, node2, weight=0.0)

        # If contigs do not share genes, they are unlikely haplotigs.
        # Non-haplotigs should be connected in phasing, so give them connected weight of one
        else:
          G.add_edge(node1, node2, weight=1.0)

  partition_contigs = community.best_partition(G)

  if len(set(partition_contigs.values())) == 2.0:
    list_1, list_2 = [], []
    # The set was neatly split into the two haplotypes
    contig_partition = [nodes for nodes in partition_contigs.keys() if partition_contigs[nodes] == 0]
    length_of_contig_partition = total_length(contig_partition, LENGTHS)
    list_1 += contig_partition
    contig_partition = [nodes for nodes in partition_contigs.keys() if partition_contigs[nodes] == 1]
    length_of_contig_partition = total_length(contig_partition, LENGTHS)
    list_2 += contig_partition

    GENE_BINS['Bin_' + str(count)] = [list_1, list_2]

  else:
    print('Partition was split into more than two haplotypes, maybe there are chimerics or gene duplications?')
    print(partition_contigs)
    print('This partition will not be included.')

  return GENE_BINS
#--------------------------------------
def read_in_BUSCOs(BUSCO_TABLE):
  """ Function: read_in_BUSCOs()
      Purpose:  This function reads in a BUSCO results table and records contigs that share BUSCOs == likely haplotigs.

                # Busco id      Status  Contig  Start   End     Score   Length
                EOG092R000C     Duplicated      tig00000366_arrow_arrow_pilon_pilon     131006  143302  2922.7  2611
                EOG092R000C     Duplicated      tig00063507_arrow_arrow_pilon_pilon     1422836 1435132 2923.4  2612
                EOG092R000I     Duplicated      tig00000071_arrow_arrow_pilon_pilon     294298  305900  1987.6  1916
                EOG092R000I     Duplicated      tig00063499_arrow_arrow_pilon_pilon     2581873 2593475 2089.2  1940

      Input:    BUSCO results table.

      Return:   A dictionary that holds contigs pairs and the number of BUSCOs they share.
  """

  f = open(BUSCO_TABLE, 'r')
  content = f.readlines()
  f.close()

  BUSCO_PAIRS = {}
  BUSCOS = {}

  for line in content:
    if line.startswith('#'):
      pass
    else:
      busco_ident = line.split('\t')[0]
      status = line.split('\t')[1]

      if status == 'Duplicated':
        contig = line.split('\t')[2]
        if busco_ident in BUSCOS:
          BUSCOS[busco_ident] = BUSCOS[busco_ident] + [contig]
        else:
          BUSCOS[busco_ident] = [contig]

  for busco, contig_list in BUSCOS.items():
    # Only use BUSCOs that appear exactly twice
    if len(contig_list) == 2.0:

      pair = (contig_list[0], contig_list[1])
      if pair in BUSCO_PAIRS:
        BUSCO_PAIRS[pair] += BUSCO_WEIGHT * 1
      else:
        BUSCO_PAIRS[pair] = BUSCO_WEIGHT * 1

      pair = (contig_list[1], contig_list[0])
      if pair in BUSCO_PAIRS:
        BUSCO_PAIRS[pair] += BUSCO_WEIGHT * 1
      else:
        BUSCO_PAIRS[pair] = BUSCO_WEIGHT * 1

  return BUSCO_PAIRS
#--------------------------------------
'''
This script requires a couple of input files:
1) a file that has the gene mapping, produced by biokanga blitz
2) a busco results table (full_table)
4) a path to the assembly fasta file
###############################
# Example for biokanga blitz file
module load biokanga/4.4.2
genome="LR1-canu-guppy4-haplotigs.contigs.fasta"
biokanga index --threads=4 -i ${genome} -o biokanga_index/genome -r gene_mapping
biokanga blitz --sensitivity=2 --mismatchscore=1 --threads=4 -o Pgt_Annotation_LR1.txt --in=Puccinia_graminis_tritici_21-0.transcripts.fa --sfx=biokanga_index/genome
###############################
'''
#--------------------------------------
'''f = open('/home/spe12g/LeafRust_LR1/A_141tigsTo18Chr_use.gff3')
content = f.readlines()
A_dic = {}
for line in content[1:]:
  chromosome = line.split()[0]
  contig = line.split()[8].replace('Name=','').strip()
  A_dic[contig] = chromosome
f.close()

f = open('/home/spe12g/LeafRust_LR1/B_135tigsTo18Chr_use.gff3')
content = f.readlines()
B_dic = {}
for line in content[1:]:
  chromosome = line.split()[0]
  contig = line.split()[8].replace('Name=','').strip()
  B_dic[contig] = chromosome
f.close()
#--------------------------------------'''
