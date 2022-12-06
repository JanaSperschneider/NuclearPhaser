#!/usr/bin/env python3
"""
    NuclearPhaser 1.1
    Copyright (C) 2022-2023 Jana Sperschneider
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    Contact: jana.sperschneider@csiro.au
"""
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
import GeneBinning
import math
import getopt

import matplotlib.pyplot as plt
# -----------------------------------------------------------------------------------------------------------
# Global variables
# -----------------------------------------------------------------------------------------------------------
#--------------------------------------
THRESHOLD_SYNTENY = 20.0
THRESHOLD_HIC_PERCENT = 80.0
THRESHOLD_HIC_COUNT = 20.0
THRESHOLD_MIN_ALIGNMENT = 50.0
MIN_LENGTH_TO_PHASE_WITH_HIC = 20000.0 # contig length in base pairs
HIC_CONTACT_THRESHOLD = 0.0
PHASE_SWITCH_CONTIG_MIN_SIZE = 1000000
# -----------------------------------------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------------------------------------
def check_contigs_are_haplotigs(HAPLOTIGS, pair):
  """ Function: check_contigs_are_haplotigs()
      Purpose:  Given two contigs, check if they are haplotigs.

      Input:    A pair of contigs and a dictionary of haplotigs.

      Return:   The % of bases that align to each other for each contig.
  """
  contig1, contig2 = pair[0], pair[1]
  contig1_aligned, contig2_aligned = 0.0, 0.0

  # Are these contigs likely haplotigs?
  if contig1 in HAPLOTIGS:
    haplotig_list_contig1 = HAPLOTIGS[contig1]
    for (bases_aligned, haplotig, alignment_blocks) in haplotig_list_contig1:
      if haplotig == contig2:
        contig1_aligned = 100.0*bases_aligned/LENGTHS[contig1]

  # Are these contigs likely haplotigs?
  if contig2 in HAPLOTIGS:
    haplotig_list_contig2 = HAPLOTIGS[contig2]
    for (bases_aligned, haplotig, alignment_blocks) in haplotig_list_contig2:
      if haplotig == contig1:
        contig2_aligned = 100.0*bases_aligned/LENGTHS[contig2]

  return contig1_aligned, contig2_aligned
#--------------------------------------
def read_in_contact_map_per_contig(CONTACT_MAP):

  CONTACTS_CIS_TRANS_PER_CONTIG, CONTACTS_PER_CONTIG = {}, {}

  with open(CONTACT_MAP ,'r') as f:
    for line in f:
      contig1 = line.split('\t')[0]
      start1 = int(line.split('\t')[1])
      end1 = int(line.split('\t')[2])
      contig2 = line.split('\t')[3]
      start2 = int(line.split('\t')[4])
      end2 = int(line.split('\t')[5])
      freq = float(line.split('\t')[6])

      if freq >= HIC_CONTACT_THRESHOLD:

        #---------------
        # Include both cis & trans contacts
        #---------------
        if contig1 in CONTACTS_CIS_TRANS_PER_CONTIG:
        	CONTACTS_CIS_TRANS_PER_CONTIG[contig1] = CONTACTS_CIS_TRANS_PER_CONTIG[contig1] + [(start1, end1, contig2, freq)]
        else:
        	CONTACTS_CIS_TRANS_PER_CONTIG[contig1] = [(start1, end1, contig2, freq)]

        if contig2 in CONTACTS_CIS_TRANS_PER_CONTIG:
        	CONTACTS_CIS_TRANS_PER_CONTIG[contig2] = CONTACTS_CIS_TRANS_PER_CONTIG[contig2] + [(start2, end2, contig1, freq)]
        else:
        	CONTACTS_CIS_TRANS_PER_CONTIG[contig2] = [(start2, end2, contig1, freq)]

        #---------------
        # Include only trans contacts
        #---------------
        if contig1 != contig2:
          #---------------
          if contig1 in CONTACTS_PER_CONTIG:
            CONTACTS_PER_CONTIG[contig1] = CONTACTS_PER_CONTIG[contig1] + [(start1, end1, contig2, freq)]
          else:
            CONTACTS_PER_CONTIG[contig1] = [(start1, end1, contig2, freq)]
          #---------------
          if contig2 in CONTACTS_PER_CONTIG:
            CONTACTS_PER_CONTIG[contig2] = CONTACTS_PER_CONTIG[contig2] + [(start2, end2, contig1, freq)]
          else:
            CONTACTS_PER_CONTIG[contig2] = [(start2, end2, contig1, freq)]


  return CONTACTS_CIS_TRANS_PER_CONTIG, CONTACTS_PER_CONTIG
#--------------------------------------
def read_in_contact_map(CONTACT_MAP):

  CONTACTS_CIS_TRANS, CONTACTS = {}, {}

  with open(CONTACT_MAP ,'r') as f:
    for line in f:
      start1 = int(line.split('\t')[1])
      end1 = int(line.split('\t')[2])
      BINSIZE = end1-start1
      break

  with open(CONTACT_MAP ,'r') as f:
    for line in f:
      contig1 = line.split('\t')[0]
      start1 = int(line.split('\t')[1])
      end1 = int(line.split('\t')[2])
      contig2 = line.split('\t')[3]
      start2 = int(line.split('\t')[4])
      end2 = int(line.split('\t')[5])
      freq = float(line.split('\t')[6])

      if freq >= HIC_CONTACT_THRESHOLD:

        #---------------
        # Include both cis & trans contacts
        #---------------
        pair = (contig1, contig2)

        if pair in CONTACTS_CIS_TRANS:
        	previous_freq = CONTACTS_CIS_TRANS[pair]
        	CONTACTS_CIS_TRANS[pair] = previous_freq + freq
        else:
        	CONTACTS_CIS_TRANS[pair] = freq

        pair = (contig2, contig1)

        if pair in CONTACTS_CIS_TRANS:
        	previous_freq = CONTACTS_CIS_TRANS[pair]
        	CONTACTS_CIS_TRANS[pair] = previous_freq + freq
        else:
        	CONTACTS_CIS_TRANS[pair] = freq

        #---------------
        # Include only trans contacts
        #---------------
        if contig1 != contig2: ## This allows only trans contacts
          #---------------
          pair = (contig1, contig2)

          if pair in CONTACTS:
            previous_freq = CONTACTS[pair]
            CONTACTS[pair] = previous_freq + freq
          else:
            CONTACTS[pair] = freq
          #---------------
          pair = (contig2, contig1)

          if pair in CONTACTS:
            previous_freq = CONTACTS[pair]
            CONTACTS[pair] = previous_freq + freq
          else:
            CONTACTS[pair] = freq          
      #---------------

  return CONTACTS_CIS_TRANS, CONTACTS, BINSIZE
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
def HiC_links_shared_per_MB(length1, length2, counts):

  hic_density = 0.0

  shorter_contig_length = min(length1, length2)

  hic_density = counts/float(shorter_contig_length) * 1000000.0

  return hic_density
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
def check_haplotype_assignment(HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT):
  # Check the haplotype assignments
  # Swap contigs if necessary

  haplotype0, haplotype1 = HAPLOTYPE_BINS['Haplotype_0'], HAPLOTYPE_BINS['Haplotype_1']
  haplotype0_updated, haplotype1_updated = [], []

  assigned_to_wrong_haplotype = 0

  for contig in sorted(haplotype0):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      if haplotype_1_contacts > haplotype_0_contacts and 100.0*haplotype_1_contacts/total_contacts > THRESHOLD_HIC_PERCENT:
        #print(contig, 'Haplotype_0 is wrongly assigned.', 100.0*haplotype_0_contacts/total_contacts, 100.0*haplotype_1_contacts/total_contacts, '(', LENGTHS[contig], 'bps)')
        assigned_to_wrong_haplotype += 1
        haplotype1_updated.append(contig)
      else:
        haplotype0_updated.append(contig)
    else:
      haplotype0_updated.append(contig)

  for contig in sorted(haplotype1):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      if haplotype_0_contacts > haplotype_1_contacts and 100.0*haplotype_0_contacts/total_contacts > THRESHOLD_HIC_PERCENT:
        #print(contig, 'Haplotype_1 is wrongly assigned', 100.0*haplotype_0_contacts/total_contacts, 100.0*haplotype_1_contacts/total_contacts, '(', LENGTHS[contig], 'bps)')
        assigned_to_wrong_haplotype += 1
        haplotype0_updated.append(contig)
      else:
        haplotype1_updated.append(contig)
    else:
      haplotype1_updated.append(contig)

  return haplotype0_updated, haplotype1_updated
#--------------------------------------
def phase_contigs_with_HiC(unphased_contigs, HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT):

  for unphased_contig in sorted(unphased_contigs):
    if LENGTHS[unphased_contig] > MIN_LENGTH_TO_PHASE_WITH_HIC:

      haplotype_0_contacts = 0.0
      for contig in HAPLOTYPE_BINS['Haplotype_0']:
        if (unphased_contig, contig) in CONTACTS:
          if (unphased_contig, contig) not in BLACKLISTED_HIC_LINKS:
            haplotype_0_contacts += CONTACTS[(unphased_contig, contig)]

      haplotype_1_contacts = 0.0
      for contig in HAPLOTYPE_BINS['Haplotype_1']:
        if (unphased_contig, contig) in CONTACTS:
          if (unphased_contig, contig) not in BLACKLISTED_HIC_LINKS:
            haplotype_1_contacts += CONTACTS[(unphased_contig, contig)]

      total_contacts = haplotype_0_contacts + haplotype_1_contacts

      if total_contacts > THRESHOLD_HIC_COUNT:

        if 100.0*haplotype_0_contacts/total_contacts > THRESHOLD_HIC_PERCENT:
          HAPLOTYPE_BINS['Haplotype_0'] = HAPLOTYPE_BINS['Haplotype_0'] + [unphased_contig]

        if 100.0*haplotype_1_contacts/total_contacts > THRESHOLD_HIC_PERCENT:
          HAPLOTYPE_BINS['Haplotype_1'] = HAPLOTYPE_BINS['Haplotype_1'] + [unphased_contig]

  return HAPLOTYPE_BINS
#--------------------------------------
def place_contigs_with_synteny(HAPLOTYPE_BINS, identifiers, sequences):

  all_contigs_phased = HAPLOTYPE_BINS['Haplotype_0'] + HAPLOTYPE_BINS['Haplotype_1']
  unphased_contigs = [contig for contig in LENGTHS if contig not in all_contigs_phased]

  print('Unplaced contigs before synteny assignment:', len(unphased_contigs), total_length(unphased_contigs, LENGTHS)/1000000.0, 'MB')

  additional_contigs_haplotype0, additional_contigs_haplotype1 = [], []

  # Now find the most likely haplotig and the bin it belongs to
  for unphased_contig in unphased_contigs:
    does_it_belong_to_haplotype0 = False
    does_it_belong_to_haplotype1 = False

    for contig0 in HAPLOTYPE_BINS['Haplotype_0']:
      unphased_contig_aligned, contig0_aligned = check_contigs_are_haplotigs(HAPLOTIGS, (unphased_contig, contig0))
      if unphased_contig_aligned > THRESHOLD_MIN_ALIGNMENT:
        # This contig aligns well with haplotype 0, so it belongs in the opposite haplotype 1
        does_it_belong_to_haplotype1 = True

    for contig1 in HAPLOTYPE_BINS['Haplotype_1']:
      unphased_contig_aligned, contig1_aligned = check_contigs_are_haplotigs(HAPLOTIGS, (unphased_contig, contig1))
      if unphased_contig_aligned > THRESHOLD_MIN_ALIGNMENT:
        # This contig aligns well with haplotype 1, so it belongs in the opposite haplotype 0
        does_it_belong_to_haplotype0 = True

    if does_it_belong_to_haplotype0 == True and does_it_belong_to_haplotype1 == False:
      additional_contigs_haplotype0.append(unphased_contig)

    if does_it_belong_to_haplotype0 == False and does_it_belong_to_haplotype1 == True:
      additional_contigs_haplotype1.append(unphased_contig)

  HAPLOTYPE_BINS['Haplotype_0'] = HAPLOTYPE_BINS['Haplotype_0'] + additional_contigs_haplotype0
  HAPLOTYPE_BINS['Haplotype_1'] = HAPLOTYPE_BINS['Haplotype_1'] + additional_contigs_haplotype1

  return HAPLOTYPE_BINS
#--------------------------------------
def print_haplotype_assignment(HAPLOTYPE_BINS, THRESHOLD_HIC_COUNT, OUTPUT_FILE_WITH_ALLELIC, OUTPUT_FILE_WITHOUT_ALLELIC, OUTPUT_FILE_WITHOUT_ALLELIC_CIS_TRANS, OUTPUT_FILE_WITH_ALLELIC_CIS_TRANS):
  header = '#Contig\t% Hi-C contacts to haplotype 0\t% Hi-C contacts to haplotype 1\tHi-C contact frequency\tContig length\n'
  #-------------------------
  # Print without allelic links and include only trans
  #-------------------------  
  output_string = ""

  # Print out the haplotype assignments for plotting
  haplotype0, haplotype1 = HAPLOTYPE_BINS['Haplotype_0'], HAPLOTYPE_BINS['Haplotype_1']

  for contig in sorted(haplotype0):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      output_string += contig + '\t' + 'Haplotype_0' + '\t' + str(100.0*haplotype_0_contacts/total_contacts) + '\t'
      output_string += str(100.0*haplotype_1_contacts/total_contacts) + '\t' + str(total_contacts) + '\t' + str(LENGTHS[contig]) + '\n'

  for contig in sorted(haplotype1):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      output_string += contig + '\t' + 'Haplotype_1' + '\t' + str(100.0*haplotype_0_contacts/total_contacts) + '\t'
      output_string += str(100.0*haplotype_1_contacts/total_contacts) + '\t' + str(total_contacts) + '\t' + str(LENGTHS[contig]) + '\n'

  #-------------------------
  if OUTPUT_FILE_WITHOUT_ALLELIC != None:
    f = open(OUTPUT_DIRECTORY_PATH + '/' + OUTPUT_FILE_WITHOUT_ALLELIC, 'w')
    f.writelines(header)    
    f.writelines(output_string)
    f.close()
  #-------------------------
  # Print without allelic links and include both cis and trans
  #------------------------- 
  output_string = ""

  # Print out the haplotype assignments for plotting
  haplotype0, haplotype1 = HAPLOTYPE_BINS['Haplotype_0'], HAPLOTYPE_BINS['Haplotype_1']

  for contig in sorted(haplotype0):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS_CIS_TRANS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS_CIS_TRANS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS_CIS_TRANS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS_CIS_TRANS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      output_string += contig + '\t' + 'Haplotype_0' + '\t' + str(100.0*haplotype_0_contacts/total_contacts) + '\t'
      output_string += str(100.0*haplotype_1_contacts/total_contacts) + '\t' + str(total_contacts) + '\t' + str(LENGTHS[contig]) + '\n'

  for contig in sorted(haplotype1):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS_CIS_TRANS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS_CIS_TRANS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS_CIS_TRANS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS_CIS_TRANS and (contig, x) not in BLACKLISTED_HIC_LINKS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      output_string += contig + '\t' + 'Haplotype_1' + '\t' + str(100.0*haplotype_0_contacts/total_contacts) + '\t'
      output_string += str(100.0*haplotype_1_contacts/total_contacts) + '\t' + str(total_contacts) + '\t' + str(LENGTHS[contig]) + '\n'

  #-------------------------
  if OUTPUT_FILE_WITHOUT_ALLELIC_CIS_TRANS != None:
    f = open(OUTPUT_DIRECTORY_PATH + '/' + OUTPUT_FILE_WITHOUT_ALLELIC_CIS_TRANS, 'w')
    f.writelines(header)        
    f.writelines(output_string)
    f.close()
  #------------------------- 
  # Print with allelic links and include only trans
  #-------------------------  
  output_string = ""
  # Print out the haplotype assignments for plotting
  haplotype0, haplotype1 = HAPLOTYPE_BINS['Haplotype_0'], HAPLOTYPE_BINS['Haplotype_1']

  for contig in sorted(haplotype0):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      output_string += contig + '\t' + 'Haplotype_0' + '\t' + str(100.0*haplotype_0_contacts/total_contacts) + '\t'
      output_string += str(100.0*haplotype_1_contacts/total_contacts) + '\t' + str(total_contacts) + '\t' + str(LENGTHS[contig]) + '\t' + str(haplotype_0_contacts) + '\t' + str(haplotype_1_contacts) + '\n'

  for contig in sorted(haplotype1):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      output_string += contig + '\t' + 'Haplotype_1' + '\t' + str(100.0*haplotype_0_contacts/total_contacts) + '\t'
      output_string += str(100.0*haplotype_1_contacts/total_contacts) + '\t' + str(total_contacts) + '\t' + str(LENGTHS[contig]) + '\t' + str(haplotype_0_contacts) + '\t' + str(haplotype_1_contacts)  + '\n'

  if OUTPUT_FILE_WITH_ALLELIC != None:
    f = open(OUTPUT_DIRECTORY_PATH + '/' + OUTPUT_FILE_WITH_ALLELIC, 'w')
    f.writelines(header)        
    f.writelines(output_string)
    f.close()
  #------------------------- 
  # Print with allelic links and include both cis and trans
  #-------------------------  
  output_string = ""
  # Print out the haplotype assignments for plotting
  haplotype0, haplotype1 = HAPLOTYPE_BINS['Haplotype_0'], HAPLOTYPE_BINS['Haplotype_1']

  for contig in sorted(haplotype0):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS_CIS_TRANS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS_CIS_TRANS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS_CIS_TRANS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS_CIS_TRANS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      output_string += contig + '\t' + 'Haplotype_0' + '\t' + str(100.0*haplotype_0_contacts/total_contacts) + '\t'
      output_string += str(100.0*haplotype_1_contacts/total_contacts) + '\t' + str(total_contacts) + '\t' + str(LENGTHS[contig]) + '\n'

  for contig in sorted(haplotype1):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS_CIS_TRANS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS_CIS_TRANS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS_CIS_TRANS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS_CIS_TRANS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      output_string += contig + '\t' + 'Haplotype_1' + '\t' + str(100.0*haplotype_0_contacts/total_contacts) + '\t'
      output_string += str(100.0*haplotype_1_contacts/total_contacts) + '\t' + str(total_contacts) + '\t' + str(LENGTHS[contig]) + '\n'

  if OUTPUT_FILE_WITH_ALLELIC_CIS_TRANS != None:
    f = open(OUTPUT_DIRECTORY_PATH + '/' + OUTPUT_FILE_WITH_ALLELIC_CIS_TRANS, 'w')
    f.writelines(header)        
    f.writelines(output_string)
    f.close()
  return
#--------------------------------------
def print_phased_contig_stats(HAPLOTYPE_BINS):

  all_contigs_phased = HAPLOTYPE_BINS['Haplotype_0'] + HAPLOTYPE_BINS['Haplotype_1']
  unphased_contigs = [contig for contig in LENGTHS if contig not in all_contigs_phased]

  haplotype0, haplotype1 = HAPLOTYPE_BINS['Haplotype_0'], HAPLOTYPE_BINS['Haplotype_1']

  haplotype_0_contacts_to_haplotype0, haplotype_0_contacts_to_haplotype1 = [], []
  haplotype_1_contacts_to_haplotype0, haplotype_1_contacts_to_haplotype1 = [], []

  for contig in sorted(haplotype0):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      haplotype_0_contacts_to_haplotype0.append(100.0*haplotype_0_contacts/total_contacts)
      haplotype_0_contacts_to_haplotype1.append(100.0*haplotype_1_contacts/total_contacts)

  for contig in sorted(haplotype1):
    # Check if this contig really belongs to its haplotype
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0

    haplotype_0_contacts = [CONTACTS[(contig, x)] for x in haplotype0 if (contig, x) in CONTACTS]
    haplotype_0_contacts = sum([contact_count for contact_count in haplotype_0_contacts])

    haplotype_1_contacts = [CONTACTS[(contig, x)] for x in haplotype1 if (contig, x) in CONTACTS]
    haplotype_1_contacts = sum([contact_count for contact_count in haplotype_1_contacts])

    total_contacts = haplotype_0_contacts + haplotype_1_contacts

    if total_contacts > THRESHOLD_HIC_COUNT:
      haplotype_1_contacts_to_haplotype0.append(100.0*haplotype_0_contacts/total_contacts)
      haplotype_1_contacts_to_haplotype1.append(100.0*haplotype_1_contacts/total_contacts)

  mean_haplotype_0_contacts_to_haplotype0 = round(sum(haplotype_0_contacts_to_haplotype0)/len(haplotype_0_contacts_to_haplotype0), 2)
  mean_haplotype_1_contacts_to_haplotype1 = round(sum(haplotype_1_contacts_to_haplotype1)/len(haplotype_1_contacts_to_haplotype1), 2)

  print()
  print('Total number of contigs:', len(LENGTHS))
  print('Number of contigs in haplotype 0:', len(HAPLOTYPE_BINS['Haplotype_0']))
  print('Number of contigs in haplotype 1:', len(HAPLOTYPE_BINS['Haplotype_1']))
  print('Number of unphased contigs:', len(unphased_contigs))
  print()
  print('Total length of haplotype 0 (bps):', total_length(HAPLOTYPE_BINS['Haplotype_0'], LENGTHS), '(', round(total_length(HAPLOTYPE_BINS['Haplotype_0'], LENGTHS)/1000000.0,3),'Mb)')
  print('Total length of haplotype 1 (bps):', total_length(HAPLOTYPE_BINS['Haplotype_1'], LENGTHS), '(', round(total_length(HAPLOTYPE_BINS['Haplotype_1'], LENGTHS)/1000000.0,3),'Mb)')
  print('Total length of unphased contigs (bps):', total_length(unphased_contigs, LENGTHS), '(', round(total_length(unphased_contigs, LENGTHS)/1000000.0,3),'Mb)')
  print()
  print('Haplotype 0: % of Hi-C trans contacts that are within the same haplotype:', str(mean_haplotype_0_contacts_to_haplotype0))
  print('Haplotype 1: % of Hi-C trans contacts that are within the same haplotype:', str(mean_haplotype_1_contacts_to_haplotype1))

  return unphased_contigs
# -----------------------------------------------------------------------------------------------------------
def usage():
    """ Function: usage()
        Purpose:  Print helpful information for the user.        
        
        Input:    None.
    
        Return:   Print options for running EffectorP 3.0.       
    """
    print('''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NuclearPhaser 1.1: Phasing of dikaryotic fungal genome assemblies
# Copyright (C) 2022-2023 Jana Sperschneider.
# Freely distributed under the GNU General Public License (GPLv3).
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ''')
    print("Usage for NuclearPhaser 1.1: ")
    print("python NuclearPhaser.py [-options] -g <gene mapping file> -b <BUSCO output table> -c <Hi-C contact map> -f <genome FASTA file> -o <output folder>")
    print()
    print("where required options are:")
    print("-g : biokanga blitz gene mapping file")    
    print("-b : BUSCO full output table")    
    print("-c : Hi-C contact map")    
    print("-f : Genome FASTA file")    
    print("-o : Output folder name")            
    print()
    print("where other options are:")   
    print("-p : p*BINSIZE is the minimum size of a phase switch region. For example, if your Hi-C contact map has binsize 100 Kb, NuclearPhaser will warn about phase switch regions >= p*100 KB (default: p=4).")    
    print("-h : show brief help on version and usage")
    print()
    print("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -")
    print()
    sys.exit()    

    return
# -----------------------------------------------------------------------------------------------------------
def scan_arguments(commandline):
    """ Function: scan_arguments()
        Purpose:  Scan the input options given to the NuclearPhaser program.        
        
        Input:    Input options given by the user.
    
        Return:   Parsed options.
    """
    try:
        opts, args = getopt.getopt(commandline, "hg:b:c:f:o:p:", ["help"])        
    except getopt.GetoptError as err:
        # print help information and exit:
        print(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    GENE_MAPPING = None
    BUSCO_TABLE = None
    CONTACT_MAP = None
    PATH_TO_GENOME_FASTA = None
    OUTPUT_DIRECTORY_PATH = None
    PHASE_SWITCH_REGION_MIN_SIZE = None

    g_count, b_count, c_count, f_count, o_count = 0, 0, 0, 0, 0
   
    for opt, arg in opts:
        if opt in ("-g"):
            GENE_MAPPING = arg
            g_count += 1
        
        elif opt in ("-b"):
            BUSCO_TABLE = arg
            b_count += 1

        elif opt in ("-c"):
            CONTACT_MAP = arg
            c_count += 1

        elif opt in ("-f"):
            PATH_TO_GENOME_FASTA = arg        
            f_count += 1    

        elif opt in ("-o"):
            OUTPUT_DIRECTORY_PATH = arg
            o_count += 1

        elif opt in ("-p"):
          PHASE_SWITCH_REGION_MIN_SIZE = arg
          if float(PHASE_SWITCH_REGION_MIN_SIZE) < 1.0:
            print("Please provide a number >= 1 for the -p parameter")
            usage()

        elif opt in ("-h", "--help"):
            usage()
        else:
            print()
            print ("Commandline option was supplied that was not recognized!")
            usage()

    if g_count == 0 or b_count == 0 or c_count == 0 or f_count == 0 or o_count == 0:
       print('One of the arguments was not provided.')

    if g_count != 1 or b_count != 1 or c_count != 1 or f_count != 1 or o_count != 1:
       usage()

    return GENE_MAPPING, BUSCO_TABLE, CONTACT_MAP, PATH_TO_GENOME_FASTA, OUTPUT_DIRECTORY_PATH, PHASE_SWITCH_REGION_MIN_SIZE
# -----------------------------------------------------------------------------------------------------------
def investigate_phaseswitch_on_contig(contig, length, header, PHASESWITCH_CANDIDATES, identifier):

  # Now get the haplotigs of this contig
  for gene_bin, (contigs_a, contigs_b) in sorted(GENE_BINS.items()):
    if contig in contigs_a:
      haplotigs = contigs_b
    if contig in contigs_b:
      haplotigs = contigs_a

  haplotype_0_contacts_contig = sum([CONTACTS[(contig, x)] for x in HAPLOTYPE_BINS['Haplotype_0'] if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS])
  haplotype_1_contacts_contig = sum([CONTACTS[(contig, x)] for x in HAPLOTYPE_BINS['Haplotype_1'] if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS])

  haplotype_0_contacts_contig_with_allelic = sum([CONTACTS[(contig, x)] for x in HAPLOTYPE_BINS['Haplotype_0'] if (contig, x) in CONTACTS])
  haplotype_1_contacts_contig_with_allelic = sum([CONTACTS[(contig, x)] for x in HAPLOTYPE_BINS['Haplotype_1'] if (contig, x) in CONTACTS])

  haplotig_contacts_contig = sum([CONTACTS[(contig, x)] for x in LENGTHS if (contig, x) in CONTACTS and (contig, x) in BLACKLISTED_HIC_LINKS])

  total_contacts = haplotype_0_contacts_contig + haplotype_1_contacts_contig   
  total_contacts_with_allelic = haplotype_0_contacts_contig_with_allelic + haplotype_1_contacts_contig_with_allelic

  if contig in CONTACTS_PER_CONTIG:
    # Print out the trans contact information
    f = open(OUTPUT_DIRECTORY_PATH + '/Plots/' + contig + '_HiC_Contacts_TransOnly.txt', 'w')     
    f.writelines(header)

    contact_list = CONTACTS_PER_CONTIG[contig]

    # First collect all trans contacts in each bin
    CONTACTS_PER_BIN = {}
    for (start, end, contig2, freq) in contact_list:
      bin_coordinates = (start, end)
      if bin_coordinates in CONTACTS_PER_BIN:
        CONTACTS_PER_BIN[bin_coordinates] = CONTACTS_PER_BIN[bin_coordinates] + [(contig2, freq)]
      else:
        CONTACTS_PER_BIN[bin_coordinates] = [(contig2, freq)]

    hap0_coordinates, hap1_coordinates = [], []

    x_axis = []
    y_axis_hap0, y_axis_hap1 = [], []

    # Then print out the % of contacts to each haplotype for each bin
    for (start, end), contact_list in sorted(CONTACTS_PER_BIN.items()):
      haplotype_0_contacts = sum([freq for (contig2, freq) in contact_list if contig2 in HAPLOTYPE_BINS['Haplotype_0']])
      haplotype_1_contacts = sum([freq for (contig2, freq) in contact_list if contig2 in HAPLOTYPE_BINS['Haplotype_1']])

      haplotig_contacts = sum([freq for (contig2, freq) in contact_list if contig2 in haplotigs])

      if haplotype_0_contacts > haplotype_1_contacts and 100.0*haplotype_0_contacts/(haplotype_0_contacts+haplotype_1_contacts) > 80.0: # A phaseswitch region should have > 80% of Hi-C links going to the other haplotype
        hap0_coordinates.append((start, end))
      if haplotype_0_contacts < haplotype_1_contacts and 100.0*haplotype_1_contacts/(haplotype_0_contacts+haplotype_1_contacts) > 80.0: # A phaseswitch region should have > 80% of Hi-C links going to the other haplotype
        hap1_coordinates.append((start, end))

      total_contacts = haplotype_0_contacts+haplotype_1_contacts

      if total_contacts > 0.0:
        output = str(start) + '\t' + str(end) + '\t'
        output += "".join(['|' for i in range(0,round(100.0*haplotype_0_contacts/total_contacts))])
        output += "".join(['*' for i in range(0,round(100.0*haplotype_1_contacts/total_contacts))])
        output += '\t'
        output += str(round(100.0*haplotype_0_contacts/total_contacts,2)) + '\t' + str(round(100.0*haplotype_1_contacts/total_contacts,2)) + '\t'
        output += str(haplotype_0_contacts) + '\t' + str(haplotype_1_contacts) + '\t'
        output += str(round(100.0*haplotig_contacts/total_contacts,2)) + '\t' + str(haplotig_contacts) + '\n'
        f.writelines(output)

        x_axis.append(start)
        y_axis_hap0.append(round(100.0*haplotype_0_contacts/total_contacts,2))
        y_axis_hap1.append(round(100.0*haplotype_1_contacts/total_contacts,2))
    f.close()

    x_axis_complete = [*range(0,length,500000)]

    try:
      plt.bar(x_axis, y_axis_hap0, width=BINSIZE*0.8)
      plt.bar(x_axis, y_axis_hap1, width=BINSIZE*0.8, bottom=y_axis_hap0)

      plt.xticks(x_axis_complete, rotation='vertical')

      plt.xlabel('Position on contig (Mb)')
      plt.ylabel('% of Hi-C trans contacts to haplotype')
      plt.title(identifier + ': ' + contig)

      #plt.show()
      plt.savefig(OUTPUT_DIRECTORY_PATH + '/Plots/' + identifier + '_' + contig + '.png')
      plt.clf()
      plt.cla()
    except:
      pass

    merged_intervals_hap0 = GeneBinning.merge_intervals(hap0_coordinates, 0)
    merged_intervals_hap1 = GeneBinning.merge_intervals(hap1_coordinates, 0)

    merged_intervals_hap0_threshold = [(x,y) for (x,y) in merged_intervals_hap0 if y-x >= PHASE_SWITCH_REGION_MIN_SIZE]
    merged_intervals_hap1_threshold = [(x,y) for (x,y) in merged_intervals_hap1 if y-x >= PHASE_SWITCH_REGION_MIN_SIZE]

    if merged_intervals_hap0_threshold != [] and merged_intervals_hap1_threshold != []:
      if round(100.0*haplotig_contacts_contig/total_contacts_with_allelic, 2) > 0.0:
        PHASESWITCH_CANDIDATES.append((contig, length, haplotype_0_contacts_contig, haplotype_1_contacts_contig, merged_intervals_hap0_threshold, merged_intervals_hap1_threshold, haplotig_contacts_contig, total_contacts_with_allelic))

  f_alignments = open(OUTPUT_DIRECTORY_PATH + '/Plots/' + contig + '_Haplotigs.txt', 'w')
  # Now print out the alignment coordinates for the haplotigs, too
  ALIGNMENTS = []
  for (bases_aligned, potential_haplotig, merged_hits) in HAPLOTIGS[contig]:
      contig_aligned, potential_haplotig_aligned = check_contigs_are_haplotigs(HAPLOTIGS, (contig, potential_haplotig))

      if potential_haplotig_aligned > 25.0:
        for (start, end) in merged_hits:
          if end-start > 20000:
            ALIGNMENTS.append((start, end, potential_haplotig))

  for (start, end, potential_haplotig) in sorted(ALIGNMENTS):
    haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0
    haplotype_0_contacts = sum([CONTACTS[(potential_haplotig, x)] for x in HAPLOTYPE_BINS['Haplotype_0'] if (potential_haplotig, x) in CONTACTS and (potential_haplotig, x) not in BLACKLISTED_HIC_LINKS])
    haplotype_1_contacts = sum([CONTACTS[(potential_haplotig, x)] for x in HAPLOTYPE_BINS['Haplotype_1'] if (potential_haplotig, x) in CONTACTS and (potential_haplotig, x) not in BLACKLISTED_HIC_LINKS])
    sum_of_contacts = haplotype_0_contacts + haplotype_1_contacts
    output = potential_haplotig + '\t' + str(start) + '\t' + str(end)
    if sum_of_contacts > 0.0:
      output += '\t' + str(round(100.0*haplotype_0_contacts/sum_of_contacts, 2)) + '\t' + str(round(100.0*haplotype_1_contacts/sum_of_contacts, 2)) + '\n'
    else:
      output += '\t' + '0.0' + '\t' + '0.0' + '\n'
    f_alignments.writelines(output)
  f_alignments.close()

  return PHASESWITCH_CANDIDATES
# -----------------------------------------------------------------------------------------------------------

#--------------------------------------
#--------------------------------------
#--------------------------------------
'''
This script requires a couple of input files:
1) a file that has the gene mapping, produced by biokanga blitz
2) a busco results table (full_table)
3) a HiC contact map in ginteractions format, produced by HiC-Pro
4) a path to the assembly fasta file
###############################
# Example for biokanga blitz file
module load biokanga/4.4.2
genome="LR1-canu-guppy4-haplotigs.contigs.fasta"
biokanga index --threads=4 -i ${genome} -o biokanga_index/genome -r gene_mapping
biokanga blitz --sensitivity=2 --mismatchscore=1 --threads=4 -o Pgt_Annotation_LR1.txt --in=Puccinia_graminis_tritici_21-0.transcripts.fa --sfx=biokanga_index/genome
###############################
'''
commandline = sys.argv[1:]
# -----------------------------------------------------------------------------------------------------------
if commandline:
    GENE_MAPPING, BUSCO_TABLE, CONTACT_MAP, PATH_TO_GENOME_FASTA, OUTPUT_DIRECTORY_PATH, PHASE_SWITCH_REGION_MIN_SIZE = scan_arguments(commandline)
else:
    usage()
#--------------------------------------
#--------------------------------------
#--------------------------------------
#--------------------------------------
#--------------------------------------
# Try to create folder where results will be stored
try:
    os.mkdir(OUTPUT_DIRECTORY_PATH)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise
try:
    os.mkdir(OUTPUT_DIRECTORY_PATH + '/Plots/')
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise        
#--------------------------------------
# Read in the contigs from the genome FASTA file
identifiers, sequences = [], []
for identifier, sequence in SimpleFastaParser(open(PATH_TO_GENOME_FASTA, 'r')):
  # Use simplified headers without whitespaces
  identifiers.append(identifier.split()[0])
  sequences.append(sequence)
#--------------------------------------
# Read in the lengths of each contig for later use
LENGTHS = {}
for identifier, sequence in zip(identifiers, sequences):
  LENGTHS[identifier] = len(sequence)
#--------------------------------------
# Read in haplotig information from minimap2 alignment file
HAPLOTIGS = GeneBinning.find_haplotigs(PATH_TO_GENOME_FASTA, OUTPUT_DIRECTORY_PATH, LENGTHS)
#--------------------------------------
print('Find contig pairs that share genes.')
BUSCO_PAIRS = GeneBinning.read_in_BUSCOs(BUSCO_TABLE)
CONTIG_PAIRS = GeneBinning.read_gene_mapping(GENE_MAPPING)
#--------------------------------------
# Now merge the gene mapping information with the BUSCO information
# to get the total numbers of genes shared per contig pair
SHARED_GENES = GeneBinning.genes_shared(BUSCO_PAIRS, CONTIG_PAIRS, LENGTHS)
#--------------------------------------
# Write the shared gene content to a file
GeneBinning.write_shared_genes_to_file(OUTPUT_DIRECTORY_PATH, SHARED_GENES, BUSCO_PAIRS, CONTIG_PAIRS, LENGTHS)
#--------------------------------------
print('Read in contact map')
CONTACTS_CIS_TRANS, CONTACTS, BINSIZE  = read_in_contact_map(CONTACT_MAP)
CONTACTS_CIS_TRANS_PER_CONTIG, CONTACTS_PER_CONTIG = read_in_contact_map_per_contig(CONTACT_MAP)
print('Hi-C contact matrix has resolution:', BINSIZE)
print('Done')
#--------------------------------------
if PHASE_SWITCH_REGION_MIN_SIZE == None:
  PHASE_SWITCH_REGION_MIN_SIZE = BINSIZE*4.0
else:
  PHASE_SWITCH_REGION_MIN_SIZE = BINSIZE*float(PHASE_SWITCH_REGION_MIN_SIZE)
#--------------------------------------
# First step is to generate bins of contig pairs that share genes == parts of ~chromosomes a/b in each bin
#--------------------------------------
#--------------------------------------
print('Construct graph from contig pairs to find well-connected groups')
G = nx.Graph()

for pair, gene_density in SHARED_GENES.items():
  pair1_aligned, pair2_aligned = check_contigs_are_haplotigs(HAPLOTIGS, pair)

  if pair1_aligned > THRESHOLD_SYNTENY or pair2_aligned > THRESHOLD_SYNTENY:
    G.add_edge(pair[0], pair[1], weight=gene_density)

print('Construct community')
partition = community.best_partition(G, randomize=False)
#--------------------------------------
GENE_BINS = GeneBinning.partition_to_genebins(partition, LENGTHS, SHARED_GENES)
if len(GENE_BINS) == 0.0:
  print('The gene binning did not work. Check that the BUSCO and biokanga blitz input files are valid and contain contigs that have genes that occur exactly twice.')
  sys.exit()

print('Use', len(GENE_BINS), 'gene bins.')
print()
#--------------------------------------
all_contigs_phased = []
for gene_bin, (contigs_a, contigs_b) in GENE_BINS.items():
  all_contigs_phased += [c for c in contigs_a]
  all_contigs_phased += [c for c in contigs_b]

print('Contigs in bins:', len(all_contigs_phased), 'with a total length of', total_length(all_contigs_phased, LENGTHS)/1000000.0, 'MB')
unphased_contigs = [contig for contig in LENGTHS if contig not in all_contigs_phased]
print('Still unphased:', len(unphased_contigs), 'contigs with a total length of', total_length(unphased_contigs, LENGTHS)/1000000.0, 'MB')
print()
#--------------------------------------
#--------------------------------------
# Now print the gene binning dictionary to output file
#--------------------------------------
#--------------------------------------
GeneBinning.gene_bins_to_file(OUTPUT_DIRECTORY_PATH, GENE_BINS, LENGTHS)
#--------------------------------------
#--------------------------------------
# Now make FASTA files for bins and their haplotypes for dot-plot alignments
for key, values in GENE_BINS.items():
  make_fasta_files(values[0], OUTPUT_DIRECTORY_PATH + '/' + key + '_0.fasta', identifiers, sequences)
  make_fasta_files(values[1], OUTPUT_DIRECTORY_PATH + '/' + key + '_1.fasta', identifiers, sequences)
#--------------------------------------
# For Hi-C phasing, do not allow Hi-C links between haplotigs
BLACKLISTED_HIC_LINKS = []
for gene_bin, (contigs_a, contigs_b) in sorted(GENE_BINS.items()):
  for x in contigs_a:
    for y in contigs_b:
      BLACKLISTED_HIC_LINKS.append((x,y))
      BLACKLISTED_HIC_LINKS.append((y,x))
#--------------------------------------
#--------------------------------------
# Now try and phase the gene bins with HiC
#--------------------------------------
#--------------------------------------
print('Construct graph to phase the gene bins with Hi-C trans contact data')
G = nx.Graph()

CONNECTIONS = {}

for gene_bin1, (contigs1_a, contigs1_b) in sorted(GENE_BINS.items()):
  for gene_bin2, (contigs2_a, contigs2_b) in sorted(GENE_BINS.items()):

    if gene_bin1 != gene_bin2:

      total_links_between_bins = 0.0
      a_to_a, a_to_b, b_to_a, b_to_b = 0.0, 0.0, 0.0, 0.0

      for x in contigs1_a:
        for y in contigs2_a:
          if (x,y) in CONTACTS:
            a_to_a += CONTACTS[(x,y)]

        for y in contigs2_b:
          if (x,y) in CONTACTS:
            a_to_b += CONTACTS[(x,y)]

      for x in contigs1_b:
        for y in contigs2_a:
          if (x,y) in CONTACTS:
            b_to_a += CONTACTS[(x,y)]

        for y in contigs2_b:
          if (x,y) in CONTACTS:
            b_to_b += CONTACTS[(x,y)]

      # Collect HiC links between the two bins
      total_links_between_bins = a_to_a + a_to_b + b_to_a + b_to_b
      same_haplotype, different_haplotype = 0.0, 0.0

      if total_links_between_bins > 0.0:

        same_haplotype = 100.0*(a_to_a + b_to_b)/total_links_between_bins
        different_haplotype = 100.0*(a_to_b + b_to_a)/total_links_between_bins

        G.add_edge(gene_bin1 + 'a', gene_bin2 + 'a', weight=same_haplotype)
        G.add_edge(gene_bin1 + 'a', gene_bin2 + 'b', weight=different_haplotype)
        G.add_edge(gene_bin1 + 'b', gene_bin2 + 'a', weight=different_haplotype)
        G.add_edge(gene_bin1 + 'b', gene_bin2 + 'b', weight=same_haplotype)

      CONNECTIONS[gene_bin1 + 'a', gene_bin2 + 'a'] = round(same_haplotype, 2)
      CONNECTIONS[gene_bin1 + 'b', gene_bin2 + 'b'] = round(same_haplotype, 2)
      CONNECTIONS[gene_bin1 + 'a', gene_bin2 + 'b'] = round(different_haplotype, 2)
      CONNECTIONS[gene_bin1 + 'b', gene_bin2 + 'a'] = round(different_haplotype, 2)

print()
#--------------------------------------
bin_numbers = sorted([int(gene_bin.replace('Bin_', '')) for gene_bin in GENE_BINS])

f = open(OUTPUT_DIRECTORY_PATH + '/gene_bins_hic_connections.txt', 'w')

header = '\t' + "\t".join(['Bin_' + str(bin_number) + 'a' + '\t' + 'Bin_' + str(bin_number) + 'b' for bin_number in bin_numbers])
f.writelines(header + '\n')

for bin_number1 in bin_numbers:
  gene_bin1a = 'Bin_' + str(bin_number1) + 'a'
  line = 'Bin_' + str(bin_number1) + 'a' + '\t'
  for bin_number2 in bin_numbers:
    if (gene_bin1a, 'Bin_' + str(bin_number2) + 'a') in CONNECTIONS:
      line += str(CONNECTIONS[(gene_bin1a,'Bin_' + str(bin_number2) + 'a')]) + '\t'
    else:
      line += '0.0' + '\t'
    if (gene_bin1a, 'Bin_' + str(bin_number2) + 'b') in CONNECTIONS:
      line += str(CONNECTIONS[(gene_bin1a,'Bin_' + str(bin_number2) + 'b')]) + '\t'
    else:
      line += '0.0' + '\t'
  f.writelines(line + '\n')

  gene_bin1b = 'Bin_' + str(bin_number1) + 'b'
  line = 'Bin_' + str(bin_number1) + 'b' + '\t'
  for bin_number2 in bin_numbers:
    if (gene_bin1b, 'Bin_' + str(bin_number2) + 'a') in CONNECTIONS:
      line += str(CONNECTIONS[(gene_bin1b,'Bin_' + str(bin_number2) + 'a')]) + '\t'
    else:
      line += '0.0' + '\t'
    if (gene_bin1b, 'Bin_' + str(bin_number2) + 'b') in CONNECTIONS:
      line += str(CONNECTIONS[(gene_bin1b,'Bin_' + str(bin_number2) + 'b')]) + '\t'
    else:
      line += '0.0' + '\t'
  f.writelines(line + '\n')

f.close()
#--------------------------------------
import community
print('Construct community to phase gene bins into the haplotypes')
partition = community.best_partition(G, randomize=False)
#--------------------------------------
print('----------')
#--------------------------------------
HAPLOTYPE_BINS = {}

for com in set(partition.values()) :
    list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]

    contigs_lengths = 0
    contigs = []

    for gene_bin in list_nodes:
      if 'a' in gene_bin:
        for item in GENE_BINS[gene_bin.replace('a','')][0]:
          contigs_lengths += LENGTHS[item]
          contigs.append(item)
      if 'b' in gene_bin:
        for item in GENE_BINS[gene_bin.replace('b','')][1]:
          contigs_lengths += LENGTHS[item]
          contigs.append(item)

    print('----------')
    print('Haplotype_' + str(com), sorted(list_nodes))
    print(contigs_lengths/1000000.0, 'MB')
    HAPLOTYPE_BINS['Haplotype_' + str(com)] = contigs
    print('----------')

if len(HAPLOTYPE_BINS.keys()) != 2.0:
  print('There are more than two haplotype bins, check what is going on. There might be contaminants or major phase switching.')
  sys.exit()
#--------------------------------------
# Make fasta files of the haplotypes, recommended to do a DGenies dot-plot alignment of these at this stage
print('----------------------------------------')
print('Recommended to do a DGenies dot-plot alignment of these two files at this stage to confirm that the gene binning & phasing went well:')
print(OUTPUT_DIRECTORY_PATH + '/Haplotype_0_genephasing.fasta')
print(OUTPUT_DIRECTORY_PATH + '/Haplotype_1_genephasing.fasta')
print('----------------------------------------')
make_fasta_files(HAPLOTYPE_BINS['Haplotype_0'], OUTPUT_DIRECTORY_PATH + '/Haplotype_0_genephasing.fasta', identifiers, sequences)
make_fasta_files(HAPLOTYPE_BINS['Haplotype_1'], OUTPUT_DIRECTORY_PATH + '/Haplotype_1_genephasing.fasta', identifiers, sequences)
#--------------------------------------
# Check the haplotype assignments, but do not swap contigs at this stage
#--------------------------------------
haplotype0_updated, haplotype1_updated = check_haplotype_assignment(HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)

OUTPUT_FILE_WITH_ALLELIC = 'Haplotype_Assignment_Contigs_in_ScaffoldBins_with_allelic_links_trans.txt'
OUTPUT_FILE_WITHOUT_ALLELIC = 'Haplotype_Assignment_Contigs_in_ScaffoldBins_without_allelic_links_trans.txt'
OUTPUT_FILE_WITH_ALLELIC_CIS_TRANS = 'Haplotype_Assignment_Contigs_in_ScaffoldBins_with_allelic_links_cis_trans.txt'
OUTPUT_FILE_WITHOUT_ALLELIC_CIS_TRANS = 'Haplotype_Assignment_Contigs_in_ScaffoldBins_without_allelic_links_cis_trans.txt'

print_haplotype_assignment(HAPLOTYPE_BINS, 0.0, OUTPUT_FILE_WITH_ALLELIC, OUTPUT_FILE_WITHOUT_ALLELIC, OUTPUT_FILE_WITHOUT_ALLELIC_CIS_TRANS, OUTPUT_FILE_WITH_ALLELIC_CIS_TRANS)
#--------------------------------------
all_contigs_phased = HAPLOTYPE_BINS['Haplotype_0'] + HAPLOTYPE_BINS['Haplotype_1']
unphased_contigs = [contig for contig in LENGTHS if contig not in all_contigs_phased]
print()
print('Haplotype 0 length', total_length(HAPLOTYPE_BINS['Haplotype_0'], LENGTHS)/1000000.0, 'MB')
print('Haplotype 1 length', total_length(HAPLOTYPE_BINS['Haplotype_1'], LENGTHS)/1000000.0, 'MB')
print('Unphased contigs', total_length(unphased_contigs, LENGTHS)/1000000.0, 'MB')
#--------------------------------------
#--------------------------------------
# Now print haplotype phasing information for each scaffold bin
for gene_bin, (contigs_a, contigs_b) in sorted(GENE_BINS.items()):

  haplotype_0_contacts, haplotype_1_contacts = 0.0, 0.0
  haplotype_0_contacts_bin_a, haplotype_1_contacts_bin_a = [], []
  haplotype_0_contacts_bin_b, haplotype_1_contacts_bin_b = [], []

  for contig in contigs_a:
    haplotype_0_contacts = sum([CONTACTS[(contig, x)] for x in HAPLOTYPE_BINS['Haplotype_0'] if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS])
    haplotype_1_contacts = sum([CONTACTS[(contig, x)] for x in HAPLOTYPE_BINS['Haplotype_1'] if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS])
    total_contacts = haplotype_0_contacts + haplotype_1_contacts
    if total_contacts > 0.0:
      haplotype_0_contacts_bin_a.append(round(100.0*haplotype_0_contacts/total_contacts,2))
      haplotype_1_contacts_bin_a.append(round(100.0*haplotype_1_contacts/total_contacts,2))
    else:
      haplotype_0_contacts_bin_a.append(0.0)
      haplotype_1_contacts_bin_a.append(0.0)

  for contig in contigs_b:
    haplotype_0_contacts = sum([CONTACTS[(contig, x)] for x in HAPLOTYPE_BINS['Haplotype_0'] if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS])
    haplotype_1_contacts = sum([CONTACTS[(contig, x)] for x in HAPLOTYPE_BINS['Haplotype_1'] if (contig, x) in CONTACTS and (contig, x) not in BLACKLISTED_HIC_LINKS])
    total_contacts = haplotype_0_contacts + haplotype_1_contacts
    if total_contacts > 0.0:
      haplotype_0_contacts_bin_b.append(round(100.0*haplotype_0_contacts/total_contacts,2))
      haplotype_1_contacts_bin_b.append(round(100.0*haplotype_1_contacts/total_contacts,2))
    else:
      haplotype_0_contacts_bin_b.append(0.0)
      haplotype_1_contacts_bin_b.append(0.0)

  print('------- Hi-C trans contacts in this bin', gene_bin, '-------')
  print(gene_bin)
  print('Haplotype 0 contigs:', contigs_a)
  print('Haplotype 1 contigs:', contigs_b)
  print('Haplotype 0 contigs - trans contact to haplotype 0:', haplotype_0_contacts_bin_a)
  print('Haplotype 1 contigs - trans contact to haplotype 0:', haplotype_0_contacts_bin_b)
  print()

print('----------------------------------------')
print('At this stage it is advised to break phase switches in the gene bins')
print('----------------------------------------')
#--------------------------------------
# Produce plots for each contig (> 1Mb) that visualize phase switches
PHASESWITCH_CANDIDATES = []
header = '#Start\tEnd\tContacts (|: haplotype 0 and *: haplotype 1)\t% Hi-C contacts to haplotype 0\t% Hi-C contacts to haplotype 1\tHi-C frequency to haplotype 0\tHi-C frequency to haplotype 1\t'
header += '%Hi-C frequency to haplotigs\tHi-C frequency to haplotigs\n'

for gene_bin, (contigs_a, contigs_b) in sorted(GENE_BINS.items()):

  for contig in contigs_a:
    length = LENGTHS[contig]
    identifier = gene_bin + '_haplotype_0'
    PHASESWITCH_CANDIDATES = investigate_phaseswitch_on_contig(contig, length, header, PHASESWITCH_CANDIDATES, identifier)
  for contig in contigs_b:
    length = LENGTHS[contig]    
    identifier = gene_bin + '_haplotype_1'    
    PHASESWITCH_CANDIDATES = investigate_phaseswitch_on_contig(contig, length, header, PHASESWITCH_CANDIDATES, identifier)

if PHASESWITCH_CANDIDATES == []:
  print('----------------------------------------')
  print('Found no phase switch regions >=', str(PHASE_SWITCH_REGION_MIN_SIZE), 'bps on contigs >', str(PHASE_SWITCH_CONTIG_MIN_SIZE/1000000.0), 'Mb.')
  print('----------------------------------------')
else:
  print('----------------------------------------')
  print('Found potential phase switch regions >=', str(PHASE_SWITCH_REGION_MIN_SIZE), 'bps on', len(PHASESWITCH_CANDIDATES), 'contigs >', str(PHASE_SWITCH_CONTIG_MIN_SIZE/1000000.0), 'Mb.')
  print('----------------------------------------')  

  for contig, length, haplotype_0_contacts_contig, haplotype_1_contacts_contig, merged_intervals_hap0_threshold, merged_intervals_hap1_threshold, haplotig_contacts_contig, total_contacts_with_allelic in PHASESWITCH_CANDIDATES:
    print('----------------------------------------')
    print('Potential phase switch in this contig >', str(PHASE_SWITCH_CONTIG_MIN_SIZE/1000000.0), 'Mb:', contig, '(', length, ' bps)')
    print('Hi-C trans contact to haplotype 0 (without allelic contacts):', str(round(100.0*haplotype_0_contacts_contig/(haplotype_0_contacts_contig+haplotype_1_contacts_contig), 2)), '%')
    print('Haplotype 0 regions on this contig:')
    print(merged_intervals_hap0_threshold)
    print('Hi-C trans contact to haplotype 1 (without allelic contacts):', str(round(100.0*haplotype_1_contacts_contig/(haplotype_0_contacts_contig+haplotype_1_contacts_contig), 2)), '%')  
    print('Haplotype 1 regions on this contig:')      
    print(merged_intervals_hap1_threshold)
    print('Allelic contacts:', str(round(100.0*haplotig_contacts_contig/(total_contacts_with_allelic), 2)), '%') 
    print('----------------------------------------')

#--------------------------------------
#--------------------------------------
# Now place the remaining contigs based on synteny
#--------------------------------------
HAPLOTYPE_BINS = place_contigs_with_synteny(HAPLOTYPE_BINS, identifiers, sequences)
#--------------------------------------
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
#--------------------------------------
# Check the haplotype assignments, but do not swap contigs at this stage
#--------------------------------------
haplotype0_updated, haplotype1_updated = check_haplotype_assignment(HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)
#--------------------------------------
# Now place the remaining contigs based on synteny
#--------------------------------------
HAPLOTYPE_BINS = place_contigs_with_synteny(HAPLOTYPE_BINS, identifiers, sequences)
#--------------------------------------
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
#--------------------------------------
# Check the haplotype assignments, but do not swap contigs at this stage
#--------------------------------------
haplotype0_updated, haplotype1_updated = check_haplotype_assignment(HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)
#--------------------------------------
#--------------------------------------
print('--------------------------------------')
print('Contigs have been placed based on synteny, now assign contigs with Hi-C')
#--------------------------------------
#--------------------------------------
print('--------------------------------------')
print('Round1')
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
HAPLOTYPE_BINS = phase_contigs_with_HiC(unphased_contigs, HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
print('--------------------------------------')
print('Round2')
HAPLOTYPE_BINS = phase_contigs_with_HiC(unphased_contigs, HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
print('--------------------------------------')
#--------------------------------------
#--------------------------------------
print('Now place the remaining contigs based on synteny')
print('--------------------------------------')
# Now place the remaining contigs based on synteny
#--------------------------------------
print('Round1')
HAPLOTYPE_BINS = place_contigs_with_synteny(HAPLOTYPE_BINS, identifiers, sequences)
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
#--------------------------------------
#--------------------------------------
print('Round2')
HAPLOTYPE_BINS = place_contigs_with_synteny(HAPLOTYPE_BINS, identifiers, sequences)
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
#--------------------------------------
#--------------------------------------
print('Contigs have been placed based on synteny, now assign contigs with Hi-C')
print('--------------------------------------')
#--------------------------------------
#--------------------------------------
print('Round1')
HAPLOTYPE_BINS = phase_contigs_with_HiC(unphased_contigs, HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
#--------------------------------------
#--------------------------------------
print('Round2')
HAPLOTYPE_BINS = phase_contigs_with_HiC(unphased_contigs, HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
print('--------------------------------------')
#--------------------------------------
#--------------------------------------
#--------------------------------------
# Now place the remaining contigs based on synteny
#--------------------------------------
HAPLOTYPE_BINS = place_contigs_with_synteny(HAPLOTYPE_BINS, identifiers, sequences)
#--------------------------------------
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
#--------------------------------------
HAPLOTYPE_BINS = place_contigs_with_synteny(HAPLOTYPE_BINS, identifiers, sequences)
#--------------------------------------
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
#--------------------------------------
all_contigs_phased = HAPLOTYPE_BINS['Haplotype_0'] + HAPLOTYPE_BINS['Haplotype_1']
unphased_contigs = [contig for contig in LENGTHS if contig not in all_contigs_phased]
print('--------------------------------------')
#--------------------------------------
#--------------------------------------
#--------------------------------------
#--------------------------------------
# Now as a last check, go through the
# HiC links for each contig and swap if unsupported
haplotype0_updated, haplotype1_updated = check_haplotype_assignment(HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)
HAPLOTYPE_BINS['Haplotype_0'] = haplotype0_updated
HAPLOTYPE_BINS['Haplotype_1'] = haplotype1_updated

# Run it again to resolve the swapped contigs in the step before
haplotype0_updated, haplotype1_updated = check_haplotype_assignment(HAPLOTYPE_BINS, THRESHOLD_HIC_PERCENT, THRESHOLD_HIC_COUNT)
HAPLOTYPE_BINS['Haplotype_0'] = haplotype0_updated
HAPLOTYPE_BINS['Haplotype_1'] = haplotype1_updated
#--------------------------------------
#--------------------------------------
#--------------------------------------
OUTPUT_FILE_WITH_ALLELIC = 'Haplotype_Assignment_Contigs_FINAL_with_allelic_links_trans.txt'
OUTPUT_FILE_WITHOUT_ALLELIC = 'Haplotype_Assignment_Contigs_FINAL_without_allelic_links_trans.txt'
OUTPUT_FILE_WITH_ALLELIC_CIS_TRANS = 'Haplotype_Assignment_Contigs_FINAL_with_allelic_links_cis_trans.txt'
OUTPUT_FILE_WITHOUT_ALLELIC_CIS_TRANS = 'Haplotype_Assignment_Contigs_FINAL_without_allelic_links_cis_trans.txt'

print_haplotype_assignment(HAPLOTYPE_BINS, 0.0, OUTPUT_FILE_WITH_ALLELIC, OUTPUT_FILE_WITHOUT_ALLELIC, OUTPUT_FILE_WITHOUT_ALLELIC_CIS_TRANS, OUTPUT_FILE_WITH_ALLELIC_CIS_TRANS)
#--------------------------------------
#--------------------------------------
make_fasta_files(HAPLOTYPE_BINS['Haplotype_0'], OUTPUT_DIRECTORY_PATH + '/Haplotype_0_final.fasta', identifiers, sequences)
make_fasta_files(HAPLOTYPE_BINS['Haplotype_1'], OUTPUT_DIRECTORY_PATH + '/Haplotype_1_final.fasta', identifiers, sequences)
make_fasta_files(unphased_contigs, OUTPUT_DIRECTORY_PATH + '/unphased_contigs_final.fasta', identifiers, sequences)
#--------------------------------------
print('--------------------------------------')
unphased_contigs = print_phased_contig_stats(HAPLOTYPE_BINS)
print
print('--------------------------------------')
print('All done.')
print('--------------------------------------')
#--------------------------------------
