#!/usr/bin/python
#import math
import sys
import os
import commands
import numpy as np
from collections import defaultdict

def fasta_iterator(fh):
  while True:
    line = fh.readline()
    if line.startswith('>'): break  
  while True:
    header = line[1:-1].rstrip()
    sequence = fh.readline().rstrip()
    while True:
      line = fh.readline()
      if not line: break
      if line.startswith('>'): break
      sequence += line.rstrip()
    yield(header, sequence)
    if not line: return

def Write_out(out, file):
  fh=open(file,"a")
  fh.write(out)
  fh.close()
  return()

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

def Check_consensus_raw(primer_tag_file_count, primer_tag_file):
  fh=open(primer_tag_file_count,"r")
  seqs1,inverse_seq1 = Tree(),{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      key,id = l[1]+"\t"+l[2],l[0]
      seqs1[key][id].value =1 
      inverse_seq1[id] = key
  fh.close()
  fh=open(primer_tag_file,"r")
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      key,id = l[3]+"\t"+l[4],l[0]
      if(key in seqs1):
        if(len(seqs1[key])>1):
          print len(seqs1[key]), l[0],l[1]
  fh.close()
  return()

def Get_consensus_raw(primer_tag_file,seq_file,consensus_counts_file):
  fh=open(primer_tag_file,"r")
  seqs2 ,inverse_seq2,cons_count = Tree(),{},{}
  for l in fh:
    if(l[0]!="#"):
      l=l.strip().split()
      consensus_count, id, seq = int(l[1]),l[0].split("__")[0],l[7]
      seqs2[seq][id].value = 1
      cons_count[id] = consensus_count
      inverse_seq2[id] = seq
  fh.close()
  done,total = [],0
  fh=open(consensus_counts_file,"w")
  fh.close()
  fh=open(seq_file, "r")
  out,ind = "id\tconsensus_count\tduplicate_count\n", 0
  for header,sequence in fasta_iterator(fh):
    total_count = map(int, header.split("|")[0].split("__")[1].split("_"))
    seqs_match,count = [],0
    for seq in seqs2:
      if(seq not in done):
        if(seq.count(sequence)==1):
          seqs_match.append(seq)
          done.append(seq)
          for i in seqs2[seq]:
            count = count+cons_count[i]
    ind = ind+1
    if(sum(total_count)!=count):
        print header, len(seqs_match), count, sum(total_count)
    if(sum(total_count)>count):
      print len(seqs_match), count, sum(total_count)
      count = sum(total_count)
    out=out+header+"\t"+str(count)+"\t"+str(sum(total_count))+"\n"
    if(ind>=100):
      Write_out(out, consensus_counts_file)
      out, ind = '',0
  fh.close()
  Write_out(out, consensus_counts_file)
  out, ind = '',0
  print len(done), len(seqs2), total, len(cons_count)
  return()


################################################
dir = sys.argv[1]
id = sys.argv[2]

################################################
primer_tag_file = dir+"ORIENTATED_SEQUENCES/TMP/Barcode_filtering_information_"+id+".txt"
primer_tag_file_count = dir+"ORIENTATED_SEQUENCES/TMP/All_barcodes_"+id+".txt"
nn_orf_filtered = dir+"ORIENTATED_SEQUENCES/Nucleotide_ORF_filtered_all_"+id+".fasta"
file_seqs = dir+"ORIENTATED_SEQUENCES/NETWORKS/Sequences_"+id+".txt"
Reduced_file = dir+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+id+".fasta"

consensus_counts_file =  dir+"ORIENTATED_SEQUENCES/NETWORKS/Consensus_counts_"+id+".txt"
################################################
Check_consensus_raw(primer_tag_file_count, primer_tag_file)
Get_consensus_raw(primer_tag_file,Reduced_file,consensus_counts_file)




