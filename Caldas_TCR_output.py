#!/usr/bin/python
#import math
import sys
import os
import commands
import numpy as np
from collections import defaultdict

class Tree(defaultdict):
  def __init__(self, value=None):
    super(Tree, self).__init__(Tree)
    self.value = value

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

def Get_sample_files(batch_file):
  id, dir = [],[]
  fh=open(batch_file,"r")
  ind = 0
  for l in fh:
    ind = ind+1
    if(l[0]!="#" and len(l)>3):
      l=l.strip().split()
      id.append(l[1])
      dir.append(l[5])
  fh.close()
  return(id, dir)

def Split_functional_non_function(id,dir): 
  for s in range(len(id)):
    print "\t\t", s, id
    sample,dir_use = id[s],dir[s]
    dir_IMGT = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT/"
    fasta_file = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+".fasta"
    fh= open(fasta_file,"r")
    ids = {}
    for header,sequence in fasta_iterator(fh):
      ids[header.split("__")[0]] = ">"+header+"\n"+sequence
    print "\r","read samples:",sample, "n sequences:",len(ids)
    fh.close()
    files = ['10_V-REGION-mutation-hotspots.txt','1_Summary.txt','2_IMGT-gapped-nt-sequences.txt','3_Nt-sequences.txt','4_IMGT-gapped-AA-sequences.txt',
        '5_AA-sequences.txt','6_Junction.txt','7_V-REGION-mutation-and-AA-change-table.txt','8_V-REGION-nt-mutation-statistics.txt','9_V-REGION-AA-change-statistics.txt']
    type_seq = Tree()
    file = dir_IMGT+"IMGT_"+sample+"_1_Summary.txt"
    fh=open(file,"r")
    for l in fh:
      l=l.strip().split("\t")
      if(l[0]!='Sequence number'):
        typ = l[2]
        if(typ=="productive (see comment)"):typ = "productive"
        if(typ=="unproductive (see comment)"):typ = "unproductive"
        if(typ in [ "productive", "unproductive"]):
          type_seq[typ][l[1].split("__")[0]].value = 1
    fh.close()
    for t in type_seq: 
      print "\t",t, len(type_seq[t])
      ### fasta file
      out, ind = '',0
      file = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+"_"+t+".fasta"
      fh=open(file, "w")
      fh.close()
      for id1 in type_seq[t]: 
        if(id1 in ids):
          out=out+ids[id1]+"\n"
          ind = ind+1
          if(ind>500):
            Write_out(out, file)
            out, ind = '',0
        else:
          print "not found",id1
      fh.close()
      Write_out(out, file)
      ### cluster file
      file_in = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+sample+".txt"
      file_out = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+sample+"_"+t+".txt"
      h=open(file_out, "w")
      fh.close()
      fh=open(file_in,"r")
      out, ind = '',0
      for l in fh:
        if(l[0]!='#'):
          l1=l
          l=l.strip().split()
          if(l[2].split("__")[0] in type_seq[t]):
            out=out+l1
            ind = ind+1
            if(ind>500):
              Write_out(out, file_out)
              out, ind = '',0
      fh.close()
      Write_out(out, file)
      out, ind = '',0
      for f in range(len(files)): 
        file_in = dir_IMGT+"IMGT_"+sample+"_"+files[f]
        file_out = dir_IMGT+"IMGT_"+sample+"_"+t+"_"+files[f]
        fh=open(file_out, "w")
        fh.close()
        fh=open(file_in,"r")
        out, ind = '',0
        for l in fh:
          l1 = l
          l=l.strip().split("\t")
          if(l[0]!='Sequence number'):
            id1 = l[1].split("__")[0]
            if(id1 in type_seq[t]):
              out = out+l1
              ind = ind+1
              if(ind>500):
                Write_out(out, file_out)
                out, ind = '',0
        fh.close()
        Write_out(out, file_out)
        out, ind = '',0
        print files[f]
  return()

def Extract_sequences_for_IMGT(id,dir):
  dir_use = dir[0]
  for s in range(len(id)):
    sample,dir_use = id[s],dir[s]
    fasta_file = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+".fasta"
    fh= open(fasta_file,"r")
    seqs =  {}
    for header,sequence in fasta_iterator(fh):
      seqs[header.split("__")[0]] = [header, sequence]
    dir_IMGT = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_ISOTYPER/"
    IMGT_file = dir_IMGT+"IMGT_"+sample.replace("PUBLIC_REDUCED_","")+"_1_Summary.txt"
    fh= open(IMGT_file,"r")
    annot1 = {}
    for l in fh:
      l=l.strip().split("\t")
      if(l[0]!='Sequence number'):
        if(l[1].split("__")[0] in seqs):
          v, j = l[3].split()[1], l[9].split()[1]
          cdr3 = l[20]
          annot1[l[1].split("__")[0]] = [v,j,cdr3]
    fh.close()
    IMGT_file = dir_IMGT+"IMGT_"+sample.replace("PUBLIC_REDUCED_","")+"_3_Nt-sequences.txt"
    fh= open(IMGT_file,"r")
    annot2 = {}
    for l in fh:
      l=l.strip().split("\t")
      if(l[0]!='Sequence number'):
        if(l[1].split("__")[0] in seqs):
          cdr3 = l[14]
          VDJ,VJ = l[6],l[7]
          if(len(VDJ)==0):VDJ = VJ
          annot2[l[1].split("__")[0]] = [VDJ, cdr3]
    fh.close()
    file_out = dir_use+"ORIENTATED_SEQUENCES/TCR_SUMMARY_FILES/TCR_sequence_annotation_"+sample+".txt"
    fh=open(file_out,"w")
    fh.close()
    out="#sample\tsequence_id\tduplicate_count\tchain\tV_gene\tJ_gene\tCDR3(nn)\tCDR3(aa)\tV(D)J_sequence\n"
    ind = 0
    for id1 in seqs: 
      freq = map(int, seqs[id1][0].split("__")[1].split("|")[0].split("_"))
      out=out+"\t".join([sample, seqs[id1][0], str(sum(freq)),annot1[id1][0] [0:3], annot1[id1][0], annot1[id1][1], annot2[id1][1],  annot1[id1][2], annot2[id1][0]])+"\n"
      ind = ind+1
      if(ind>100):
        Write_out(out, file_out)
        out, ind = '',0
    Write_out(out, file_out)
    out, ind = '',0
    print sample, "\t",len(seqs),"\t", len(annot1),"\t", len(annot2)
  return()

def Write_out(out, file):
  fh=open(file,"a")
  fh.write(out)
  fh.close()
  return()

################################################
batch_file = "/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/Samples_Caldas_BREAST_TCR1.txt"
id, dir = Get_sample_files(batch_file) # samples.post file as argument 

##### extract sequences from batched files for IMGT
Extract_sequences_for_IMGT(id,dir)


