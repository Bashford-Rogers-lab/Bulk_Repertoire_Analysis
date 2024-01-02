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
      id.append(l[0])
      dir.append(l[7])
  fh.close()
  return(id, dir)

def Split_functional_non_function(id,dir): 
  for s in range(len(id)):
    print "\t\t", s, id
    sample,dir_use = id[s],dir[s]
    dir_IMGT = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT/"
    fasta_file = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+".fasta"
    if(os.path.exists(fasta_file)):
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
      file_in1 = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+sample+".txt"
      file_out1 = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Cluster_identities_"+sample+"_"+t+".txt"
      h=open(file_out1, "w")
      fh.close()
      fh=open(file_in1,"r")
      out, ind = '',0
      for l in fh:
        if(l[0]!='#'):
          l1=l
          l=l.strip().split()
          if(l[2].split("__")[0] in type_seq[t]):
            out=out+l1
            ind = ind+1
            if(ind>500):
              Write_out(out, file_out1)
              out, ind = '',0
      fh.close()
      Write_out(out, file_out1)
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

def Extract_sequences_for_IMGT(id,dir,batch_name):
  dir_use = dir[0]
  extract = True
  #extract = False
  if(extract):
    for batch in batch_name: 
      IMGT_file_compressed=dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/"+batch+".txz"
      command = "tar Jxvf "+IMGT_file_compressed+" -C "+dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/"+batch+"/"
      os.system("mkdir "+dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/"+batch+"/")
      os.system(command)
  ids = {}
  samples_count = Tree()
#  id = ["BCR_2_UK02XX0029_T5"] ##### remove
  for s in range(len(id)):
    sample,dir_use = id[s],dir[s]
    fasta_file = dir_use+"ORIENTATED_SEQUENCES/NETWORKS/Fully_reduced_"+sample+".fasta"
    if(os.path.exists(fasta_file)):
      fh= open(fasta_file,"r")
      for header,sequence in fasta_iterator(fh):
        ids[header.split("__")[0]] = sample
        samples_count[sample][header.split("__")[0]].value = 1
      print "\r","read samples:",s, "n sequences:",len(ids),
    fh.close()
  for s in samples_count:
    print s, "\t",len(samples_count[s]), "unique sequences"
  dir_IMGT = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_SPLIT/"
  ## Added by LEO as caused an error
  if(os.path.exists(dir_IMGT)==False):
    os.system("mkdir "+dir_IMGT)
  ##################################
  files = ['10_V-REGION-mutation-hotspots.txt',
  '1_Summary.txt',
  '2_IMGT-gapped-nt-sequences.txt',
  '3_Nt-sequences.txt',
  '4_IMGT-gapped-AA-sequences.txt',
  '5_AA-sequences.txt',
  '6_Junction.txt',
  '7_V-REGION-mutation-and-AA-change-table.txt',
  '8_V-REGION-nt-mutation-statistics.txt',
  '9_V-REGION-AA-change-statistics.txt']
  # files = ['9_V-REGION-AA-change-statistics.txt'] #by Sakina
  for f in range(len(files)):
    for s in range(len(id)):
       out_file = dir_IMGT+"IMGT_"+id[s]+"_"+files[f]
       fh=open(out_file, "w")
       fh.close()
    info = {}
    for batch in batch_name:
      input_file = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/"+batch+"/"+files[f]
      print input_file
      fh=open(input_file,"r")
      for l in fh:
        if(l[0]!="S"):
          l1 = l.strip()
          l=l1.split()
          if(l[1].split("__")[0] in ids): 
            sam = ids[l[1].split("__")[0]]
            if(sam in info): info[sam] = info[sam]+l1+"\n"
            else:info[sam]=l1+"\n"
            if(len(info[sam])>10000):
              Write_out(info[sam], dir_IMGT+"IMGT_"+sam+"_"+files[f])
              info[sam] = ''
      fh.close()
      for sam in info: 
        Write_out(info[sam], dir_IMGT+"IMGT_"+sam+"_"+files[f])
        info[sam] = ''
    print f, files[f]
  return()

def Write_out(out, file):
  fh=open(file,"a")
  fh.write(out)
  fh.close()
  return()

################################################
args=sys.argv
batch_file = args[1]
id, dir = Get_sample_files(batch_file) # samples.post file as argument 

# Edited by LEO to run on command line without need to edit. 
output_dir = args[2]
path_to_IMGTfiles = output_dir + 'ORIENTATED_SEQUENCES/ANNOTATIONS/IMGT_RAW/'
batch_name = os.listdir(path_to_IMGTfiles) 
batch_name1 = []
for i in range(len(batch_name)): 
  if(batch_name[i].count('.txz')!=0):
    batch_name1 = batch_name1+[batch_name[i].split('.txz')[0]]

batch_name = list(set(batch_name1))
#batch_name = ["IMGT_BCR1/","IMGT_BCR2/"]
#batch_name = ["IMGT_BCR3/"]
print batch_name
##### extract sequences from batched files for IMGT
Extract_sequences_for_IMGT(id,dir,batch_name)

#### split functional and non-functional sequences + IMGT files for IsoTyper
#Split_functional_non_function(id,dir)


