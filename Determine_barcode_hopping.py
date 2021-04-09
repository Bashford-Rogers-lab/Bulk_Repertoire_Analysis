#!/usr/bin/python
#import math
import sys
import os
import commands
import numpy as np

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
  id, dir, batch = [],[],{}
  fh=open(batch_file,"r")
  ind = 0
  for l in fh:
    ind = ind+1
    if(l[0]!="#" and len(l)>3):
      l=l.strip().split()
      id.append(l[0])
      dir.append(l[7])
      bat = l[1]
      if(bat in batch):batch[bat] = batch[bat]+[[l[0],l[7]]]
      else:batch[bat]=[[l[0],l[7]]]
  fh.close()
  return(id, dir,batch)

def Get_ID_overlap(batch,batch_file):
  for b  in batch: 
    for b1 in batch[b]: 
      dir = b1[1]
      break
    break
  fh=open(dir+"ORIENTATED_SEQUENCES/Barcode_hopping_"+batch_file.replace("Samples_",""),"w")
  fh.close()
  for s in batch:
    seqs = {}
    bat = batch[s]
    for b in bat: 
      sample,dir_use = b[0],b[1]
      fasta_file = dir_use+"ORIENTATED_SEQUENCES/TMP/Untrimmed_"+sample+".fasta"
      fh= open(fasta_file,"r")
      for header,sequence in fasta_iterator(fh):
        if(header in seqs):seqs[header] = seqs[header] +[sample]
        else:seqs[header] =[sample]
      fh.close()
    overlap = {}
    totals = {}
    for id in seqs:
      if(len(seqs[id])>1):
        samples = []
        for sample in seqs[id]:
          samples.append(sample)
        name = ",".join(samples)
        if(name in overlap): overlap[name] =  overlap[name] +1
        else: overlap[name]  = 1
    #for name in  overlap: 
    #  print bat, name,  overlap[name] 

  #Write_out(out, file_out_pre+str(batch)+".fasta")
  #out, ind = '',0
  #print "scp -p mfj169@rescomp1.well.ox.ac.uk:"+file_out_pre+"* ./"
  return()

def Write_out(out, file):
  fh=open(file,"a")
  fh.write(out)
  fh.close()
  return()

def Get_isotype_depth(id,dir):
  isotypes_uniq,isotypes_total = {},{}
  for s in range(len(id)):
    sample,dir_use = id[s],dir[s]
    isotype_file = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/IsoTyper_chain_repertoire_statistics_file_"+sample+".txt"
    fh=open(isotype_file,"r")
    for l in fh:
      if(l[0]!="#"):
        l=l.strip().split()
        iso, umis, uniqs = l[1],int(l[2]), int(l[3])
        if(iso in isotypes_uniq):
          isotypes_uniq[iso], isotypes_total[iso] = isotypes_uniq[iso]+[uniqs], isotypes_total[iso] +[umis]
        else:isotypes_uniq[iso], isotypes_total[iso] = [uniqs],[umis]
    fh.close()
  out = "#isotype\ttype\tmin\t5th percentile\t10th percentile\t20th percentile\n"
  for iso in isotypes_uniq:
    quantiles = [np.percentile(isotypes_uniq[iso],5), np.percentile(isotypes_uniq[iso],10),np.percentile(isotypes_uniq[iso],20)]
    out = out+"\t".join(map(str, [iso,"UNIQ",min(isotypes_uniq[iso])]+quantiles))+"\n"
  for iso in isotypes_total:
    quantiles = [np.percentile(isotypes_total[iso],5), np.percentile(isotypes_total[iso],10),np.percentile(isotypes_total[iso],20)]
    out = out+"\t".join(map(str, [iso,"TOTAL",min(isotypes_total[iso])]+quantiles))+"\n"
  out_file = dir_use+"ORIENTATED_SEQUENCES/ANNOTATIONS/Sampling_depth_per_isotype_"+batch_file.replace(".txt","")+".txt"
  print out_file
  fh=open(out_file, "w")
  fh.write(out)
  fh.close()
  return()

################################################
#args=sys.argv
#batch_file = args[1]
batch_file = "Samples_COMBAT_BCR_post1.txt"

id, dir,batch = Get_sample_files(batch_file)

Get_ID_overlap(batch,batch_file)




