# Subsampling_IDs_V.V.py developed by Rachael Bashford-Rogers (2020)
# at the University of Oxford and Univeristy of Cambridge
# E-mail: rbr1@well.ox.ac.uk

# If you find the methods in ImmuneReceptor_PROCESSING_PIPELINE, please cite the following reference:
# Bashford-Rogers, R. et al. Nature 2019 (https://www.nature.com/articles/s41586-019-1595-3.pdf)

# Copyright (C) 2022  Dr Rachael Bashford-Rogers

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#!/usr/bin/python
import math
import sys
from collections import defaultdict
import os
import commands
import sys
from operator import itemgetter
import networkx as nx
import re
import numpy as np

def Generate_network(edge_file,file_vertex):
  G=nx.Graph()
  ids = {}
  G.rtt={}
  fh=open(file_vertex,"r")
  for l in fh:
    id = l.strip().split()[0].replace(">","")
    G.add_node(id)
    G.rtt[id] = 1
    ids[id] = 1
  fh.close()
  fh=open(edge_file,"r")
  for l in fh:
    l=l.strip().split()
    if (l[0] in ids and l[1] in ids):
      G.add_edge(l[0],l[1])
  fh.close()
  return(G)

def Unique(lst):
  uniq = Uniq(lst)
  count=[]
  for i in range(0,len(uniq)):
    count.append(len(uniq[i]))
  return(uniq,count)

def Uniq(v):
  C=set(v)
  return list(C)

def Subsampling_CC(G, cluster,ids,sample_depth,csizes,type,cluster_all,sample):
  subsample_file = dir+"NETWORKS/Subsample_IDs_"+type+"_"+sample+".txt"
  subsample_stats = dir+"NETWORKS/Subsample_stats_"+type+"_"+sample+".txt"
  sample_name = type
  print "Done clustering"
  repeats = 2
  l1 = len(cluster)
  out="#sample\tmethod\tsample depth\tmax cluster %\n"
  metrix_real = max(csizes)*100.0/sum(csizes)
  out=out+"\t".join(map(str, [sample_name, "REAL",len(csizes), max(csizes)*100.0/sum(csizes)]))+"\n"
  X, metric = [],[]
  for r in range(0,repeats):
    print r
    rand = np.random.choice(xrange(l1), l1,replace = False)
    for i in range(sample_depth,len(rand)):
      l = len(np.unique(cluster[rand[0:i]]))
      if(l>=sample_depth):
        break
    cluster_sub = cluster[rand[0:i]]
    ids_sub = ids[rand[0:i]]
    unique, counts = Unique(cluster_sub)
    print unique, counts
    #unique, counts = np.unique(cluster_sub, return_counts=True)
    nz = [i for i in range(len(counts)) if counts[i]>2]
    for i1 in nz:
      w = np.where(cluster_sub == unique[i1])[0]
      cluster_ids =ids_sub[w]
      ids_original_cluster = np.array(cluster_all[unique[i1]])
      if(len(cluster_ids)<len(ids_original_cluster)):
        for i in cluster_ids:
          if(i not in ids_original_cluster):
            print "ERROR",unique[i1], i
        degrees = list(G.degree(ids_original_cluster))
        degs,new_ids = np.array([d for n, d in degrees]), [n for n, d in degrees]
        new_ids = new_ids[np.where(degs == max(degs))[0][0]]
        nn_possible = []
        sample_id = new_ids
        for i in range(1,len(cluster_ids)):
          nn_possible = np.append(nn_possible,list(G.neighbors(sample_id)))
          nn_possible1 = [x for x in nn_possible if x not in set(new_ids)]
          nn_possible = nn_possible1
          if(len(nn_possible)>0):
            sample_id = np.random.choice(nn_possible, 1,replace = False)[0]
            new_ids  = np.append(new_ids, sample_id)
            if(len(new_ids)>=len(cluster_ids)):
              break
          else:
            nn_possible = np.setdiff1d(ids_original_cluster, new_ids)
            print "ERROR",unique[i1], nn_possible
        ids_sub[w] = new_ids
    H = G.subgraph(ids_sub)
    csizes_sub = []
    con= nx.connected_components(H)
    for i in con:
      cs = 0
      for j in i:
        cs = cs+G.rtt[j]
      csizes_sub.append(cs)
    metrix =  max(csizes_sub)*100.0/sum(csizes_sub)
    out = out+"\t".join(map(str, [sample_name,type,sample_depth,metrix]))+"\n"
    X.append(H)
    metric.append(metrix)
    print r
  Print_out_metrix(out, subsample_stats, metric, metrix_real,subsample_file, type, X)
  return()


def Subsampling_CLONE_DEPTH(G, cluster,ids,sample_depth,subsample_file, subsample_stats,sample_name,csizes,type,cluster_all):
  print "Done clustering CLONE DEPTH"
  repeats = 50
  l1 = len(cluster)
  out="#sample\tmethod\tsample depth\tmax cluster %\n"
  metrix_real = max(csizes)*100.0/sum(csizes)
  X, metric = [],[]
  clusters_unique = []
  max_clone,max_size,t= '',0,0
  for c in cluster_all:
    clusters_unique.append(c)
    t = t+len(cluster_all[c])
    if(len(cluster_all[c])>max_size):
      max_clone,max_size = c, len(cluster_all[c])
  metrix_real = max_size*100.0/t
  out=out+"\t".join(map(str, [sample_name, "REAL",len(csizes), max_size*100.0/t]))+"\n"
  print max_clone,max_size, len(cluster_all[c])
  for r in range(0,repeats):
    rand = np.random.choice(clusters_unique,sample_depth,replace = False)
    if(max_clone not in rand):
      rand = np.append(rand,max_clone)
    ids_sub = []
    for c in rand:
      ids_sub = np.append(ids_sub, cluster_all[c])
    H = G.subgraph(ids_sub)
    csizes_sub = []
    con= nx.connected_components(H)
    for i in con:
      csizes_sub.append(len(i))
    metrix =  max(csizes_sub)*100.0/sum(csizes_sub)
    out = out+"\t".join(map(str, [sample_name, type,sample_depth,metrix]))+"\n"
    X.append(H)
    metric.append(metrix)
    print r
  metrix_real = metrix
  Print_out_metrix(out, subsample_stats, metric, metrix_real,subsample_file, type, X)
  return()

def Subsampling_reads(G, cluster,ids,sample_depth,csizes,type,cluster_all,sample):
  subsample_file = dir+"NETWORKS/Subsample_IDs_"+type+"_"+sample+".txt"
  subsample_stats = dir+"NETWORKS/Subsample_stats_"+type+"_"+sample+".txt"
  print "Done clustering"
  repeats = 20
  l1 = len(cluster)
  out="#sample\tmethod\tsample depth\tmax cluster %\n"
  metrix_real = max(csizes)*100.0/sum(csizes)
  out=out+"\t".join(map(str, [sample, "REAL",len(csizes), max(csizes)*100.0/sum(csizes)]))+"\n"
  X, metric = [],[]
  for r in range(0,repeats):
    rand = np.random.choice(xrange(l1), l1,replace = False)
    prev = sample_depth
    ids_sub= ids[rand[0:sample_depth]]
    uniq = Uniq(ids_sub)
    while(len(uniq)<sample_depth):
      diff = sample_depth - len(uniq)
      ld = prev+diff
      ids_sub= ids[rand[0:ld]]
      uniq = Uniq(ids_sub)
      prev = ld
      if(len(uniq)==sample_depth):
        break
      if(len(uniq)>sample_depth):
        diff =  len(uniq) - sample_depth
        break
    H = G.subgraph(ids_sub)
    n_con= nx.number_connected_components(H)
    if(n_con>=sample_depth):
      break
    csizes_sub = []
    con= nx.connected_components(H)
    for i in con:
      csizes_sub.append(len(i))
    metrix =  max(csizes_sub)*100.0/sum(csizes_sub)
    out = out+"\t".join(map(str, [sample,type,sample_depth,metrix]))+"\n"
    X.append(H)
    metric.append(metrix)
    print r
  Print_out_metrix(out, subsample_stats, metric, metrix_real,subsample_file, type, X)
  return()

def Subsampling_CC_fast(G, cluster,ids,sample_depth,csizes,type,cluster_all,sample):
  subsample_file = dir+"NETWORKS/Subsample_IDs_"+type+"_"+sample+".txt"
  subsample_stats = dir+"NETWORKS/Subsample_stats_"+type+"_"+sample+".txt"
  sample_name = type
  print "Done clustering"
  repeats = 50
  l1 = len(cluster)
  out="#sample\tmethod\tsample depth\tmax cluster %\n"
  metrix_real = max(csizes)*100.0/sum(csizes)
  out=out+"\t".join(map(str, [sample_name, "REAL",sum(csizes), max(csizes)*100.0/sum(csizes)]))+"\n"
  X, metric = [],[]
  N_ids = len(ids)
  cluster, ids = np.array(cluster), np.array(ids)
  for r in range(0,repeats):
    print r
    #rand = [int(N_ids* random.random()) for i in range(sample_depth)]
    rand = np.random.choice(xrange(N_ids), sample_depth,replace = False)
    #for i in range(sample_depth,len(rand)):
    #  l = len(np.unique(cluster[rand[0:i]]))
    #  if(l>=sample_depth):
    #    break
    cluster_sub = cluster[rand]
    ids_sub = ids[rand]
    unique, counts = Unique(cluster_sub)
    nz = [i for i in range(len(counts)) if counts[i]>3]
    from_ids, to_ids = [],[]
    for i1 in nz:
      w = np.where(cluster == unique[i1])[0]
      if(len(w)>counts[i1]): #not all ids in cluster sampled
        cluster_ids =ids[w]
        degrees = list(G.degree(cluster_ids))
        degs = [[n,int(d)] for n, d in degrees]
        degs.sort(key=lambda x: x[1], reverse = True)
        if(len(w)>counts[i1]*0.5): # faster by subtracting lower degree IDs
          n_remove = len(w)-counts[i1]
          nsub = len(w)
          id_remove = [degs[nsub-i-1][0] for i in range(n_remove)]
          ids_keep = list(set(cluster_ids) - set(id_remove))
          from_ids = from_ids+id_remove
          to_ids = to_ids+ids_keep
        else:
          id_keep = [degs[i][0] for i in range(counts[i1])]
          id_remove = list(set(cluster_ids) - set(id_keep))
          from_ids = from_ids+id_remove
          to_ids = to_ids+ids_keep
    #### change old IDs to new IDs
    ids_sub1 = ids_sub.tolist()+to_ids
    ids_sub1= list(set(ids_sub1) - set(from_ids))
    print len(ids_sub1), len(ids_sub)
    H = G.subgraph(ids_sub1)
    csizes_sub = []
    con= nx.connected_components(H)
    for i in con:
      cs = [G.rtt[j] for j in i]
      csizes_sub.append(sum(cs))
    metrix =  max(csizes_sub)*100.0/sum(csizes_sub)
    out = out+"\t".join(map(str, [sample_name,type,sample_depth,metrix]))+"\n"
    X.append(H)
    metric.append(metrix)
    print r
  Print_out_metrix(out, subsample_stats, metric, metrix_real,subsample_file, type, X)
  print out
  return()

def Subsample_network(file_edges, file_vertex, sample,sample_depth):
  G = Generate_network(file_edges, file_vertex) 
  con= nx.connected_components(G)
  print "Generated network"
  cluster,ids,cluster_all, c,csizes = [],[],{},0,[]
  cluster_sizes = {}
  for i in con:
    c = c+1
    c1 = "C"+str(c)
    cluster_sizes[c1] = len(i)
    ids1 = [j for j in i]
    ids = ids+ids1
    cluster = cluster+[c1]*len(ids1)
    cs = [G.rtt[j] for j in i]
    csizes.append(sum(cs))
    for j in i:
      freq = map(int, j.split("__")[1].split("|")[0].split("_"))
      freq = sum(freq)
      cluster, ids = cluster+[c1]*freq, ids+[j]*freq
      if(c1 in cluster_all):cluster_all[c1] = cluster_all[c1] +[j]
      else:cluster_all[c1] = [j]
  cluster = np.array(cluster)
  ids=np.array(ids)
  print len(csizes),sample_depth, len(ids)
  if(len(ids)>sample_depth):
    print len(csizes), len(ids)
    #Subsampling_CC(G, cluster,ids,sample_depth,csizes,"CC",cluster_all,sample)
    Subsampling_CC_fast(G, cluster,ids,sample_depth,csizes,"CCFAST",cluster_all,sample)
    #Subsampling_reads(G, cluster,ids,sample_depth,csizes,"READS",cluster_all,sample)
    #Subsampling_CLONE_DEPTH(G, cluster,ids,sample_depth,subsample_file, subsample_stats,sample_name,csizes,"CLONES",cluster_all)
    #Subsampling_weighted_CC(G, cluster,ids,sample_depth,subsample_file, cluster_size_file,sample_name,csizes,"WEIGHTED_CC",cluster_all,cluster_sizes)
  return()

def Print_out_metrix(out, subsample_stats, metric, metrix_real,subsample_file, type, X):
  print out
  fh=open(subsample_stats,"w")
  fh.write(out)
  fh.close()
  metric = ((np.array(metric) - metrix_real)**2)**0.5
  min_value, min_index = 0,-1
  min_value1, min_index1 = 0,-1
  for i in range(0,len(metric)):
    if(min_value<=metric[i]):
      min_value1, min_index1 = metric[i], i
      if(metric[i]<=50 ):
        min_value, min_index = metric[i], i
  w = min_index
  if(w ==-1):w = min_value1
  #w1 = np.where(((metric-metrix_real)**2)**0.5 <50)[0]
  #w2 =np.where(metric*1.2>metrix_real)[0]
  #w2 = w1
  #y = np.median(metric[np.intersect1d(w1,w2)])
  #diff = (metric-y)**2 #np.absolute(metric-y)
  #w = np.where(diff==min(diff))[0][0]
  print "TRUE METRIC",metrix_real, "SUBSAMPLED METRIC",metric[w]
  H = X[w]
  out="#Id:rep"+str(w)+"\n"
  for id in H:
    out = out+id+"\n"
  fh=open(subsample_file.replace("XXXTYPEXXX",type),"w")
  fh.write(out)
  fh.close()
  return()


###########################
if(len(sys.argv)!=4):
  print "\nUsage:\npython Subsample_IDs.py <output directory> <sample ID> <input file> <Subsample depth for plotting> "
  print "\nNote:\n\tPlease check line XXX location of directory with reference libraries\nTry 8,000 depth in first instance\n\n"
  quit()

dir = sys.argv[1]
id = sys.argv[2]
sample_depth = int(sys.argv[3]) # default try 8000
########################### Files for network analysis
sequences_raw = dir+"Raw_sequences_"+id+".txt"
sequence_file = dir+"Original_sequences_"+id+".fasta"
att_file = dir+"NETWORKS/Vertex_relations_"+id+".txt"
file_vertex = dir+"NETWORKS/Att_"+id+".txt"
file_edges = dir+"NETWORKS/Edges_"+id+".txt"
cluster_file = dir+"NETWORKS/Cluster_identities_"+id+".txt"
Reduced_file = dir+"NETWORKS/Fully_reduced_"+id+".fasta"
checked_edges = dir+"NETWORKS/Checked_edges_"+id+".txt"
plot_ids_file = dir+"NETWORKS/Plot_ids_"+id+".txt"
file_seqs = dir+"NETWORKS/Sequences_"+id+".txt"
tmp_file0 = dir+"NETWORKS/Decon_0_"+id+".txt"
tmp_reduced_sequences = dir+"NETWORKS/Primary_reduced_sequences"+id+".fasta"
tmp_pre = dir+"NETWORKS/Pre_tmp_"+id
tmp_file1 = dir+"NETWORKS/NN_Tmp_cluster_"+id+".1"
edge_lengths=0.85
tmp_file = dir+"NETWORKS/NN_Tmp_cluster_"+id+"."
read_number_division = "__"
######################### Commands

### Subsampling BCs for plotting
Subsample_network(file_edges, file_vertex, id,sample_depth)
  

