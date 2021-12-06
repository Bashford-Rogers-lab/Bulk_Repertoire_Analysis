
IDS=$(awk -F '\t' "{ print \$1 }" $SAMPLES_FILE_POST)
ID=$(awk -F '\t' "{if (NR==$SGE_TASK_ID) print \$1}" $SAMPLES_FILE_POST) 

file=/gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/COMMANDLOGS/job_LEO_SEPSIS_BCR_nocombat_post.txt_2.txt

if fgrep –qx "$ID" word_list



list=$IDS 
x=$ID
for item in $list
do 
	echo $item
done

for item in $list 
do
    case "$x" in
      item1|item2)
        echo "In the list"
        ;;
      not_an_item)
        echo "Error" >&2
        exit 1
        ;;
    esac
done


WORD="$1"
for w in dog cat horse
do
  if [ "$w" == "$WORD" ] 
  then
      yes=1
      break
  fi
done;
[ "$yes" == "1" ] && echo "$WORD is in the list" || 
      
	  


if fgrep –qx "$ID" /gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/COMMANDLOGS/job_LEO_SEPSIS_BCR_nocombat_post.txt_2.txt
	then
		echo "YES"
	else 
		FAILED="echo ${ID} >> COMMANDLOGS/job_${SAMPLES_FILE_POST}_${TASK}_FAILED_SAMPLES.txt"
		eval "${FAILED}"
fi 


if fgrep –qx "$ID" /gpfs2/well/immune-rep/shared/CODE/BCR_TCR_PROCESSING_PIPELINE/COMMANDLOGS/job_LEO_SEPSIS_BCR_nocombat_post.txt_2.txt


python 

import math
import sys
from collections import defaultdict
import os
import commands
import sys
from operator import itemgetter
import networkx as nx
import re


source='/well/immune-rep/shared/MISEQ/SEPSIS_DATA_LEO/LANE_1/FASTQ/WTCHG_847696_D702501/WTCHG_847696_D702501_'


pre = source
reads1 = source+"R1_001.fastq"
reads2 = source+"R2_001.fastq"
#reads1 = source
#reads2 = source.replace("1.fq","2.fq")
reads1 = source+"1.fastq"
reads2 = source+"2.fastq"
  #reads1 = source.replace("XXXXX","1")
  #reads2 = source.replace("XXXXX","2")
 # os.system("gunzip "+reads1+".gz")
 # os.system("gunzip "+reads2+".gz")

reads1 = source+"R1_001.fastq"
reads2 = source+"R2_001.fastq"

if os.path.isfile(source+"R1_001.fastq"):
	reads1 = source+"R1_001.fastq"
	reads2 = source+"R2_001.fastq"
elif os.path.isfile(source+"R1_001.fastq.gz"):
	reads1 = source+"R1_001.fastq"
	reads2 = source+"R2_001.fastq"
	os.system("gunzip "+reads1+".gz")
	os.system("gunzip "+reads2+".gz")
elif os.path.isfile(source+"1.fastq"):
	reads1 = source+"1.fastq"
	reads2 = source+"2.fastq"
elif os.path.isfile(source+"1.fastq.gz"):
	reads1 = source+"1.fastq"
	reads2 = source+"2.fastq"
	os.system("gunzip "+reads1+".gz")
	os.system("gunzip "+reads2+".gz")
elif os.path.isfile(source.replace("XXXXX","1")):
	reads1 = source.replace("XXXXX","1")
	reads2 = source.replace("XXXXX","2")
elif os.path.isfile(source.replace("XXXXX","1")+".gz"):
	reads1 = source.replace("XXXXX","1")
	reads2 = source.replace("XXXXX","2")
	os.system("gunzip "+source.replace("XXXXX","1")+".gz")
	os.system("gunzip "+source.replace("XXXXX","2")+".gz")
elif os.path.isfile(source):
	reads1 = source
	reads2 = source.replace("1.fq","2.fq")
elif os.path.isfile(source+".gz"):
	reads1 = source
	reads2 = source.replace("1.fq","2.fq")
	os.system("gunzip "+reads1+".gz")
	os.system("gunzip "+reads2+".gz")