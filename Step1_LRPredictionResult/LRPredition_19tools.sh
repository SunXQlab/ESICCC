#!/bin/bash
tools="CellPhoneDB2 CellPhoneDB3 CellTalker Connectome NATMI ICELLNET scConnect CellChat SingleCellSignalR CellCall scSeqComm NicheNet Domino PyMINEr iTALK cell2cell scMLnet"

sampleID='CID44971 CID4465 CK357 CK358 CK161 CK165 CK361 CK362 CK162 CK368 Slide14 pbmc4k pbmc6k pbmc8k"
  
for i in $sampleID;
do
if [[ "$i" =~ 'CID' ]]
then
fpath='./Step1_LRPredictionResult/NG_BC_'$i"/"
elif [[ "$i" =~ 'CK' ]]
then
fpath='./Step1_LRPredictionResult/N_MI_'$i"/"
elif [[ "$i" =~ 'Slide' ]]
then
fpath='./Step1_LRPredictionResult/MouseEmbryo_'$i"/"
else
fpath='./Step1_LRPredictionResult/GSE106487_'$i"/"
m='ST'
n='SSC'
fi

if [ ! -d $fpath ]
then
mkdir $fpath
fi
  for j in $tools;
  do
  fpathout=$fpath"/"$j"/"
  if [ ! -d $fpathout ]
  then
  mkdir $fpathout
  fi
  
  if [[ "$i" =~ 'Slide'  && !($j == 'CytoTalk') ]]
  then
  /usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /usr/bin/Rscript ../Script/Step1_LRPredictionResult/RunScript.R $j $i mouse;
  elif [[ "$i" =~ 'Slide'  && ($j == 'CytoTalk') ]]
  then
  /usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /home/ljx/software/R-4.1.0/bin/Rscript ../Script/Step1_LRPredictionResult/RunScript.R $j $i mouse;
  elif [[ ($j == 'CytoTalk') ]]
  then
  /usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /home/ljx/software/R-4.1.0/bin/Rscript ../Script/Step1_LRPredictionResult/RunScript.R $j $i human;
  else
  /usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /usr/bin/Rscript ../Script/Step1_LRPredictionResult/RunScript.R $j $i human;
  fi

  done
done
