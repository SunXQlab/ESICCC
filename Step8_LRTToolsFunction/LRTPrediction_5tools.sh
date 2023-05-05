#!/bin/bash
tools="NicheNet MISTy stMLnet HoloNet CytoTalk"

sampleID='CID44971 CID4465 control_P7 control_P8 UKF243_T_ST UKF260_T_ST UKF266_T_ST UKF334_T_ST'

for i in $sampleID;
do
fpath='./Step6_LRTPredictionResult/'$i"/"

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

if [[ !($j == 'CytoTalk') ]]
then
/usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /usr/bin/Rscript ../Script/Step8_LRTToolsFunction/RunSTScript.R $j $i human;
else
/usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /home/ljx/software/R-4.1.0/bin/Rscript ../Script/Step8_LRTToolsFunction/RunSTScript.R $j $i human;
fi

done
done
