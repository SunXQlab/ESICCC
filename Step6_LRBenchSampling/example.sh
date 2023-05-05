#!/bin/bash
tools="CellPhoneDB2 CellPhoneDB3 CellTalker Connectome NATMI ICELLNET scConnect CellChat SingleCellSignalR CellCall scSeqComm NicheNet Domino PyMINEr iTALK scMLnet"

sampleID='CID44971 CID4465 CK357 CK358 CK161 CK165 CK361 CK362 CK162 CK368 Slide14 pbmc4k pbmc6k pbmc8k'
ratios='90 80 70 60 50'

for i in $sampleID;
do
for ratio in $ratios;
do

fpath='./Step10_LRBenchSamplingResult/'$i"_"$ratio"/"

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
/usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /usr/bin/Rscript ../Script/Step10_LRBenchSampling/RunScript.R $j $i mouse $ratio;
elif [[ "$i" =~ 'Slide'  && ($j == 'CytoTalk') ]]
then
/usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /home/ljx/software/R-4.1.0/bin/Rscript ../Script/Step10_LRBenchSampling/RunScript.R $j $i mouse $ratio;
elif [[ ($j == 'CytoTalk') ]]
then
/usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /home/ljx/software/R-4.1.0/bin/Rscript ../Script/Step10_LRBenchSampling/RunScript.R $j $i human $ratio;
else
/usr/bin/time -v -o $fpathout"TimeMemRecord.txt" /usr/bin/Rscript ../Script/Step10_LRBenchSampling/RunScript.R $j $i human $ratio;
fi

done
done
done