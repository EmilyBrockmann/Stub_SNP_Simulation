#!/bin/bash

#path to reads
PATH_DAT=/media/rna/WHAKAHAO/Emily/Data/simulatedData/readSimulation/Try2

THREADS=6

INDV1=""
INDV2=""
for i in {50..99}; do
INDV1=$(echo "${INDV1} outputNGSNGS_Try2_Indv${i}_R1.fq")
INDV2=$(echo "${INDV2} outputNGSNGS_Try2_Indv${i}_R2.fq")

done 

cd $PATH_DAT
fastqc --threads $THREADS $INDV1 $INDV2
