#!/bin/bash 
start=$(date +%s)
PATH_REF=/media/rna/WHAKAHAO/Emily/Data/ReferenceGenomes
#PATH_REF=~/Emily/Data/test_NGSNGS/reduced_Data
PATH_SIM=/media/rna/WHAKAHAO/Emily/Data/simulatedData
#PATH_SIM=~/Emily/Data/test_NGSNGS/simulated
PATH_TOOLS=/media/rna/WHAKAHAO/Emily/Tools

LABEL=Try2

mkdir $PATH_SIM/readSimulation/${LABEL}
for indiv in {0..99}; do 
	$PATH_TOOLS/NGSNGS/ngsngs -i $PATH_REF/Genome_for_sim_data4VCF_${LABEL}.fa \
				-c 60 \
				-ld Norm,400,50 \
				-seq PE -f fq -o $PATH_SIM/readSimulation/${LABEL}/outputNGSNGS_${LABEL}_Indv${indiv} \
				-q1 $PATH_TOOLS/NGSNGS/Test_Examples/AccFreqL150R1.txt \
				-q2 $PATH_TOOLS/NGSNGS/Test_Examples/AccFreqL150R2.txt \
				-vcf $PATH_SIM/vcfSimulation/simulatedVCF_${LABEL}.vcf -chr Chr5 \
				-id ${indiv} \
				-s 1234567 \
				-ll 35 \
				-t 6
				#-cl 150 \
				# -q1 ../Test_ExamplesNGSNGS/AccFreqL150R1.txt -q2 ../Test_ExamplesNGSNGSTest_Examples/AccFreqL150R2.txt 
done;

end=$(date +%s)
runtime=$((end - start))
echo "Program took $runtime seconds"
echo -e "generateReads.sh\t${runtime}" >> $PATH_SIM/timingSimulation_${LABEL}.txt