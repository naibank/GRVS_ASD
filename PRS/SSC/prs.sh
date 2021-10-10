#!/bin/bash -x

#PBS -N ssc.prs
#PBS -l nodes=1:ppn=1
#PBS -l vmem=32g
#PBS -d /hpf/largeprojects/tcagstor/users/worrawat/NFLD/SSCPRS
#PBS -l walltime=24:00:00
#PBS -e /hpf/largeprojects/tcagstor/users/worrawat/NFLD/SSCPRS/prs.error
#PBS -o /hpf/largeprojects/tcagstor/users/worrawat/NFLD/SSCPRS/prs.stdout

../AdaPRS/PRSice_linux/PRSice_linux --A1 A1 \
    --A2 A2 \
    --bar-levels 0.001,0.05,0.1 \
    --base "iPSYCH-PGC_ASD_Nov2017_P0.1_INFO0.9_SNPs_clean_hg38_for_SSC.tsv"  \
    --binary-target T \
    --bp BP \
    --chr CHR \
    --clump-kb 250 \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --interval 0.10 \
    --lower 0.000 \
    --model add \
    --out SSC.PRS \
    --pvalue P \
    --seed 335560616 \
    --snp SNP \
    --stat OR \
    --target "plinkfiles/SSC_chr1_22_Filter_Pass_iPSYCH_SNPs_7261_Subjects" \
    --thread 1 \
    --upper 0.500000 \
    --all-score
