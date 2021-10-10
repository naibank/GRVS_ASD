#!/bin/bash -x

#PBS -N PRS
#PBS -l nodes=1:ppn=1
#PBS -l vmem=32g
#PBS -d /hpf/largeprojects/tcagstor/users/worrawat/AdaPRS
#PBS -l walltime=24:00:00
#PBS -e /hpf/largeprojects/tcagstor/users/worrawat/AdaPRS/prs.error
#PBS -o /hpf/largeprojects/tcagstor/users/worrawat/AdaPRS/prs.stdout

./PRSice_linux/PRSice_linux --A1 A1 \
    --A2 A2 \
    --bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
    --base "iPSYCH-ASD_May2021_P0.1_INFO0.9_SNPs_clean_hg19_for_NFLD.tsv"  \
    --binary-target T \
    --bp BP \
    --chr CHR \
    --clump-kb 250 \
    --clump-p 1.000000 \
    --clump-r2 0.100000 \
    --base-info INFO,0.9 \
    --interval 0.10 \
    --lower 0.000 \
    --model add \
    --out ASD.PRS.2021 \
    --pvalue P \
    --seed 335560616 \
    --snp SNP \
    --stat OR \
    --target "bed/NFLD" \
    --thread 1 \
    --upper 0.500000 \
    --all-score