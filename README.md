# GRVS_ASD
Contains scripts used in the study of Genome-wide rare variant score associates with morphological subtypes of autism spectrum disorder. Most of the scripts are stand-alone and sequential but a few of them are functions containing scripts, which were called by other scripts.

#Scripts included
1. Under GS_NFLD folder (to analyse the gold-standard morphological data)\n
  1.1 coding/cnvBurdenAnalysis.R = to perform rare coding CNV (>10kb) global and gene-set burden analysis\n
  1.2 coding/cnvBurdenAnalysisSmaller.R = to perform rare coding CNV (<= 10kb) global and gene-set burden analysis\n
  1.3 coding/getSNVGenesetMatrix.R = to process coding SNVs data and put them in a gene-set matrix format for a burden analysis\n
  1.4 coding/snvBurdenAnalysis.R = to perform rare coding SNV burden analysis\n
  1.5 coding/snvBurdenAnalysisDenovo.R = to perform de novo coding SNV burden analysis\n
  1.6 noncoding/burdenRequireFunctions.R = a library script contains function required for analysis\n
  1.7 noncoding/burdenTestCNV.R = burden analysis of rare noncoding CNVs (>10kb)\n
  1.8 noncoding/burdenTestSmallerCNV.R = burden analysis of rare noncoding CNVs (<= 10kb)\n
  1.9 noncoding/burdenTestSNVDenovo.R = burden analysis of de novo noncoding SNVs\n
  1.10 noncoding/burdenTestSNVRare.R = burden analysis of rare noncoding SNVs\n
  1.11 noncoding/getCNV.R = to obtain rare noncoding CNV data for burden analysis\n
  1.12 noncoding/getDataSNV.R = to obtain rare noncoding SNV data for burden analysis\n
  1.13 noncoding/getDenovoSNV.R = to obtain de novo noncoding SNV data for burden analysis\n
2. 
