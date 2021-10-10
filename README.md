# GRVS_ASD
Contains scripts used in the study of Genome-wide rare variant score associates with morphological subtypes of autism spectrum disorder. Most of the scripts are stand-alone and sequential but a few of them are functions containing scripts, which were called by other scripts.

#Scripts included
1. Under GS_NFLD folder (to analyse the gold-standard morphological data)<br/>
  1.1 coding/cnvBurdenAnalysis.R = to perform rare coding CNV (>10kb) global and gene-set burden analysis<br/>
  1.2 coding/cnvBurdenAnalysisSmaller.R = to perform rare coding CNV (<= 10kb) global and gene-set burden analysis<br/>
  1.3 coding/getSNVGenesetMatrix.R = to process coding SNVs data and put them in a gene-set matrix format for a burden analysis<br/>
  1.4 coding/snvBurdenAnalysis.R = to perform rare coding SNV burden analysis<br/>
  1.5 coding/snvBurdenAnalysisDenovo.R = to perform de novo coding SNV burden analysis<br/>
  1.6 noncoding/burdenRequireFunctions.R = a library script contains function required for analysis<br/>
  1.7 noncoding/burdenTestCNV.R = burden analysis of rare noncoding CNVs (>10kb)<br/>
  1.8 noncoding/burdenTestSmallerCNV.R = burden analysis of rare noncoding CNVs (<= 10kb)<br/>
  1.9 noncoding/burdenTestSNVDenovo.R = burden analysis of de novo noncoding SNVs<br/>
  1.10 noncoding/burdenTestSNVRare.R = burden analysis of rare noncoding SNVs<br/>
  1.11 noncoding/getCNV.R = to obtain rare noncoding CNV data for burden analysis<br/>
  1.12 noncoding/getDataSNV.R = to obtain rare noncoding SNV data for burden analysis<br/>
  1.13 noncoding/getDenovoSNV.R = to obtain de novo noncoding SNV data for burden analysis<br/>
  1.14 getCrossValidatePredictiveScore.R = to get and analyse GRVS <br/>
  1.15 getCrossValidatePredictiveScoreWithCRV.R = get GRVS for CRVs <br/><br/>
2. Under ADM_NFLD folder (to analyse the ADM morphological data)<br/>
  2.1 getGenesetCoeffADM.R = to get coefficient of variants to be used in calculating GRVS in an indenpendent cohort (SSC)<br/>
  2.2 getCrossValidatePredictiveADM.R = to get and analyse GRVS based on ADM classification data<br/><br/>
3. Under ADM_SSC folder (to analyse the ADM morphological data of SSC cohort)<br/>
  3.1 filterSNVs.R = to obtain and process rare SNVs data of SSC<br/>
  3.2 getCNVGeneset.R = to transform rare CNVs (>10kb) data in a format suitable for burden analysis<br/>
  3.3 getDenovoSNVGS.R = to transform de novo SNVs data in a format suitable for burden analysis<br/>
  3.4 getNC.R = to obtain noncoding features for SNVs data<br/>
  3.5 getSmallCNVGeneset.R = to transform rare CNVs (<10kb) data in a format suitable for burden analysis<br/>
  3.6 getSNVGeneset.R =  to transform rare SNVs data in a format suitable for burden analysis<br/>
  3.7 
