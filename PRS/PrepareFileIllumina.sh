
../plink_linux_x86_64/plink --const-fid --make-bed --out bgi --real-ref-alleles --recode vcf-iid --snps-only --vcf /hpf/largeprojects/tcagstor/users/btg1/dropbox/ada/wgs_ada/extract_merge_snps/BGI/merge.all.geno.id.norm.vcf.gz --vcf-half-call m

../plink_linux_x86_64/plink --const-fid --make-bed --out macrogen --real-ref-alleles --recode vcf-iid --snps-only --vcf /hpf/largeprojects/tcagstor/users/btg1/dropbox/ada/wgs_ada/extract_merge_snps/Macrogen/merge.all.geno.id.norm.vcf.gz --vcf-half-call m

#../plink_linux_x86_64/plink --const-fid --make-bed --out TCAG --real-ref-alleles --recode vcf-iid --snps-only --vcf /hpf/largeprojects/tcagstor/users/btg1/dropbox/ada/wgs_ada/extract_merge_snps/TCAG/merge.all.geno.id.norm.vcf.gz --vcf-half-call m

../plink_linux_x86_64/plink --make-bed --out BGI.0.9 --bfile bgi --geno 0.1 --real-ref-alleles
../plink_linux_x86_64/plink --make-bed --out Macrogen.0.9 --bfile macrogen --geno 0.1 --real-ref-alleles
../plink_linux_x86_64/plink --make-bed --out TCAG.0.9 --bfile TCAG --geno 0.1 --real-ref-alleles

Rscript getBedLocation.R

../plink_linux_x86_64/plink --bfile BGI.0.9 --extract common.snps.txt --make-bed --out BGI.filtered --real-ref-alleles --recode vcf
../plink_linux_x86_64/plink --bfile Macrogen.0.9 --extract common.snps.txt --make-bed --out Macrogen.filtered --real-ref-alleles --recode vcf
../plink_linux_x86_64/plink --bfile TCAG.0.9 --extract common.snps.txt --make-bed --out TCAG.filtered --real-ref-alleles --recode vcf

bgzip Macrogen.filtered.vcf -c > Macrogen.filtered.vcf.gz 
bgzip BGI.filtered.vcf -c > BGI.filtered.vcf.gz 
bgzip TCAG.filtered.vcf -c > TCAG.filtered.vcf.gz 

tabix -p vcf Macrogen.filtered.vcf.gz 
tabix -p vcf BGI.filtered.vcf.gz 
tabix -p vcf TCAG.filtered.vcf.gz 

vcfFiles="*filtered.vcf.gz"
bcftools merge $vcfFiles --output "merge.nfld.vcf"

../plink_linux_x86_64/plink --const-fid --make-bed --out NFLD --real-ref-alleles --recode vcf-iid --snps-only --vcf merge.nfld.vcf --vcf-half-call m --maf 0.05

rm TCAG.*
rm Macrogen.*
rm BGI.*
rm bgi.*
rm macrogen.*
rm merge.nfld.vcf  

Rscript fixFam.R  
