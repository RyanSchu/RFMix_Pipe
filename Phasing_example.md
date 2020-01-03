Here We will format data for local ancestry estimation and run RFMix


1. Process NAT genotypes from martin et al and merge them with the 1000G from PEL that have > 90% NAT ancestry

```
#first download all the nat data from alicia martin et al
wget https://personal.broadinstitute.org/armartin/tgp_admixture/nat_ref/nam_fwd_cleaned_hg19_ref.bed
wget https://personal.broadinstitute.org/armartin/tgp_admixture/nat_ref/nam_fwd_cleaned_hg19_ref.bim
wget https://personal.broadinstitute.org/armartin/tgp_admixture/nat_ref/nam_fwd_cleaned_hg19_ref.fam
wget https://personal.broadinstitute.org/armartin/tgp_admixture/AMR_lai.txt
```
2. Select inividuals from AMR_lai.txt with > 90% NAT ancestry for use as reference
```
awk '{if ( $4 >= 0.9) print $1}' AMR_lai.txt > NAT_90_perc.txt
## Exract those individulas from the 1000G population
bcftools view -S NAT_90_perc.txt 1000G.vcf -O z -o 1000G_NAT_90_perc
##can also extract them into individual chromosomes which is what I have done
for i in {1..22}
do
  bcftools view -S NAT_90_perc.txt chr.${i}.recode.vcf -O z -o 1000G_NAT_90_perc_chr${i}
done
```
3. Next merge the files with plink 

process for merging one file at a time where one is a vcf and one is bed/bim/fam
```
for i in {1..22}
do
  plink --vcf 1000G_NAT_90_perc_chr${i}.recode.bcf.gz --bmerge nat_clean_no_dups --chr ${i} --make-bed --out chr${i}.NAT.initial.merge
  # Usually there will be snps that cause the merge to fail so we must remove those
  if [ -e chr${i}.NAT.initial.merge-merge.missnp ]
  then
  
    plink --exclude chr${i}.NAT.initial.merge-merge.missnp \
      --vcf 1000G_NAT_90_perc_chr${i}.recode.bcf.gz \
      --make-bed \
      --out 1000G_exclude_missing_chr${i}
      
    plink --exclude chr${i}.NAT.initial.merge-merge.missnp \
      --bfile nat_clean_no_dups \
      --make-bed \
      --out NAT_exclude_missing_chr${i}
      
    plink --bfile NAT_exclude_missing_chr${i} \
      --bmerge 1000G_exclude_missing_chr${i} \
      --make-bed \
      --out merged_NAT_chr${i}
     
    plink --bfile merged_NAT_chr${i} \
      --geno 0.01 \
      --make-bed \
      --out merged_NAT_chr${i}_geno0.01
  else
    plink --bfile chr${i}.NAT.initial.merge \
      --geno 0.01 \
      --make-bed \
      --out merged_NAT_chr${i}_geno0.01
  fi
done
```
Once those are merged phase with shapeit

```
/home/ryan/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit \
--input-bed /home/ryan/nat_reference_martin/merged_NAT_pops/merged_NAT_chr${chr}_geno0.01.bed \
/home/ryan/nat_reference_martin/merged_NAT_pops/merged_NAT_chr${chr}_geno0.01.bim \
/home/ryan/nat_reference_martin/merged_NAT_pops/merged_NAT_chr${chr}_geno0.01.fam \
--input-map /home/ryan/1000-genomes-genetic-maps/genetic_map/genetic_map_GRCh37_chr${chr}_shapeit.txt \
--output-max /home/ryan/nat_reference_martin/merged_NAT_pops/merged_NAT_chr${chr}.phased.haps \
/home/ryan/nat_reference_martin/merged_NAT_pops/merged_NAT_chr${chr}.phased.sample

/home/ryan/software/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit \
-convert \
--input-haps /home/ryan/nat_reference_martin/merged_NAT_pops/phased/merged_NAT_chr${chr}.phased \
--output-vcf /home/ryan/nat_reference_martin/merged_NAT_pops/phased/merged_NAT_chr${chr}.phased.vcf
```
Next we will have to intersect all reference data and the genetic map

```
```
```
