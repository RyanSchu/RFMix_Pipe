#!/bin/bash
#Clean Local Ancestry input to contain concordant snps


while :
do
    case "$1" in
      --ref) #reference population for LAI, all in one vcf
                ref="$2"
                shift 2
                ;;
      --query) #Experimental population with unknown LA to be inferred
                query="$2"
                shift 2
                ;;
      --pop) #Reference population ids. First column sample ids. Second column population names.
                pop="$2"
                shift 2
                ;;
      --map) #population code file. First column sample ids. Second column population codes
                map="$2"
                shift 2
                ;;
      --out) #as in plink prefix of output including optional path specification
                outDir="$2"
                shift 2
                ;;
      -*) #unknown
                echo "Error: Unknown option: $1" >&2
                echo "./vcf_lift.sh --help or ./vcf_lift.sh -h for option help"
                exit 1
                ;;
      *)  # No more options
                shift
                break
                ;;
     esac
done

rm -f ${outDir}*

if [ ! -f ${ref}.tbi ]
then
	echo "Indexing Reference"
	/usr/local/bin/tabix ${ref}
fi
if [ ! -f ${query}.tbi ]
then
        echo "Indexing Query"
	/usr/local/bin/tabix ${query}
fi 

echo "finding snp intersection subset"
/usr/local/bin/bcftools query -f '%ID %CHROM %POS\n' ${query} -o ${outDir}_query.snps
/usr/local/bin/bcftools query -f '%ID %CHROM %POS\n' ${ref} -o ${outDir}_ref.snps

Rscript /home/ryan/software/RFMix_Pipe/Clean_input_intersect.R --query ${outDir}_query.snps --ref ${outDir}_ref.snps --map ${map} --out ${outDir}
awk '{print $4}' ${outDir}_intersect_snp.list.txt > ${outDir}.pos

echo "Subsetting Files"
/usr/local/bin/bcftools view -R <( cut -f 2,3 -d$'\t' ${outDir}_intersect_snp.list.txt ) -v snps ${ref} -o ${outDir}_ref.intersect.vcf.gz -O z
/usr/local/bin/bcftools view -R <( cut -f 2,3 -d$'\t' ${outDir}_intersect_snp.list.txt ) -v snps ${query} -o ${outDir}_query.intersect.vcf.gz -O z
echo "Sorting"
/usr/local/bin/bcftools sort ${outDir}_ref.intersect.vcf.gz -o ${outDir}_ref.intersect_sorted.vcf.gz -O z
/usr/local/bin/bcftools sort ${outDir}_query.intersect.vcf.gz -o ${outDir}_query.intersect_sorted.vcf.gz -O z
/usr/local/bin/tabix ${outDir}_ref.intersect_sorted.vcf.gz
/usr/local/bin/tabix ${outDir}_query.intersect_sorted.vcf.gz
echo "Merging"
/usr/local/bin/bcftools merge ${outDir}_ref.intersect_sorted.vcf.gz ${outDir}_query.intersect_sorted.vcf.gz -o ${outDir}_merged.vcf.gz -O z
/usr/local/bin/bcftools view ${outDir}_merged.vcf.gz -U -e 'GT[*] = "mis"' -o ${outDir}_merged_exclude_uncalled.vcf.gz -O z
echo "making haplotype file"
/usr/local/bin/bcftools convert --hapsample --vcf-ids ${outDir}_merged_exclude_uncalled.vcf.gz -o ${outDir}_merged.haps
zcat  ${outDir}_merged.haps.hap.gz | awk '{ $1=""; $2=""; $3=""; $4=""; $5=""; print}' | sed 's/\s//g' | sed s/\\*//g > ${outDir}_merged.haps

echo "making class file"
/usr/local/bin/bcftools query -l ${outDir}_merged.vcf.gz > ${outDir}_sample_list.txt
Rscript /home/ryan/software/RFMix_Pipe/Clean_input_make_classes.R --pop ${pop} --samples ${outDir}_sample_list.txt --out ${outDir}.classes

