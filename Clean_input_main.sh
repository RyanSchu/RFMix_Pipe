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

bcftools query -f '%ID %CHROM %POS\n' '${query}' -o ${outDir}_query.snps
bcftools query -f '%ID %CHROM %POS\n' '${ref} -o ${outDir}_ref.snps

Rscript Clean_input_intersect.R --query ${outDir}_query.snps --ref ${outDir}_ref.snps --map '${map}' --out ${outDir}



