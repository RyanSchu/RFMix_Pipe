# RFMix_Pipe
Local ancestry pipeline for running RFMix

Contained in this repository is a handful of scripts to create the files necessary to run RFMix from VCF.

`Clean_input_main.sh` is the main script while `Clean_input_intersect.R` and `Clean_input_make_classes.R` are helper scripts.
All are designed to be run from the command line, an example of which can be found in the file labelled wrap `wrap`. This wrapper file can be modified with the users files to easily create the necessary inputs

```
nano wrap ###modify with your files, see inputs below
./wrap
```

### Example Data
Reference files were sourced from 1000 Genomes project YRI and CEU populations. 

Example query data can be taken from the 1000 Genomes ASW population

Map files were sourced from: https://github.com/joepickrell/1000-genomes-genetic-maps

### Inputs

There are five main inputs. 

```
--query <query.vcf> VCF file containing the admixed experimental population on which local ancestry inferrence is to be conducted.
--ref <reference.vcf> VCF file that contains all samples to be used as reference populations
--pop <pop_codes.txt> A two column tab delimeted file that contains the reference samples in column 1 and their respective populations in column 2. Should not contain the admixed sample IDs as those will be inferred from <query.vcf>.
--map <genetic_map.txt> A three column text file similar to the genetic map used by SHAPEIT. Column 1 is SNP ID, column 2 is base pair position, and column 3 is centiMorgan position.
--out <outdir/output_prefix> As in plink, out specifies the both the output directory as well as the output prefix for files to be output as.
```

All inputs should be separated by chromosome (barring `pop` which is the same across chromosomes). 
Make sure all inputs are in the same build and that their SNP ID formats are concordant. 

### Workflow

0. Phase data if not phased. See Phasing_example.md

1.Run `Clean_input_main.sh`
* Remove snps containing duplicate positions/IDs
* Find the subset of snps that exist across all files and make pos file
* Subset files
* Merge vcfs and create haplotype files
* Make class file for RFMix input

2.Run RFMix (See Below)

The script should take a few minutes to execute depending on the number of snps/samples being processed. 
After it has executed users can take the output and run RFMix as normal

**Example: Running RFMix**
```
 python RunRFMix.py -e 2 -w 0.2 --num-threads <n> --forward-backward PopPhased <example_merged.haps> <example.classes> <example.pos> -o <output.results>
```
