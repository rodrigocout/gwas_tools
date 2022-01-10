You will probably need the following conversion from vcf (HRC files) to BIMBAM.
BIMBAM format consists of three files, a mean genotype file, a phenotype file, and an optional SNP annotation file.
The first column of the genotype file is SNP id, the second and third columns are allele types with minor allele first, and the remaining columns are the posterior/imputed mean genotypes of different individuals numbered between 0 and 2.

An example of a mean genotype file with two SNPs and three individuals is as follows:
```
rs1, A, T, 0.02, 0.80, 1.50
rs2, G, C, 0.98, 0.04, 1.00

qctool_v2.0-rc9 -g gatk_exome.vcf -ofiletype bimbam_dosage -og abc.txt

-g $dose.vcf.gz -ofiletype bimbam_dosage -oq $chr.txt
```

•	plink2:
```
plink2 --vcf ${vcf} dosage=DS --export oxford --out genfile
cat genfile | awk -v s=[number of samples/individuals] '{ printf $2 "," $4 "," $5; for(i=1; i<=s; i++) printf "," $(i*3+3)2+$(i3+4); printf "\n" }' > bimbam
bcftools query -f '%ID %AF %REF %ALT [%DS ]\n' ${vcf} |perl -lane '{if($F[1]>0.5){for(my $i=4; $i<=$#F;$i++){$F[$i]=2-$F[$i];}}else{my$tmp=$F[2];$F[2]=$F[3];$F[3]=$tmp;}print("$F[0], ", join(", ", @F[2..$#F]))}' > bimbam
```

To run the association it’s quite straight forward (pages 17-18 in the manual).

For running multiple phenotypes (in chunks which makes it faster), have a look at the attached.
Usage:

```
./ GEMMA_assoc_forGARP_cov.bash      
```

You need to create some input files before setting it run such as covariates & phenotypes file, matrix etc. You will find the commands in the manual.
Also, you will need to change the locations of your files and software in the above command line, the  (eg KNEE_OA) and the .
If you find the attached wrapper script too confusing then just prepare your files and run it by using the assoc command line (line 42).

### 1. Generate the Relatedness Matrix
### plink

```
sed -i 's/-9/1/' ${output_plink}/chrALL.${cohort}.maf005.hwe1e-5.for-rel-mat.ldpruned.missnp.fam
~ag15/local_programs/gsub 10g -I ${gemma} -bfile ${output_plink}/chrALL.${cohort}.maf005.hwe1e-5.for-rel-mat.ldpruned.missnp -gk 1 -n 1 -o chrALL.${cohort}.maf005.hwe1e-5.for-rel-mat.ldpruned.missnp
```

### bimbam

```
~ag15/local_programs/gsub 10g -I ${gemma} -g ${bimbam}/chrALL_${cohort}_bimbam_geno_maf005.hwe1e-5.ldpruned.for-rel-mat  -p ${bimbam}/bimbam-pheno.txt -gk 1 -n 1 -o bimbam-matrix
```

#### 2. Creating annotation file from the vcf imputed files

```
for i in {1..22}; do bcftools query -f '%ID %CHROM  %POS \n' ${input}/chr"$i".dose.vcf.gz | awk '{print $1","$3","$2}' > ${bimbam}/ids_pos_chr"$i" ;
done
 ```
 
#### 3. Spliting the bimbam genotype file (input for assoc)

```
for i in {1..22}; do split -l 100000 ${bimbam}/ATHENS_bimbam_geno_chr"$i" ${input_bimbamfiles}/${cohort}_bimbam_geno_chr"$i"_chunk; 
done
```
