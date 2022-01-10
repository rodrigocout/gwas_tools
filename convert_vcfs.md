## Convert to plink, and filter by call rate. __

```
for i in *vcf;do
    plink1.9 --allow-extra-chr --allow-no-sex \
        --vcf ${i}  \
        --make-bed \
        --geno 0 \
        --chr chr2L, chr2R, chr3L, chr3R \
        --keep-allele-order \
        --out ${i%.vcf}
                done;
```

## Convert back to vcf (I couldn't get the input file names right so they're written in full) . ___ 

```
for i in hclone_lhm_only_dbSNP_biSNP dgrp2_dm6_dbSNP; do
    plink1.9 --allow-extra-chr --allow-no-sex \
    --bfile ${i} \
    --recode vcf \
    --keep-allele-order \
        --out ${i}.bas
            done;
```

## Combine the two vcfs. ___

```
java -jar ~/SOFTWARE/GenomeAnalysisTK-3.4-46/GenomeAnalysisTK.jar \
        -R ${refseq} \
        -T CombineVariants \
            -V ${vcf1%.vcf}.bas.vcf \
            -V ${vcf2%.vcf}.bas.vcf \
                --unsafe LENIENT_VCF_PROCESSING \
                --genotypemergeoption UNSORTED \
                    -o comb.vcf
```
 
