## 1. Download the files as VCF.gz (and tab-indices)
```
prefix="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr" ;

suffix=".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz" ;

for chr in {1..22}; do
    wget "${prefix}""${chr}""${suffix}" "${prefix}""${chr}""${suffix}".tbi ;
done
```

## 2. Download 1000 Genomes PED file
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped ;
```

## 3. Download the GRCh37 / hg19 reference genome
```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz ;

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.fai ;

gunzip human_g1k_v37.fasta.gz ;
```

# NB - if wget is not working, try curl:
```
curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz -O human_g1k_v37.fasta.gz
```

## 4. Convert the 1000 Genomes files to BCF

Ensure that multi-allelic calls are split and that indels are left-aligned compared to reference genome (1st pipe)
Sets the ID field to a unique value: CHROM:POS:REF:ALT (2nd pipe)
Removes duplicates (3rd pipe)
```
-I +'%CHROM:%POS:%REF:%ALT' means that unset IDs will be set to CHROM:POS:REF:ALT

-x ID -I +'%CHROM:%POS:%REF:%ALT' first erases the current ID and then sets it to CHROM:POS:REF:ALT

for chr in {1..22}; do
    bcftools norm -m-any --check-ref w -f human_g1k_v37.fasta \
      ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz | \
      bcftools annotate -x ID -I +'%CHROM:%POS:%REF:%ALT' | \
        bcftools norm -Ob --rm-dup both \
          > ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf ;

    bcftools index ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf ;
done
```

## 5. Convert the BCF files to PLINK format

```
for chr in {1..22}; do
    plink --noweb \
      --bcf ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.bcf \
      --keep-allele-order \
      --vcf-idspace-to _ \
      --const-fid \
      --allow-extra-chr 0 \
      --split-x b37 no-fail \
      --make-bed \
      --out ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;
done
```

## 6. Exclude variants not on the coding strand
NB - This step is only for microarray studies where the probes may only target one strand or the other (sense or non-sense)

## 7. Prune variants from each chromosome

```
--maf 0.10, only retain SNPs with MAF greater than 10%
--indep [window size] [step size/variant count)] [Variance inflation factor (VIF) threshold]

e.g. indep 50 5 1.5, Generates a list of markers in approx. linkage equilibrium - takes 50 SNPs at a time and then shifts by 5 for the window. VIF (1/(1-r^2)) is the cut-off for linkage disequilibrium

mkdir Pruned ;

for chr in {1..22}; do
    plink --noweb \
      --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
      --maf 0.10 --indep 50 5 1.5 \
      --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;

    plink --noweb \
      --bfile ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes \
      --extract Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.prune.in \
      --make-bed \
      --out Pruned/ALL.chr"${chr}".phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes ;
done
```

## 8. Get a list of all PLINK files

```
find . -name "*.bim" | grep -e "Pruned" > ForMerge.list ;

sed -i 's/.bim//g' ForMerge.list ;
```

## 9. Merge all projects into a single PLINK file

```
plink --merge-list ForMerge.list --out Merge ;
```

NB - if you have your own data that you want to merge with 1000 Genomes
Process the 1000 Genomes data as per this tutorial from Steps 1-9. In this way, you will have already identified the population-specific variants / markers. Then, after Step 9, do

find common variants between your dataset and the merged 1000 Genomes dataset (and filter both for these common variants)
merge the 1000 Genomes data with your own data
proceed to Step 10
Depending on its size, your own dataset may be divided by chromosome; so, you may have to do some pre-processing before aligning to 1000 Genomes. Either way, the population specific markers will be defined by just the 1000 Genomes dataset (Step 7). If your dataset is microarray, you'll have to pre-filter it for coding (plus / +) strand variants (Step 6).

## 10. Perform PCA

```
plink --bfile Merge --pca
```

## 11. Generate plots in R
```
R

options(scipen=100, digits=3)

# read in the eigenvectors, produced in PLINK
eigenvec <- read.table('plink.eigenvec', header = FALSE, skip=0, sep = ' ')
rownames(eigenvec) <- eigenvec[,2]
eigenvec <- eigenvec[,3:ncol(eigenvec)]
colnames(eigenvec) <- paste('Principal Component ', c(1:20), sep = '')

# read in the PED data
PED <- read.table('20130606_g1k.ped', header = TRUE, skip = 0, sep = '\t')
PED <- PED[which(PED$Individual.ID %in% rownames(eigenvec)), ]
PED <- PED[match(rownames(eigenvec), PED$Individual.ID),]
all(PED$Individual.ID == rownames(eigenvec)) == TRUE
[1] TRUE

# set colours
require('RColorBrewer')

# from: http://www.internationalgenome.org/category/population/
PED$Population <- factor(PED$Population, levels=c(
  "ACB","ASW","ESN","GWD","LWK","MSL","YRI",
  "CLM","MXL","PEL","PUR",
  "CDX","CHB","CHS","JPT","KHV",
  "CEU","FIN","GBR","IBS","TSI",
  "BEB","GIH","ITU","PJL","STU"))

col <- colorRampPalette(c(
  "yellow","yellow","yellow","yellow","yellow","yellow","yellow",
  "forestgreen","forestgreen","forestgreen","forestgreen",
  "grey","grey","grey","grey","grey",
  "royalblue","royalblue","royalblue","royalblue","royalblue",
  "black","black","black","black","black"))(length(unique(PED$Population)))[factor(PED$Population)]

# generate PCA bi-plots
project.pca <- eigenvec
summary(project.pca)


par(mar = c(5,5,5,5), cex = 2.0,
  cex.main = 7, cex.axis = 2.75, cex.lab = 2.75, mfrow = c(1,2))

plot(project.pca[,1], project.pca[,2],
  type = 'n',
  main = 'A',
  adj = 0.5,
  xlab = 'First component',
  ylab = 'Second component',
  font = 2,
  font.lab = 2)
points(project.pca[,1], project.pca[,2], col = col, pch = 20, cex = 2.25)
legend('bottomright',
  bty = 'n',
  cex = 3.0,
  title = '',
  c('Population 1', 'Population 2', 'Population 3',
    'Population 4', 'Population 5'),
  fill = c('yellow', 'forestgreen', 'grey', 'royalblue', 'black'))

plot(project.pca[,1], project.pca[,3],
  type="n",
  main="B",
  adj=0.5,
  xlab="First component",
  ylab="Third component",
  font=2,
  font.lab=2)
points(project.pca[,1], project.pca[,3], col=col, pch=20, cex=2.25)
```
