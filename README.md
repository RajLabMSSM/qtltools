# QTLtools: a tool set for molecular QTL discovery and analysis

## Description

Paper: [Delaneau et. al. A complete tool set for molecular QTL discovery and analysis. Nature Communications (2017)](https://www.nature.com/articles/ncomms15452)

git: https://github.com/qtltools/qtltools

web: https://qtltools.github.io/qtltools/

developers:
* Olivier Delaneau (olivier.delaneau@gmail.com)
* Halit Ongen (halit.ongen@unige.ch)

## Tutorial for running on Minerva

### Preparing your input files

Please follow instructions on [GTEx pipeline](https://github.com/RajLabMSSM/gtex-pipeline/tree/master/qtl) to prepare your files properly. 

In the end you will need the following files:

#### VCF (or BCF) file containing the genotype data.
  + The vcf must be compressed and indexed. Use the following commands to do this.
~~~
module load tabix
bgzip genotypes.vcf && tabix -p vcf genotypes.vcf.gz
~~~
  + All QTLtools functionalities can use either the GT or DS (i.e. genotype dosage) fields.
  + To validate your VCF/BCF file, use bcftools as follows:
~~~
bcftools view myGenotypes.vcf.gz | less -S
~~~
  + Below there is an example of a VCF file with 3 variants and 4 individuals
~~~
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT UNR1	UNR2	UNR3	UNR4
chr7	123	SNP1	A	G	100	PASS	.	GT:DS	0/0:0.001	0/0:0.000	0/1:0.999	1/1:1.999
chr7	456	SNP2	T	C	100	PASS	.	GT:DS	0/0:0.001	0/0:0.000	0/1:1.100	0/0:0.100
chr7	789	SNP3	A	T	100	PASS	.	GT:DS	1/1:2.000	0/1:1.001	0/0:0.010	0/1:0.890
~~~

#### BED file containing the phenotype data.
  + Here we will usually have the expression data. So in addition to the standard expression information (e.g. a table with rows as genes/transcripts and columns as samples), we'll need the genomic position of each gene or transcript. This will follow the standard BED file formart. This file is TAB delimited. Each line corresponds to a single molecular phenotype. The first 6 columns are:
    1. Chromosome ID [string]
    2. Start genomic position of the phenotype (here the TSS of gene1) [integer, 0-based]
    3. End genomic position of the phenotype (here the TSS of gene1) [integer, 1-based]
    4. Phenotype ID (here the exon IDs) [string].
    5. Phenotype group ID (here the gene IDs, multiple exons belong to the same gene) [string]
    6. Strand orientation [+/-]
  + The other columns will contain the expression for each sample. See an example below with 3 genes/transcripts and 4 samples:
~~~ 
#Chr	start	end	pid	gid	strand	UNR1	UNR2	UNR3	UNR4 chr1	99999	100000	pheno1	pheno1	+	-0.50	0.82	-0.71 0.83
chr1	199999	201000	pheno2	pheno2	+	1.18	-2.84	1.34	-1.56
chr1	299999	300000	exon1	gene1	+	-1.13	1.18	-0.03	0.11
chr1	299999	300000	exon2	gene1	+	-1.18	1.32	-0.36	1.26
~~~
  + The BED files needs to be sorted, compressed and indexed as follows:
~~~
module load tabix
bgzip myPhenotypes.bed && tabix -p bed phenotypes.bed.gz
~~~

#### Optional: A covariate data (TXT)
  + The COV file contains the covariate data in simple TXT format. Hereafter an example of 4 covariates for 4 samples.
~~~
id	UNR1	UNR2	UNR3	UNR4
PC1	-0.02	0.14	0.16	-0.02
PC2	0.01	0.11	0.10	0.01
PC3	0.03	0.05	0.08	0.07
BIN	1	0	0	1
~~~

### Running a QTL in *cis* - nominal pass

Here we are gonna run a cis QTL discovery on Minerva. The data used here is available in the **examples** folder.

We will use 3 input files:

* The phenotype data matrix for chr22 on 358 samples: BED / index
* The genotype data matrix for chr22 on 358 samples: VCF / index
* The covariate data matrix on 358 samples: TXT

*files were downloaded from [here](https://qtltools.github.io/qtltools/)*

We will run it by splitting the analysis in several jobs. 

For this example we will split in 10 chunks (one for each job). 

Just remember to set the number of chunks >= the number of chromosomes that you have.

*edit the script bellow as necessary*
~~~
module load qtltools
numChunks=10
for j in $(seq 1 ${numChunks}); do
     echo "QTLtools cis --vcf genotypes.chr22.vcf.gz --bed genes.50percent.chr22.bed.gz --cov genes.covariates.pc50.txt.gz --nominal 1 --chunk $j $numChunks --out nominals_${j}_${numChunks}.txt.gz" \
     | bsub \
      -n 1 \
      -R "rusage[mem=4000]" \
      -W 144:00 \
      -oo qtltools.chunk.${j}.out \
      -eo qtltools.chunk.${j}.err \
      -P acc_ad-omics \
      -q premium \
      -J qtltools.chunk.${j}
done
~~~

After finishing running all jobs, you can merge the chunks using the command bellow
~~~
zcat nominals_*.txt.gz | gzip -c > mergedChunks_nominals.gz
~~~

* The output will look like:
~~~
ENSG00000237438.2 chr22 17517460 17517460 + 5803 -202619 rs74762450 chr22 17314841 17314841 0.00957345 -0.347492 0
ENSG00000237438.2 chr22 17517460 17517460 + 5803 -165953 rs5748687 chr22 17351507 17351507 0.00227542 -0.225077 0
ENSG00000237438.2 chr22 17517460 17517460 + 5803 -165363 rs4006343 chr22 17352097 17352097 0.0022876 -0.226253 0
ENSG00000237438.2 chr22 17517460 17517460 + 5803 -164231 rs3875996 chr22 17353229 17353229 0.00181073 -0.233649 0
...
ENSG00000099954.14 chr22 17840837 17840837 + 6200 469530 rs5992121 chr22 18310367 18310367 0.00653855 -0.00519037 0
ENSG00000099954.14 chr22 17840837 17840837 + 6200 473554 rs367922 chr22 18314391 18314391 0.00125305 -0.00612444 0
ENSG00000099954.14 chr22 17840837 17840837 + 6200 476586 rs11705197 chr22 18317423 18317423 0.00961104 -0.00498535 0
ENSG00000099954.14 chr22 17840837 17840837 + 6200 754091 rs362027 chr22 18594928 18594928 0.00673601 0.0100518 0
ENSG00000099954.14 chr22 17840837 17840837 + 6200 754094 rs361893 chr22 18594931 18594931 0.00673601 0.0100518 0
~~~

* The columns are:
1. The phenotype ID
2. The chromosome ID of the phenotype
3. The start position of the phenotype
4. The end position of the phenotype
5. The strand orientation of the phenotype
6. The total number of variants tested in cis
7. The distance between the phenotype and the tested variant (accounting for strand orientation)
8. The ID of the tested variant
9. The chromosome ID of the variant
10. The start position of the variant
11. The end position of the variant
12. The nominal P-value of association between the variant and the phenotype
13. The corresponding regression slope
14. A binary flag equal to 1 is the variant is the top variant in cis



