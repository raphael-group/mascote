# mascote

This module generates a human-diploid genome, generates a phylogenetic tree describing the evolution of a collection of tumor clones, and generates the genome of these clones by introducing different kinds of CNAs and WGDs. 

## Input

This modules requires 2 input:

| Name | Description | Usage |
|------|-------------|-------|
| `REF` | A reference human genome | The reference human is used as a base to simulate a human diploid genome |
| `-n`, `--numclones` | Number of tumor clones | The number of tumor clones to simulate in addition to the normal diploid clone  |

## Output

This modules generates three kind of output:

| Name | Description | Usage |
|------|-------------|-------|
| <ul><li>`normal_maternal`, `clone0_maternal.fa`, ..., `cloneN-1_maternal.fa`</li><li>`normal_paternal`, `clone0_paternal.fa`, ..., `cloneN-1_paternal.fa`</li></ul> | Two haplotype-specific FASTA genomes every clone (normal diploi and tumor clones) | Every haplotype-specific FASTA genome contains all the maternal or paternal copies of each chromosome for the haplotype of the corresponding clone |
| `copynumbers.csv` | A tab-separated file describing the allele and clone-specific copy-number profiles | The fields of the file are <ul><li>`#CHR`: A name of a simulated chromosome</li><li>`START`: the genomic position representing the start of a genomic segment</li><li>`END`: the genomic position representing the end of a genomic segment</li><li>`clone0`: the allele-specific copy numbers of `clone0` in the genomic segment `(START, END)`, given in the format `<code>A&#124;B</code>` where `A` and `B` are the corresponding allele-specific opy numbers</li><li>...</li><li>`cloneN-1`: the allele-specific copy numbers of `cloneN-1` in the genomic segment `(START, END)`, given in the format `<code>A&#124;B</code>` where `A` and `B` are the corresponding allele-specific opy numbers</li></ul> |
| `tumor.dot` | A phylogenetic tree describing the tumor evolution with the corresponding CNAs and WGDs | The tree is given in the `DOT` format and the command `dot` can be used to transform it into the corresponding PDF figure as `dot -Tpdf tumor.dot -o tumor.pdf` |

## Main parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-l`, `--snplist` | Path to SNP list | Heterozygous and homozygous SNP positions are generated considering all given positions | Required |
| `-g`, `--ignore` | Path to ignore list file | Only the chromosomes of the reference genome NOT included in this list are considered for the simulation | All contigs in the reference genome are considred |
| `-e`, `--hehoratio` | Fraction of heterozygous SNPs to generate over homozygous SNPs | This fraction is used to simulate the SNPs and fix the fraction of heterozygous/homozygous SNPs. A fraction higher than `0.5` results in most of the SNPs being heterozygous, while a fraction lower than `0.5` results in most of the SNPs being homozygous | `0.67` |
| `-r`, `--adratio` | Fraction of amplifications over deletions | This fraction is used to choose between duplications and deletions when simulating focal CNAs or chromosomal's arms aberrations and fix the fraction of duplications/deletions. A fraction higher than `0.5` results in most of the aberrations being duplications, while a fraction lower than `0.5` results in most of the aberrations being deletions | `0.65` |
| `-cwgd`, `--clonalwgd` | Number of clonal WGDs | Number of WGDs simulated in the trunk of the phylogeny and affecting the genome of all tumor clones | 0 |
| `-cwcl`, `--clonalwcl` | Number of clonal chromosomal losses | Number of chromosomal losses simulated in the trunk of the phylogeny and affecting the genome of all tumor clones. A chromosomal loss is the loss of a haplotype specific copy of a chromosome | 0 |
| `-ccam`, `--clonalcam` | Number of clonal chromosomal's arm aberrations | Number of chromosomal's arm aberrations simulated in the trunk of the phylogeny and affecting the genome of all tumor clones. A chromosomal's arm aberration is either a duplication or deletion of a haplotype specific arm of a chromosome | 0 |
| `-ccna`, `--clonalcna` | Number of clonal focal CNAs | Number of focal CNAs simulated in the trunk of the phylogeny and affecting the genome of all tumor clones. A focal CNAs is either a duplication or deletion (chosen by corresponding fraction) of a region in a haplotype specific copy of a chromosome (chosen depending on chromosome's current length). Different kinds of clonal focal CNAs are specified and characterized in a white-separated list within apices such that every kind of CNAs is characterized in the following format `SZ:QT` where `SZ` is the average size of the CNAs (the size can be specified as `Mb` or `kb`) and `QT` is the number of these CNAs. Moreover, the variance of the size distribution can be added by using the format `SZ:VR:QT` where `VR` is the variance of a normal distribution whose mean is `SZ`. | None, an example is `"20Mb:5 10Mb:10 3Mb:20 1Mb:10000:30"` where 5 focal CNAs are simulated with size of `20Mb`, `10` of `10Mb`, `20` of `3Mb`, and `30` with the size drawn from a normal distribution around `1Mb` with a variance of `10000` |
| `-swgd`, `--subclonalwgd` | Number of subclonal WGDs | Number of WGDs simulated in a random branch of the phylogeny and affecting the genome of some tumor clones | 0 |
| `-swcl`, `--subclonalwcl` | Number of subclonal chromosomal losses | Number of chromosomal losses simulated in a random branch of the phylogeny and affecting the genome of some tumor clones. A chromosomal loss is the loss of a haplotype specific copy of a chromosome | 0 |
| `-scam`, `--subclonalcam` | Number of subclonal chromosomal's arm aberrations | Number of chromosomal's arm aberrations simulated in a random branch of the phylogeny and affecting the genome of some tumor clones. A chromosomal's arm aberration is either a duplication or deletion of a haplotype specific arm of a chromosome | 0 |
| `-scna`, `--subclonalcna` | Number of subclonal focal CNAs | Number of focal CNAs simulated in a random branch of the phylogeny and affecting the genome of some tumor clones. A focal CNAs is either a duplication or deletion (chosen by corresponding fraction) of a region in a haplotype specific copy of a chromosome (chosen depending on chromosome's current length). Different kinds of clonal focal CNAs are specified and characterized in a white-separated list within apices such that every kind of CNAs is characterized in the following format `SZ:QT` where `SZ` is the average size of the CNAs (the size can be specified as `Mb` or `kb`) and `QT` is the number of these CNAs. Moreover, the variance of the size distribution can be added by using the format `SZ:VR:QT` where `VR` is the variance of a normal distribution whose mean is `SZ`. | None, an example is `"20Mb:5 10Mb:10 3Mb:20 1Mb:10000:30"` where 5 focal CNAs are simulated with size of `20Mb`, `10` of `10Mb`, `20` of `3Mb`, and `30` with the size drawn from a normal distribution around `1Mb` with a variance of `10000` |
| `-s`, `--rndseed` | Random seed | This random seed is necessary for replicating an exection | None, non-deterministic execution |

## Optional parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-j`, `--jobs` | Number of parallele jobs | Chromosomes are executed on parallel | 1, the suggested value is the number of simulated chromosomes when possible |
| `-x`, `--runningdirectory` | Running directory | Running directory where all output files are generated | Current directory |
| `-p`, `--snpratio` | Fraction of SNPs | The fraction of the SNPs to simulate along whole genome. This ratio is only used when randomly generating the SNPs along the entire genome because a list of SNPs is not provided | None |
| `-b`, `--binsize` | Resolution for breakpoints | Resolution used for selecting the breakpoint when simulating chromosomal's arms aberrations and focal CNAs | `10kb` |
| `-v`, `--noverbose` | Activate non-verbose log | Decrease the verbosity of the generated log | Verbose log |
