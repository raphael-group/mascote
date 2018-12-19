# Tutorial

This tutorial illustrates how to use the complete pipeline which is encoded in the main script [runMASCoTE.sh](../script/runMASCoTE.sh) in `script` folder.
The tutorial is subdivided into some subsections and each of these describes sequential parts of the full pipeline:
1. [Requirements](#requirements)
2. [Parameters](#parameters)
3. [Set up](#setup)
4. [Simulate genomes](#simulategenomes)
5. [Map reads](#mapreads)
6. [Process reads](#processreads)
7. [Simulate mixed samples](#simulatemixedsamples)

We suggest to make a copy of the script, place the script into the designated running directory, and follow the tutorial.

***
## Requirements
<a name="requirements"></a>

```shell
MASCOTE_HOME="/path/to/MASCoTE/"
MASCOTE='python2 ${MASCOTE_HOME}src/mascotte.py'
MIX='python2 ${MASCOTE_HOME}src/MixBAMs.py'
ART='/path/to/ART/bin/art_illumina'
BWA='/path/to/BWA/bwa'
SAM='/path/to/samtools/bin/samtools'

REF='/path/to/ref.fa'
SNP='/path/to/dbsnp.tsv'
IGNORE='/path/to/ignore.txt'
```

First, MASCoTE requirements must be specified. The full path to the home of this repository is specified thorugh the variable `${MASCOTE_HOME}` and the path to the two main components are atuomatically obtain by `${MASCOTE}` and `${MIX}`. Note that in these variables to full command in specified and hence `python` is preceeding the path. In addition, the full path to the 3 required method need to be specified by the variables `${ART}`, `${BWA}`, and `${SAM}`. Note that all these paths are pointing to the corresponding binaries.

MASCoTE requires 3 kind of data. A human reference genome in FASTA format pointed by the variable `${REF}`, a list of common gemrline SNPs pointed by the variable `${SNP}` (note this is optional but reccomended), and a list of contigs in the reference genome to ignore pointed by the variable `${IGNORE}`.

***
## Parameters
<a name="parameters"></a>

```shell
DIR="/path/to/rundir"
N=2
COV=30
J=4

CWGD=1
CWCL=1
CCAM=10
CFOC="20Mb:5 10Mb:10 3Mb:20 1Mb:10000:30"

SWGD=0
SWCL=1
SCAM=14
SFOC="20Mb:5 10Mb:10 3Mb:20 1Mb:10000:30"

ADRATIO=0.67

PROPS="0.1:0.9: 0.2::0.8"
```

Second, several parameters needs to be specified for the simulation.
- The running directory `${DIR}` where all the files produced by the simulation are produced.
- The total number of tumor clones `N`. Note that the total number of all clones is `N + 1` as it includes the normal diploid clone.
- The physical coverage `COV` simulated from every clone. On average, every nucleotide in the genome of a clone is covered by this amount of reads when simulating the reads from each individual clone.
- The maximum number of parallel threads `J` used when working on each clone.
- The fraction of amplifications over the deletions which are simulated when generating focal or chromosomal arm's CNAs (i.e. values higher than `0.5` result in a higher fraction of amplifications and values lower than `0.5` result in a lower fraction of deletions when generating these type of mutations).

Other parameters specify the kind, number, and characteristics of CNAs and WGDs to simulate for the genome of tumor clones. More specifically, these mutations are divided into:
- _clonal_, which are the mutations occurring in the trunk of the phylogenetic tree modelling tumor evolution and affecting the genome of all tumor clones. The different kind of mutations are the following:
    - The number of clonal WGDs is specified by `CWGD`
    - The number of clonal whole-chromosome losses is specified by `CWCL`
    - The number of clonal chromosomal arm amplifications and deletions is specified by `CCAM`
    - Different kinds of clonal focal CNAs are specified and characterized in a white-separated list within apices such that every kind of CNAs is characterized in the following format `SZ:QT` where `SZ` is the average size of the CNAs (the size can be specified as `Mb` or `kb`) and `QT` is the number of these CNAs. Moreover, the variance of the size distribution can be added by using the format `SZ:VR:QT` where `VR` is the variance of a normal distribution whose mean is `SZ`.
- _subclonal_, which are the mutations occurring the other branches of the phylogenetic tree modelling tumor evolution and affecting the genome of all tumor clones. The different kind of mutations are the following:
    - The number of subclonal WGDs is specified by `SWGD`
    - The number of subclonal whole-chromosome losses is specified by `SWCL`
    - The number of subclonal chromosomal arm amplifications and deletions is specified by `SCAM`
    - Different kinds of subclonal focal CNAs given in the formats specified above.

Last, the clone proportions for the mixed samples to simulate are specified in the variable `PROPS` which is a white-space separated list within apices with an element for each sample. Each element specifies the clone proportions in the corresponding sample with the following format `U:U_0:U_1: ... :U_N` such that: `U` is the proportion of normal diploid clone and `U_X` is the proporion of tumor clone `X` in the corresponding sample when the clone is present, otherwise is an emprty string. Remember there is a total of `N` total tumor clones and all clone proportions have to sum up to `1.0` in every sample. For example `0.1::0.9` is a sample comprising a mixture of `10%` of normal cells and `90%` of tumor cells from tumor clone 1, and `0.2:0.8:` is a sample comprising a mixture of `20%` of normal cells and `80%` of tumor cells from tumor clone 0.

***
## Set up
<a name="setup"></a>

```shell
set -e
set -o xtrace
PS4='\''[\t]'\'


echo -e "\033[1m\033[95m### Simulating datasets with number of clones '${CLONES}' \033[0m"

SEED=${RANDOM}
echo -e "\033[1m\033[95m## Selecting random seed: '${SEED}' \033[0m"

echo -e "\033[1m\033[95m## Setting up folders \033[0m"
FASTA=${DIR}'fasta/'
FASTQ=${DIR}'fastq/'
BAM=${DIR}'bam/'
BULK=${DIR}'bulk/'
mkdir -p ${DIR}
mkdir -p ${FASTA}
mkdir -p ${FASTQ}
mkdir -p ${BAM}
mkdir -p ${BULK}
```

Third, the script sets up the simulation enviroment. Three commands activate the log trace for the script which terminates in case of error and add time stamps to this. The `echo` commands are simply used for organizing the log message. A random seed is selected arbitrarily using the `${RANDOM}` enviromental variable of the OS; the random seed is used as a seed for the entire simulation pipeline. This value can also be fixed by the used instead of drawing a random value and is crucial for being able of replicating the entire simulation. Last, the script prepares the directories that are used to organize all the outputs from this run of MASCoTE. To avoid condlicts, user should make sure the running directory is an empty directory, expect when re-executing the same run.

***
## Simulate genomes
<a name="simulategenomes"></a>

```shell
echo -e "\033[1m\033[95m## Simulating genomes \033[0m"
\time -v ${MASCOTE} ${REF} -n ${N} -s ${SEED} -g ${IGNORE} -l ${SNP} \
                    -x ${FASTA} -b 1kb -j 20 -r ${ADRATIO} -cwgd ${CWGD} \
                    -cwcl ${CWCL} -ccam ${CCAM} -ccna "${CFOC}" -swgd ${SWGD} \
                    -swcl ${SWCL} -scam ${SCAM} -scna "${SFOC}" \
                    |& tee ${FASTA}mascotte.log


echo -e "\033[1m\033[95m## Merging the haplotypes FASTA in unique diploid FASTA \033[0m"
sed '/^>chr/ s/$/-A/' ${FASTA}human.maternal.fa > ${FASTA}normal.fa && sed '/^>chr/ s/$/-B/' ${FASTA}human.paternal.fa >> ${FASTA}normal.fa &
TCLONES=""
for (( i=0; i<${N}; i++ ))
do
    CLONE='clone'${i}
    TCLONES="${TCLONES} ${BAM}${CLONE}.bam"
    sed '/^>chr/ s/$/-A/' ${FASTA}${CLONE}.maternal.fa > ${FASTA}${CLONE}.fa && sed '/^>chr/ s/$/-B/' ${FASTA}${CLONE}.paternal.fa >> ${FASTA}${CLONE}.fa &
done
wait


echo -e "\033[1m\033[95m## Remove haplotype specific FASTA \033[0m"
rm -rf ${FASTA}*.maternal.fa
rm -rf ${FASTA}*.paternal.fa
```

Fourth, MASCoTE performs the first step of the simulation framework. This step aims to simulate via the module `mascote` a haplotype-specific diploid human genome and the haplotype-specific genome of every tumor clone characterized by the specified CNAs and WGDs. This steps produces:
- the haplotype-specific genome of every clone in [FASTA](https://en.wikipedia.org/wiki/FASTA) format. The two haplotypes of each clone are subsequently merged in the same FASTA format to be sequenced jointly. Only the combined genome with the two haplotypes for each clone is kept and the haplotype specific FASTA files are removed for simplicity. The merging is performed in parallele for every clone, so please consider the availalble number of processes and in case of need remove the ending `&` to avoid the parallel execution.
- a phylogenetic tree describin the evolution of the corresponding clones in [DOT](https://en.wikipedia.org/wiki/DOT) format. Note that this phylogenetic tree can be drawn in the corresponding PDF figure using the `dot` program.
- A tab-separated file describing the copy-number profile of every tumor clone (whose total number is `N` and the names are `clone0`, ..., `cloneN-1`) with the following fields:

| Field | Comment |
|-------|---------|
| `#CHR` | A name of a simulated chromosome |
| `START` | The genomic position representing the start of a genomic segment |
| `END` | The genomic position representing the end of a genomic segment |
| `clone0` | The allele-specific copy numbers of `clone0` in the genomic segment `(START, END)`, given in the format `A|B` where `A` and `B` are the corresponding allele-specific opy numbers |
| ... | ... | ... |
| `cloneN-1` | The allele-specific copy numbers of `cloneN-1` in the genomic segment `(START, END)`, given in the format `A|B` where `A` and `B` are the corresponding allele-specific opy numbers |

Note that all these output files are generated in the `fasta` directory in the running directory.

***
## Sequence reads
<a name="sequencereads"></a>

```shell
echo -e "\033[1m\033[95m## Simulating sequencing reads\033[0m"
NORMAL=${FASTQ}normal/
mkdir -p ${NORMAL}
\time -v ${ART} -i ${FASTA}normal.fa -p -ss HS25 -f ${COV} -na -l 150 -m 200 -s 10 -o ${NORMAL}normal &> ${NORMAL}normal.log &
for (( i=0; i<${N}; i++ ))
do
    CLONE=clone${i}
    DIRCLONE=${FASTQ}${CLONE}/
    mkdir -p ${DIRCLONE}
    \time -v ${ART} -i ${FASTA}${CLONE}.fa -p -ss HS25 -f ${COV} -na -l 150 -m 200 -s 10 -o ${DIRCLONE}${CLONE} &> ${DIRCLONE}${CLONE}.log &
done
wait
```

Fifth, MASCoTE simulates the sequencing reads from the genome of each individual clone using [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm). There are `N + 1` clones including a normal diploid clone and tumor clones with names from `clone0` to `cloneN-1`. The sequencing is performed using the standard profile `HS25` for Illumina techonolgies. The user can consider different profiles from ART documentation. The simulated sequecing of all clones is performed in parallel, so please consider the availalble number of processes and in case of need remove the ending `&` to avoid the parallele execution.

***
## Map reads
<a name="mapreads"></a>

```shell
echo -e "\033[1m\033[95m## Mapping sequencing reads to reference-human genome \033[0m"
NORMAL=${FASTQ}'normal/'
\time -v ${BWA} mem ${REF} ${NORMAL}normal1.fq ${NORMAL}normal2.fq -t ${J} 1> ${BAM}normal.sam 2> ${BAM}normal.sam.log &
for (( i=0; i<${N}; i++ ))
do
    CLONE=clone${i}
    DIRCLONE=${FASTQ}${CLONE}'/'
    \time -v ${BWA} mem ${REF} ${DIRCLONE}${CLONE}1.fq ${DIRCLONE}${CLONE}2.fq -t ${J} 1> ${BAM}${CLONE}.sam 2> ${BAM}${CLONE}.sam.log &
done
wait
```

Sixth, MASCoTE maps the simulated sequencing reads to the reference human genome using [BWA](http://bio-bwa.sourceforge.net/). The mapping of the reads from all clones is performed in parallel, so please consider the availalble number of processes and in case of need remove the ending `&` to avoid the parallele execution.

***
## Process reads
<a name="processreads"></a>

```shell
echo -e "\033[1m\033[95m## Sorting BAMs \033[0m"
mkdir ${BAM}tmp_normal/
\time -v ${SAM} sort -O bam -o ${BAM}normal.bam -T ${BAM}tmp_normal/ ${BAM}normal.sam -@ ${J} &> ${BAM}normal.bam.log &
for (( i=0; i<${N}; i++ ))
do
    CLONE=clone${i}
    mkdir -p ${BAM}tmp_${CLONE}/
    \time -v ${SAM} sort -O bam -o ${BAM}${CLONE}.bam -T ${BAM}tmp_${CLONE}/ ${BAM}${CLONE}.sam -@ ${J} &> ${BAM}${CLONE}.bam.log &
done
wait


echo -e "\033[1m\033[95m## Cleaning SAM files and mapping temporary files \033[0m"
rm -rf ${BAM}*.sam
rm -rf ${BAM}tmp_*


echo -e "\033[1m\033[95m## Indexing BAMs \033[0m"
\time -v ${SAM} index ${BAM}normal.bam &> ${BAM}normal.bam.bai.log &
for (( i=0; i<${N}; i++ ))
do
    CLONE=clone${i}
    \time -v ${SAM} index ${BAM}${CLONE}.bam &> ${BAM}${CLONE}.bam.bai.log &
done
wait
```

Seventh, MASCoTE processes the mapped sequencing reads. More specifically, MASCoTE sorts the mapped sequencing reads, remove temporary and unmapped read files for simplicity, and indexes the mapped sequencing reads. If one wants to perform more processing steps (e.g. marking and removing duplicates), these can be added at this point of the pipeline. These processes are executed in parallel for all clones, so please consider the availalble number of processes and in case of need remove the ending `&` to avoid the parallele execution.

***
## Simulate mixed samples
<a name="simulatemixedsamples"></a>

```shell
cd ${BULK}
\time -v ${MIX} -n ${BAM}normal.bam -t ${TCLONES} -p "${PROPS}" \
                -c ${FASTA}copynumbers.csv -st $(dirname ${SAM}) \
                -T ${DIR}tmp -o bulk -j ${J} |& tee mix.log
```

Last, MASCoTE simulate the mixed samples. More specifically, MASCoTE simulates each mixed sample by sampling the reads from the ones sequenced from each individual clone present in the sample and mix these reads according to the corrected proportions. In particular, the corrected proportions are computed by taking into account the genome lengths of all the clones computed from the clone copy-number profiles and the given clone proportions.
