#!/usr/bin/bash


MASCOTE_HOME="/path/to/MASCoTE/"
MASCOTE="python2 ${MASCOTE_HOME}src/mascotte.py"
MIX="python2 ${MASCOTE_HOME}src/MixBAMs.py"
ART='/path/to/ART/bin/art_illumina'
BWA='/path/to/BWA/bwa'
SAM='/path/to/samtools/bin/samtools'

REF='/path/to/ref.fa'
SNP='/path/to/dbsnp.tsv'
IGNORE='/path/to/ignore.txt'

DIR="/path/to/rundir"
N=2
COV=30
J=22

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


echo -e "\033[1m\033[95m## Simulating genomes \033[0m"
\time -v ${MASCOTE} ${REF} -n ${N} -s ${SEED} -g ${IGNORE} -l ${SNP} -x ${FASTA} -b 1kb -j ${J} -r ${ADRATIO} -cwgd ${CWGD} -cwcl ${CWCL} -ccam ${CCAM} -ccna "${CFOC}" -swgd ${SWGD} -swcl ${SWCL} -scam ${SCAM} -scna "${SFOC}" |& tee ${FASTA}mascotte.log


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


echo -e "\033[1m\033[95m## Mapping sequencing reads to reference-human genome and sorting results\033[0m"
NORMAL=${FASTQ}'normal/'
mkdir ${BAM}tmp_normal/
(\time -v ${BWA} mem ${REF} ${NORMAL}normal1.fq ${NORMAL}normal2.fq -t ${J} | ${SAM} sort - -O bam -o ${BAM}normal.bam -T ${BAM}tmp_normal/ -@ ${J} 2> ${BAM}normal.bam.log) &
for (( i=0; i<${N}; i++ ))
do
    CLONE=clone${i}
    DIRCLONE=${FASTQ}${CLONE}'/'
    mkdir -p ${BAM}tmp_${CLONE}/
    (\time -v ${BWA} mem ${REF} ${DIRCLONE}${CLONE}1.fq ${DIRCLONE}${CLONE}2.fq -t ${J} | ${SAM} sort - -O bam -o ${BAM}${CLONE}.bam -T ${BAM}tmp_${CLONE}/ -@ ${J} 2> ${BAM}${CLONE}.bam.log) &
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


cd ${BULK}
\time -v ${MIX} -n ${BAM}normal.bam -t ${TCLONES} -p "${PROPS}" -c ${FASTA}copynumbers.csv -st $(dirname ${SAM}) -T ${DIR}tmp -o bulk -j ${J} |& tee mix.log
