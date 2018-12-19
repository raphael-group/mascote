# MixBAMs

This module considers mapped sequencing reads obtained from the genome of individual clones and produces multiple mixed samples according to given clone proportions and with appropriate corrections for the different genome lengths of the clones. More specifically, MASCoTE simulates each mixed sample by sampling the reads from the ones sequenced from each individual clone present in the sample and mix these reads according to the corrected proportions. In particular, the corrected proportions are computed by taking into account the genome lengths of all the clones computed from the clone copy-number profiles and the given clone proportions.

## Input

| Name | Description | Usage |
|------|-------------|-------|
| `-n`, `--normal` | Path to indexed-sorted BAM for matched-normal sample | Sequencing reads from normal diploid clone |
| `-t`, `--tumors` | List of paths to indexed-sorted BAM for tumor samples sample | The white-separated list has to be given between apices |
| `-p`,`--proportions` | Clone proportions | The clone proportions are given as a white-space separated list within apices with an element for each sample. Each element specifies the clone proportions in the corresponding sample with the following format `U:U_0:U_1: ... :U_N` such that: `U` is the proportion of normal diploid clone and `U_X` is the proporion of tumor clone `X` in the corresponding sample when the clone is present, otherwise is an emprty string. Remember there is a total of `N` total tumor clones and all clone proportions have to sum up to `1.0` in every sample. For example `0.1::0.9` is a sample comprising a mixture of `10%` of normal cells and `90%` of tumor cells from tumor clone 1, and `0.2:0.8:` is a sample comprising a mixture of `20%` of normal cells and `80%` of tumor cells from tumor clone 0 |
| `-c`, `--copynumbers` | Path to files with clone and allele-specific copy-number profiles | The profiles are used to compute the genome length of each clone |

This module also requires to execute SAMtools whose binary should be specified using the flag `-st`, `--samtools`.

## Output

For each target mixed sample specified with the given clone proportions, this module produces a sorted-and-index BAM file of the corresponding sample. To avoid issues with sampling, the total number of reads for each sample is fixed to be equal to the minimum total numbers of reads obtained individually from any present clone.

## Parameters

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-u`, `--totalcounts` | Total read counts from every clone | The total numbers can be specified instead of computed to improve the performance. They should be specified in a tab-separated file with the same of a clone sample and the corresponding number | None, total counts are computed directly |
| `-s`,`--names` | Name of the samples for each clone used in totalcounts file | The name are used to map the names in totalcounts file and the actual clones | None, required only when totalcounts are given |
| `-T`, `--temp` | Temporary directory | The temporary directoy is used to store several large temporary files and is removed at the end of the process | `./tmp` |
| `-o`, `--output` | Prefix for output files |   | Automatically selected |
| `-q`, `--qualityreads` | Minimum read quality | Only reads with at least this quality are considered and sampled | `20` |
| `-j`, `--jobs` | Number of threads | Number of parallel threads to use | `1` |
| `-rs`, `--seed` | Randomg seed  |   | `12` |
