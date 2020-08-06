import sys
import os.path
import argparse
import subprocess
import textwrap

import Support as sp



def parse_mascotte_arguments():
    """
    Parse command line arguments
    Returns:
    """
    description = "Simulate mixtures of distinct tumor clones characterized by copy-number aberrations (CNAs) resulted from an evolutionary process."
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('REFERENCE', type=str, help='A human-reference genome in FASTA format used to simulate a human genome and the descending tumor clones')
    parser.add_argument('-n', '--numclones', type=int, required=True, help='The number of clones present in the mixture to simulate')
    parser.add_argument('-s', '--rndseed', type=int, required=False, default=None, help='The number of clones present in the mixture to simulate')
    parser.add_argument('-g', '--ignore', type=str, required=False, default=None, help='File-name containing a line-per-line list with chromosome names to ignore in the reference')
    parser.add_argument('-l', '--snplist', type=str, required=False, default=None, help=textwrap.dedent('File-name containing a SNP positions to add in the simulate human genome in the\nformat "#CHR POSITION REF_ALLELES ALT_ALLELES" with first line as headline (default: SNPs are placed randomly)'))
    parser.add_argument('-p', '--snpratio', type=float, required=False, default=None, help='Ratio of SNPs to place randomly when a snplist is not given (default: None, snpratio is requried only whether SNPLIST is not provided)')
    parser.add_argument('-e', '--hehoratio', type=float, required=False, default=0.67, help='Ratio of heterozygous SNPs compared to homozygous ones (default: 0.67)')
    parser.add_argument('-x', '--runningdirectory', type=str, required=False, default='./', help='Running directory where resulting files and logs are created (default: current directory)')
    parser.add_argument('-b', '--binsize', type=str, required=False, default='10kb', help='Size of bins to simulate tumor-clone genomes and corresponding CNAs (default: 10kb)')
    parser.add_argument('-r', '--adratio', type=float, required=False, default=0.65, help='Proportion of amplification-deletion in the simulated events (default: 0.65)')
    parser.add_argument('-cwgd', '--clonalwgd', type=int, required=False, default=0, help='Number of clonal whole-genome duplications (WGDs) to introduce in the ancestor of tumor evolution (default: 0)')
    parser.add_argument('-cwcl', '--clonalwcl', type=int, required=False, default=0, help='Number of clonal whole-chromosome losses (WCLs) to introduce in the ancestor of tumor evolution (default: 0)')
    parser.add_argument('-ccam', '--clonalcam', type=int, required=False, default=0, help='Number of clonal chromosomal-arm changes (CAMs) to introduce in the ancestor of tumor evolution (default: 0)')
    parser.add_argument('-ccna', '--clonalcna', type=str, required=False, default=None, help=textwrap.dedent("A list of different types of clonal focal copy-number aberrations to introduce in the ancestor of tumor\nevolution where different types of CNAs are separated by white-spaces and each is given in the format 'MEAN_LENGTH:QUANTITY' or 'MEAN_LENGTH:STD_DEVIATION:QUANTITY' (default: 0)"))
    parser.add_argument('-swgd', '--subclonalwgd', type=int, required=False, default=0, help='Number of clonal whole-genome duplications (WGDs) to introduce in subclonal branches of tumor evolution (default: 0)')
    parser.add_argument('-swcl', '--subclonalwcl', type=int, required=False, default=0, help='Number of clonal whole-chromosome losses (WCLs) to introduce in subclonal branches of tumor evolution (default: 0)')
    parser.add_argument('-scam', '--subclonalcam', type=int, required=False, default=0, help='Number of clonal chromosomal-arm changes (CAMs) to introduce in subclonal branches of tumor evolution (default: 0)')
    parser.add_argument('-scna', '--subclonalcna', type=str, required=False, default=None, help=textwrap.dedent("A list of different types of subclonal focal copy-number aberrations to introduce in the ancestor of tumor\nevolution in the format 'MEAN_LENGTH:STD_DEVIATION:QUANTITY [MEAN_LENGTH:STD_DEVIATION:QUANTITY] where standard deviation can be omitted and is computed as 20%s of mean' (default: None)" % '%%'))
    parser.add_argument('-j', '--jobs', type=int, required=False, default=1, help='The number of parallel jobs to use (default: 1)')
    parser.add_argument("-v", "--noverbose", action='store_false', default=True, required=False, help="Silence verbose log messages")
    args = parser.parse_args()

    if not os.path.isfile(args.REFERENCE):
        raise ValueError(sp.error("The specified human reference-genome file does not exist!"))
    if args.numclones < 0:
        raise ValueError(sp.error("The number of clones must be a positive integer or zero, when only a matched-normal sample should be generated!"))
    if args.rndseed != None and args.rndseed < 0:
        raise ValueError(sp.error("The random seed must be a positive integer!"))
    if args.ignore != None and not os.path.isfile(args.ignore):
        raise ValueError(sp.error("The specified ignore-list file does not exist!"))
    if args.snplist != None and not os.path.isfile(args.snplist):
        raise ValueError(sp.error("The specified SNP-list file does not exist!"))
    if args.snpratio != None and args.snpratio < 0.0 and args.snpratio > 1.0:
        raise ValueError(sp.error("The SNP ratio must be in [0.0, 1.0]!"))
    if args.hehoratio != None and args.hehoratio < 0.0 and args.snpratio > 1.0:
        raise ValueError(sp.error("The heterozygous-homozygous ratio must be in [0.0, 1.0]!"))
    if args.adratio < 0.0 and args.adratio > 1.0:
        raise ValueError(sp.error("The SNP ratio must be in [0.0, 1.0]!"))
    if not os.path.isdir(args.runningdirectory):
        raise ValueError(sp.error("The specified running directory does not exist!"))
    if args.clonalwgd < 0:
        raise ValueError(sp.error("The number of clonal WGD must be a positive integer!"))
    if args.clonalwcl < 0:
        raise ValueError(sp.error("The number of clonal WCL must be a positive integer!"))
    if args.clonalcam < 0:
        raise ValueError(sp.error("The number of clonal CAM must be a positive integer!"))
    if args.subclonalwgd < 0:
        raise ValueError(sp.error("The number of subclonal WGD must be a positive integer!"))
    if args.subclonalwcl < 0:
        raise ValueError(sp.error("The number of subclonal WCL must be a positive integer!"))
    if args.subclonalcam < 0:
        raise ValueError(sp.error("The number of subclonal CAM must be a positive integer!"))
    if args.jobs <= 0:
        raise ValueError(sp.error("The number of jobs must be a non-zero positive integer!"))

    if args.clonalcna != None and sum(len(event.split(':')) != 2 and len(event.split(':')) != 3 for event in args.clonalcna.split()) > 0:
        raise ValueError(sp.error('The clonal focal CNAs are given in wrong format!'))
    if args.subclonalcna != None and sum(len(event.split(':')) != 2 and len(event.split(':')) != 3 for event in args.subclonalcna.split()) > 0:
        raise ValueError(sp.error('The subclonal focal CNAs are given in wrong format!'))

    mutations = {}
    mutations['ratioAD'] = args.adratio

    mutations['clonalwgd'] = args.clonalwgd
    mutations['clonalwcl'] = args.clonalwcl
    mutations['clonalcam'] = args.clonalcam
    mutations['clonalfocal'] = {}
    if args.clonalcna != None:
        mutations['clonalfocal'] = {(sp.basesize(event.split(':')[0]), sp.basesize(event.split(':')[1])) if len(event.split(':')) == 3 else (sp.basesize(event.split(':')[0]), int(sp.basesize(event.split(':')[0]) * 0.2)) : int(event.split(':')[2]) if len(event.split(':')) == 3 else int(event.split(':')[1]) for event in args.clonalcna.strip().split()}

    mutations['subclonalwgd'] = args.subclonalwgd
    mutations['subclonalwcl'] = args.subclonalwcl
    mutations['subclonalcam'] = args.subclonalcam
    mutations['subclonalfocal'] = {}
    if args.subclonalcna != None:
        mutations['subclonalfocal'] = {(sp.basesize(event.split(':')[0]), sp.basesize(event.split(':')[1])) if len(event.split(':')) == 3 else (sp.basesize(event.split(':')[0]), int(sp.basesize(event.split(':')[0]) * 0.2)) : int(event.split(':')[2]) if len(event.split(':')) == 3 else int(event.split(':')[1]) for event in args.subclonalcna.strip().split()}

    if args.snplist == None and args.snpratio == None:
        raise ValueError(sp.error('One among SNP list and SNP ratio must be provided!'))

    ignorelist = []
    if args.ignore != None:
        with open(args.ignore) as f:
            for line in f:
                if line != '':
                    ignorelist.append(line.strip())

    binsize = sp.basesize(args.binsize)

    return {'reference' : args.REFERENCE,
            'numclones' : args.numclones,
            'ignore' : ignorelist,
            'snplist' : args.snplist,
            'snpratio' : args.snpratio,
            'HEHOratio' : args.hehoratio,
            'adratio' : args.adratio,
            'xdir' : args.runningdirectory,
            'rndseed' : args.rndseed,
            'mutations' : mutations,
            'binsize' : binsize,
            'jobs' : args.jobs,
            'noverbose' : args.noverbose}
