#!/usr/bin/python2

import sys
import os.path
import argparse
import copy
import shlex
import subprocess
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value
import random

import Genomics
import Evolution
import Mutation
import Support
import Argparser
import Builder


def main():
    Support.log(msg="# Parsing and checking the input arguments\n", level="STEP")
    args = Argparser.parse_mascotte_arguments()
    logArgs(args)
    if args['rndseed'] != None:
        random.seed(args['rndseed'])

    Support.log(msg="# Setting up for simulating human diploid genome\n", level="STEP")
    human = Genomics.HumanGenome(reference=args['reference'], snplist=args['snplist'], snpratio=args['snpratio'], HEHOratio=args['HEHOratio'], ignorelist=args['ignore'])
    maternalhuman = os.path.join(args['xdir'], 'human.maternal.fa')
    paternalhuman = os.path.join(args['xdir'], 'human.paternal.fa')
    Support.log(msg="# Simulating human diploid genome\n", level="STEP")
    human.buildGenome(maternalout=maternalhuman, paternalout=paternalhuman)
    Support.log('Chromosomes: {}\n'.format(', '.join(human.chromosomes)), level='INFO')
    Support.log('Number of simulated SNPs: {}\n'.format(human.numsnps), level='INFO')
    Support.log('Number of heterozygous SNPs: {}\n'.format(human.hetsnps), level='INFO')
    Support.log('Maternal chromosome of human genome written in {}\n'.format(maternalhuman), level='INFO')
    Support.log('Maternal chromosome of human genome written in {}\n'.format(paternalhuman), level='INFO')

    if args['numclones'] > 0:
        Support.log(msg="# Simulating tumor clones and their evolution through specified CNAs\n", level="STEP")
        tumor = Mutation.simulateEvolution(numclones=args['numclones'], humanGenome=human, binsize=args['binsize'], mutations=args['mutations'])
        Support.log('Simulated tumor clones: {}\n'.format(', '.join([clone.label for clone in tumor.clones])), level='INFO')
        Support.log('Founder tumor clone: {}\n'.format(tumor.root.label), level='INFO')
        with open(os.path.join(args['xdir'], 'tumor.dot'), 'w') as o: o.write("{}\n".format(tumor.draw()))
        Support.log('The resulting tumor evolution of clones and related CNAs have been drawn in {} as dot format\n'.format(os.path.join(args['xdir'], 'tumor.dot')), level='INFO')
        Support.log('Genome length of the various tumor clones:\n\t{}\n'.format('\n\t'.join(['{}: {}'.format(clone.label, clone.genomeLength()) for clone in tumor.clones])), level='INFO')
        Support.log('Computing and segmenting the copy-number profiles jointly for all tumor clones\n', level='INFO')
        segments = segmentation(evolution=tumor)
        Support.log('Total number of resulting segments= {}\n'.format(sum(len(segments[chro]) for chro in tumor.human.chromosomes)), level='INFO')
        segout = os.path.join(args['xdir'], 'copynumbers.csv')
        with open(segout, 'w') as o:
            o.write('\t'.join(['#CHR', 'START', 'END'] + [clone.label for clone in tumor.clones]) + '\n')
            o.write('\n'.join(['\t'.join(map(str, [chro, seg[0], seg[1]]+['{}|{}'.format(segments[chro][seg][clone.idx]['m'], segments[chro][seg][clone.idx]['p']) for clone in tumor.clones])) for chro in tumor.human.chromosomes for seg in sorted(segments[chro], key=(lambda x : x[0]))]))
            o.write('\n')
        Support.log('The allele-specific copy number profiles for every tumor clone has been written in {}\n'.format(segout), level='INFO')
        Support.log('Writing the FASTA-format genomes of tumor clones\n', level='INFO')
        if args['jobs'] == 1:
            for clone in tumor.clones:
                maternalout = os.path.join(args['xdir'], '{}.maternal.fa'.format(clone.label))
                paternalout = os.path.join(args['xdir'], '{}.paternal.fa'.format(clone.label))
                clone.buildGenome(maternalout, paternalout)
        else:
            builder = Builder.CloneGenomeBuilder(tumor, args['xdir'])
            builder.parallelbuild(args['jobs'])
        Support.log('Tumor-clone genomes wrote in:\n{}\n'.format('\n'.join(['\t{}: maternal > {} and paternal > {}'.format(clone.label, os.path.join(args['xdir'], '{}.maternal.fa'.format(clone.label)), os.path.join(args['xdir'], '{}.paternal.fa'.format(clone.label))) for clone in tumor.clones])), level='INFO')
    else:
        Support.log(msg="# No tumor clones will be generated as input tumor clones is 0\n", level="INFO")
    Support.log('KTHXBY!\n', level='STEP')


def segmentation(evolution):
    profiles = {clone.idx : clone.copyNumberProfile() for clone in evolution.clones}
    merge = {chro : {seg : {clone.idx : {'m' : profiles[clone.idx][chro]['m'][seg] if seg in profiles[clone.idx][chro]['m'] else 0, 'p' : profiles[clone.idx][chro]['p'][seg] if seg in profiles[clone.idx][chro]['p'] else 0} for clone in evolution.clones} for seg in evolution.root.genome[chro].reference} for chro in evolution.human.chromosomes}
    segments = {chro : {} for chro in evolution.human.chromosomes}
    for chro in evolution.human.chromosomes:
        prev = sorted(evolution.root.genome[chro].reference, key=(lambda x : x[0]))[:-1]
        succ = sorted(evolution.root.genome[chro].reference, key=(lambda x : x[0]))[1:]
        start = prev[0][0]
        for p, s in zip(prev, succ):
            flag = True
            for clone in evolution.clones:
                if merge[chro][p][clone.idx]['m'] != merge[chro][s][clone.idx]['m'] or merge[chro][p][clone.idx]['p'] != merge[chro][s][clone.idx]['p']:
                    flag = False
                    break
            if not flag:
                segments[chro][(start, p[1])] = {clone.idx : {'m' : merge[chro][p][clone.idx]['m'], 'p' : merge[chro][p][clone.idx]['p']} for clone in evolution.clones}
                start = s[0]
        segments[chro][(start, succ[-1][1])] = {clone.idx : {'m' : merge[chro][succ[-1]][clone.idx]['m'], 'p' : merge[chro][succ[-1]][clone.idx]['p']} for clone in evolution.clones}
    return segments


def logArgs(args):
    Support.log('\n'.join(["Arguments:"]+['\t{}:\t{}'.format(arg, args[arg]) for arg in args if arg != 'mutations'] + [""]), level="INFO")
    Support.log('\n'.join(["Mutations:"]+['\t{}:\t{}'.format(mut, args['mutations'][mut]) for mut in args['mutations'] if mut != 'clonalfocal' and mut != 'subclonalfocal'] + [""]), level="INFO")
    if 'clonalfocal' in args['mutations']:
        Support.log('\n'.join(['Clonal focal CNAs:']+['\tmean_len= {}, standard_deviation= {}, quantity= {}'.format(key[0], key[1], args['mutations']['clonalfocal'][key]) for key in args['mutations']['clonalfocal']] + [""]), level="INFO")
    if 'subclonalfocal' in args['mutations']:
        Support.log('\n'.join(['Subclonal focal CNAs:']+['\tmean_len= {}, standard_deviation= {}, quantity= {}'.format(key[0], key[1], args['mutations']['subclonalfocal'][key]) for key in args['mutations']['subclonalfocal']] + [""]), level="INFO")



if __name__ == '__main__':
        main()
