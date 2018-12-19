import math
import copy
import random
from collections import Counter

import Support



class HumanGenome:

    def __init__(self, reference, snplist, snpratio, HEHOratio, ignorelist):
        self.reference = reference
        self.snplist = Support.parseSNPList(snplist)
        self.snpratio = snpratio
        self.HEHOratio = HEHOratio
        self.ignorelist = ignorelist
        self.chromosomes = []
        self.lengths = {}
        self.maternalfa = None
        self.paternalfa = None
        self.numsnps = 0
        self.hetsnps = 0

    def buildGenome(self, maternalout, paternalout):
        self.maternalfa = maternalout
        self.paternalfa = paternalout
        with open(self.maternalfa, 'w') as mafa:
            with open(self.paternalfa, 'w') as pafa:
                name = ''
                sequence = []
                with open(self.reference, 'r') as ref:
                    for line in ref:
                        if line != '':
                            if line[0] == '>':
                                if name != '' and not name in self.ignorelist:
                                    mafa.write(">{}\n".format(name))
                                    pafa.write(">{}\n".format(name))
                                    maternalhap, paternalhap, length = self.buildHaplotypes(chromosome=name, sequence=sequence)
                                    mafa.write("{}\n".format(''.join(maternalhap)))
                                    pafa.write("{}\n".format(''.join(paternalhap)))
                                    self.chromosomes.append(name)
                                    self.lengths[name] = length
                                name = line.strip()[1:]
                                sequence = []
                            else:
                                if not name is self.ignorelist:
                                    sequence += list(line.strip())
                    if len(sequence) > 0:
                        if name != '' and not name in self.ignorelist:
                            mafa.write(">{}\n".format(name))
                            pafa.write(">{}\n".format(name))
                            maternalhap, paternalhap, length = self.buildHaplotypes(chromosome=name, sequence=sequence)
                            mafa.write("{}\n".format(''.join(maternalhap)))
                            pafa.write("{}\n".format(''.join(paternalhap)))
                            self.chromosomes.append(name)
                            self.lengths[name] = length


    def buildHaplotypes(self, chromosome, sequence):
        mhap = copy.deepcopy(sequence)# if s != 'N' and s != 'n' and s != '*']
        phap = copy.deepcopy(sequence)
        lenchr = len(mhap)
        if self.snplist == None:
            assert(0 <= self.snpratio <= 1.0)
            snplist = {snp : set(random.sample(['A','T','C','G'], 2)) for snp in random.sample(xrange(len(mhap)), int(round(len(mhap) * self.snpratio)))}
        else:
            snplist = self.snplist[chromosome]

        for snp in snplist:
            if snp < lenchr:
                self.numsnps += 1
                if random.random() < self.HEHOratio:
                    # Heterozygous SNP
                    alleles = random.sample(snplist[snp], 2)
                    mhap[snp] = alleles[0].lower() if mhap[snp].islower() else alleles[0].upper()
                    phap[snp] = alleles[1].lower() if phap[snp].islower() else alleles[1].upper()
                    self.hetsnps += 1
                else:
                    # Homozygous SNP
                    allele = random.choice(list(snplist[snp]))
                    mhap[snp] = allele.lower() if mhap[snp].islower() else allele.upper()
                    phap[snp] = allele.lower() if phap[snp].islower() else allele.upper()

        return mhap, phap, lenchr


class Clone:

    def __init__(self, idx, humanGenome, binsize, label=None):
        self.idx = idx
        if label == None:
            self.label = str(idx)
        else:
            self.label = label
        assert(humanGenome.maternalfa != None and humanGenome.paternalfa != None)
        self.humanGenome = humanGenome
        self.chromosomes = humanGenome.chromosomes
        self.genome = {c : Chromosome(name=c, length=humanGenome.lengths[c], binsize=binsize) for c in self.chromosomes}
        self.parent = None
        self.children = []
        self.mutationLabels = []

    def inherit(self, parent):
        assert(isinstance(parent, Clone))
        assert(self.humanGenome is parent.humanGenome)
        self.genome = {copy.deepcopy(c) : copy.deepcopy(parent.genome[c]) for c in parent.genome}

    def genomeLength(self):
        return sum(self.genome[c].maternalHaplotypeLength + self.genome[c].paternalHaplotypeLength for c in self.genome)

    def reference(self):
        return {chro : copy.deepcopy(self.genome[chro].reference) for chro in self.chromosomes}

    def copyNumberProfile(self):
        return {chro : {'m' : dict(Counter(self.genome[chro].maternalHaplotype)), 'p' : dict(Counter(self.genome[chro].paternalHaplotype))} for chro in self.chromosomes}

    def wgd(self):
        for c in self.genome:
            if self.genome[c].maternalHaplotypeLength > 0:
                self.genome[c].tandemDuplicate(start=0, size=self.genome[c].maternalHaplotypeLength, allele='m')
            if self.genome[c].paternalHaplotypeLength > 0:
                self.genome[c].tandemDuplicate(start=0, size=self.genome[c].paternalHaplotypeLength, allele='p')
        self.mutationLabels.append("WGD")

    def wcl(self, chromosome, allele):
        size = self.genome[chromosome].maternalHaplotypeLength if allele == 'm' else self.genome[chromosome].paternalHaplotypeLength
        self.genome[chromosome].delete(start=0, size=size, allele=allele)
        self.mutationLabels.append("{}-{} loss".format(allele.upper(), chromosome))

    def tandemDuplicate(self, chromosome, start, size, allele, arm=False):
        self.genome[chromosome].tandemDuplicate(start=start, size=size, allele=allele)
        if not arm:
            self.mutationLabels.append("({},{}) tdup in {}-{}".format(start, start+size, allele.upper(), chromosome))
        else:
            self.mutationLabels.append("({},{}) dup of {}-{} arm".format(start, start+size, allele.upper(), chromosome))

    def delete(self, chromosome, start, size, allele, arm=False):
        self.genome[chromosome].delete(start=start, size=size, allele=allele)
        if not arm:
            self.mutationLabels.append("({},{}) del in {}-{}".format(start, start+size, allele.upper(), chromosome))
        else:
            self.mutationLabels.append("({},{}) del of {}-{} arm".format(start, start+size, allele.upper(), chromosome))

    def buildGenome(self, maternaloutput, paternaloutput):

        def buildChromosome(sequence, haplotype, reference):
            referencesequence = {bi : sequence[bi[0]:bi[1]] for bi in reference}
            result = ''
            for bi in haplotype:
                result += ''.join(referencesequence[bi])
            return result

        with open(self.humanGenome.maternalfa, 'r') as mafa:
            with open(maternaloutput, 'w') as maout:
                name = ''
                sequence = []
                for line in mafa:
                    if line != '':
                        if line[0] == '>':
                            if name != '' and name in self.chromosomes and self.genome[name].maternalHaplotypeLength > 0:
                                maout.write(">{}\n".format(name))
                                assert(len(sequence) == self.genome[name].length)
                                maout.write("{}\n".format(buildChromosome(sequence, self.genome[name].maternalHaplotype, self.genome[name].reference)))
                            name = line.strip()[1:]
                            sequence = []
                        else:
                            if name in self.chromosomes and self.genome[name].maternalHaplotypeLength > 0:
                                sequence += list(line.strip())
                if len(sequence) > 0:
                    if name != '' and name in self.chromosomes:
                        maout.write(">{}\n".format(name))
                        assert(len(sequence) == self.genome[name].length)
                        maout.write("{}\n".format(buildChromosome(sequence, self.genome[name].maternalHaplotype, self.genome[name].reference)))

        with open(self.humanGenome.paternalfa, 'r') as pafa:
            with open(paternaloutput, 'w') as paout:
                name = ''
                sequence = []
                for line in pafa:
                    if line != '':
                        if line[0] == '>':
                            if name != '' and name in self.chromosomes and self.genome[name].paternalHaplotypeLength > 0:
                                paout.write(">{}\n".format(name))
                                assert(len(sequence) == self.genome[name].length)
                                paout.write("{}\n".format(buildChromosome(sequence, self.genome[name].paternalHaplotype, self.genome[name].reference)))
                            name = line.strip()[1:]
                            sequence = []
                        else:
                            if name in self.chromosomes:
                                sequence += list(line.strip())
                if len(sequence) > 0:
                    if name != '' and name in self.chromosomes and self.genome[name].paternalHaplotypeLength > 0:
                        paout.write(">{}\n".format(name))
                        assert(len(sequence) == self.genome[name].length)
                        paout.write("{}\n".format(buildChromosome(sequence, self.genome[name].paternalHaplotype, self.genome[name].reference)))



class Chromosome:

    def __init__(self, name, length, binsize):
        self.name = name
        self.length = length
        self.binsize = binsize
        self.reference = Support.bins(size=self.length, step=self.binsize)
        self.maternalHaplotype = copy.deepcopy(self.reference)
        self.maternalHaplotypeLength = copy.deepcopy(self.length)
        self.paternalHaplotype = copy.deepcopy(self.reference)
        self.paternalHaplotypeLength = copy.deepcopy(self.length)

    def tandemDuplicate(self, start, size, allele):
        if allele.lower() == "m":
            if start > self.maternalHaplotypeLength:
                raise ValueError("The starting point of a mutation cannot be over maternal-haplotype-chromosome length!")
            else:
                startbin, endbin = Support.startendbins(bins=self.maternalHaplotype, start=start, size=size)
                self.maternalHaplotype = self.maternalHaplotype[0:startbin] + self.maternalHaplotype[startbin:endbin] + self.maternalHaplotype[startbin:endbin] + self.maternalHaplotype[endbin:]
                self.maternalHaplotypeLength = sum(seg[1] - seg[0] for seg in self.maternalHaplotype)
        elif allele.lower() == "p":
            if start > self.paternalHaplotypeLength:
                raise ValueError("The starting point of a mutation cannot be over paternal-haplotype-chromosome length!")
            else:
                startbin, endbin = Support.startendbins(bins=self.paternalHaplotype, start=start, size=size)
                self.paternalHaplotype = self.paternalHaplotype[0:startbin] + self.paternalHaplotype[startbin:endbin] + self.paternalHaplotype[startbin:endbin] + self.paternalHaplotype[endbin:]
                self.paternalHaplotypeLength = sum(seg[1] - seg[0] for seg in self.paternalHaplotype)
        else:
            raise ValueError("The specified allele should be equal to either M or P (non-case sensitive)!")

    def delete(self, start, size, allele):
        if allele.lower() == "m":
            if start > self.maternalHaplotypeLength:
                raise ValueError("The starting point of a mutation cannot be over maternal-haplotype-chromosome length!")
            else:
                startbin, endbin = Support.startendbins(bins=self.maternalHaplotype, start=start, size=size)
                self.maternalHaplotype = self.maternalHaplotype[0:startbin] + self.maternalHaplotype[endbin:]
                self.maternalHaplotypeLength = sum(seg[1] - seg[0] for seg in self.maternalHaplotype)
        elif allele.lower() == "p":
            if start > self.paternalHaplotypeLength:
                raise ValueError("The starting point of a mutation cannot be over paternal-haplotype-chromosome length!")
            else:
                startbin, endbin = Support.startendbins(bins=self.paternalHaplotype, start=start, size=size)
                self.paternalHaplotype = self.paternalHaplotype[0:startbin] + self.paternalHaplotype[endbin:]
                self.paternalHaplotypeLength = sum(seg[1] - seg[0] for seg in self.paternalHaplotype)
        else:
            raise ValueError("The specified allele should be equal to either M or P (non-case sensitive)!")
