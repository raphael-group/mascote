import random

import Genomics
import Evolution
import Support


def simulateEvolution(numclones, humanGenome, binsize, mutations):
    # clonalwgd, clonalwcl, clonalcam, clonalfocal, subclonalwgd, subclonalwcl, subclonalcam, subclonalfocal = mutations
    evolution = Evolution.RandomTree(n=numclones, humanGenome=humanGenome, binsize=binsize)
    sublocations = locateSubclonal(clones=evolution.clones, mutations=mutations)
    mutate(clone=evolution.root, mutations=mutations, sublocations=sublocations)
    return evolution


def locateSubclonal(clones, mutations):
    cidx = [clone.idx for clone in clones if clone.parent != None]
    sublocations = {}
    if len(cidx) > 0:
        sublocations["subclonalwgd"] = list(random.sample(cidx, mutations["subclonalwgd"]))
        sublocations["subclonalwcl"] = [random.choice(cidx) for i in range(mutations["subclonalwcl"])]
        sublocations["subclonalcam"] = [random.choice(cidx) for i in range(mutations["subclonalcam"])]
    return sublocations


def mutate(clone, mutations, sublocations):
    if clone.parent == None:
        mutateClonal(clone=clone, wgd=mutations["clonalwgd"], wcl=mutations["clonalwcl"], cam=mutations["clonalcam"], focal=mutations["clonalfocal"], ratioAD=mutations["ratioAD"])
    else:
        clone.inherit(clone.parent)
        mutateSubclonal(clone=clone, sublocations=sublocations, focal=mutations["subclonalfocal"], ratioAD=mutations["ratioAD"])
    for child in clone.children:
        assert(child.parent is clone)
        # child.inherit(clone)
        mutate(clone=child, mutations=mutations, sublocations=sublocations)


def mutateClonal(clone, wgd, wcl, cam, focal, ratioAD):
    timing = [(0,) for i in xrange(wgd)] + [(1,) for i in xrange(wcl)] + [(2,) for i in xrange(cam)]
    for size in focal:
        timing += [(3, size[0], size[1]) for i in xrange(focal[size])]
    random.shuffle(timing)
    for event in timing:
        mutation(clone=clone, event=event, ratioAD=ratioAD)


def mutateSubclonal(clone, sublocations, focal, ratioAD):
    timing = [(0,) for i in sublocations["subclonalwgd"] if clone.idx == i]
    timing += [(1,) for i in sublocations["subclonalwcl"] if clone.idx == i]
    timing += [(2,) for i in sublocations["subclonalcam"] if clone.idx == i]
    for size in focal:
        timing += [(3, size[0], size[1]) for i in range(focal[size])]
    random.shuffle(timing)
    for event in timing:
        mutation(clone=clone, event=event, ratioAD=ratioAD)


def mutation(clone, event, ratioAD):
    if event[0] == 0:
        ### Whole-genome duplication (WGD)
        clone.wgd()
    else:
        actives = [(clone.genome[c].name, 'm', clone.genome[c].maternalHaplotypeLength) for c in clone.genome if clone.genome[c].maternalHaplotypeLength > 0]
        actives += [(clone.genome[c].name, 'p', clone.genome[c].paternalHaplotypeLength) for c in clone.genome if clone.genome[c].paternalHaplotypeLength > 0]
        if len(actives) > 0:
            if event[0] == 1:
                ### Whole-chromosome loss (WCL)
                target = random.choice([(e[0], e[1]) for e in actives])
                clone.wcl(chromosome=target[0], allele=target[1])
            elif event[0] == 2:
                ### Chromosome-arm mutation (CAM)
                target = random.choice([(e[0], e[1]) for e in actives])
                breakpoint = int(clone.genome[target[0]].maternalHaplotypeLength / 2.0) if target[1] == 'm' else int(clone.genome[target[0]].paternalHaplotypeLength / 2.0)
                if random.random() < 0.5:
                    # Left arm
                    start = 0
                    size = breakpoint
                else:
                    # Right arm
                    size = clone.genome[target[0]].maternalHaplotypeLength if target[1] == 'm' else clone.genome[target[0]].paternalHaplotypeLength
                    size = max(breakpoint, size-breakpoint)
                    start = breakpoint
                if random.random() < ratioAD:
                    # Arm duplication
                    clone.tandemDuplicate(chromosome=target[0], start=start, size=size, allele=target[1], arm=True)
                else:
                    # Arm deletion
                    clone.delete(chromosome=target[0], start=start, size=size, allele=target[1], arm=True)
            elif event[0] == 3:
                ### Focal copy-number aberration
                target = Support.wchoice([(e[0], e[1]) for e in actives], [e[2] for e in actives])
                size = max(int(round(random.gauss(event[1], event[2]))), 0)
                size = min(size, clone.genome[target[0]].maternalHaplotypeLength if target[1] == 'm' else clone.genome[target[0]].paternalHaplotypeLength)
                breakpoint = random.randint(0, max((clone.genome[target[0]].maternalHaplotypeLength-1-size if target[1] == 'm' else clone.genome[target[0]].paternalHaplotypeLength-1-size), 0))
                assert(breakpoint < clone.genome[target[0]].maternalHaplotypeLength if target[1] == 'm' else clone.genome[target[0]].paternalHaplotypeLength)
                if random.random() < ratioAD:
                    # Focal duplication
                    clone.tandemDuplicate(chromosome=target[0], start=breakpoint, size=size, allele=target[1])
                else:
                    # Focal deletion
                    clone.delete(chromosome=target[0], start=breakpoint, size=size, allele=target[1])
            else:
                assert(False)
