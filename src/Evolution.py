
import random

from Genomics import Clone


class RandomTree:

    def __init__(self, n, humanGenome, binsize, labels=None):
        self.n = n
        self.human = humanGenome
        if self.n <= 0:
            raise ValueError("A tree cannot be built with 0 or negative number of nodes!")
        if labels != None:
            self.labels = labels
            if n != len(labels):
                raise ValueError("A labels should be given for every node of the tree!")
        else:
            self.labels = ["clone" + str(i) for i in range(self.n)]
        self.clones =[Clone(idx=i, humanGenome=humanGenome, binsize=binsize, label=self.labels[i]) for i in range(self.n)]
        self.mapclone = {clone.idx : clone for clone in self.clones}
        self.root = None
        self.buildRandom()
        noparent = set(clone for clone in self.clones if clone.parent == None)
        assert(len(noparent) == 1)
        assert(self.root != None)

    def buildRandom(self):
        '''
        Based on Prufer sequence
        '''
        if self.n == 1:
            self.sequence = []
            self.root = self.clones[0]
        elif self.n == 2:
            self.sequence = []
            self.root = self.clones[0]
            self.clones[0].children.append(self.clones[1])
            self.clones[1].parent = self.clones[0]
        elif self.n == 3:
            self.sequence = []
            self.root = self.clones[0]
            if random.random() < 0.5:
                self.clones[0].children.append(self.clones[1])
                self.clones[1].parent = self.clones[0]
                self.clones[1].children.append(self.clones[2])
                self.clones[2].parent = self.clones[1]
            else:
                self.clones[0].children.append(self.clones[1])
                self.clones[1].parent = self.clones[0]
                self.clones[0].children.append(self.clones[2])
                self.clones[2].parent = self.clones[0]
        else:
            self.sequence = [random.choice([clone.idx for clone in self.clones]) for i in range(self.n - 2)]
            #self.sequence.append(0)
            adj = {clone.idx : [] for clone in self.clones}
            degrees = {idx : 1 + self.sequence.count(idx) for idx in [clone.idx for clone in self.clones]}

            for i in self.sequence:
                for clone in self.clones:
                    if degrees[clone.idx] == 1:
                        adj[i].append(clone.idx)
                        adj[clone.idx].append(i)
                        degrees[i] -= 1
                        degrees[clone.idx] -= 1
                        break
            assert(sum(degrees[idx] == 1 for idx in degrees) == 2)

            u = -1
            v = -1
            for clone in self.clones:
                if degrees[clone.idx] == 1:
                    if u == -1:
                        u = clone
                    else:
                        v = clone
                        break
            adj[u.idx].append(v.idx)
            adj[v.idx].append(u.idx)
            degrees[u.idx] -= 1
            degrees[v.idx] -= 1
            assert(len(set(degrees[idx] == 0 for idx in range(self.n))) == 1)

            def orient(clone, adj, placed, mapclone):
                random.shuffle(adj[clone.idx])
                for child in adj[clone.idx]:
                    if not child in placed:
                        clone.children.append(mapclone[child])
                        mapclone[child].parent = clone
                        placed.add(child)
                        orient(mapclone[child], adj, placed, mapclone)

            placed = set()
            #self.root = random.choice(self.clones)
            self.root = self.clones[0]
            assert(self.root != None)
            assert(self.root.parent == None)
            placed.add(self.root.idx)
            orient(self.root, adj, placed, self.mapclone)



    def draw(self):
        colors = ["red", "blue", "purple", "green", "brown", "cadetblue", "chartreuse","cyan", "pink", "Grey", "orange"]

        s = "digraph EvolutionaryCloneTree {\n"
        s += 'splines=true;\nsep="+25,25";\noverlap=scalexy;\nnodesep=0.6;\n'

        s += "\tsubgraph T {\n"
        rank = []
        s += "\t\tN[label=<<B>Normal</B>>,color=black]\n"
        for clone in self.clones:
            s += "\t\t{}[label=<<B>clone</B><SUB>{}</SUB>>,color={}]\n".format(clone.idx, clone.idx, colors[clone.idx%11])
            if len(clone.children) == 0:
                rank.append(clone.idx)

        s += "\t{rank = same; "+ "; ".join(map(str, rank)) + "}\n"
        s += "\t}\n"

        s += '\tN -> {} [label="{}", fontsize=5, fixedsize=true]\n'.format(self.root.idx, "\n".join(self.root.mutationLabels))
        for clone in self.clones:
            for child in clone.children:
                s += '\t{} -> {} [label="{}", fontsize=5, fixedsize=true]\n'.format(clone.idx, child.idx, "\n".join(child.mutationLabels))

        s += "}\n"

        return s
        #proc = subprocess.call("dot -Tpdf " + dotfile + " -o " + outname,shell=True)
