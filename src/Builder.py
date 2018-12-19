import os
import copy
from multiprocessing import Process, Queue, JoinableQueue, Lock, Value


class CloneGenomeBuilder:

    def __init__(self, tumor, xdir):
        self.tumor = tumor
        self.xdir = xdir

    def parallelbuild(self, numworkers):

        class Worker(Process):

            def __init__(self, task_queue, results, chromosomes, reference, lengths):
                Process.__init__(self)
                self.task_queue = task_queue
                self.result_queue = results
                self.chromosomes = chromosomes
                self.reference = reference
                self.lengths = lengths

            def run(self):
                while True:
                    next_task = self.task_queue.get()
                    if next_task is None:
                        # Poison pill means shutdown
                        self.task_queue.task_done()
                        break
                    self.buildHaplotype(next_task)
                    self.task_queue.task_done()
                    self.result_queue.put(next_task[-1])
                return

            def buildChromosome(self, sequence, haplotype, reference):
                referencesequence = {bi : sequence[bi[0]:bi[1]] for bi in reference}
                result = ''
                for bi in haplotype:
                    result += ''.join(referencesequence[bi])
                return result

            def buildHaplotype(self, task):
                originalfa, haplotypes, haplengths, output = task
                with open(originalfa, 'r') as fa:
                    with open(output, 'w') as out:
                        name = ''
                        sequence = []
                        for line in fa:
                            if line != '':
                                if line[0] == '>':
                                    if name != '' and name in self.chromosomes and haplengths[name] > 0:
                                        out.write(">{}\n".format(name))
                                        assert(len(sequence) == self.lengths[name])
                                        out.write("{}\n".format(self.buildChromosome(sequence, haplotypes[name], self.reference[name])))
                                    name = line.strip()[1:]
                                    sequence = []
                                else:
                                    if name in self.chromosomes and haplengths[name] > 0:
                                        sequence += list(line.strip())
                        if len(sequence) > 0:
                            if name != '' and name in self.chromosomes:
                                out.write(">{}\n".format(name))
                                assert(len(sequence) == self.lengths[name])
                                out.write("{}\n".format(self.buildChromosome(sequence, haplotypes[name], self.reference[name])))

        # Establish communication queues
        tasks = JoinableQueue()
        results = Queue()

        # Enqueue jobs
        c = copy.deepcopy
        jobs_count = 0
        for clone in self.tumor.clones:
            originalfa = c(clone.humanGenome.maternalfa)
            haplotypes = {c(chro) : c(clone.genome[chro].maternalHaplotype) for chro in clone.humanGenome.chromosomes}
            haplengths = {c(chro) : c(clone.genome[chro].maternalHaplotypeLength) for chro in clone.humanGenome.chromosomes}
            output = os.path.join(self.xdir, '{}.maternal.fa'.format(c(clone.label)))
            tasks.put((originalfa, haplotypes, haplengths, output))
            jobs_count += 1

            originalfa = c(clone.humanGenome.paternalfa)
            haplotypes = {c(chro) : c(clone.genome[chro].paternalHaplotype) for chro in clone.humanGenome.chromosomes}
            haplengths = {c(chro) : c(clone.genome[chro].paternalHaplotypeLength) for chro in clone.humanGenome.chromosomes}
            output = os.path.join(self.xdir, '{}.paternal.fa'.format(c(clone.label)))
            tasks.put((originalfa, haplotypes, haplengths, output))
            jobs_count += 1

        # Setting up the workers
        workers = [Worker(task_queue=tasks, results=results, chromosomes=c(self.tumor.human.chromosomes), reference={c(chro) : c(self.tumor.root.genome[chro].reference) for chro in self.tumor.human.chromosomes}, lengths={c(chro) : c(self.tumor.root.genome[chro].length) for chro in self.tumor.human.chromosomes}) for i in range(min(numworkers, jobs_count))]

        # Add a poison pill for each worker
        for i in range(len(workers)):
            tasks.put(None)

        # Start the workers
        for w in workers:
            w.start()

        # Wait for all of the tasks to finish
        tasks.join()

        # Collect results
        collect = []
        for i in range(jobs_count):
            collect.append(results.get())

        # Close Queues
        tasks.close()
        results.close()

        # Ensure each worker terminates
        for w in workers:
            w.terminate()
            w.join()

        return collect
