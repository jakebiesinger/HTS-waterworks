"""
Utility functions for sampling
"""

import re
import optparse
import sys
import os
import cPickle as pickle
import scipy
import random
import itertools
from pygr import worldbase

from hts_waterworks.utils.common import (nSeq, repeatSeq, reverseComplement,
                                     pygrSeqToBed)
import hts_waterworks.utils.sge as sge

def main(argv=None):
    """
    Sample from a given genome or annotationDB
    """
    usage = "%prog [options] output.fasta \n" + main.__doc__
    parser = optparse.OptionParser(usage)
    parser.add_option('--genome', '-g', dest='sample_genome', type='string', default=None,
                      help="""sample from the given genome""")
    parser.add_option('--sample_resource', '-r', dest='sample_resource', type='string', default=None,
                      help='sample from the given resource or bed file')
    parser.add_option('--sample_length', '-l', dest='sample_length', type='int', default=500,
                      help='size of sequence samples, default=%default')
    parser.add_option('--num_samples', '-n', dest='num_samples', type='int', default=10000,
                      help='number of samples to generate')
    parser.add_option('--output_bed', '-b', dest='out_bed_file', type='string',  default='',
                      help='Generate a BED file with the genomic coordinates of sampled regions')
    parser.add_option('--no_fasta', dest='no_fasta', action='store_true',
                      help='Forego generating a fasta file for the samples')
    parser.add_option('--parallel_jobs', '-j', dest='num_jobs', type='int', default=1,
                      help='Use num_jobs to generate the sample, concatenating the sequences at the end')
    parser.add_option('--no_repeats', dest='no_repeats', action='store_true',
                      help='Exclude any repeat sequence (lower case nucleotides) from samples.')
    if argv is None:
        argv = sys.argv[1:]
    opts, args = parser.parse_args(argv)
    if len(args) != 1 or not (opts.sample_genome or opts.sample_resource):
        parser.print_help()
        print 'Please specify an output fasta file!'
        sys.exit(-1)
    
    outfileDir, outfileName = os.path.split(args[0])
    codeDir = os.path.abspath(os.path.dirname(sys.argv[0]))
    
    if opts.num_jobs > 1:
        samplesPerJob = opts.num_samples / opts.num_jobs
        print 'Submitting %s sampling jobs of %s samples each...' % (opts.num_jobs, samplesPerJob)
        cmd = '%s %s/sampling.py %s.$SGE_TASK_ID ' % (sge.python_cmd, codeDir, args[0])
        cmd += '--sample_length=%s ' % opts.sample_length
        if opts.sample_genome:
            cmd += '--sample_genome=%s ' % opts.sample_genome
        else:
            cmd += '--sample_resource=%s ' % opts.sample_resource
        if opts.no_repeats:
            cmd += '--no_repeats '
        if opts.no_fasta:
            cmd += '--no_fasta '
        cmd += '--num_samples=$num_samples '
        sampleSizes = [str(samplesPerJob)] * opts.num_jobs + [str(opts.num_samples - samplesPerJob * opts.num_samples)]
        sampleJobs = sge.JobGroup('sample_for_%s' % outfileName, cmd,
                                  arguments={'num_samples':sampleSizes})
        concatJob = sge.Job('sample_for_%s_concat' % outfileName,
                            'cat %s.* > %s' % (args[0], args[0]))
        concatJob.addDependency(sampleJobs)
        sge.build_submission(outfileDir, [sampleJobs, concatJob])
        concatJob.wait()
        
    else:
        
        if opts.sample_genome:
            genome = worldbase(opts.sample_genome)
            sample_gen = sample_genome(genome, [opts.sample_length], sampleSize=opts.num_samples, excludeRepeat=opts.no_repeats)
        else:  # opts.sample_resource:
            res1Map = worldbase(res1Name)
            sample_gen = sample_resource(annotDB, [opts.sample_length], sampleSize=opts.num_samples, excludeRepeat=opts.no_repeats)
            
            
        print '# Generating sequence samples and writing to disk...'
        if not opts.no_fasta:
            outfile = open(args[0], 'w')
        if opts.out_bed_file != '':
            bedOutfile = open(opts.out_bed_file, 'w')
        for index, seq in enumerate(sample_gen):
            if not opts.no_fasta:
                outfile.write('>sample_%s\n%s\n' % (index, seq))
            if opts.out_bed_file != '':
                bedOutfile.write(pygrSeqToBed(seq, name='sample_%s'%index) + '\n')    
        if opts.out_bed_file != '':
            bedOutfile.close()
        if not opts.no_fasta:
            outfile.close()
        
        print '# Sampling complete!'



def sample_genome(genome, sizeDistribution, sampleSize=100000, excludeRepeat=True,
                  excludeN=True, ignoreCharacters='', weighted=True):
    """ Sample the given genome, selecting sample sizes from size distribution
        Can optionally ignore repeat sequences (lower-case) and N's.
        weighted indicates if sampling should be weighted by chromosome length"""
    if weighted:
        # Sample chromosomes according to their length
        weightedChroms = [(len(genome[chrom]), chrom) for chrom in genome ]
                                    #if not re.search(ignoreCharacters, chrom)]
        chroms = weighted_sample(weightedChroms, sampleSize)
    else:
        # equal weights on all chromosomes
        chroms = [chrom for chrom in genome.keys() ]
                        #if not re.search(ignoreCharacters, chrom)]
        chroms = sample_with_replacement(chroms, sampleSize)
    samplesLeft = sampleSize
    for sampleIndex in xrange(sampleSize):
        #if sampleIndex % 100 == 0:
        #    print 'sample', sampleIndex
        sampleLength = random.sample(sizeDistribution,1)[0]
        chrom = chroms.next()
        badSample = True
        t = 0
        while badSample:
            #if t % 100 == 0:
            #    print 'try', t
            t += 1
            #start = random.sample(genome[chrom],1)[0].start
            biggest = len(genome[chrom]) - sampleLength
            if biggest < 2:
                print 'bug: len(chrom):%s, sampleLength=%s' % (len(genome[chrom]), sampleLength)
                sampleLength = random.sample(sizeDistribution,1)[0]
                chrom = chroms.next()
                continue
            start = random.randrange(biggest)
            stop = start + sampleLength
            region = genome[chrom][start:stop]
            seq = str(region)
            if len(seq) == sampleLength \
                and (not excludeRepeat or not re.search(repeatSeq, seq)) \
                and (not excludeN or not re.search(nSeq, seq)):
                    badSample = False
                    yield region

def sample_middles(genome, sizeDistribution, bedLinesForMiddle, sampleSize=10000):
    """sample the genome using bedLinesForMiddle as the middle of the samples"""
    middles = []
    bedLinesForMiddle.seek(0)
    lineCount = sum(1 for l in bedLinesForMiddle)
    bedLinesForMiddle.seek(0)
    def getMiddles(bedlines):
        for line in bedlines:
            fields = line.strip().split('\t')
            start, stop = map(int, fields[1:3])
            yield (fields[0], (start + stop) / 2)  # save chrom and middle
    print 'done getting middles, now sampling'
    sys.stdout.flush()
    # randomly sample middles from bedLines and sizes from sizeDistribution
    sizes = sample_with_replacement(sizeDistribution, sampleSize)
    for chrom, mid in sample_with_replacement(getMiddles(bedLinesForMiddle), sampleSize, lineCount):
        nextSize = sizes.next()
        newStart = max(0, mid - nextSize)
        newStop = min(len(genome[chrom]), mid + nextSize)
        yield genome[chrom][newStart:newStop]

def sample_resource(annotDB, sizeDistribution, sampleSize=100000,
                    excludeRepeat=False, excludeN=False, weighted=False):
    """ Sample segments of annotDB, """
    genome = annotDB.itervalues().next().sequence.db
    if weighted:
        # Sample indexes according to their length
        weightedIndexes = [(len(annotDB[id]), id) for id in annotDB]
        keys = weighted_sample(weightedIndexes, sampleSize)
    else:
        # equal weights on all indexes
        keys = sample_with_replacement(annotDB.keys(), sampleSize)
    allSampleLengths = sample_with_replacement(sizeDistribution, sampleSize)
    for sampleIndex in xrange(sampleSize):
        sampleLength = allSampleLengths.next()
        badSample = True
        while badSample:
            try:
                resourceRegionId = keys.next()
                resourceRegion = annotDB[resourceRegionId].sequence
                oldOrientation = resourceRegion.orientation
                if resourceRegion.orientation == -1:
                    # random range is easier if assume its on the positive strand
                    resourceRegion  = -resourceRegion
                # sample must overlap region by at least 1 bp, so sampleStart could be
                # before region.start and sampleStop could be after region.stop
                sampleStart = random.randrange(max(0,abs(resourceRegion.start) - sampleLength), resourceRegion.stop)
                sampleSeq = genome[resourceRegion.id][sampleStart:sampleStart+sampleLength]
                if oldOrientation == -1:
                    sampleSeq = -sampleSeq
                if len(sampleSeq) == sampleLength \
                    and (not excludeRepeat or not re.search(repeatSeq, seq)) \
                    and (not excludeN or not re.search(nSeq, seq)):
                        badSample = False
                        yield sampleSeq
            except StopIteration, e:
                # ran out of index samples
                if weighted:
                    keys = weighted_sample(weightedIndexes, sampleSize - sampleIndex)
                else:
                    keys = sample_with_replacement(annotDB.keys(), sampleSize - sampleIndex)


#def calc_pwm_bg_dist(matrix, genome, keepScores=False, sampleSize=100000, samples=None, **kwargs):
def calc_pwm_bg_dist(matrix, genome, keepScores=False, sampleSize=100000, samples=None, **kwargs):
    """ Sample the given genome to approximate the background distribution for the given matrix
            returns the (mean,stdev) of the distribution
        The function takes about 40 seconds for 100,000 samples or 6 minutes for 1 million

        matrix should be a motility matrix
                (import motility; matrix=motility.PWM([[1,2,3,4],[5,6,7,8]])
        genome should be a pygr resource genome
                (from pygr import worldbase; genome=worldbase.Bio.Seq.Genome.MOUSE.mm9())
        returns: mu, sd [, allScores]
    """
    if not samples:
        sizeDistribution = [len(matrix)]
        allScores = [0.0] * sampleSize
        if len(kwargs) > 0:
            samples = sample_genome(genome, sizeDistribution, sampleSize, **kwargs)
        else:
            samples = sample_genome(genome, sizeDistribution, sampleSize)
    else:
        allScores = [0.0] * len(samples)
    for index, seq in enumerate(samples):
        seq = str(seq)
        if len(seq) < len(matrix):
            raise RuntimeError("seq is too short to be scored! seq: %s, matrix: %sbp %s" % (seq, len(matrix), matrix.matrix))
        revcompSeq = reverseComplement(seq)
        allScores[index] = max(matrix.calc_score(seq), matrix.calc_score(revcompSeq))
    mu = scipy.mean(allScores)
    sd = scipy.std(allScores)
    if keepScores:
        return mu, sd, allScores
    else:
        return mu,sd


def sample_with_replacement(population, k, popsize=None):
    "Chooses k random elements (with replacement) from a population"
    _random, _int = random.random, int  # speed hack
    if popsize is not None and k < popsize:
        # select indices first, then sample them treating population as an iterator
        indices = set([_int(_random() * popsize) for i in itertools.repeat(None, k)])
        return itertools.compress(population, (i in indices for i in xrange(popsize)))
    else:
        population = list(population)
        popsize = len(population)
        return (population[_int(_random() * popsize)] for i in itertools.repeat(None, k))


def weighted_sample(items, n):
    """ weighted sample with replacement from stack overflow
    input: list of (weight, item) pairs
http://stackoverflow.com/questions/2140787/select-random-k-elements-from-a-list-whose-elements-have-weights/2149533
    """
    items = list(items)
    total = float(sum(w for w, v in items))
    i = 0
    w, v = items[0]
    while n:
        x = total * (1 - random.random() ** (1.0 / n))
        total -= x
        while x > w:
            x -= w
            i += 1
            w, v = items[i]
        w -= x
        yield v
        n -= 1


class Node:
    # Each node in the heap has a weight, value, and total weight.
    # The total weight, self.tw, is self.w plus the weight of any children.
    __slots__ = ['w', 'v', 'tw']
    def __init__(self, w, v, tw):
        self.w, self.v, self.tw = w, v, tw

def rws_heap(items):
    # h is the heap. It's like a binary tree that lives in an array.
    # It has a Node for each pair in `items`. h[1] is the root. Each
    # other Node h[i] has a parent at h[i>>1]. Each node has up to 2
    # children, h[i<<1] and h[(i<<1)+1].  To get this nice simple
    # arithmetic, we have to leave h[0] vacant.
    h = [None]                          # leave h[0] vacant
    for w, v in items:
        h.append(Node(w, v, w))
    for i in range(len(h) - 1, 1, -1):  # total up the tws
        h[i>>1].tw += h[i].tw           # add h[i]'s total to its parent
    return h

def rws_heap_pop(h):
    gas = h[1].tw * random.random()     # start with a random amount of gas

    i = 1                     # start driving at the root
    while gas > h[i].w:       # while we have enough gas to get past node i:
        gas -= h[i].w         #   drive past node i
        i <<= 1               #   move to first child
        if gas > h[i].tw:     #   if we have enough gas:
            gas -= h[i].tw    #     drive past first child and descendants
            i += 1            #     move to second child
    w = h[i].w                # out of gas! h[i] is the selected node.
    v = h[i].v

    h[i].w = 0                # make sure this node isn't chosen again
    while i:                  # fix up total weights
        h[i].tw -= w
        i >>= 1
    return v

def weighted_sample_no_replacement(items, n):
    """Random sampling without replacement:
http://stackoverflow.com/questions/2140787/select-random-k-elements-from-a-list-whose-elements-have-weights/2149533"""
    heap = rws_heap(items)              # just make a heap...
    for i in range(n):
        yield rws_heap_pop(heap)        # and pop n items off it.



if __name__ == '__main__':
    main()
